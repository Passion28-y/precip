#include "ellipsoid.h"
#include <cmath>
#include <algorithm>

static inline void matmul3x3(const double* M, const double* v, double* y) {
    y[0]=M[0]*v[0]+M[1]*v[1]+M[2]*v[2];
    y[1]=M[3]*v[0]+M[4]*v[1]+M[5]*v[2];
    y[2]=M[6]*v[0]+M[7]*v[1]+M[8]*v[2];
}

void build_matrices(Ellipsoid& e){
    // U assumed orthonormal. Build A = U Σ² U^T
    double S2[9]={ 1.0/(e.r[0]*e.r[0]),0,0, 0,1.0/(e.r[1]*e.r[1]),0, 0,0,1.0/(e.r[2]*e.r[2]) };
    double US2[9];
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
        US2[3*i+j] = e.U[3*i+0]*S2[3*0+j] + e.U[3*i+1]*S2[3*1+j] + e.U[3*i+2]*S2[3*2+j];
    }
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
        e.A[3*i+j] = US2[3*i+0]*e.U[3*j+0] + US2[3*i+1]*e.U[3*j+1] + US2[3*i+2]*e.U[3*j+2];
    }

    // Conservative AABB by sampling 6 principal directions
    double dirs[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    double xmin=1e300,xmax=-1e300,ymin=1e300,ymax=-1e300,zmin=1e300,zmax=-1e300;
    for(auto& d:dirs){
        double yp[3]={ d[0]*e.r[0], d[1]*e.r[1], d[2]*e.r[2] };
        double xw[3];
        // x = U y + c
        xw[0] = e.U[0]*yp[0] + e.U[1]*yp[1] + e.U[2]*yp[2] + e.c[0];
        xw[1] = e.U[3]*yp[0] + e.U[4]*yp[1] + e.U[5]*yp[2] + e.c[1];
        xw[2] = e.U[6]*yp[0] + e.U[7]*yp[1] + e.U[8]*yp[2] + e.c[2];
        xmin = std::min(xmin, xw[0]); xmax = std::max(xmax, xw[0]);
        ymin = std::min(ymin, xw[1]); ymax = std::max(ymax, xw[1]);
        zmin = std::min(zmin, xw[2]); zmax = std::max(zmax, xw[2]);
    }
    e.aabb = { xmin,xmax, ymin,ymax, zmin,zmax };
}

double implicitF(const Ellipsoid& e, const std::array<double,3>& x, std::array<double,3>* grad){
    double dx[3]={x[0]-e.c[0], x[1]-e.c[1], x[2]-e.c[2]};
    double Axd[3]; matmul3x3(e.A.data(), dx, Axd);
    double F = dx[0]*Axd[0]+dx[1]*Axd[1]+dx[2]*Axd[2] - 1.0;
    if(grad){
        (*grad) = { 2*Axd[0], 2*Axd[1], 2*Axd[2] };
    }
    return F;
}