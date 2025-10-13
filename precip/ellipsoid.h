#pragma once
#include <array>
#include <vector>

struct Ellipsoid {
    int id = -1;
    std::array<double,3> c{};     // center
    std::array<double,3> r{};     // radii r1,r2,r3
    std::array<double,9> U{};     // rotation, row-major (world <- principal)
    std::array<double,9> A{};     // 3x3 world A = U Σ² U^T  with Σ=diag(1/r1,1/r2,1/r3)
    std::array<double,6> aabb{};  // axis-aligned bbox for quick cull [xmin,xmax, ymin,ymax, zmin,zmax]
    std::array<double,3> eigenstrain_voigt{}; // principal diag eigenstrain (e11,e22,e33)
    bool   allow_cut = true;
    double tau_cut = 0.0; // [Pa]
};

// Build A, AABB from c,r,U
void build_matrices(Ellipsoid& e);

// Level set F(x)= (x-c)^T A (x-c) - 1 ; grad = 2A(x-c)
double implicitF(const Ellipsoid& e, const std::array<double,3>& x, std::array<double,3>* grad=nullptr);

// Coordinate transforms: world <-> principal
inline void toPrincipal(const Ellipsoid& e, const std::array<double,3>& xw, std::array<double,3>& yp) {
    double dx[3] = { xw[0]-e.c[0], xw[1]-e.c[1], xw[2]-e.c[2] };
    // y = U^T (x-c)
    yp[0] = e.U[0]*dx[0] + e.U[3]*dx[1] + e.U[6]*dx[2];
    yp[1] = e.U[1]*dx[0] + e.U[4]*dx[1] + e.U[7]*dx[2];
    yp[2] = e.U[2]*dx[0] + e.U[5]*dx[1] + e.U[8]*dx[2];
}
inline void toWorld(const Ellipsoid& e, const std::array<double,3>& yp, std::array<double,3>& xw) {
    // x = U y + c
    xw[0] = e.U[0]*yp[0] + e.U[1]*yp[1] + e.U[2]*yp[2] + e.c[0];
    xw[1] = e.U[3]*yp[0] + e.U[4]*yp[1] + e.U[5]*yp[2] + e.c[1];
    xw[2] = e.U[6]*yp[0] + e.U[7]*yp[1] + e.U[8]*yp[2] + e.c[2];
}