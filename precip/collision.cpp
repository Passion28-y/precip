#include "collision.h"
#include <cmath>
#include <algorithm>

static inline double dot3(const double* a,const double* b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static inline void sub3(const double* a,const double* b,double* c){ c[0]=a[0]-b[0]; c[1]=a[1]-b[1]; c[2]=a[2]-b[2]; }
static inline double clamp(double x,double lo,double hi){ return x<lo?lo:(x>hi?hi:x); }

// segment-AABB distance squared in principal frame
static double segAABBDist2(const double bmin[3], const double bmax[3], const double p0[3], const double p1[3]) {
    double best=1e300;
    for(double t: {0.0,0.25,0.5,0.75,1.0}){
        double q[3] = { p0[0] + t*(p1[0]-p0[0]),
                        p0[1] + t*(p1[1]-p0[1]),
                        p0[2] + t*(p1[2]-p0[2]) };
        double d2=0.0;
        for(int k=0;k<3;k++){
            double v = (q[k]<bmin[k])?(bmin[k]-q[k]):((q[k]>bmax[k])?(q[k]-bmax[k]):0.0);
            d2 += v*v;
        }
        best = std::min(best,d2);
    }
    return best;
}

bool segmentOBBMayCollide(const Ellipsoid& e, const Segment& seg, double inflate){
    std::array<double,3> y1,y2;
    toPrincipal(e, seg.x1, y1);
    toPrincipal(e, seg.x2, y2);
    double bmin[3] = { -e.r[0]-inflate, -e.r[1]-inflate, -e.r[2]-inflate };
    double bmax[3] = { +e.r[0]+inflate, +e.r[1]+inflate, +e.r[2]+inflate };
    double d2 = segAABBDist2(bmin,bmax, y1.data(), y2.data());
    return d2 <= inflate*inflate + 1e-30;
}

bool segmentEllipsoidIntersectAlpha(const Ellipsoid& e, const Segment& seg, double& alpha_hit){
    double d1[3] = { seg.x1[0]-e.c[0], seg.x1[1]-e.c[1], seg.x1[2]-e.c[2] };
    double d2v[3]= { seg.x2[0]-e.c[0], seg.x2[1]-e.c[1], seg.x2[2]-e.c[2] };
    double v[3] = { d2v[0]-d1[0], d2v[1]-d1[1], d2v[2]-d1[2] };
    double Ad1[3] = { e.A[0]*d1[0]+e.A[1]*d1[1]+e.A[2]*d1[2],
                      e.A[3]*d1[0]+e.A[4]*d1[1]+e.A[5]*d1[2],
                      e.A[6]*d1[0]+e.A[7]*d1[1]+e.A[8]*d1[2] };
    double Av[3]  = { e.A[0]*v[0]+e.A[1]*v[1]+e.A[2]*v[2],
                      e.A[3]*v[0]+e.A[4]*v[1]+e.A[5]*v[2],
                      e.A[6]*v[0]+e.A[7]*v[1]+e.A[8]*v[2] };
    double c0 = d1[0]*Ad1[0]+d1[1]*Ad1[1]+d1[2]*Ad1[2] - 1.0;
    double b  = 2.0*(v[0]*Ad1[0]+v[1]*Ad1[1]+v[2]*Ad1[2]);
    double a  = v[0]*Av[0]+v[1]*Av[1]+v[2]*Av[2];
    double disc = b*b - 4*a*c0;
    if (disc < 0) return false;
    double sdisc = std::sqrt(std::max(0.0, disc));
    double r1 = (-b - sdisc) / (2*a);
    double r2 = (-b + sdisc) / (2*a);
    bool ok=false;
    double best=1e300;
    for(double r: {r1,r2}){
        if (r>=0.0 && r<=1.0){
            if (r < best){ best=r; ok=true; }
        }
    }
    if (ok){ alpha_hit = best; return true; }
    return false;
}

// Projection of point onto ellipsoid using principal-frame Newton on λ (outside/inside)
static bool projectPointToEllipsoidPrincipal(const double q[3], const double r[3], double s[3], double& signed_dist, int maxIt, double tol){
    double q2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    if (q2==0.0){
        s[0]=r[0]; s[1]=0.0; s[2]=0.0;
        signed_dist = -std::sqrt(r[0]*r[0]);
        return true;
    }
    double phi = (q[0]*q[0])/(r[0]*r[0]) + (q[1]*q[1])/(r[1]*r[1]) + (q[2]*q[2])/(r[2]*r[2]);
    bool outside = (phi >= 1.0);
    double lam = 0.0;
    if (outside) lam = 0.0; else {
        // avoid initializer_list overload of std::min
        double minri2 = std::min(r[0]*r[0], std::min(r[1]*r[1], r[2]*r[2]));
        lam = -0.5*minri2;
    }
    for(int it=0; it<maxIt; ++it){
        double f=0.0, df=0.0;
        for(int i=0;i<3;i++){
            double Ri2 = r[i]*r[i];
            double d   = Ri2 + lam;
            double gi  = q[i]*q[i]*Ri2/(d*d);
            f  += gi;
            df += -2.0*q[i]*q[i]*Ri2/(d*d*d);
        }
        f -= 1.0;
        if (std::abs(f) < tol) break;
        double step = -f/df;
        lam += step;
        if (std::abs(step) < tol) break;
    }
    for(int i=0;i<3;i++){
        double Ri2 = r[i]*r[i];
        s[i] = q[i] * (Ri2/(Ri2+lam));
    }
    double dvec[3] = { s[0]-q[0], s[1]-q[1], s[2]-q[2] };
    double dist = std::sqrt(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
    signed_dist = outside ? dist : -dist;
    return true;
}

bool projectPointToEllipsoid(const Ellipsoid& e, const std::array<double,3>& p,
                             std::array<double,3>& s_out, double& signed_dist,
                             int maxIt, double tol){
    std::array<double,3> q;
    toPrincipal(e, p, q);
    double sp[3];
    if (!projectPointToEllipsoidPrincipal(q.data(), e.r.data(), sp, signed_dist, maxIt, tol)) return false;
    std::array<double,3> sw;
    toWorld(e, {sp[0],sp[1],sp[2]}, sw);
    s_out = sw;
    return true;
}

MinDistResult minDistEllipsoid(const Ellipsoid& e, const Segment& seg){
    MinDistResult out;
    out.intersects=false;
    double alpha_hit;
    if (segmentEllipsoidIntersectAlpha(e, seg, alpha_hit)){
        out.intersects = true;
        out.alpha = alpha_hit;
        out.x = { seg.x1[0] + alpha_hit*(seg.x2[0]-seg.x1[0]),
                  seg.x1[1] + alpha_hit*(seg.x2[1]-seg.x1[1]),
                  seg.x1[2] + alpha_hit*(seg.x2[2]-seg.x1[2]) };
        out.s = out.x;
        out.d = 0.0;
        return out;
    }
    auto dist_eval = [&](double a, std::array<double,3>& s, std::array<double,3>& x)->double{
        x = { seg.x1[0] + a*(seg.x2[0]-seg.x1[0]),
              seg.x1[1] + a*(seg.x2[1]-seg.x1[1]),
              seg.x1[2] + a*(seg.x2[2]-seg.x1[2]) };
        double sd;
        projectPointToEllipsoid(e, x, s, sd);
        return std::abs(sd);
    };
    double lo=0.0, hi=1.0;
    const double gr = (std::sqrt(5.0)-1.0)/2.0;
    double c = hi - gr*(hi-lo);
    double d = lo + gr*(hi-lo);
    std::array<double,3> sC,sD,xC,xD;
    double fC = dist_eval(c, sC, xC);
    double fD = dist_eval(d, sD, xD);
    for(int it=0; it<40; ++it){
        if (fC > fD){
            lo = c; c = d; fC = fD; sC=sD; xC=xD;
            d = lo + gr*(hi-lo);
            fD = dist_eval(d, sD, xD);
        } else {
            hi = d; d = c; fD = fC; sD=sC; xD=xC;
            c = hi - gr*(hi-lo);
            fC = dist_eval(c, sC, xC);
        }
    }
    if (fC < fD){
        out.alpha = c; out.d = fC; out.s = sC; out.x = xC;
    } else {
        out.alpha = d; out.d = fD; out.s = sD; out.x = xD;
    }
    return out;
}