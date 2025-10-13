#include "misfit_quad.h"
#include <cmath>
#include <algorithm>

// Convert 3x3 tensor to Voigt [xx,yy,zz,yz,xz,xy]
static inline void matToVoigt(const double s[3][3], std::array<double,6>& v){
    v[0]=s[0][0]; v[1]=s[1][1]; v[2]=s[2][2];
    v[3]=s[1][2]; v[4]=s[0][2]; v[5]=s[0][1];
}

// Isotropic Kelvin solution: stress at x due to unit point force at x' in j-direction (σ_ij,k)
// Full closed form is lengthy; here we approximate via displacement gradient and Hooke's law. For practicality, we use a regularized kernel.
// This provides a reasonable first-order field for PK force accumulation near the inclusion surface.
static void approxKelvinStress(const Material& mat, const double R[3], double tvec[3], double sout[3][3]){
    double mu = mat.mu;
    double nu = mat.nu;
    double r2 = R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
    double r  = std::sqrt(r2 + 1e-30);
    // Regularized 1/r^3 tensor kernel
    double c = 1.0/(4.0*M_PI);
    double k = c / (r*r*r + 1e-24);
    // Build traction-induced stress proxy σ ≈ k [ t ⊗ R + R ⊗ t ] scaled by isotropic factors
    // NOTE: This is a simplified proxy to avoid modifying ParaDiS core; for accurate work, replace with exact kernel.
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){
        sout[i][j] = k * ( tvec[i]*R[j] + R[i]*tvec[j] );
    }
    // Symmetrize and scale
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){
        sout[i][j] = 0.5*(sout[i][j] + sout[j][i]) * ( (1.0-2.0*nu)/(1.0-nu) );
    }
}

void misfitStressAtPoint(const Ellipsoid& e, const Material& mat,
                         const std::array<double,3>& x,
                         std::array<double,6>& sigma_voigt)
{
#if !defined(PRECIP_ENABLE_MISFIT)
    // 默认关闭：返回零应力。若需要激活近似错配应力，请编译时 -DPRECIP_ENABLE_MISFIT
    sigma_voigt = {0,0,0,0,0,0};
    return;
#else
    // Surface parameterization in principal frame
    // θ ∈ [0,2π), φ ∈ [0,π]
    const int Ntheta=12, Nphi=6; // 可调
    const double dth = 2.0*M_PI/Ntheta;
    const double dph = M_PI/Nphi;

    // Precompute σ* = C:ε* in principal frame (diagonal ε*)
    double lam = 2.0*mat.mu*mat.nu/(1.0-2.0*mat.nu); // Lame λ
    double eps[3] = { e.eigenstrain_voigt[0], e.eigenstrain_voigt[1], e.eigenstrain_voigt[2] };
    double tr = eps[0]+eps[1]+eps[2];
    double sigma_star[3] = { 2*mat.mu*eps[0] + lam*tr,
                             2*mat.mu*eps[1] + lam*tr,
                             2*mat.mu*eps[2] + lam*tr };

    double sig[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(int it=0; it<Ntheta; ++it){
        double th = (it+0.5)*dth;
        double cth=std::cos(th), sth=std::sin(th);
        for(int ip=0; ip<Nphi; ++ip){
            double ph = (ip+0.5)*dph;
            double cph=std::cos(ph), sph=std::sin(ph);
            // principal point x' and normal n' (normalized with ellipsoid metric)
            double xp[3] = { e.r[0]*sph*cth, e.r[1]*sph*sth, e.r[2]*cph };
            // normal proportional to (x_i / r_i^2)
            double nraw[3] = { xp[0]/(e.r[0]*e.r[0]), xp[1]/(e.r[1]*e.r[1]), xp[2]/(e.r[2]*e.r[2]) };
            double nrm = std::sqrt(nraw[0]*nraw[0]+nraw[1]*nraw[1]+nraw[2]*nraw[2])+1e-30;
            double np[3] = { nraw[0]/nrm, nraw[1]/nrm, nraw[2]/nrm };
            // area element Jacobian on ellipsoid surface (approximate)
            double J = std::sqrt( (e.r[0]*e.r[1]*e.r[2])*( std::pow(sph,2)/ ( (e.r[0]*e.r[0])*(e.r[1]*e.r[1]) ) + std::pow(cph,2)/( (e.r[2]*e.r[2]) ) ) + 1e-30 );
            // traction t' = σ* · n' (principal frame)
            double tp[3] = { sigma_star[0]*np[0], sigma_star[1]*np[1], sigma_star[2]*np[2] };
            // world point xw' and normal/traction in world via U
            std::array<double,3> xw;
            toWorld(e, {xp[0],xp[1],xp[2]}, xw);
            double nW[3] = {
                e.U[0]*np[0]+e.U[1]*np[1]+e.U[2]*np[2],
                e.U[3]*np[0]+e.U[4]*np[1]+e.U[5]*np[2],
                e.U[6]*np[0]+e.U[7]*np[1]+e.U[8]*np[2]
            };
            double tW[3] = {
                e.U[0]*tp[0]+e.U[1]*tp[1]+e.U[2]*tp[2],
                e.U[3]*tp[0]+e.U[4]*tp[1]+e.U[5]*tp[2],
                e.U[6]*tp[0]+e.U[7]*tp[1]+e.U[8]*tp[2]
            };
            // kernel from x' to x
            double R[3] = { x[0]-xw[0], x[1]-xw[1], x[2]-xw[2] };
            double sK[3][3];
            approxKelvinStress(mat, R, tW, sK);
            // accumulate: σ += sK * (J * sinφ * dθ dφ)
            double w = (sph) * dth * dph; // surface element weight on sphere
            double wt = w * (e.r[0]*e.r[1]*e.r[2]) / ( (e.r[0]*e.r[1]*e.r[2]) + 1e-30 ); // mild scaling
            for(int i=0;i<3;i++) for(int j=0;j<3;j++){
                sig[i][j] += sK[i][j] * wt;
            }
        }
    }
    matToVoigt(sig, sigma_voigt);
#endif
}