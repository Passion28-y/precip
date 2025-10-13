#include "misfit_aniso.h"
#include "c_tensor.h"
#include "Home.h"
#include "FM.h"
#include "ellipsoid.h"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

// 旋转：out = U * v（U 指针指向 9 个连续元素）
static inline void mat33_vec3(const double* U, const double v[3], double out[3]) {
    out[0]=U[0]*v[0]+U[1]*v[1]+U[2]*v[2];
    out[1]=U[3]*v[0]+U[4]*v[1]+U[5]*v[2];
    out[2]=U[6]*v[0]+U[7]*v[1]+U[8]*v[2];
}
static inline void mat33T_vec3(const double* U, const double v[3], double out[3]) {
    out[0]=U[0]*v[0]+U[3]*v[1]+U[6]*v[2];
    out[1]=U[1]*v[0]+U[4]*v[1]+U[7]*v[2];
    out[2]=U[2]*v[0]+U[5]*v[1]+U[8]*v[2];
}

static inline double ellipsoid_volume(const Ellipsoid& e) {
    return 4.0*M_PI/3.0 * e.r[0]*e.r[1]*e.r[2];
}

// 将 [-1,1]^3 的节点映射到椭球体内（简化中点法），并检查是否在体内
static inline bool map_unitcube_to_ellipsoid(const Ellipsoid& e, double u, double v, double w, double Xq[3])
{
    double q[3] = { 0.5*u, 0.5*v, 0.5*w };
    double r2 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    if (r2>0.25) {
        double s = 0.5/std::sqrt(r2);
        q[0]*=s; q[1]*=s; q[2]*=s;
    }
    double qa[3] = { q[0]*e.r[0], q[1]*e.r[1], q[2]*e.r[2] };
    double qo[3]; mat33_vec3(e.U.data(), qa, qo);
    Xq[0] = e.c[0] + qo[0];
    Xq[1] = e.c[1] + qo[1];
    Xq[2] = e.c[2] + qo[2];

    double qc[3] = { Xq[0]-e.c[0], Xq[1]-e.c[1], Xq[2]-e.c[2] };
    double qa2[3]; mat33T_vec3(e.U.data(), qc, qa2);
    double el = (qa2[0]/e.r[0])*(qa2[0]/e.r[0]) + (qa2[1]/e.r[1])*(qa2[1]/e.r[1]) + (qa2[2]/e.r[2])*(qa2[2]/e.r[2]);
    return (el<=1.0);
}

// 调 ParaDiS 的一阶 Green 导数（使用 8 参数变体：含 qMax）
static inline void get_dGdx(Home_t* home, const double R[3], double dGdx[3][6]) {
    int qMax = 10;
    if (home && home->param && home->param->anisoHarmonicsNumTermsBase > 0) {
        qMax = home->param->anisoHarmonicsNumTermsBase;
    }
    DerivativeOfGreensFunction(home, 0.0, 0.0, 0.0, R[0], R[1], R[2], qMax, dGdx);
}

// 用中心差分近似二阶导数：D2[a][vin][b] = ∂/∂x_b ( dGdx[a][vin] )
static inline void approx_d2G(Home_t* home, const double R[3], double h, double D2[3][6][3]) {
    double dP[3][6], dM[3][6];
    for (int b=0;b<3;b++) {
        double Rp[3] = { R[0], R[1], R[2] };
        double Rm[3] = { R[0], R[1], R[2] };
        Rp[b]+=h; Rm[b]-=h;
        get_dGdx(home, Rp, dP);
        get_dGdx(home, Rm, dM);
        for (int a=0;a<3;a++) for (int vin=0;vin<6;vin++) {
            D2[a][vin][b] = (dP[a][vin]-dM[a][vin])/(2.0*h);
        }
    }
}

// 显式收缩：σ_ij = C_ijkl ∑_{m,n} ∂^2 G_km / ∂x_l∂x_n ε*_mn
static inline void accumulate_sigma_with_C(
    const std::array<std::array<std::array<std::array<double,3>,3>,3>,3>& C4,
    const double D2[3][6][3],
    const std::array<double,6>& eigenVoigt,
    double weight,
    std::array<double,6>& sigmaVoigt)
{
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        double sij = 0.0;
        for (int k=0;k<3;k++) for (int l=0;l<3;l++) {
            double Cijkl = C4[i][j][k][l];
            double Skl = 0.0;
            for (int vin=0; vin<6; vin++) {
                int vLM = (vin<=2) ? ij_voigt(l, vin) : ij_voigt(l, (vin==3?2:(vin==4?0:1)));
                for (int n=0;n<3;n++) {
                    int baseM = (vin<=2? vin : (vin==3?1:(vin==4?2:0)));
                    int vMN = ij_voigt(baseM, n);
                    double d2 = D2[k][vLM][n];
                    Skl += d2 * eigenVoigt[vMN];
                }
            }
            sij += Cijkl * Skl;
        }
        sigmaVoigt[ij_voigt(i,j)] += weight * sij;
    }
}

bool misfitStressAtPoint_Aniso_C(_home* h, const Ellipsoid& e,
                                 const std::array<double,3>& Xp,
                                 const std::array<double,6>& eigenVoigt,
                                 std::array<double,6>& sigmaVoigt,
                                 int quadN, double fd_eps)
{
    Home_t* home = reinterpret_cast<Home_t*>(h);
    sigmaVoigt = {0,0,0,0,0,0};

    // 构建 C_ijkl（来自 Param->elasticConstantMatrix 或 C11/C12/C44 组装）
    std::array<std::array<std::array<std::array<double,3>,3>,3>,3> C4{};
    build_C4_from_C66(home, C4);

    // 积分节点
    std::vector<double> nodes(quadN);
    for (int i=0;i<quadN;i++) nodes[i] = -1.0 + 2.0*(i+0.5)/quadN;
    const double vol = ellipsoid_volume(e);
    const double baseW = vol / (quadN*quadN*quadN);

    for (int iu=0; iu<quadN; iu++)
    for (int iv=0; iv<quadN; iv++)
    for (int iw=0; iw<quadN; iw++) {
        double u=nodes[iu], v=nodes[iv], w=nodes[iw];
        double Xq[3];
        if (!map_unitcube_to_ellipsoid(e, u,v,w, Xq)) continue;

        double R[3] = { Xp[0]-Xq[0], Xp[1]-Xq[1], Xp[2]-Xq[2] };
        double Rn = std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
        if (Rn<1e-10) continue;

        double hstep = std::max(fd_eps*Rn, 1e-6);
        double D2[3][6][3];
        approx_d2G(home, R, hstep, D2);

        accumulate_sigma_with_C(C4, D2, eigenVoigt, baseW, sigmaVoigt);
    }
    return true;
}