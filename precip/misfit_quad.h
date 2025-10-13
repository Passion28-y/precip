#pragma once
#include <array>
#include <cmath>
#include "ellipsoid.h"
#include "precip_params.h"

// 各向同性错配应力 σ(X) 的求值（实现见 misfit_quad.cpp，不要重复定义）
void misfitStressAtPoint(const Ellipsoid& e,
                         const Material& mat,
                         const std::array<double,3>& X,
                         std::array<double,6>& sigmaVoigt);

// Voigt -> 3x3
inline void voigt_to_mat(const std::array<double,6>& v, double M[3][3]) {
    M[0][0]=v[0]; M[1][1]=v[1]; M[2][2]=v[2];
    M[1][2]=M[2][1]=v[3];
    M[2][0]=M[0][2]=v[4];
    M[0][1]=M[1][0]=v[5];
}

// 线积分将错配应力场转为等效端点力（模板定义必须在头文件）
template <class SigmaAtFunc>
inline void nodeForceFromMisfitOnSegment(const std::array<double,3>& X1,
                                         const std::array<double,3>& X2,
                                         const std::array<double,3>& b,
                                         const std::array<double,3>& xi,   // 需为单位切向
                                         SigmaAtFunc sigma_at,             // 回调：sigma_at(X, sigmaVoigt)
                                         std::array<double,3>& f1_out)     // 返回：分配到节点1的力（节点2取相反数）
{
    // 3 点 Gauss-Legendre（[-1,1] 的点权映射到 [0,1]）
    constexpr double sqrt35 = 0.7745966692414834;     // sqrt(3/5)
    const double tp[3] = { (1.0 - sqrt35)/2.0, 0.5, (1.0 + sqrt35)/2.0 };
    const double tw[3] = { 5.0/18.0, 4.0/9.0, 5.0/18.0 };

    // 几何
    const double dx[3] = { X2[0]-X1[0], X2[1]-X1[1], X2[2]-X1[2] };
    const double L = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (L <= 0.0) {
        f1_out = {0.0,0.0,0.0};
        return;
    }

    // 累积总力（线力密度沿长积分）
    double Ftot[3] = {0.0, 0.0, 0.0};

    for (int i=0; i<3; ++i) {
        const double t = tp[i];
        const double w = tw[i];

        std::array<double,3> X = { X1[0] + t*dx[0],
                                   X1[1] + t*dx[1],
                                   X1[2] + t*dx[2] };

        // 应力（Voigt）→ 矩阵
        std::array<double,6> sv{};
        sigma_at(X, sv);

        double S[3][3];
        voigt_to_mat(sv, S);

        // s_b = σ · b
        double sb[3] = {
            S[0][0]*b[0] + S[0][1]*b[1] + S[0][2]*b[2],
            S[1][0]*b[0] + S[1][1]*b[1] + S[1][2]*b[2],
            S[2][0]*b[0] + S[2][1]*b[1] + S[2][2]*b[2]
        };

        // PK 力密度 f = (σ·b) × ξ
        double f[3] = {
            sb[1]*xi[2] - sb[2]*xi[1],
            sb[2]*xi[0] - sb[0]*xi[2],
            sb[0]*xi[1] - sb[1]*xi[0]
        };

        // 加权积分：∫ f ds ≈ Σ f(t_i) * w_i * L
        Ftot[0] += f[0] * w * L;
        Ftot[1] += f[1] * w * L;
        Ftot[2] += f[2] * w * L;
    }

    // 对称分配到两端点（节点2取相反数，包裹处会设置 f2 = -f1）
    f1_out = { 0.5*Ftot[0], 0.5*Ftot[1], 0.5*Ftot[2] };
}