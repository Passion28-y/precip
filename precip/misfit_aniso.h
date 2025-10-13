#pragma once
#include <array>
struct _home;
struct Ellipsoid;

// 各向异性错配应力（显式 C_ijkl 收缩 + 体积分 + 二阶导数有限差分）
// 注意：参数名用 Xp 避免与 ParaDiS 的 X/Y/Z 宏冲突
bool misfitStressAtPoint_Aniso_C(_home* home,
                                 const Ellipsoid& e,
                                 const std::array<double,3>& Xp,
                                 const std::array<double,6>& eigenVoigt,
                                 std::array<double,6>& sigmaVoigt,
                                 int quadN = 3,
                                 double fd_eps = 1e-3);