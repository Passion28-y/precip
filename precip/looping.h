#pragma once
#include "ellipsoid.h"
#include <array>
#include <algorithm>

// l_loop = 0.5 * ( r_min^2 / r_max )
inline double loopingLength(const Ellipsoid& e){
    // avoid initializer_list overloads to be compatible with older stdlib/flags
    double rmin = std::min(e.r[0], std::min(e.r[1], e.r[2]));
    double rmax = std::max(e.r[0], std::max(e.r[1], e.r[2]));
    return 0.5*(rmin*rmin/rmax);
}

// 给定铰接段上距铰接点 l_loop 的点在 t 与 t+Δt 的位置，直线轨迹求 β 使其落在椭球面
bool loopingCollisionBeta(const Ellipsoid& e,
                          const std::array<double,3>& x_old,
                          const std::array<double,3>& x_new,
                          double& beta);