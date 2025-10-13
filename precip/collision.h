#pragma once
#include "ellipsoid.h"
#include <array>

struct Segment { std::array<double,3> x1, x2; };

struct MinDistResult {
    double d=0.0;      // signed dist (>0 outside; 0 on)
    double alpha=0.0;  // location on segment
    std::array<double,3> x{}; // closest on segment
    std::array<double,3> s{}; // closest on surface
    bool intersects=false;     // pierce?
};

// OBB预检（此处使用椭球主轴对齐的AABB）：将段变换到主轴系，计算段-盒距离；inflate 为安全缓冲
bool segmentOBBMayCollide(const Ellipsoid& e, const Segment& seg, double inflate);

// 段-椭球穿越（Case 3）：解 (x(α)-c)^T A (x(α)-c) = 1，返回 α∈[0,1] 是否存在（取最小正根）
bool segmentEllipsoidIntersectAlpha(const Ellipsoid& e, const Segment& seg, double& alpha_hit);

// 点投影到椭球面（主轴系牛顿/闭式混合）：返回表面点 s_out
bool projectPointToEllipsoid(const Ellipsoid& e, const std::array<double,3>& p,
                             std::array<double,3>& s_out, double& signed_dist,
                             int maxIt=50, double tol=1e-12);

// 椭球一般情形：若不穿越 → 沿段用黄金分割搜索最小距离（每点调用 point->ellipsoid 投影）
MinDistResult minDistEllipsoid(const Ellipsoid& e, const Segment& seg);