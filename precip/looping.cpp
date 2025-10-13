#include "looping.h"
#include "collision.h"

bool loopingCollisionBeta(const Ellipsoid& e,
                          const std::array<double,3>& x_old,
                          const std::array<double,3>& x_new,
                          double& beta){
    Segment seg{ x_old, x_new };
    return segmentEllipsoidIntersectAlpha(e, seg, beta);
}