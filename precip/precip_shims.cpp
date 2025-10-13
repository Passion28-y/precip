// 将析出相包裹与 ParaDiS 核心打通的两个 shim：
// - Precip_ForEachSegment：遍历本地域所有位错段，回调给包裹
// - Precip_AddForceToSegmentEndpoints：把等效端点力加到节点
//
// 依据（LLNL/ParaDiS 官方仓库）：
// - 节点/段的表达：include/Node.h（Node_t: x,y,z,numNbrs,burgX/Y/Z[] 等）
// - 段遍历与 PBC 最小镜像：src/ParadisStep.cc（获取邻居、ZImage 后得到 x2）
// - 力累加 API：include/Force.h（AddtoNodeForce 原型），src/NodeForce.cc（实现，内部带锁）

#include <cmath>
#include <cfloat>
#include <cstdlib>   // getenv
#include <cstdio>    // fprintf

#include "Home.h"
#include "Node.h"
#include "Util.h"    // GetNeighborNode, ZImage
#include "Force.h"   // AddtoNodeForce

// 在本地域节点表中按坐标匹配查找节点（容差 tol，单位=b）
static Node_t* findLocalNodeByPos(Home_t* home, double x, double y, double z, double tol=1e-9)
{
    if (!home || !home->nodeKeys) return nullptr;
    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t* node = home->nodeKeys[i];
        if (!node) continue;
        double dx = node->x - x;
        double dy = node->y - y;
        double dz = node->z - z;
        if (std::fabs(dx) <= tol && std::fabs(dy) <= tol && std::fabs(dz) <= tol) {
            return node;
        }
    }
    return nullptr;
}

extern "C" {

// 遍历本地域的所有段并回调 (x1,x2,b,xi)
void Precip_ForEachSegment(
    Home_t* home,
    void (*cb)(const double x1[3], const double x2[3],
               const double b [3], const double xi[3],
               void* user),
    void* user)
{
    if (!home || !cb) return;

    Param_t* param = home->param;
    const bool dbg = (std::getenv("PRECIP_DEBUG_SHIM") != nullptr);
    long segCount = 0;

    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t* node = home->nodeKeys[i];
        if (!node) continue;

        const double x1 = node->x;
        const double y1 = node->y;
        const double z1 = node->z;

        for (int arm = 0; arm < node->numNbrs; arm++) {
            Node_t* nbr = GetNeighborNode(home, node, arm);
            if (!nbr) continue;

            // 该段的 Burgers（以本端臂为准）
            double b[3] = { node->burgX[arm], node->burgY[arm], node->burgZ[arm] };

            // 二端点差向量，做最小镜像
            double dx = nbr->x - x1;
            double dy = nbr->y - y1;
            double dz = nbr->z - z1;
            ZImage(param, &dx, &dy, &dz);

            double x2 = x1 + dx;
            double y2 = y1 + dy;
            double z2 = z1 + dz;

            double L = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (L < DBL_EPSILON) continue;

            double xi[3] = { dx/L, dy/L, dz/L };
            const double X1[3] = { x1, y1, z1 };
            const double X2[3] = { x2, y2, z2 };

            cb(X1, X2, b, xi, user);
            segCount++;
        }
    }

    if (dbg) {
        std::fprintf(stderr, "[precip] shim: enumerated %ld segments in this scan\n", segCount);
    }
}

// 把等效端点力累加到节点（始终加到本地域节点；邻居若也在本地域则一并加）
void Precip_AddForceToSegmentEndpoints(
    Home_t* home,
    const double x1[3], const double x2[3],
    const double f1[3], const double f2[3])
{
    if (!home || !x1 || !x2 || !f1 || !f2) return;

    const bool dbg = (std::getenv("PRECIP_DEBUG_SHIM") != nullptr);
    double l1 = 0.0;

    if (Node_t* n1 = findLocalNodeByPos(home, x1[0], x1[1], x1[2])) {
        double f[3] = { f1[0], f1[1], f1[2] };
        AddtoNodeForce(n1, f);
        l1 += std::fabs(f1[0]) + std::fabs(f1[1]) + std::fabs(f1[2]);
    }
    if (Node_t* n2 = findLocalNodeByPos(home, x2[0], x2[1], x2[2])) {
        double f[3] = { f2[0], f2[1], f2[2] };
        AddtoNodeForce(n2, f);
        l1 += std::fabs(f2[0]) + std::fabs(f2[1]) + std::fabs(f2[2]);
    }

    if (dbg) {
        std::fprintf(stderr, "[precip] shim: add force L1=%.3e to endpoints\n", l1);
    }
}

} // extern "C"