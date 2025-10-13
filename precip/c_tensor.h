// 从 ParaDiS 的 Param->elasticConstantMatrix[6][6] 构建 C_ijkl（3x3x3x3）
// 并提供 Voigt <-> 张量索引工具。
#pragma once
#include "Home.h"
#include <array>
#include <cassert>

inline int voigt_ij(int v, int& i, int& j) {
    static const int map[6][2] = {{0,0},{1,1},{2,2},{1,2},{2,0},{0,1}};
    i = map[v][0]; j = map[v][1];
    return 0;
}
inline int ij_voigt(int i, int j) {
    if (i==j) return i; // 0,1,2 -> xx,yy,zz
    if ((i==1 && j==2) || (i==2 && j==1)) return 3; // yz
    if ((i==2 && j==0) || (i==0 && j==2)) return 4; // zx
    return 5; // xy
}

inline void build_C4_from_C66(Home_t* home, std::array<std::array<std::array<std::array<double,3>,3>,3>,3>& C4)
{
    Param_t* p = home->param;
    // 若 elasticConstantMatrix 未设置（全 0），尝试由 C11/C12/C44 等自动构建
    bool allZero = true;
    for (int a=0;a<6;a++) for (int b=0;b<6;b++) if (p->elasticConstantMatrix[a][b]!=0.0) { allZero=false; break; }

    double C66[6][6] = {};
    if (!allZero) {
        for (int a=0;a<6;a++) for (int b=0;b<6;b++) C66[a][b] = p->elasticConstantMatrix[a][b];
    } else {
        // 立方晶体：用 SetElasticConstantMatrix 的同样布局（这里手工拼一份）
        C66[0][0] = p->C11; C66[0][1] = p->C12; C66[0][2] = p->C12;
        C66[1][0] = p->C12; C66[1][1] = p->C11; C66[1][2] = p->C12;
        C66[2][0] = p->C12; C66[2][1] = p->C12; C66[2][2] = p->C11;
        C66[3][3] = p->C44; C66[4][4] = p->C44; C66[5][5] = p->C44;
    }

    // 将 6x6 映射到 3x3x3x3
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
    for (int k=0;k<3;k++) for (int l=0;l<3;l++) {
        int a = ij_voigt(i,j);
        int b = ij_voigt(k,l);
        C4[i][j][k][l] = C66[a][b];
    }
}