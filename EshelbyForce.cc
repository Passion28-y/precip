#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Home.h"
#include "EshelbyTable.h"

// 5点高斯积分权重和位置
static const double gauss_x[5] = {0.0, 0.538469, -0.538469, 0.90618, -0.90618};
static const double gauss_w[5] = {0.568889, 0.478629, 0.478629, 0.236927, 0.236927};

void EshelbyForce(Home_t *home, real8 EPos[3], real8 radius[3], 
                  real8 EStrain[6], real8 pos1[3], real8 pos2[3],
                  real8 burg[3], real8 f1[3], real8 f2[3])
{
    // 1. 初始化表 (支持环境变量)
    if (!g_EshelbyTable) {
        const char* env_path = getenv("ESHELBY_TABLE_PATH");
        if (env_path) {
            printf("Loading Eshelby Table: %s\n", env_path);
            InitEshelbyTable(env_path);
        } else {
            printf("Loading default table: EshelbyStressTable.bin\n");
            InitEshelbyTable("EshelbyStressTable.bin");
        }
    }
    
    f1[0] = 0.0; f1[1] = 0.0; f1[2] = 0.0;
    f2[0] = 0.0; f2[1] = 0.0; f2[2] = 0.0;

    // --- 1. 自动识别变体 (基于半径形状) ---
    // 假设: 最小的半径对应的轴就是 c 轴 (错配轴)
    int variant = 3; // 默认为 Z 轴变体
    if (radius[0] < radius[1] && radius[0] < radius[2]) {
        variant = 1; // X 轴变体 (短轴在 X)
    } else if (radius[1] < radius[0] && radius[1] < radius[2]) {
        variant = 2; // Y 轴变体 (短轴在 Y)
    } else {
        variant = 3; // Z 轴变体 (短轴在 Z)
    }

    // 2. 几何准备
    real8 diff[3], t[3];
    for (int i=0; i<3; i++) diff[i] = pos2[i] - pos1[i];
    real8 L = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    if (L < 1e-12) return;
    for (int i=0; i<3; i++) t[i] = diff[i] / L; 

    // 3. 高斯积分
    for (int k = 0; k < 5; k++) {
        double xi = gauss_x[k];
        double weight = gauss_w[k];
        double s = (L * 0.5) * (xi + 1.0);
        double integration_factor = weight * (L * 0.5);

        double p_integ[3];
        for (int i=0; i<3; i++) p_integ[i] = pos1[i] + t[i] * s;

        // 计算相对于粒子中心的物理坐标
        double dx = p_integ[0] - EPos[0];
        double dy = p_integ[1] - EPos[1];
        double dz = p_integ[2] - EPos[2];

        // --- 核心：坐标映射 (Physical -> Table) ---
        // 我们假设 Table 是按照 Z 轴为短轴生成的。
        // 我们需要把当前变体的短轴“旋转”到 Z 轴去查表。
        double tx, ty, tz; // Table coordinates

        if (variant == 1) { // X is short -> Map X to Table's Z
            // Mapping: Phys(Y, Z, X) -> Table(X, Y, Z)
            tx = dy; ty = dz; tz = dx;
        } else if (variant == 2) { // Y is short -> Map Y to Table's Z
            // Mapping: Phys(Z, X, Y) -> Table(X, Y, Z)
            tx = dz; ty = dx; tz = dy;
        } else { // Z is short -> Direct map
            tx = dx; ty = dy; tz = dz;
        }

        // 查表 (得到的是 Table 坐标系下的应力)
        double sigma_table[3][3];
        GetStressFromTable(tx, ty, tz, sigma_table);
		// 关键：将表中的物理应力 (Pa) 转换为 ParaDiS 内部的归一化应力
		// 这样就抵消了 $10^9$ 的量级错误，同时保留了材料的模量比例
		real8 invMu = 1.0 / home->param->shearModulus;
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				sigma_table[i][j] *= invMu;
			}
		}

        // --- 核心：应力映射回去 (Table -> Physical) ---
        double sigma[3][3];

        if (variant == 1) {
            // Table(xx, yy, zz) -> Phys(yy, zz, xx)
            // Table(xy, xz, yz) -> Phys(yz, yx, zx)
            sigma[0][0] = sigma_table[2][2]; // s_xx = t_zz
            sigma[1][1] = sigma_table[0][0]; // s_yy = t_xx
            sigma[2][2] = sigma_table[1][1]; // s_zz = t_yy

            sigma[0][1] = sigma_table[2][0]; // s_xy = t_zx
            sigma[0][2] = sigma_table[2][1]; // s_xz = t_zy
            sigma[1][2] = sigma_table[0][1]; // s_yz = t_xy
        } else if (variant == 2) {
            // Table(xx, yy, zz) -> Phys(zz, xx, yy)
            sigma[0][0] = sigma_table[1][1]; // s_xx = t_yy
            sigma[1][1] = sigma_table[2][2]; // s_yy = t_zz
            sigma[2][2] = sigma_table[0][0]; // s_zz = t_xx

            sigma[0][1] = sigma_table[1][2]; // s_xy = t_yz
            sigma[0][2] = sigma_table[1][0]; // s_xz = t_yx
            sigma[1][2] = sigma_table[2][0]; // s_yz = t_zx
        } else {
            // No rotation
            for(int i=0;i<3;i++) for(int j=0;j<3;j++) sigma[i][j] = sigma_table[i][j];
        }
        // 补全对称分量
        sigma[1][0] = sigma[0][1];
        sigma[2][0] = sigma[0][2];
        sigma[2][1] = sigma[1][2];

        // --- 计算力 (常规步骤) ---
        double sigma_b[3] = {0, 0, 0};
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                sigma_b[i] += sigma[i][j] * burg[j];
            }
        }
        
        double f_density[3];
        f_density[0] = sigma_b[1]*t[2] - sigma_b[2]*t[1];
        f_density[1] = sigma_b[2]*t[0] - sigma_b[0]*t[2];
        f_density[2] = sigma_b[0]*t[1] - sigma_b[1]*t[0];
        
        double N2 = s / L;
        double N1 = 1.0 - N2;
        
        for (int i=0; i<3; i++) {
            f1[i] += f_density[i] * N1 * integration_factor;
            f2[i] += f_density[i] * N2 * integration_factor;
        }
    }
}