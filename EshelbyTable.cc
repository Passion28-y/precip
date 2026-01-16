#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm> // for std::lower_bound
#include "EshelbyTable.h"

// 全局静态变量，确保只加载一次
EshelbyTable_t *g_EshelbyTable = NULL;

void InitEshelbyTable(const char* filename) {
    if (g_EshelbyTable != NULL) return; // 已经加载过

    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        printf("Error: Cannot open Eshelby table file %s\n", filename);
        exit(1);
    }

    g_EshelbyTable = (EshelbyTable_t*)calloc(1, sizeof(EshelbyTable_t));
    
    // 1. 读取 Header (Nx, Ny, Nz)
    fread(&g_EshelbyTable->Nx, sizeof(int), 1, fp);
    fread(&g_EshelbyTable->Ny, sizeof(int), 1, fp);
    fread(&g_EshelbyTable->Nz, sizeof(int), 1, fp);
    
    int Nx = g_EshelbyTable->Nx;
    int Ny = g_EshelbyTable->Ny;
    int Nz = g_EshelbyTable->Nz;
    
    // 2. 分配并读取网格坐标
    g_EshelbyTable->gridX = (double*)malloc(Nx * sizeof(double));
    g_EshelbyTable->gridY = (double*)malloc(Ny * sizeof(double));
    g_EshelbyTable->gridZ = (double*)malloc(Nz * sizeof(double));
    
    fread(g_EshelbyTable->gridX, sizeof(double), Nx, fp);
    fread(g_EshelbyTable->gridY, sizeof(double), Ny, fp);
    fread(g_EshelbyTable->gridZ, sizeof(double), Nz, fp);
    
    // 3. 读取应力数据
    long totalSize = (long)Nx * Ny * Nz * 6;
    g_EshelbyTable->stressData = (double*)malloc(totalSize * sizeof(double));
    fread(g_EshelbyTable->stressData, sizeof(double), totalSize, fp);
    
    g_EshelbyTable->isLoaded = true;
    fclose(fp);
    printf("Eshelby Table Loaded: %dx%dx%d grid.\n", Nx, Ny, Nz);
}

// 辅助函数：在非均匀网格中找到索引
int FindIndex(double val, double *grid, int N) {
    // 使用二分查找找到第一个 >= val 的位置
    double *ptr = std::lower_bound(grid, grid + N, val);
    int idx = (int)(ptr - grid);
    
    if (idx == 0) return 0;
    if (idx >= N) return N - 2;
    return idx - 1; // 返回 val 所在的区间左侧索引 [idx, idx+1]
}

// 核心插值函数
void GetStressFromTable(double x, double y, double z, double sigma[3][3]) {
    if (!g_EshelbyTable || !g_EshelbyTable->isLoaded) {
        // 实际使用时应处理此错误或自动加载
        return; 
    }
    
    EshelbyTable_t *T = g_EshelbyTable;
    
    // 1. 确定 x, y, z 在网格中的位置索引
    int ix = FindIndex(x, T->gridX, T->Nx);
    int iy = FindIndex(y, T->gridY, T->Ny);
    int iz = FindIndex(z, T->gridZ, T->Nz);
    
    // 2. 计算局部归一化坐标 u, v, w (0到1之间)
    double x0 = T->gridX[ix], x1 = T->gridX[ix+1];
    double y0 = T->gridY[iy], y1 = T->gridY[iy+1];
    double z0 = T->gridZ[iz], z1 = T->gridZ[iz+1];
    
    double u = (x - x0) / (x1 - x0);
    double v = (y - y0) / (y1 - y0);
    double w = (z - z0) / (z1 - z0);
    
    // 3. 三线性插值 (Trilinear Interpolation)
    // 我们需要对 6 个分量分别插值
    // Data indexing: index = (ix * Ny * Nz + iy * Nz + iz) * 6
    
    int NyNz = T->Ny * T->Nz;
    int Nz = T->Nz;
    
    // 8个顶点的偏移量
    int offsets[8];
    offsets[0] = ((ix  ) * NyNz + (iy  ) * Nz + (iz  )) * 6; // 000
    offsets[1] = ((ix  ) * NyNz + (iy  ) * Nz + (iz+1)) * 6; // 001
    offsets[2] = ((ix  ) * NyNz + (iy+1) * Nz + (iz  )) * 6; // 010
    offsets[3] = ((ix  ) * NyNz + (iy+1) * Nz + (iz+1)) * 6; // 011
    offsets[4] = ((ix+1) * NyNz + (iy  ) * Nz + (iz  )) * 6; // 100
    offsets[5] = ((ix+1) * NyNz + (iy  ) * Nz + (iz+1)) * 6; // 101
    offsets[6] = ((ix+1) * NyNz + (iy+1) * Nz + (iz  )) * 6; // 110
    offsets[7] = ((ix+1) * NyNz + (iy+1) * Nz + (iz+1)) * 6; // 111
    
    double s[6]; // 暂存插值后的6个分量
    
    for (int comp = 0; comp < 6; comp++) {
        double c00 = T->stressData[offsets[0] + comp] * (1-u) + T->stressData[offsets[4] + comp] * u;
        double c01 = T->stressData[offsets[1] + comp] * (1-u) + T->stressData[offsets[5] + comp] * u;
        double c10 = T->stressData[offsets[2] + comp] * (1-u) + T->stressData[offsets[6] + comp] * u;
        double c11 = T->stressData[offsets[3] + comp] * (1-u) + T->stressData[offsets[7] + comp] * u;
        
        double c0 = c00 * (1-v) + c10 * v;
        double c1 = c01 * (1-v) + c11 * v;
        
        s[comp] = c0 * (1-w) + c1 * w;
    }
    
    // 4. 填充 3x3 张量 (对称)
    // python 保存顺序: xx, xy, xz, yy, yz, zz
    sigma[0][0] = s[0]; sigma[0][1] = s[1]; sigma[0][2] = s[2];
    sigma[1][0] = s[1]; sigma[1][1] = s[3]; sigma[1][2] = s[4];
    sigma[2][0] = s[2]; sigma[2][1] = s[4]; sigma[2][2] = s[5];
}