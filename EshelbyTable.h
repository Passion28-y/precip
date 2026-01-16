#ifndef _ESHELBY_TABLE_H
#define _ESHELBY_TABLE_H

#include "Home.h" // 包含 ParaDiS 的基础定义

// 定义存储应力场数据的结构体
struct EshelbyTable_t {
    int Nx, Ny, Nz;
    double *gridX; // X轴网格坐标 (长度 Nx)
    double *gridY; // Y轴网格坐标 (长度 Ny)
    double *gridZ; // Z轴网格坐标 (长度 Nz)
    
    // 应力数据: 扁平化的一维数组
    // 大小: Nx * Ny * Nz * 6
    // 存储顺序: x -> y -> z -> component(xx, xy, xz, yy, yz, zz)
    double *stressData; 
    
    bool isLoaded;
};

extern EshelbyTable_t *g_EshelbyTable;
// 全局函数声明
void InitEshelbyTable(const char* filename);
//void FreeEshelbyTable();
void GetStressFromTable(double x, double y, double z, double sigma[3][3]);

#endif