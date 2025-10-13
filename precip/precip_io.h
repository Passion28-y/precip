#pragma once
#include <vector>
#include <string>
#include "precip_params.h"
#include "ellipsoid.h"

// 读取 ctrl，填充 prm。返回：true=成功打开并解析，false=文件不存在或打开失败（字段保持默认值）
bool load_precip_global_params(PrecipGlobalParams& prm);

// 读取析出体数据文件（简单格式：每行至少 6 个数）
// 格式（空白分隔，# 注释）：
//   cx cy cz  r1 r2 r3  [alphaZ betaY gammaX deg]  [exx eyy ezz]
// 角度为可选，单位度；若省略，用单位旋转。应变三分量可选，省略则为 0。
bool load_precipitates_dat(const char* path, std::vector<Ellipsoid>& out);