#pragma once
#include <string>
#include <vector>

struct Ellipsoid;

// 控制文件默认路径
#ifndef DEFAULT_PRECIP_CTRL
#define DEFAULT_PRECIP_CTRL "precip/precip.ctrl"
#endif

// 各向同性材料参数（唯一来源）
struct Material {
    double mu{0.0};
    double nu{0.3};
};

struct PrecipGlobalParams {
    // 文件路径
    std::string ctrlFile{DEFAULT_PRECIP_CTRL};
    std::string inputFile;

    // 开关
    bool   enabled{false};
    bool   enable_misfit{false};

    // 各向同性参数
    double mu{0.0};
    double nu{0.3};

    // 碰撞粗判膨胀
    double obb_inflate{0.0};

    // 各向异性错配控制（由 precip.ctrl 与/或环境变量驱动）
    int    aniso{0};
    int    aniso_quad_n{3};
    double aniso_fd_eps{1e-3};
    int    force_on_hit_only{0};
};

// 读取 precip.ctrl → 填充 PrecipGlobalParams
bool load_precip_global_params(PrecipGlobalParams& gp);

// 读取沉淀数据文件（每行一个椭球），填充 out
// 返回 true 表示至少读取到一个有效椭球
bool load_precipitates_dat(const char* path, std::vector<Ellipsoid>& out);