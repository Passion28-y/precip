#include "precip_params.h"
#include "ellipsoid.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#if defined(_WIN32) || defined(_WIN64)
#define strtok_r strtok_s
#endif

static inline void trim(std::string& s) {
    size_t i=0, j=s.size();
    while (i<j && std::isspace((unsigned char)s[i])) i++;
    while (j>i && std::isspace((unsigned char)s[j-1])) j--;
    s = s.substr(i, j-i);
}
static inline std::string lower(std::string s) {
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}
static bool parse_kv(const std::string& line, std::string& key, std::string& val) {
    std::string s = line;
    auto pos_hash = s.find('#');
    if (pos_hash != std::string::npos) s = s.substr(0, pos_hash);
    trim(s);
    if (s.empty()) return false;
    auto pos_eq = s.find('=');
    if (pos_eq == std::string::npos) return false;
    key = s.substr(0, pos_eq);
    val = s.substr(pos_eq+1);
    trim(key); trim(val);
    return !key.empty();
}
static int to_int(const std::string& v, int defv) {
    if (v.empty()) return defv;
    return std::atoi(v.c_str());
}
static double to_double(const std::string& v, double defv) {
    if (v.empty()) return defv;
    return std::atof(v.c_str());
}
static bool to_bool(const std::string& v, bool defv) {
    std::string s = lower(v);
    if (s=="1" || s=="true" || s=="yes" || s=="on")  return true;
    if (s=="0" || s=="false"|| s=="no"  || s=="off") return false;
    char* end=nullptr;
    long x = std::strtol(v.c_str(), &end, 10);
    if (end && *end==0) return (x!=0);
    return defv;
}

bool load_precip_global_params(PrecipGlobalParams& gp)
{
    PrecipGlobalParams def;
    gp.enabled            = def.enabled;
    gp.enable_misfit      = def.enable_misfit;
    gp.mu                 = def.mu;
    gp.nu                 = def.nu;
    gp.obb_inflate        = def.obb_inflate;
    gp.aniso              = def.aniso;
    gp.aniso_quad_n       = def.aniso_quad_n;
    gp.aniso_fd_eps       = def.aniso_fd_eps;
    gp.force_on_hit_only  = def.force_on_hit_only;

    // 支持环境变量 PRECIP_CTRL 优先
    const char* envCtrl = std::getenv("PRECIP_CTRL");
    std::string ctrlPath = envCtrl && *envCtrl
                         ? std::string(envCtrl)
                         : (!gp.ctrlFile.empty() ? gp.ctrlFile : std::string(DEFAULT_PRECIP_CTRL));
    gp.ctrlFile = ctrlPath;

    FILE* fp = std::fopen(ctrlPath.c_str(), "r");
    if (!fp) {
        std::fprintf(stderr, "[precip] Warning: cannot open ctrl file: %s (using defaults; enabled=0)\n",
                     ctrlPath.c_str());
        return false;
    }

    char buf[4096];
    while (std::fgets(buf, sizeof(buf), fp)) {
        std::string key, val;
        if (!parse_kv(buf, key, val)) continue;

        if (key.rfind("precip.", 0) != 0) continue;
        std::string k = lower(key);

        if (k == "precip.enabled")            { gp.enabled       = to_bool(val, gp.enabled); continue; }
        if (k == "precip.enable_misfit")      { gp.enable_misfit = to_bool(val, gp.enable_misfit); continue; }
        if (k == "precip.mu")                 { gp.mu            = to_double(val, gp.mu); continue; }
        if (k == "precip.nu")                 { gp.nu            = to_double(val, gp.nu); continue; }
        if (k == "precip.obb_inflate")        { gp.obb_inflate   = to_double(val, gp.obb_inflate); continue; }
        if (k == "precip.input_file")         { gp.inputFile     = val; continue; }

        if (k == "precip.aniso")              { gp.aniso         = to_bool(val, gp.aniso) ? 1 : 0; continue; }
        if (k == "precip.aniso_quad_n")       { gp.aniso_quad_n  = std::max(2, to_int(val, gp.aniso_quad_n)); continue; }
        if (k == "precip.aniso_fd_eps")       { gp.aniso_fd_eps  = std::max(1e-6, to_double(val, gp.aniso_fd_eps)); continue; }
        if (k == "precip.force_on_hit_only")  { gp.force_on_hit_only = to_bool(val, gp.force_on_hit_only) ? 1 : 0; continue; }
    }

    std::fclose(fp);
    return true;
}

static void eulerZYX_deg(double yaw, double pitch, double roll, double U[9])
{
    const double pi = std::acos(-1.0);
    double az = yaw   * pi/180.0;
    double ay = pitch * pi/180.0;
    double ax = roll  * pi/180.0;

    double cz = std::cos(az), sz = std::sin(az);
    double cy = std::cos(ay), sy = std::sin(ay);
    double cx = std::cos(ax), sx = std::sin(ax);

    double Rz[9] = { cz, -sz, 0,  sz,  cz, 0,   0,  0, 1 };
    double Ry[9] = { cy,  0, sy,  0,  1, 0,  -sy, 0, cy };
    double Rx[9] = { 1,   0,  0,  0,  cx, -sx,  0, sx, cx };

    auto mm = [](const double A[9], const double B[9], double C[9]){
        C[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
        C[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
        C[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
        C[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
        C[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
        C[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
        C[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
        C[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
        C[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
    };

    double Rzy[9]; mm(Rz,Ry,Rzy);
    mm(Rzy,Rx,U);
}

bool load_precipitates_dat(const char* path, std::vector<Ellipsoid>& out)
{
    out.clear();
    if (!path || !*path) return false;

    FILE* fp = std::fopen(path, "r");
    if (!fp) {
        std::fprintf(stderr, "[precip] Warning: cannot open precipitates file: %s\n", path);
        return false;
    }

    char line[4096];
    int lineNo = 0;
    while (std::fgets(line, sizeof(line), fp)) {
        lineNo++;
        char* h = std::strchr(line, '#');
        if (h) *h = '\0';

        std::vector<double> tok;
        char* save = nullptr;
        for (char* p = ::strtok_r(line, " \t\r\n", &save); p; p = ::strtok_r(nullptr, " \t\r\n", &save)) {
            char* end=nullptr;
            double v = std::strtod(p, &end);
            if (end && (*end==0)) tok.push_back(v);
            else { tok.clear(); break; }
        }
        if (tok.empty()) continue;

        Ellipsoid e;
        if ((int)tok.size() < 6) {
            std::fprintf(stderr, "[precip] skip line %d: need at least 6 numbers (cx cy cz r1 r2 r3)\n", lineNo);
            continue;
        }
        e.c = { tok[0], tok[1], tok[2] };
        e.r = { tok[3], tok[4], tok[5] };
        e.U = { 1,0,0, 0,1,0, 0,0,1 };

        if ((int)tok.size() == 9) {
            e.eigenstrain_voigt = { tok[6], tok[7], tok[8] };
        } else if ((int)tok.size() >= 12) {
            eulerZYX_deg(tok[6], tok[7], tok[8], e.U.data());
            e.eigenstrain_voigt = { tok[9], tok[10], tok[11] };
        } else {
            e.eigenstrain_voigt = { 0.0, 0.0, 0.0 };
        }

        build_matrices(e);
        out.push_back(e);
    }

    std::fclose(fp);

    if (out.empty()) {
        std::fprintf(stderr, "[precip] Warning: no valid precipitates parsed from: %s\n", path);
        return false;
    }

    std::fprintf(stderr, "[precip] info: parsed %zu precipitates from %s\n", out.size(), path);
    return true;
}