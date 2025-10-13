#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <climits>  // LONG_MIN

#include "precip_params.h"
#include "precip_io.h"
#include "ellipsoid.h"
#include "collision.h"
#include "looping.h"
#include "misfit_quad.h"
#include "misfit_aniso.h"
#include "Home.h"  // 读取 Home_t->cycle

struct _home;

// 声明两种真实 NodeForce（不同签名），用弱引用避免缺失时报链接错误
extern "C" void __real_NodeForce2(_home*, int)   asm("__real__Z9NodeForceP5_homei") __attribute__((weak));
extern "C" void __real_NodeForce1(_home*)        asm("__real__Z9NodeForceP5_home")  __attribute__((weak));

// 为两种签名分别提供 wrap（名称需与 --wrap 选项一致）
extern "C" void __wrap_NodeForce2(_home*, int)   asm("__wrap__Z9NodeForceP5_homei");
extern "C" void __wrap_NodeForce1(_home*)        asm("__wrap__Z9NodeForceP5_home");

// 其他已有 wrap（保持不变）
extern "C" void __real_HandleCollisions(_home*)  asm("__real__Z16HandleCollisionsP5_home");
extern "C" void __real_Remesh(_home*)            asm("__real__Z6RemeshP5_home");

extern "C" void __wrap_HandleCollisions(_home*)  asm("__wrap__Z16HandleCollisionsP5_home");
extern "C" void __wrap_Remesh(_home*)            asm("__wrap__Z6RemeshP5_home");

// 段遍历/注力 shims
extern "C" {
void Precip_ForEachSegment(_home*,
    void (*cb)(const double x1[3], const double x2[3],
               const double b[3],  const double xi[3],
               void* user),
    void* user) __attribute__((weak));

void Precip_AddForceToSegmentEndpoints(_home*,
    const double x1[3], const double x2[3],
    const double f1[3], const double f2[3]) __attribute__((weak));
}

static inline bool str_true_env(const char* e) {
    return (e && (e[0]=='1'||e[0]=='t'||e[0]=='T'||e[0]=='y'||e[0]=='Y'));
}
static bool is_abs_path(const std::string& p) {
#if defined(_WIN32) || defined(_WIN64)
    if (p.size() > 1 && std::isalpha((unsigned char)p[0]) && p[1] == ':') return true;
    if (!p.empty() && (p[0] == '\\' || p[0] == '/')) return true;
    return false;
#else
    return (!p.empty() && p[0] == '/');
#endif
}
static std::string dirname_of(const std::string& p) {
#if defined(_WIN32) || defined(_WIN64)
    const char* seps = "/\\";
#else
    const char* seps = "/";
#endif
    size_t pos = p.find_last_of(seps);
    return (pos == std::string::npos) ? std::string(".") : p.substr(0, pos);
}

static void set3(std::array<double,3>& a, double x, double y, double z) { a[0]=x; a[1]=y; a[2]=z; }

static inline void mat_to_voigt(const double M[3][3], std::array<double,6>& v) {
    v[0]=M[0][0]; v[1]=M[1][1]; v[2]=M[2][2];
    v[3]=M[1][2]; v[4]=M[2][0]; v[5]=M[0][1];
}
static inline void mat33_mul(const double* U, const double E[3][3], double UE[3][3]) {
    double A[3][3]={{U[0],U[1],U[2]},{U[3],U[4],U[5]},{U[6],U[7],U[8]}};
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        UE[i][j]=A[i][0]*E[0][j]+A[i][1]*E[1][j]+A[i][2]*E[2][j];
}
static inline void mat33_mulT(const double A[3][3], const double* UT, double C[3][3]) {
    double U[3][3]={{UT[0],UT[3],UT[6]},{UT[1],UT[4],UT[7]},{UT[2],UT[5],UT[8]}};
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        C[i][j]=A[i][0]*U[0][j]+A[i][1]*U[1][j]+A[i][2]*U[2][j];
}
static inline void rotate_voigt_by_U(const double* U,
                                     const std::array<double,6>& v_local,
                                     std::array<double,6>& v_global)
{
    double ELoc[3][3]; voigt_to_mat(v_local, ELoc); // misfit_quad.h
    double tmp[3][3];  mat33_mul(U, ELoc, tmp);
    double EGlob[3][3]; mat33_mulT(tmp, U, EGlob);
    mat_to_voigt(EGlob, v_global);
}

namespace {

struct PrecipRuntime {
    bool inited = false;
    bool enabled = false;
    bool enableMisfit = false;

    bool   aniso = false;
    int    aniso_quad_n = 3;
    double aniso_fd_eps = 1e-3;
    bool   forceOnHitOnly = false;

    std::string ctrlFile = DEFAULT_PRECIP_CTRL;
    std::string inputFile;

    Material mat{0.0, 0.3};
    std::vector<Ellipsoid> precips;

    double obbInflate = 0.0;
    int    verboseEvery = 100;

    double lastInjectedL1 = 0.0;

    bool debugWrap = false;
} G;

static void lazy_init() {
    if (G.inited) return;
    G.inited = true;

    PrecipGlobalParams gp{};
    (void)load_precip_global_params(gp);

    if (!gp.ctrlFile.empty()) G.ctrlFile = gp.ctrlFile;
    G.enabled      = gp.enabled;
    G.enableMisfit = gp.enable_misfit;
    G.mat.mu       = gp.mu;
    G.mat.nu       = gp.nu;
    G.obbInflate   = gp.obb_inflate;

    G.aniso          = (gp.aniso != 0);
    G.aniso_quad_n   = (gp.aniso_quad_n > 1 ? gp.aniso_quad_n : 3);
    G.aniso_fd_eps   = (gp.aniso_fd_eps > 0 ? gp.aniso_fd_eps : 1e-3);
    G.forceOnHitOnly = (gp.force_on_hit_only != 0);

    if (const char* e = std::getenv("PRECIP_ANISO"))              G.aniso = str_true_env(e);
    if (const char* e = std::getenv("PRECIP_ANISO_QUAD_N"))       { int v=std::atoi(e); if (v>1) G.aniso_quad_n=v; }
    if (const char* e = std::getenv("PRECIP_ANISO_FD_EPS"))       { double v=std::atof(e); if (v>0) G.aniso_fd_eps=v; }
    if (const char* e = std::getenv("PRECIP_FORCE_ON_HIT_ONLY"))  G.forceOnHitOnly = str_true_env(e);
    if (const char* e = std::getenv("PRECIP_VERBOSE_EVERY"))      { int v=std::atoi(e); if (v>0) G.verboseEvery=v; }
    if (const char* e = std::getenv("PRECIP_DEBUG_WRAP"))         G.debugWrap = str_true_env(e);

    // 解析 inputFile 相对 ctrl 目录
    std::string inFile = gp.inputFile;
    if (!is_abs_path(inFile) && !G.ctrlFile.empty()) {
        std::string base = dirname_of(G.ctrlFile);
        if (!base.empty() && base != ".") {
            if (inFile.size() >= 2 && (inFile.compare(0, 2, "./") == 0 || inFile.compare(0, 2, ".\\") == 0)) {
                inFile = base + "/" + inFile.substr(2);
            } else if (!inFile.empty()) {
                inFile = base + "/" + inFile;
            }
        }
    }
    G.inputFile = inFile;

    if (G.enabled) {
        if (!load_precipitates_dat(G.inputFile.c_str(), G.precips)) {
            std::fprintf(stderr, "[precip] Warning: cannot open precip data file: %s (ctrl=%s)\n",
                         G.inputFile.c_str(), G.ctrlFile.c_str());
            return;
        }
        std::fprintf(stderr,
            "[precip] info: loaded %zu precipitates from %s (mu=%.3g, nu=%.3g, misfit=%s)\n",
            G.precips.size(), G.inputFile.c_str(), G.mat.mu, G.mat.nu, G.enableMisfit ? "on" : "off");
        std::fprintf(stderr, "[precip] info: ctrl = %s\n", G.ctrlFile.c_str());
        std::fprintf(stderr,
            "[precip] info: aniso=%s, quadN=%d, fd_eps=%.2e, force_on_hit_only=%s, obb_inflate=%.2f\n",
            (G.aniso ? "on":"off"), G.aniso_quad_n, G.aniso_fd_eps, (G.forceOnHitOnly ? "on":"off"), G.obbInflate);
    } else {
        std::fprintf(stderr, "[precip] info: precip.enabled=0 or ctrl not loaded (ctrl=%s)\n", G.ctrlFile.c_str());
    }
}

struct SegEvalCtx {
    _home* home=nullptr;
    size_t near_hits=0,true_hits=0,gauss_forces=0;
    double injectedL1=0.0;
};

static void eval_one_segment_against_precips(
    const double x1[3], const double x2[3],
    const double b[3],  const double xi[3],
    void* user)
{
    SegEvalCtx* C = reinterpret_cast<SegEvalCtx*>(user);

    Segment seg;
    for (int i=0;i<3;i++) { seg.x1[i]=x1[i]; seg.x2[i]=x2[i]; }

    for (const auto& e : G.precips) {
        if (!segmentOBBMayCollide(e, seg, G.obbInflate)) continue;
        C->near_hits++;

        double alpha=0.0;
        bool hit = segmentEllipsoidIntersectAlpha(e, seg, alpha);
        if (hit) C->true_hits++;

        if (!G.enableMisfit) continue;
        if (G.forceOnHitOnly && !hit) continue;

        auto sigma_at = [&](const std::array<double,3>& Xp, std::array<double,6>& sv) {
            const double eabs = std::fabs(e.eigenstrain_voigt[0])
                              + std::fabs(e.eigenstrain_voigt[1])
                              + std::fabs(e.eigenstrain_voigt[2]);
            if (eabs < 1e-20) { sv = {0,0,0,0,0,0}; return; }

            if (G.aniso) {
                std::array<double,6> e_local = { e.eigenstrain_voigt[0], e.eigenstrain_voigt[1], e.eigenstrain_voigt[2], 0.0, 0.0, 0.0 };
                std::array<double,6> e_global;
                rotate_voigt_by_U(e.U.data(), e_local, e_global);
                misfitStressAtPoint_Aniso_C(C->home, e, Xp, e_global, sv,
                                            std::max(2, G.aniso_quad_n),
                                            std::max(1e-6, G.aniso_fd_eps));
            } else {
                misfitStressAtPoint(e, G.mat, Xp, sv);
            }
        };

        std::array<double,3> b3{b[0],b[1],b[2]};
        std::array<double,3> xi3{xi[0],xi[1],xi[2]};
        std::array<double,3> X1{x1[0],x1[1],x1[2]};
        std::array<double,3> X2{x2[0],x2[1],x2[2]};

        std::array<double,3> f1{};
        nodeForceFromMisfitOnSegment(X1, X2, b3, xi3, sigma_at, f1);
        double f2[3] = {-f1[0], -f1[1], -f1[2]};

        if (Precip_AddForceToSegmentEndpoints) {
            Precip_AddForceToSegmentEndpoints(C->home, x1, x2, f1.data(), f2);
        }
        C->gauss_forces++;
        C->injectedL1 += std::fabs(f1[0])+std::fabs(f1[1])+std::fabs(f1[2])
                       + std::fabs(f2[0])+std::fabs(f2[1])+std::fabs(f2[2]);
    }
}

static void maybe_scan_segments(_home* home, const char* where) {
    if (!G.enabled || G.precips.empty()) return;
    if (!Precip_ForEachSegment) {
        static bool told=false;
        if (!told) { std::fprintf(stderr,
            "[precip] note: skipping segment scan (no shim Precip_ForEachSegment provided)\n"); told=true; }
        return;
    }
    SegEvalCtx ctx; ctx.home=home;
    Precip_ForEachSegment(home, &eval_one_segment_against_precips, &ctx);
    G.lastInjectedL1 = ctx.injectedL1;
    std::fprintf(stderr,
        "[precip] NodeForce: near=%zu, hit=%zu, misfit-forced=%zu, injected|f|_L1=%.3e (precips=%zu)\n",
        ctx.near_hits, ctx.true_hits, ctx.gauss_forces, ctx.injectedL1, G.precips.size());

    static unsigned long tick=0; if ((++tick % (unsigned long)G.verboseEvery)==1) {
        std::fprintf(stderr, "[precip] %s: near=%zu, hit=%zu, misfit-forced=%zu, injected|f|_L1=%.3e (precips=%zu)\n",
            where, ctx.near_hits, ctx.true_hits, ctx.gauss_forces, ctx.injectedL1, G.precips.size());
    }
}

// 每步只注力一次：依据 home->cycle 判定
static void scan_once_per_step(_home* home, const char* where) {
    static long lastCycle = LONG_MIN;
    Home_t* H = reinterpret_cast<Home_t*>(home);
    long cyc = (H ? (long)H->cycle : LONG_MIN+1);

    if (G.debugWrap) {
        std::fprintf(stderr, "[precip] wrap: NodeForce called, cycle=%ld (%s)\n", cyc, where);
    }

    if (cyc != lastCycle) {
        maybe_scan_segments(home, "NodeForce");
        lastCycle = cyc;
    }
}

} // namespace

// 二参 NodeForce 的 wrap：先调真实，再注力（若一参存在则兜底）
extern "C" void __wrap_NodeForce2(_home* home, int node) {
    lazy_init();
    if (__real_NodeForce2)      __real_NodeForce2(home, node);
    else if (__real_NodeForce1) __real_NodeForce1(home);
    else                        std::fprintf(stderr, "[precip] warn: no real NodeForce symbol found (2-arg wrap)\n");
    scan_once_per_step(home, "sig=(_home*,int)");
}

// 一参 NodeForce 的 wrap：先调真实，再注力（若二参存在则兜底）
extern "C" void __wrap_NodeForce1(_home* home) {
    lazy_init();
    if (__real_NodeForce1)      __real_NodeForce1(home);
    else if (__real_NodeForce2) __real_NodeForce2(home, 0);
    else                        std::fprintf(stderr, "[precip] warn: no real NodeForce symbol found (1-arg wrap)\n");
    scan_once_per_step(home, "sig=(_home*)");
}

extern "C" void __wrap_HandleCollisions(_home* home) {
    lazy_init();
    __real_HandleCollisions(home);
}
extern "C" void __wrap_Remesh(_home* home) {
    lazy_init();
    __real_Remesh(home);
}

extern "C" {
size_t Precip_Count() {
    lazy_init();
    return (G.enabled ? G.precips.size() : 0);
}
int Precip_GetCenter(size_t i, double out_c[3]) {
    lazy_init();
    if (!G.enabled || i>=G.precips.size() || !out_c) return 0;
    out_c[0]=G.precips[i].c[0]; out_c[1]=G.precips[i].c[1]; out_c[2]=G.precips[i].c[2];
    return 1;
}
int Precip_GetAxes(size_t i, double out_r[3], double out_U[9]) {
    lazy_init();
    if (!G.enabled || i>=G.precips.size()) return 0;
    if (out_r) { out_r[0]=G.precips[i].r[0]; out_r[1]=G.precips[i].r[1]; out_r[2]=G.precips[i].r[2]; }
    if (out_U) { for (int k=0;k<9;k++) out_U[k]=G.precips[i].U[k]; }
    return 1;
}
}