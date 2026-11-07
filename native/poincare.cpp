// Native DirectX 12 renderer for the automorphic function
//     f(z) = P_{-1,-3}(z) / P_{0,-4}(z),   P_{k,l}(z) = sum_{g in G} (az+b)^k (cz+d)^l
// for the Kleinian group G = <g1, g2>.
//
// The Poincaré series are summed per pixel in a compute shader (fp32 or fp64),
// accumulated progressively in a double-precision buffer, and displayed with
// domain colouring. No external dependencies: Win32 + DirectX 12 + the Windows SDK.
//
// Controls
//   drag            pan                 wheel        zoom about cursor
//   + / -           double / halve the number of group elements
//   D               toggle fp32 / fp64 summation
//   M               cycle plot style    F            cycle function (f, P0,-4, P-1,-3)
//   [ / ]           contour density     R            reset view
//   I               invariance check (console)
//   P               save PNG of the current window
//   Esc             quit

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <queue>
#include <string>
#include <vector>

// Ask hybrid-graphics laptops to run this process on the discrete GPU.
extern "C" {
__declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}

static void fatal(const char* msg) {
    std::fprintf(stderr, "\nERROR: %s\n", msg);
    std::fflush(nullptr);
    // std::exit would run static destructors that release D3D objects while
    // GPU work may still be stuck, which can hang the process on the way out.
    TerminateProcess(GetCurrentProcess(), 1);
}

#include "dx12.h"
static DxRenderer renderer;

// ---------------------------------------------------------------------------
// Group: generators and enumeration of reduced words by matrix norm
// ---------------------------------------------------------------------------
using cd = std::complex<double>;
struct Mat { cd a, b, c, d; };

static Mat mul(const Mat& p, const Mat& q) {
    return {p.a * q.a + p.b * q.c, p.a * q.b + p.b * q.d, p.c * q.a + p.d * q.c, p.c * q.b + p.d * q.d};
}
static Mat inv(const Mat& m) { return {m.d, -m.b, -m.c, m.a}; }  // det = 1
static double norm2(const Mat& m) { return std::norm(m.a) + std::norm(m.b) + std::norm(m.c) + std::norm(m.d); }
static cd mobius(const Mat& m, cd z) { return (m.a * z + m.b) / (m.c * z + m.d); }

static const cd B(-0.10997580175045825, -0.0011537600196144);
static const cd Cc(0.07184238006772357, -1.8193355782014326);
static const Mat G1 = {cd(1, 0.1), B, Cc, cd(1, 0.1)};
static const Mat G2 = {cd(1.1, -1), cd(1.1, 0), cd(1.1, 0), cd(1.1, 1)};
static const Mat GENS[4] = {G1, inv(G1), G2, inv(G2)};  // g1, g1^-1, g2, g2^-1
static const int INVERSE[4] = {1, 0, 3, 2};

// All reduced words, smallest Frobenius norm first (G is free on g1, g2, so
// each reduced word is a distinct group element).
static std::vector<Mat> enumerateGroup(size_t N) {
    struct Node {
        double n2;
        Mat m;
        int last;
        bool operator<(const Node& o) const { return n2 > o.n2; }  // min-heap
    };
    std::priority_queue<Node> pq;
    pq.push({2.0, {1, 0, 0, 1}, -1});
    std::vector<Mat> out;
    out.reserve(N);
    while (out.size() < N && !pq.empty()) {
        Node e = pq.top();
        pq.pop();
        out.push_back(e.m);
        for (int k = 0; k < 4; k++) {
            if (e.last >= 0 && k == INVERSE[e.last]) continue;
            Mat m = mul(GENS[k], e.m);
            pq.push({norm2(m), m, k});
        }
    }
    return out;
}

struct Series { cd P1, P2, f; };
static Series evalSeries(const std::vector<Mat>& mats, cd z) {
    cd P1 = 0, P2 = 0;
    for (const Mat& m : mats) {
        cd w1 = m.a * z + m.b, w2 = m.c * z + m.d;
        cd r = 1.0 / (w1 * w2 * w2 * w2 * w2);
        P1 += w2 * r;  // 1 / ((az+b)(cz+d)^3)
        P2 += w1 * r;  // 1 / (cz+d)^4
    }
    return {P1, P2, P1 / P2};
}

// ---------------------------------------------------------------------------
// PNG writer (uncompressed deflate, no dependencies)
// ---------------------------------------------------------------------------
static uint32_t crc32(const uint8_t* d, size_t n, uint32_t c = 0xFFFFFFFFu) {
    static uint32_t T[256];
    static bool init = false;
    if (!init) {
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t r = i;
            for (int k = 0; k < 8; k++) r = (r & 1) ? 0xEDB88320u ^ (r >> 1) : r >> 1;
            T[i] = r;
        }
        init = true;
    }
    for (size_t i = 0; i < n; i++) c = T[(c ^ d[i]) & 0xFF] ^ (c >> 8);
    return c;
}
static bool writePNG(const char* path, int w, int h, const std::vector<uint8_t>& rgbTopDown) {
    std::vector<uint8_t> raw;
    raw.reserve(size_t(h) * (1 + 3 * w));
    for (int y = 0; y < h; y++) {
        raw.push_back(0);
        raw.insert(raw.end(), rgbTopDown.begin() + size_t(y) * 3 * w, rgbTopDown.begin() + size_t(y + 1) * 3 * w);
    }
    std::vector<uint8_t> z = {0x78, 0x01};
    uint32_t a = 1, b = 0;
    for (uint8_t v : raw) { a = (a + v) % 65521; b = (b + a) % 65521; }
    for (size_t off = 0; off < raw.size() || off == 0; off += 65535) {
        size_t n = std::min<size_t>(65535, raw.size() - off);
        bool last = off + n >= raw.size();
        z.push_back(last ? 1 : 0);
        z.push_back(uint8_t(n & 0xFF)); z.push_back(uint8_t(n >> 8));
        z.push_back(uint8_t(~n & 0xFF)); z.push_back(uint8_t((~n >> 8) & 0xFF));
        z.insert(z.end(), raw.begin() + off, raw.begin() + off + n);
        if (last) break;
    }
    uint32_t adler = (b << 16) | a;
    for (int s = 24; s >= 0; s -= 8) z.push_back(uint8_t(adler >> s));

    FILE* f = std::fopen(path, "wb");
    if (!f) return false;
    auto be32 = [&](uint32_t v) { uint8_t t[4] = {uint8_t(v >> 24), uint8_t(v >> 16), uint8_t(v >> 8), uint8_t(v)}; std::fwrite(t, 1, 4, f); };
    auto chunk = [&](const char* type, const std::vector<uint8_t>& data) {
        be32(uint32_t(data.size()));
        std::vector<uint8_t> td(type, type + 4);
        td.insert(td.end(), data.begin(), data.end());
        std::fwrite(td.data(), 1, td.size(), f);
        be32(crc32(td.data(), td.size()) ^ 0xFFFFFFFFu);
    };
    const uint8_t sig[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    std::fwrite(sig, 1, 8, f);
    std::vector<uint8_t> ihdr = {uint8_t(w >> 24), uint8_t(w >> 16), uint8_t(w >> 8), uint8_t(w),
                                 uint8_t(h >> 24), uint8_t(h >> 16), uint8_t(h >> 8), uint8_t(h), 8, 2, 0, 0, 0};
    chunk("IHDR", ihdr);
    chunk("IDAT", z);
    chunk("IEND", {});
    std::fclose(f);
    return true;
}

// ---------------------------------------------------------------------------
// Application state
// ---------------------------------------------------------------------------
using Clock = std::chrono::steady_clock;
static double msSince(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

struct App {
    HWND hwnd = nullptr;
    int W = 0, H = 0;
    double cx = 0.0, cy = -0.35, span = 5.0;  // span = width of the shorter screen side in z-units
    size_t N = size_t(1) << 17;
    bool fp64 = false;
    int mode = 0, func = 0;
    float dens = 1.0f;

    std::vector<Mat> mats;
    int count = 0, done = 0, chunk = 256;
    Clock::time_point startT;
    double renderMs = 0;

    bool dragging = false;
    int mx = 0, my = 0;
    std::string gpu;
    Clock::time_point titleT = Clock::now();
    bool saveRequested = false;

    double scale() const { return span / std::max(1, std::min(W, H)); }
    cd toZ(int x, int y) const {  // window pixel (top-left origin) -> z
        double s = scale();
        return {cx + (x + 0.5 - 0.5 * W) * s, cy + ((H - 1 - y) + 0.5 - 0.5 * H) * s};
    }
};
static App app;

static void restart() {
    app.done = 0;
    app.startT = Clock::now();
}

static void buildGroup() {
    auto t0 = Clock::now();
    app.mats = enumerateGroup(app.N);
    app.count = int(app.mats.size());
    std::vector<float> f(app.mats.size() * 8);
    std::vector<double> d(app.mats.size() * 8);
    for (size_t i = 0; i < app.mats.size(); i++) {
        const Mat& m = app.mats[i];
        double v[8] = {m.a.real(), m.a.imag(), m.b.real(), m.b.imag(), m.c.real(), m.c.imag(), m.d.real(), m.d.imag()};
        for (int k = 0; k < 8; k++) { d[8 * i + k] = v[k]; f[8 * i + k] = float(v[k]); }
    }
    renderer.matrices(f, d);
    std::printf("group: %d elements, max ||g|| = %.1f  (enumerated in %.0f ms)\n", app.count,
                std::sqrt(norm2(app.mats.back())), msSince(t0));
    restart();
}

static void resize(int w, int h) {
    if (w <= 0 || h <= 0 || (w == app.W && h == app.H)) return;
    app.W = w;
    app.H = h;
    renderer.resize(w, h);
    restart();
}

static void invarianceCheck() {
    cd z0(0.3, 0.6);
    cd fz = evalSeries(app.mats, z0).f;
    std::printf("\ninvariance check with %d elements, z0 = 0.3+0.6i\n  f(z0)        = %+.8f %+.8fi\n", app.count, fz.real(), fz.imag());
    const char* names[4] = {"g1", "g1^-1", "g2", "g2^-1"};
    for (int k = 0; k < 4; k++) {
        cd fg = evalSeries(app.mats, mobius(GENS[k], z0)).f;
        std::printf("  f(%-5s(z0)) = %+.8f %+.8fi   rel.err %.1e\n", names[k], fg.real(), fg.imag(), std::abs(fg - fz) / std::abs(fz));
    }
}

static void updateTitle() {
    if (msSince(app.titleT) < 120) return;
    app.titleT = Clock::now();
    // read the accumulated sums under the cursor straight from the GPU buffer
    double a[4] = {0, 0, 0, 0};
    int px = std::clamp(app.mx, 0, app.W - 1), py = std::clamp(app.my, 0, app.H - 1);
    renderer.readPixel(px, py, a);
    cd P1(a[0], a[1]), P2(a[2], a[3]);
    cd val = app.func == 0 ? P1 / P2 : (app.func == 1 ? P2 : P1);
    cd z = app.toZ(app.mx, app.my);
    const char* fn[3] = {"f", "P0,-4", "P-1,-3"};
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "Poincaré series  |  %s  |  %s  |  terms %d/%d %s |  z = %.6f%+.6fi   %s(z) = %.6g%+.6gi  |  width %.3g",
                  app.gpu.c_str(), app.fp64 ? "fp64" : "fp32", app.done, app.count,
                  app.done < app.count ? "" : (std::string("(") + std::to_string(int(app.renderMs)) + " ms) ").c_str(),
                  z.real(), z.imag(), fn[app.func], val.real(), val.imag(), app.scale() * app.W);
    wchar_t wbuf[512];
    MultiByteToWideChar(CP_UTF8, 0, buf, -1, wbuf, 512);
    SetWindowTextW(app.hwnd, wbuf);
}

static std::string snapshotPath;  // --snapshot <file.png>: render to completion, save, exit
static bool verifyGPU = false;

static void verifyGpuValues() {
    const int points[][2] = {{app.W/2,app.H/2}, {app.W/3,app.H/4}, {3*app.W/4,2*app.H/3}};
    double worst = 0;
    for (const auto& p : points) {
        double sums[4];
        renderer.readPixel(p[0],p[1],sums);
        auto expected = evalSeries(app.mats,app.toZ(p[0],p[1]));
        cd actual[] = {{sums[0],sums[1]},{sums[2],sums[3]}};
        cd reference[] = {expected.P1,expected.P2};
        for (int i=0;i<2;++i) {
            double error=std::abs(actual[i]-reference[i])/std::max(1e-30,std::abs(reference[i]));
            std::printf("GPU check (%d,%d) P%d: relative error %.3e\n",p[0],p[1],i+1,error);
            if (!std::isfinite(error)) fatal("Non-finite GPU verification result");
            worst=std::max(worst,error);
        }
    }
    if (worst > (app.fp64 ? 1e-9 : 5e-3)) fatal("GPU sums do not match CPU reference");
    std::printf("GPU verification passed (%s)\n",app.fp64 ? "fp64" : "fp32");
}

static void savePNG(const char* path = nullptr) {
    auto rgb = renderer.pixels();
    SYSTEMTIME t;
    GetLocalTime(&t);
    char name[128];
    std::snprintf(name, sizeof name, "poincare_%04d%02d%02d_%02d%02d%02d.png", t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond);
    const char* output = path ? path : name;
    if (!writePNG(output, app.W, app.H, rgb)) fatal("Could not write PNG snapshot");
    std::printf("saved %s (%dx%d, %d/%d terms)\n", output, app.W, app.H, app.done, app.count);
}

static ViewConstants viewConstants(int start, int end) {
    return {UINT(start), UINT(end), UINT(app.W), UINT(app.H), app.cx, app.cy, app.scale(),
            UINT(app.mode), UINT(app.func), app.dens};
}

static void frame() {
    if (app.done < app.count) {
        auto t0 = Clock::now();
        while (app.done < app.count && msSince(t0) < 30.0) {
            int end = int(std::min<int64_t>(app.count, int64_t(app.done) + app.chunk));
            auto td = Clock::now();
            renderer.compute(viewConstants(app.done, end), app.fp64);
            double dt = msSince(td);
            app.done = end;
            if (dt < 4.0) app.chunk = std::min(app.chunk * 2, 1 << 16);
            else if (dt > 12.0) app.chunk = std::max(app.chunk / 2, 16);
        }
        if (app.done >= app.count) {
            app.renderMs = msSince(app.startT);
            std::printf("render complete: %d terms, %dx%d, %s, %.0f ms\n", app.count, app.W, app.H,
                        app.fp64 ? "fp64" : "fp32", app.renderMs);
        }
    }
    renderer.draw(viewConstants(0, app.done));
    if (app.saveRequested) { savePNG(); app.saveRequested = false; }
    if (!snapshotPath.empty() && app.done >= app.count) {
        if (verifyGPU) verifyGpuValues();
        invarianceCheck();
        savePNG(snapshotPath.c_str());
        PostQuitMessage(0);
        snapshotPath.clear();
    }
    renderer.present(app.done >= app.count);
    updateTitle();
}

// ---------------------------------------------------------------------------
// Window
// ---------------------------------------------------------------------------
static void zoomAt(int x, int y, double k) {
    cd z0 = app.toZ(x, y);
    app.span *= k;
    cd z1 = app.toZ(x, y);
    app.cx += z0.real() - z1.real();
    app.cy += z0.imag() - z1.imag();
    restart();
}

static LRESULT CALLBACK wndProc(HWND h, UINT msg, WPARAM wp, LPARAM lp) {
    switch (msg) {
    case WM_SIZE:
        if (renderer.ready()) resize(LOWORD(lp), HIWORD(lp));
        return 0;
    case WM_LBUTTONDOWN:
        app.dragging = true;
        app.mx = short(LOWORD(lp)), app.my = short(HIWORD(lp));
        SetCapture(h);
        return 0;
    case WM_LBUTTONUP:
        app.dragging = false;
        ReleaseCapture();
        return 0;
    case WM_MOUSEMOVE: {
        int x = short(LOWORD(lp)), y = short(HIWORD(lp));
        if (app.dragging) {
            double s = app.scale();
            app.cx -= (x - app.mx) * s;
            app.cy += (y - app.my) * s;
            restart();
        }
        app.mx = x, app.my = y;
        return 0;
    }
    case WM_MOUSEWHEEL: {
        POINT p = {short(LOWORD(lp)), short(HIWORD(lp))};
        ScreenToClient(h, &p);
        zoomAt(p.x, p.y, std::exp(-GET_WHEEL_DELTA_WPARAM(wp) / 120.0 * 0.18));
        return 0;
    }
    case WM_CHAR:
        switch (wp) {
        case '+': case '=': app.N = std::min<size_t>(app.N * 2, size_t(1) << 22); buildGroup(); break;
        case '-': case '_': app.N = std::max<size_t>(app.N / 2, 256); buildGroup(); break;
        case 'd': case 'D': app.fp64 = !app.fp64; app.chunk = 64; restart(); break;
        case 'm': case 'M': app.mode = (app.mode + 1) % 4; break;
        case 'f': case 'F': app.func = (app.func + 1) % 3; break;
        case '[': app.dens = std::max(0.25f, app.dens / 1.25f); break;
        case ']': app.dens = std::min(8.0f, app.dens * 1.25f); break;
        case 'r': case 'R': app.cx = 0; app.cy = -0.35; app.span = 5.0; restart(); break;
        case 'i': case 'I': invarianceCheck(); break;
        case 'p': case 'P': app.saveRequested = true; break;
        case 27: PostQuitMessage(0); break;
        }
        return 0;
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcW(h, msg, wp, lp);
}

int main(int argc, char** argv) {
    SetProcessDPIAware();
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    auto readNumber = [](const char* text, double& value) {
        char* end = nullptr;
        errno = 0;
        value = std::strtod(text, &end);
        return end != text && *end == '\0' && errno != ERANGE && std::isfinite(value);
    };
    bool help = false;
    for (int i = 1; i < argc; i++) {
        const std::string option = argv[i];
        if (option == "--help" || option == "-h") {
            help = true;
        } else if (option == "--fp64") {
            app.fp64 = true;
        } else if (option == "--verify-gpu") {
            verifyGPU = true;
        } else if (option == "--center") {
            if (i + 2 >= argc || !readNumber(argv[i + 1], app.cx) || !readNumber(argv[i + 2], app.cy)) {
                std::fprintf(stderr, "--center requires two finite numbers: <real> <imag>\n");
                return 1;
            }
            i += 2;
        } else if (option == "--span") {
            if (i + 1 >= argc || !readNumber(argv[i + 1], app.span) || app.span <= 0) {
                std::fprintf(stderr, "--span requires a finite number greater than zero\n");
                return 1;
            }
            ++i;
        } else if (option == "--snapshot") {
            if (i + 1 >= argc || argv[i + 1][0] == '\0' || std::string(argv[i + 1]).rfind("--", 0) == 0) {
                std::fprintf(stderr, "--snapshot requires a file path\n");
                return 1;
            }
            snapshotPath = argv[++i];
        } else if (option == "--terms") {
            char* end = nullptr;
            errno = 0;
            const char* text = i + 1 < argc ? argv[i + 1] : "";
            const auto count = std::strtoull(text, &end, 10);
            if (*text < '0' || *text > '9' || *end != '\0' || errno == ERANGE ||
                count == 0 || count > (std::numeric_limits<int>::max)()) {
                std::fprintf(stderr, "--terms requires a positive integer no greater than INT_MAX\n");
                return 1;
            }
            app.N = static_cast<size_t>(count);
            ++i;
        } else {
            std::fprintf(stderr, "Unknown option: %s (use --help)\n", option.c_str());
            return 1;
        }
    }
    if (help) {
        std::printf(
            "Usage: poincare.exe [options]\n"
            "  --center <real> <imag>  Plot center (default: 0 -0.35)\n"
            "  --span <value>          Shorter-side extent in complex units (default: 5)\n"
            "  --terms <count>         Number of group elements (default: 131072)\n"
            "  --fp64                  Use double-precision summation\n"
            "  --snapshot <file.png>   Render, save, and exit\n"
            "  --verify-gpu            Check GPU sums against CPU when saving a snapshot\n"
            "  --help, -h              Show this help\n");
        return 0;
    }
    std::printf(
        "Poincaré series  f(z) = P_{-1,-3}(z) / P_{0,-4}(z),  G = <g1, g2>\n"
        "  drag: pan   wheel: zoom   +/-: elements x2 / /2   D: fp32/fp64   M: plot style\n"
        "  F: function   [ ]: contour density   R: reset   I: invariance check   P: save PNG   Esc: quit\n\n");

    HINSTANCE hi = GetModuleHandleW(nullptr);
    WNDCLASSW wc = {};
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = wndProc;
    wc.hInstance = hi;
    wc.hCursor = LoadCursor(nullptr, IDC_CROSS);
    wc.lpszClassName = L"PoincareSeries";
    RegisterClassW(&wc);
    RECT r = {0, 0, 1400, 1000};
    AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE);
    app.hwnd = CreateWindowW(wc.lpszClassName, L"Poincaré series", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
                             r.right - r.left, r.bottom - r.top, nullptr, nullptr, hi, nullptr);
    RECT cr;
    GetClientRect(app.hwnd, &cr);
    app.gpu = renderer.initialize(app.hwnd, cr.right, cr.bottom);
    std::printf("DirectX 12\nGPU    %s\n\n", app.gpu.c_str());
    resize(cr.right, cr.bottom);
    buildGroup();
    ShowWindow(app.hwnd, SW_SHOW);

    MSG msg;
    for (;;) {
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) return 0;
            TranslateMessage(&msg);
            DispatchMessageW(&msg);
        }
        if (IsIconic(app.hwnd)) { WaitMessage(); continue; }
        frame();
    }
}
