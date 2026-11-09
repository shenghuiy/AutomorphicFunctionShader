#pragma once
// Shared native application: Win32 window, DirectX 12 renderer, progressive
// GPU summation, domain colouring, PNG export and command-line handling.
//
// A program describes its group with a `Group` and calls runApp(group, argc, argv).
// Each pixel accumulates sum_i term(w1_i, w2_i, i) over the uploaded matrices into
// a double4 (two complex sums); the group then says what value to colour.
//
// Controls
//   drag            pan                 wheel        zoom about cursor
//   + / -           double / halve the number of group elements (infinite groups)
//   D               toggle fp32 / fp64 summation
//   M               cycle plot style    F            cycle displayed function
//   [ / ]           contour density     R            reset view
//   S               plane / Riemann sphere view (groups with Group::sphere);
//                   in the sphere view, drag rotates and the wheel zooms
//   I               invariance check (console)
//   P               save PNG of the current window
//   Esc             quit

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
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

#include "mobius.h"
#include "png.h"
#include "dx12.h"
#include "shaders.h"

// Everything a program supplies about its group.
struct Group {
    const char* title;          // window title and console header, e.g. "Poincaré series"
    const char* formula;        // one-line description printed at start-up
    const char* exeName;        // for --help
    const char* pngPrefix;      // default PNG file name prefix
    const char* funcNames[3];   // names of the three displayed functions (F key)
    const char* sumNames[2];    // names of the two accumulated sums (GPU check output)
    double cx, cy, span;        // default view: centre and shorter-side extent
    size_t defaultTerms;        // 0: finite group, fixed size; otherwise --terms and +/-
    const char* hlsl;           // term / plotValue / FADE_POLES, see shaders.h

    // Build the group with (up to) `terms` elements, print a summary, and
    // return the matrices to upload. Their entries define w1 and w2 in the shader.
    std::vector<Mat> (*build)(size_t terms);
    // CPU reference for the two accumulated sums at z (GPU verification).
    void (*reference)(cd z, cd out[2]);
    // Displayed value from the two accumulated sums; func selects funcNames[func].
    cd (*value)(cd sum0, cd sum1, int func, cd z);
    // Print an invariance check to the console.
    void (*invariance)();
    // Offer the Riemann sphere view (S key, --sphere).
    bool sphere = false;
};

using Clock = std::chrono::steady_clock;
static double msSince(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

struct App {
    const Group* group = nullptr;
    HWND hwnd = nullptr;
    int W = 0, H = 0;
    double cx = 0.0, cy = 0.0, span = 5.0;  // span = width of the shorter screen side in z-units
    // Riemann sphere view: rot[j] is the world-space image of view axis j
    // (x right, y up, z toward the viewer); zoom 1 fits the sphere in the window.
    bool sphere = false;
    double rot[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double sphereZoom = 1.0;
    size_t N = 0;
    bool fp64 = false;
    int mode = 0, func = 0;
    float dens = 1.0f;

    int count = 0, done = 0, chunk = 256;
    Clock::time_point startT;
    double renderMs = 0;

    bool dragging = false;
    int mx = 0, my = 0;
    std::string gpu;
    Clock::time_point titleT = Clock::now();
    bool saveRequested = false;

    double scale() const { return span / std::max(1, std::min(W, H)); }
    double radius() const { return 0.45 * std::max(1, std::min(W, H)) * sphereZoom; }
    // Window pixel (top-left origin) -> z; matches pixelZ in shaders.h. NaN off the sphere.
    cd toZ(int x, int y) const {
        double ox = x + 0.5 - 0.5 * W, oy = (H - 1 - y) + 0.5 - 0.5 * H;
        if (!sphere) return {cx + ox * scale(), cy + oy * scale()};
        double R = radius(), R2 = R * R, o2 = ox * ox + oy * oy;
        if (o2 >= R2) return {NAN, NAN};
        double S = std::sqrt(R2 - o2), V[3];
        for (int i = 0; i < 3; i++) V[i] = ox * rot[0][i] + oy * rot[1][i] + S * rot[2][i];
        return sphereToZ(V, R);
    }
    // Stereographic projection from the south pole of the point V on the sphere
    // of radius R: 0 at the north pole, infinity at the south pole.
    static cd sphereToZ(const double V[3], double R = 1) {
        if (V[2] >= 0) return cd(V[0], V[1]) / (R + V[2]);
        return cd(V[0], V[1]) * ((R - V[2]) / std::max(V[0] * V[0] + V[1] * V[1], 1e-300));
    }
};
static App app;
static DxRenderer renderer;
static std::string snapshotPath;  // --snapshot <file.png>: render to completion, save, exit
static bool verifyGPU = false;

static void restart() {
    app.done = 0;
    app.startT = Clock::now();
}

static void normalize3(double v[3]) {
    double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    for (int i = 0; i < 3; i++) v[i] /= n;
}
static void cross3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}
// Re-orthonormalise the rotation (Gram-Schmidt) so repeated drags do not drift.
static void orthonormalize() {
    double *x = app.rot[0], *y = app.rot[1], *z = app.rot[2];
    normalize3(z);
    double d = x[0] * z[0] + x[1] * z[1] + x[2] * z[2];
    for (int i = 0; i < 3; i++) x[i] -= d * z[i];
    normalize3(x);
    cross3(z, x, y);
}

// Turn the sphere so the point z faces the viewer with Re z to the right and
// Im z up. The projection is conformal and orientation preserving, so the
// neighbourhood of z looks like the plane view around the same centre.
static void faceSphereToward(cd z) {
    double x = z.real(), y = z.imag(), r2 = std::norm(z);
    if (!std::isfinite(r2)) return;
    double *right = app.rot[0], *front = app.rot[2];
    front[0] = 2 * x / (1 + r2); front[1] = 2 * y / (1 + r2); front[2] = (1 - r2) / (1 + r2);
    // d/dx of the inverse projection, times (1 + r2)^2 / 2, divided by max(1, r2) against overflow
    double s = std::max(1.0, r2);
    right[0] = (1 - x * x + y * y) / s; right[1] = -2 * x * y / s; right[2] = -2 * x / s;
    orthonormalize();
}

// Rotate the sphere so the surface follows a drag of (dx, dy) pixels (y up).
static void dragSphere(double dx, double dy) {
    double len = std::hypot(dx, dy);
    if (len == 0) return;
    // In view space the drag is a rotation A about a = (-dy, dx, 0)/len by len/radius.
    // World points stay fixed, so rot becomes rot * A^-1 = rot * A(a, -angle).
    double angle = -len / app.radius(), c = std::cos(angle), sn = std::sin(angle);
    double a[3] = {-dy / len, dx / len, 0};
    double A[3][3];  // Rodrigues' formula; A[j] = image of view axis j
    for (int j = 0; j < 3; j++) {
        double e[3] = {0, 0, 0}, axe[3];
        e[j] = 1;
        cross3(a, e, axe);
        for (int i = 0; i < 3; i++) A[j][i] = e[i] * c + axe[i] * sn + a[i] * a[j] * (1 - c);
    }
    double R[3][3];
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) R[j][i] = app.rot[0][i] * A[j][0] + app.rot[1][i] * A[j][1] + app.rot[2][i] * A[j][2];
    std::memcpy(app.rot, R, sizeof R);
    orthonormalize();
}

static void resetView() {
    app.cx = app.group->cx;
    app.cy = app.group->cy;
    app.span = app.group->span;
    app.sphereZoom = 1.0;
    faceSphereToward(cd(app.cx, app.cy));
}

// S key: switch views, keeping the point at the centre of the window.
static void toggleSphere() {
    if (app.sphere) {
        cd z = App::sphereToZ(app.rot[2]);
        if (std::abs(z) < 1e12) { app.cx = z.real(); app.cy = z.imag(); }
    } else {
        faceSphereToward(cd(app.cx, app.cy));
    }
    app.sphere = !app.sphere;
    std::printf("view: %s\n", app.sphere ? "Riemann sphere (drag rotates, wheel zooms)" : "complex plane");
}

static void buildGroup() {
    std::vector<Mat> mats = app.group->build(app.N);
    app.count = int(mats.size());
    std::vector<float> f(mats.size() * 8);
    std::vector<double> d(mats.size() * 8);
    for (size_t i = 0; i < mats.size(); i++) {
        const Mat& m = mats[i];
        double v[8] = {m.a.real(), m.a.imag(), m.b.real(), m.b.imag(), m.c.real(), m.c.imag(), m.d.real(), m.d.imag()};
        for (int k = 0; k < 8; k++) { d[8 * i + k] = v[k]; f[8 * i + k] = float(v[k]); }
    }
    renderer.matrices(f, d);
    restart();
}

static void resize(int w, int h) {
    if (w <= 0 || h <= 0 || (w == app.W && h == app.H)) return;
    app.W = w;
    app.H = h;
    renderer.resize(w, h);
    restart();
}

static void updateTitle() {
    if (msSince(app.titleT) < 120) return;
    app.titleT = Clock::now();
    // read the accumulated sums under the cursor straight from the GPU buffer
    double a[4] = {0, 0, 0, 0};
    int px = std::clamp(app.mx, 0, app.W - 1), py = std::clamp(app.my, 0, app.H - 1);
    renderer.readPixel(px, py, a);
    cd z = app.toZ(app.mx, app.my);
    char where[160], extent[64];
    if (std::isnan(z.real())) {
        std::snprintf(where, sizeof where, "off the sphere");
    } else {
        cd val = app.group->value(cd(a[0], a[1]), cd(a[2], a[3]), app.func, z);
        std::snprintf(where, sizeof where, "z = %.6f%+.6fi   %s(z) = %.6g%+.6gi", z.real(), z.imag(),
                      app.group->funcNames[app.func], val.real(), val.imag());
    }
    if (app.sphere) std::snprintf(extent, sizeof extent, "sphere zoom %.3g", app.sphereZoom);
    else std::snprintf(extent, sizeof extent, "width %.3g", app.scale() * app.W);
    char buf[512];
    std::snprintf(buf, sizeof buf, "%s  |  %s  |  %s  |  terms %d/%d %s |  %s  |  %s",
                  app.group->title, app.gpu.c_str(), app.fp64 ? "fp64" : "fp32", app.done, app.count,
                  app.done < app.count ? "" : (std::string("(") + std::to_string(int(app.renderMs)) + " ms) ").c_str(),
                  where, extent);
    wchar_t wbuf[512];
    MultiByteToWideChar(CP_UTF8, 0, buf, -1, wbuf, 512);
    SetWindowTextW(app.hwnd, wbuf);
}

static void verifyGpuValues() {
    const int points[][2] = {{app.W/2,app.H/2}, {app.W/3,app.H/4}, {3*app.W/4,2*app.H/3}};
    double worst = 0;
    for (const auto& p : points) {
        double sums[4];
        renderer.readPixel(p[0],p[1],sums);
        cd actual[] = {{sums[0],sums[1]},{sums[2],sums[3]}};
        cd z = app.toZ(p[0],p[1]), reference[2];
        if (std::isnan(z.real())) { std::printf("GPU check (%d,%d): off the sphere, skipped\n",p[0],p[1]); continue; }
        app.group->reference(z, reference);
        for (int i=0;i<2;++i) {
            double error=std::abs(actual[i]-reference[i])/std::max(1e-30,std::abs(reference[i]));
            std::printf("GPU check (%d,%d) %s: relative error %.3e\n",p[0],p[1],app.group->sumNames[i],error);            if (!std::isfinite(error)) fatal("Non-finite GPU verification result");
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
    std::snprintf(name, sizeof name, "%s_%04d%02d%02d_%02d%02d%02d.png", app.group->pngPrefix,
                  t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond);
    const char* output = path ? path : name;
    if (!writePNG(output, app.W, app.H, rgb)) fatal("Could not write PNG snapshot");
    std::printf("saved %s (%dx%d, %d/%d terms)\n", output, app.W, app.H, app.done, app.count);
}

static ViewConstants viewConstants(int start, int end) {
    ViewConstants v = {UINT(start), UINT(end), UINT(app.W), UINT(app.H), app.cx, app.cy, app.scale(),
                       UINT(app.mode), UINT(app.func), app.dens, UINT(app.sphere), app.radius()};
    std::memcpy(v.rot, app.rot, sizeof v.rot);
    return v;
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
        app.group->invariance();
        savePNG(snapshotPath.c_str());
        PostQuitMessage(0);
        snapshotPath.clear();
    }
    renderer.present(app.done >= app.count);
    updateTitle();
}

static void zoomAt(int x, int y, double k) {
    if (app.sphere) {  // the sphere stays centred in the window
        app.sphereZoom = std::clamp(app.sphereZoom / k, 0.05, 1e4);
        restart();
        return;
    }
    cd z0 = app.toZ(x, y);
    app.span *= k;
    cd z1 = app.toZ(x, y);
    app.cx += z0.real() - z1.real();
    app.cy += z0.imag() - z1.imag();
    restart();
}

static LRESULT CALLBACK wndProc(HWND h, UINT msg, WPARAM wp, LPARAM lp) {
    const bool variableTerms = app.group->defaultTerms > 0;
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
            if (app.sphere) {
                dragSphere(x - app.mx, app.my - y);
            } else {
                double s = app.scale();
                app.cx -= (x - app.mx) * s;
                app.cy += (y - app.my) * s;
            }
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
        case '+': case '=': if (variableTerms) { app.N = std::min<size_t>(app.N * 2, size_t(1) << 22); buildGroup(); } break;
        case '-': case '_': if (variableTerms) { app.N = std::max<size_t>(app.N / 2, 256); buildGroup(); } break;
        case 'd': case 'D': app.fp64 = !app.fp64; app.chunk = 64; restart(); break;
        case 'm': case 'M': app.mode = (app.mode + 1) % 4; break;
        case 'f': case 'F': app.func = (app.func + 1) % 3; break;
        case '[': app.dens = std::max(0.25f, app.dens / 1.25f); break;
        case ']': app.dens = std::min(8.0f, app.dens * 1.25f); break;
        case 'r': case 'R': resetView(); restart(); break;
        case 's': case 'S': if (app.group->sphere) { toggleSphere(); restart(); } break;
        case 'i': case 'I': app.group->invariance(); break;
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

static int runApp(const Group& group, int argc, char** argv) {
    SetProcessDPIAware();
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    app.group = &group;
    app.N = group.defaultTerms;
    resetView();
    const bool variableTerms = group.defaultTerms > 0;
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
        } else if (option == "--sphere" && group.sphere) {
            app.sphere = true;
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
        } else if (option == "--terms" && variableTerms) {
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
        std::printf("Usage: %s [options]\n", group.exeName);
        std::printf("  --center <real> <imag>  Plot center (default: %g %g)\n", group.cx, group.cy);
        std::printf("  --span <value>          Shorter-side extent in complex units (default: %g)\n", group.span);
        if (variableTerms)
            std::printf("  --terms <count>         Number of group elements (default: %zu)\n", group.defaultTerms);
        if (group.sphere)
            std::printf("  --sphere                Start in the Riemann sphere view, centre facing the viewer\n");
        std::printf(
            "  --fp64                  Use double-precision summation\n"
            "  --snapshot <file.png>   Render, save, and exit\n"
            "  --verify-gpu            Check GPU sums against CPU when saving a snapshot\n"
            "  --help, -h              Show this help\n");
        return 0;
    }
    std::printf("%s  %s\n", group.title, group.formula);
    std::printf("  drag: pan   wheel: zoom   %sD: fp32/fp64   M: plot style\n",
                variableTerms ? "+/-: elements x2 / /2   " : "");
    std::printf("  F: function (%s, %s, %s)   [ ]: contour density   R: reset   I: invariance check   P: save PNG   Esc: quit\n",
                group.funcNames[0], group.funcNames[1], group.funcNames[2]);
    if (group.sphere) std::printf("  S: plane / Riemann sphere (in the sphere view, drag rotates and the wheel zooms)\n");
    std::printf("\n");
    faceSphereToward(cd(app.cx, app.cy));  // --center also picks the face the sphere view starts on

    HINSTANCE hi = GetModuleHandleW(nullptr);
    WNDCLASSW wc = {};
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = wndProc;
    wc.hInstance = hi;
    wc.hCursor = LoadCursor(nullptr, IDC_CROSS);
    wc.lpszClassName = L"AutomorphicPlot";
    RegisterClassW(&wc);
    RECT r = {0, 0, 1400, 1000};
    AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE);
    wchar_t title[256];
    MultiByteToWideChar(CP_UTF8, 0, group.title, -1, title, 256);
    app.hwnd = CreateWindowW(wc.lpszClassName, title, WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
                             r.right - r.left, r.bottom - r.top, nullptr, nullptr, hi, nullptr);
    RECT cr;
    GetClientRect(app.hwnd, &cr);
    app.gpu = renderer.initialize(app.hwnd, cr.right, cr.bottom, shaderSource(group.hlsl));
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
