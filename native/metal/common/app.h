#pragma once
// Shared native application for macOS: Cocoa window, Metal renderer,
// progressive GPU summation, domain colouring, PNG export and command-line
// handling. The Metal counterpart of ../../dx12/common/app.h.
//
// A program describes its group with a `Group` and calls runApp(group, argc, argv).
// Each pixel accumulates sum_i term(w1_i, w2_i, i) over the uploaded matrices into
// two complex sums kept as double-floats; the group then says what value to colour.
//
// Controls
//   drag            pan                 wheel / pinch  zoom about cursor
//   + / -           double / halve the number of group elements (infinite groups)
//   D               toggle fp32 / fp64 (double-float emulated) summation and display
//   M               cycle plot style    F            cycle displayed function
//   C               phase colours: hue wheel / Mathematica's GreenPinkTones
//   [ / ]           contour density     R            reset view
//   S               plane / Riemann sphere view (groups with Group::sphere);
//                   in the sphere view, drag rotates and the wheel zooms
//   I               invariance check (console)
//   P               save PNG of the current window
//   Esc             quit

#import <Cocoa/Cocoa.h>
#import <QuartzCore/CAMetalLayer.h>

#include <sys/ioctl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <limits>
#include <string>
#include <vector>

static void fatal(const char* msg) {
    std::fprintf(stderr, "\nERROR: %s\n", msg);
    std::fflush(nullptr);
    _exit(1);
}

// Console output of actions: render times, saved files, view changes, checks.
// In an interactive terminal only the last Console::LINES lines are shown, in a
// region below the start-up banner that is redrawn in place (ANSI escapes) as
// each new line arrives. Piped output and --snapshot runs print every line.
struct Console {
    static constexpr size_t LINES = 5;
    bool rolling = false;
    std::string partial;            // text after the last newline, until the line is complete
    std::deque<std::string> lines;  // the last LINES complete lines
    size_t shown = 0;               // lines currently drawn in the region

    void write(const char* text) {
        if (!rolling) { std::fputs(text, stdout); return; }
        partial += text;
        bool added = false;
        for (size_t nl; (nl = partial.find('\n')) != std::string::npos;) {
            std::string line = partial.substr(0, nl);
            partial.erase(0, nl + 1);
            if (line.empty()) continue;  // blank separators would waste the region
            lines.push_back(line);
            if (lines.size() > LINES) lines.pop_front();
            added = true;
        }
        if (added) redraw();
    }
    // Shorten to `width` columns (one per UTF-8 code point) so no line wraps,
    // which would break the cursor arithmetic in redraw().
    static std::string fit(const std::string& line, size_t width) {
        size_t columns = 0;
        for (size_t i = 0; i < line.size(); i++) {
            if ((uint8_t(line[i]) & 0xC0) == 0x80) continue;  // UTF-8 continuation byte
            if (++columns > width) return line.substr(0, i);
        }
        return line;
    }
    void redraw() {
        winsize ws = {};
        size_t width = ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 1 ? ws.ws_col - 1 : 79;
        std::string out;
        if (shown > 0) out += "\r\x1b[" + std::to_string(shown) + "A";  // back to the region's first line
        for (const std::string& line : lines) out += "\x1b[2K" + fit(line, width) + "\n";
        out += "\x1b[J";  // clear anything below
        std::fwrite(out.data(), 1, out.size(), stdout);
        std::fflush(stdout);
        shown = lines.size();
    }
};
static Console console;

// printf for actions; see Console.
__attribute__((format(printf, 1, 2))) static void consoleLog(const char* format, ...) {
    va_list args, copy;
    va_start(args, format);
    va_copy(copy, args);
    int n = std::vsnprintf(nullptr, 0, format, copy);
    va_end(copy);
    std::vector<char> text(size_t(std::max(n, 0)) + 1);
    std::vsnprintf(text.data(), text.size(), format, args);
    va_end(args);
    console.write(text.data());
}

#include "mobius.h"  // ../../dx12/common, portable
#include "png.h"     // ../../dx12/common, portable
#include "metal.h"
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
    const char* msl;            // term / plotValue / FADE_POLES, see shaders.h

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
    // 0: any partial sum can be shown (a series that converges term by term).
    // Otherwise the items are ordered so that the first m = minPartial, 2 minPartial,
    // 4 minPartial, ... of them sum to m/count of the total, and only those
    // partial sums are meaningful: the app stops at them and scales by count/m.
    int minPartial = 0;
};

using Clock = std::chrono::steady_clock;
static double msSince(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

struct App {
    const Group* group = nullptr;
    NSWindow* window = nil;
    NSView* view = nil;
    CAMetalLayer* layer = nil;  // the view's backing layer; frame() sizes its drawables
    int W = 0, H = 0;
    double cx = 0.0, cy = 0.0, span = 5.0;  // span = width of the shorter screen side in z-units
    // Riemann sphere view: rot[j] is the world-space image of view axis j
    // (x right, y up, z toward the viewer); zoom 1 fits the sphere in the window.
    bool sphere = false;
    double rot[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double sphereZoom = 1.0;
    size_t N = 0;
    bool fp64 = false;
    int mode = 0, func = 0, palette = 0;
    float dens = 1.0f;

    int count = 0, done = 0, chunk = 256;
    Clock::time_point startT;
    double renderMs = 0;
    bool dirty = true;  // the window needs redrawing

    bool dragging = false;
    int mx = 0, my = 0;  // cursor in pixels, top-left origin
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
static MetalRenderer renderer;
static std::string snapshotPath;  // --snapshot <file.png>: render to completion, save, exit
static bool verifyGPU = false;

static const char* precisionName() { return app.fp64 ? "fp64 (double-float)" : "fp32"; }

static void restart() {
    app.done = 0;
    app.startT = Clock::now();
    app.dirty = true;
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
    consoleLog("view: %s\n", app.sphere ? "Riemann sphere (drag rotates, wheel zooms)" : "complex plane");
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
    // read the accumulated sums under the cursor straight from the shared GPU buffer
    double a[4] = {0, 0, 0, 0};
    int px = std::clamp(app.mx, 0, app.W - 1), py = std::clamp(app.my, 0, app.H - 1);
    renderer.readPixel(px, py, a);
    cd z = app.toZ(px, py);
    char where[160], extent[64];
    if (std::isnan(z.real())) {
        std::snprintf(where, sizeof where, "off the sphere");
    } else {
        double k = app.group->minPartial && app.done ? double(app.count) / app.done : 1.0;
        cd val = app.group->value(k * cd(a[0], a[1]), k * cd(a[2], a[3]), app.func, z);
        std::snprintf(where, sizeof where, "z = %.6f%+.6fi   %s(z) = %.6g%+.6gi", z.real(), z.imag(),
                      app.group->funcNames[app.func], val.real(), val.imag());
    }
    if (app.sphere) std::snprintf(extent, sizeof extent, "sphere zoom %.3g", app.sphereZoom);
    else std::snprintf(extent, sizeof extent, "width %.3g", app.scale() * app.W);
    char buf[512];
    std::snprintf(buf, sizeof buf, "%s  |  %s  |  %s  |  terms %d/%d %s |  %s  |  %s",
                  app.group->title, app.gpu.c_str(), precisionName(), app.done, app.count,
                  app.done < app.count ? "" : (std::string("(") + std::to_string(int(app.renderMs)) + " ms) ").c_str(),
                  where, extent);
    app.window.title = [NSString stringWithUTF8String:buf];
}

static void verifyGpuValues() {
    const int points[][2] = {{app.W / 2, app.H / 2}, {app.W / 3, app.H / 4}, {3 * app.W / 4, 2 * app.H / 3}};
    double worst = 0;
    for (const auto& p : points) {
        double sums[4];
        renderer.readPixel(p[0], p[1], sums);
        cd actual[] = {{sums[0], sums[1]}, {sums[2], sums[3]}};
        cd z = app.toZ(p[0], p[1]), reference[2];
        if (std::isnan(z.real())) { consoleLog("GPU check (%d,%d): off the sphere, skipped\n", p[0], p[1]); continue; }
        app.group->reference(z, reference);
        for (int i = 0; i < 2; ++i) {
            double error = std::abs(actual[i] - reference[i]) / std::max(1e-30, std::abs(reference[i]));
            consoleLog("GPU check (%d,%d) %s: relative error %.3e\n", p[0], p[1], app.group->sumNames[i], error);
            if (!std::isfinite(error)) fatal("Non-finite GPU verification result");
            worst = std::max(worst, error);
        }
    }
    if (worst > (app.fp64 ? 1e-9 : 5e-3)) fatal("GPU sums do not match CPU reference");
    consoleLog("GPU verification passed (%s)\n", precisionName());
}

static void savePNG(const char* path = nullptr) {
    const auto& rgb = renderer.pixels();
    if (rgb.size() != size_t(app.W) * app.H * 3) fatal("No captured frame to save");
    std::time_t now = std::time(nullptr);
    std::tm t = *std::localtime(&now);
    char name[128];
    std::snprintf(name, sizeof name, "%s_%04d%02d%02d_%02d%02d%02d.png", app.group->pngPrefix,
                  t.tm_year + 1900, t.tm_mon + 1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec);
    const char* output = path ? path : name;
    if (!writePNG(output, app.W, app.H, rgb)) fatal("Could not write PNG snapshot");
    consoleLog("saved %s (%dx%d, %d/%d terms)\n", output, app.W, app.H, app.done, app.count);
}

static ViewConstants viewConstants(int start, int end) {
    ViewConstants v = {};
    v.start = uint32_t(start); v.end = uint32_t(end);
    v.width = uint32_t(app.W); v.height = uint32_t(app.H);
    splitDouble(app.cx, v.cx);
    splitDouble(app.cy, v.cy);
    splitDouble(app.scale(), v.scale);
    v.mode = uint32_t(app.mode); v.function = uint32_t(app.func);
    v.density = app.dens;
    v.sphere = app.sphere;
    splitDouble(app.radius(), v.radius);
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) splitDouble(app.rot[j][i], v.rot[3 * j + i]);
    v.palette = uint32_t(app.palette);
    v.count = uint32_t(app.count);
    return v;
}

// One frame: sum more terms for up to ~30 ms, then redraw if anything changed.
static void frame() {
    // Size the drawables to the view in pixels. CAMetalLayer does not follow its
    // view's size on its own: a stale drawableSize gives drawables of the old size,
    // which frame() cannot draw into, and the window would freeze after a resize.
    CGSize size = [app.view convertRectToBacking:app.view.bounds].size;
    if (size.width >= 1 && size.height >= 1 && !CGSizeEqualToSize(size, app.layer.drawableSize)) {
        app.layer.contentsScale = app.window.backingScaleFactor;
        app.layer.drawableSize = size;
    }
    resize(int(size.width), int(size.height));
    if (app.W <= 0 || app.H <= 0 || app.count == 0) return;
    if (app.done < app.count) {
        auto t0 = Clock::now();
        while (app.done < app.count && msSince(t0) < 30.0) {
            int64_t step = app.chunk;
            if (int m = app.group->minPartial) step = app.done ? app.done : m;  // stop at m, 2m, 4m, ...
            int end = int(std::min<int64_t>(app.count, int64_t(app.done) + step));
            auto td = Clock::now();
            renderer.compute(viewConstants(app.done, end), app.fp64);
            double dt = msSince(td);
            app.done = end;
            if (dt < 4.0) app.chunk = std::min(app.chunk * 2, 1 << 16);
            else if (dt > 12.0) app.chunk = std::max(app.chunk / 2, 16);
        }
        app.dirty = true;
        if (app.done >= app.count) {
            app.renderMs = msSince(app.startT);
            consoleLog("render complete: %d terms, %dx%d, %s, %.0f ms\n", app.count, app.W, app.H,
                        precisionName(), app.renderMs);
        }
    }
    bool snapshotNow = !snapshotPath.empty() && app.done >= app.count;
    bool capture = app.saveRequested || snapshotNow;
    if (app.dirty || capture) {
        id<CAMetalDrawable> drawable = [app.layer nextDrawable];
        if (drawable && drawable.texture.width == size_t(app.W) && drawable.texture.height == size_t(app.H)) {
            renderer.draw(viewConstants(0, app.done), app.fp64, drawable, capture);
            app.dirty = false;
            if (app.saveRequested) { savePNG(); app.saveRequested = false; }
            if (snapshotNow) {
                if (verifyGPU) verifyGpuValues();
                app.group->invariance();
                savePNG(snapshotPath.c_str());
                snapshotPath.clear();
                [NSApp terminate:nil];
            }
        }
    }
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

static void keyPressed(unichar key) {
    const bool variableTerms = app.group->defaultTerms > 0;
    switch (key) {
    case '+': case '=': if (variableTerms) { app.N = std::min<size_t>(app.N * 2, size_t(1) << 22); buildGroup(); } break;
    case '-': case '_': if (variableTerms) { app.N = std::max<size_t>(app.N / 2, 256); buildGroup(); } break;
    case 'd': case 'D': app.fp64 = !app.fp64; app.chunk = 64; restart(); break;
    case 'm': case 'M': app.mode = (app.mode + 1) % 4; break;
    case 'f': case 'F': app.func = (app.func + 1) % 3; break;
    case 'c': case 'C': app.palette = (app.palette + 1) % 2; break;
    case '[': app.dens = std::max(0.25f, app.dens / 1.25f); break;
    case ']': app.dens = std::min(8.0f, app.dens * 1.25f); break;
    case 'r': case 'R': resetView(); restart(); break;
    case 's': case 'S': if (app.group->sphere) { toggleSphere(); restart(); } break;
    case 'i': case 'I': app.group->invariance(); break;
    case 'p': case 'P': app.saveRequested = true; break;
    case 27: [NSApp terminate:nil]; break;
    }
    app.dirty = true;
}

// The Metal view: forwards input to the app and draws a frame per display refresh.
@interface PlotView : NSView
@end

@implementation PlotView
- (BOOL)wantsUpdateLayer { return YES; }
- (CALayer*)makeBackingLayer { return [CAMetalLayer layer]; }
- (BOOL)acceptsFirstResponder { return YES; }
- (BOOL)acceptsFirstMouse:(NSEvent*)event { return YES; }
// Cursor position in drawable pixels, top-left origin.
- (void)pixelOf:(NSEvent*)event x:(int*)x y:(int*)y {
    NSPoint p = [self convertPoint:event.locationInWindow fromView:nil];
    NSSize b = self.bounds.size, d = app.layer.drawableSize;
    *x = int(std::floor(p.x * d.width / std::max(1.0, b.width)));
    *y = int(std::floor((b.height - p.y) * d.height / std::max(1.0, b.height)));
}
- (void)mouseDown:(NSEvent*)event {
    app.dragging = true;
    [self pixelOf:event x:&app.mx y:&app.my];
}
- (void)mouseUp:(NSEvent*)event { app.dragging = false; }
- (void)mouseMoved:(NSEvent*)event { [self pixelOf:event x:&app.mx y:&app.my]; }
- (void)mouseDragged:(NSEvent*)event {
    int x, y;
    [self pixelOf:event x:&x y:&y];
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
}
- (void)scrollWheel:(NSEvent*)event {
    // A mouse wheel notch is about 1 line; trackpads report precise deltas in points.
    double notches = event.hasPreciseScrollingDeltas ? event.scrollingDeltaY / 10.0 : event.scrollingDeltaY;
    if (notches == 0) return;
    int x, y;
    [self pixelOf:event x:&x y:&y];
    zoomAt(x, y, std::exp(-notches * 0.18));
}
- (void)magnifyWithEvent:(NSEvent*)event {
    int x, y;
    [self pixelOf:event x:&x y:&y];
    zoomAt(x, y, 1.0 / std::max(0.1, 1.0 + event.magnification));
}
- (void)keyDown:(NSEvent*)event {
    NSString* chars = event.characters;
    if (chars.length > 0) keyPressed([chars characterAtIndex:0]);
}
@end

@interface AppDelegate : NSObject <NSApplicationDelegate, NSWindowDelegate>
@end
@implementation AppDelegate
- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication*)sender { return YES; }
@end

static int runApp(const Group& group, int argc, char** argv) {
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
        } else if (option == "--colors") {
            const std::string name = i + 1 < argc ? argv[i + 1] : "";
            if (name == "hue") app.palette = 0;
            else if (name == "greenpink") app.palette = 1;
            else {
                std::fprintf(stderr, "--colors requires hue or greenpink\n");
                return 1;
            }
            ++i;
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
        } else if (option == "--function") {
            double index = -1;
            if (i + 1 >= argc || !readNumber(argv[i + 1], index) || (index != 0 && index != 1 && index != 2)) {
                std::fprintf(stderr, "--function requires 0, 1 or 2\n");
                return 1;
            }
            app.func = int(index);
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
        std::printf("  --function <0|1|2>      Displayed function: %s, %s or %s (default: 0)\n",
                    group.funcNames[0], group.funcNames[1], group.funcNames[2]);
        std::printf(
            "  --colors <hue|greenpink> Phase colours: hue wheel or GreenPinkTones (default: hue)\n"
            "  --fp64                  Use emulated double precision (double-float) on the GPU\n"
            "  --snapshot <file.png>   Render, save, and exit\n"
            "  --verify-gpu            Check GPU sums against CPU when saving a snapshot\n"
            "  --help, -h              Show this help\n");
        return 0;
    }
    std::printf("%s  %s\n", group.title, group.formula);
    std::printf("  drag: pan   wheel/pinch: zoom   %sD: fp32/fp64   M: plot style\n",
                variableTerms ? "+/-: elements x2 / /2   " : "");
    std::printf("  F: function (%s, %s, %s)   C: colours (hue wheel / GreenPinkTones)\n  [ ]: contour density   R: reset   I: invariance check   P: save PNG   Esc: quit\n",
                group.funcNames[0], group.funcNames[1], group.funcNames[2]);
    if (group.sphere) std::printf("  S: plane / Riemann sphere (in the sphere view, drag rotates and the wheel zooms)\n");
    std::printf("\n");
    faceSphereToward(cd(app.cx, app.cy));  // --center also picks the face the sphere view starts on

    @autoreleasepool {
        [NSApplication sharedApplication];
        // Without an app bundle, a terminal-launched process needs this to get a
        // Dock icon, a menu bar and keyboard focus.
        [NSApp setActivationPolicy:NSApplicationActivationPolicyRegular];
        AppDelegate* delegate = [AppDelegate new];
        NSApp.delegate = delegate;

        NSMenu* menu = [NSMenu new];
        NSMenuItem* appItem = [NSMenuItem new];
        [menu addItem:appItem];
        NSMenu* appMenu = [NSMenu new];
        [appMenu addItemWithTitle:@"Quit" action:@selector(terminate:) keyEquivalent:@"q"];
        appItem.submenu = appMenu;
        NSApp.mainMenu = menu;

        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) fatal("Metal is not available on this Mac");

        // A 1400 x 1000 pixel drawable, as on Windows, so snapshots compare directly.
        NSScreen* screen = NSScreen.mainScreen;
        double backing = screen ? screen.backingScaleFactor : 1.0;
        NSRect visible = screen ? screen.visibleFrame : NSMakeRect(0, 0, 1400, 1000);
        double w = std::min(1400.0 / backing, visible.size.width - 40), h = std::min(1000.0 / backing, visible.size.height - 60);
        NSRect frameRect = NSMakeRect(0, 0, w, h);
        app.window = [[NSWindow alloc] initWithContentRect:frameRect
                                                 styleMask:NSWindowStyleMaskTitled | NSWindowStyleMaskClosable |
                                                           NSWindowStyleMaskMiniaturizable | NSWindowStyleMaskResizable
                                                   backing:NSBackingStoreBuffered
                                                     defer:NO];
        app.window.delegate = delegate;
        app.window.acceptsMouseMovedEvents = YES;
        app.window.title = [NSString stringWithUTF8String:group.title];
        [app.window center];

        PlotView* view = [[PlotView alloc] initWithFrame:frameRect];
        view.wantsLayer = YES;
        app.view = view;
        app.layer = (CAMetalLayer*)view.layer;
        app.layer.device = device;
        app.layer.pixelFormat = MTLPixelFormatBGRA8Unorm;
        app.layer.framebufferOnly = NO;  // P and --snapshot read the drawable back
        CGColorSpaceRef srgb = CGColorSpaceCreateWithName(kCGColorSpaceSRGB);
        app.layer.colorspace = srgb;
        CGColorSpaceRelease(srgb);
        app.window.contentView = view;
        [app.window makeFirstResponder:view];  // keys and mouse-moved events go to the first responder

        std::string sources[2] = {shaderSource(group.msl, false), shaderSource(group.msl, true)};
        app.gpu = renderer.initialize(device, app.layer.pixelFormat, sources);
        std::printf("Metal\nGPU    %s\n\n", app.gpu.c_str());
        console.rolling = isatty(STDOUT_FILENO) && snapshotPath.empty();
        buildGroup();  // frame() sizes the drawables and the accumulator

        [app.window makeKeyAndOrderFront:nil];
        [NSApp activateIgnoringOtherApps:YES];
        // Drive frames from a timer in the common run-loop modes, so the plot keeps
        // rendering while the window is being resized or a menu is open.
        NSTimer* timer = [NSTimer timerWithTimeInterval:1.0 / 60.0 repeats:YES block:^(NSTimer*) {
            // Release each frame's drawable promptly; holding them stalls nextDrawable.
            @autoreleasepool { frame(); }
        }];
        [NSRunLoop.mainRunLoop addTimer:timer forMode:NSRunLoopCommonModes];
        [NSApp run];
    }
    return 0;
}
