// Native OpenGL 4.3 renderer for the automorphic function
//     f(z) = P_{-1,-3}(z) / P_{0,-4}(z),   P_{k,l}(z) = sum_{g in G} (az+b)^k (cz+d)^l
// for the Kleinian group G = <g1, g2>.
//
// The Poincaré series are summed per pixel in a compute shader (fp32 or fp64),
// accumulated progressively in a double-precision buffer, and displayed with
// domain colouring. No external dependencies: Win32 + WGL + opengl32.lib.
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
#include <GL/gl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <queue>
#include <string>
#include <vector>

// Ask hybrid-graphics laptops to run this process on the discrete GPU.
extern "C" {
__declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}

// ---------------------------------------------------------------------------
// Minimal OpenGL 4.3 loader
// ---------------------------------------------------------------------------
typedef char GLchar;
typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

#define GL_FRAGMENT_SHADER 0x8B30
#define GL_VERTEX_SHADER 0x8B31
#define GL_COMPUTE_SHADER 0x91B9
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82
#define GL_INFO_LOG_LENGTH 0x8B84
#define GL_SHADER_STORAGE_BUFFER 0x90D2
#define GL_SHADER_STORAGE_BARRIER_BIT 0x00002000
#define GL_STATIC_DRAW 0x88E4
#define GL_DYNAMIC_COPY 0x88EA
#define GL_R32UI 0x8236
#define GL_RED_INTEGER 0x8D94

#define WGL_CONTEXT_MAJOR_VERSION_ARB 0x2091
#define WGL_CONTEXT_MINOR_VERSION_ARB 0x2092
#define WGL_CONTEXT_PROFILE_MASK_ARB 0x9126
#define WGL_CONTEXT_CORE_PROFILE_BIT_ARB 0x00000001

#define GL_FUNCS(X)                                                                         \
    X(GLuint, glCreateShader, (GLenum))                                                     \
    X(void, glShaderSource, (GLuint, GLsizei, const GLchar* const*, const GLint*))          \
    X(void, glCompileShader, (GLuint))                                                      \
    X(void, glGetShaderiv, (GLuint, GLenum, GLint*))                                        \
    X(void, glGetShaderInfoLog, (GLuint, GLsizei, GLsizei*, GLchar*))                       \
    X(GLuint, glCreateProgram, (void))                                                      \
    X(void, glAttachShader, (GLuint, GLuint))                                               \
    X(void, glLinkProgram, (GLuint))                                                        \
    X(void, glGetProgramiv, (GLuint, GLenum, GLint*))                                       \
    X(void, glGetProgramInfoLog, (GLuint, GLsizei, GLsizei*, GLchar*))                      \
    X(void, glDeleteShader, (GLuint))                                                       \
    X(void, glDeleteProgram, (GLuint))                                                      \
    X(void, glUseProgram, (GLuint))                                                         \
    X(GLint, glGetUniformLocation, (GLuint, const GLchar*))                                 \
    X(void, glUniform1i, (GLint, GLint))                                                    \
    X(void, glUniform1f, (GLint, GLfloat))                                                  \
    X(void, glUniform1d, (GLint, GLdouble))                                                 \
    X(void, glUniform2d, (GLint, GLdouble, GLdouble))                                       \
    X(void, glGenBuffers, (GLsizei, GLuint*))                                               \
    X(void, glDeleteBuffers, (GLsizei, const GLuint*))                                      \
    X(void, glBindBuffer, (GLenum, GLuint))                                                 \
    X(void, glBufferData, (GLenum, GLsizeiptr, const void*, GLenum))                        \
    X(void, glBindBufferBase, (GLenum, GLuint, GLuint))                                     \
    X(void, glClearBufferData, (GLenum, GLenum, GLenum, GLenum, const void*))               \
    X(void, glGetBufferSubData, (GLenum, GLintptr, GLsizeiptr, void*))                      \
    X(void, glDispatchCompute, (GLuint, GLuint, GLuint))                                    \
    X(void, glMemoryBarrier, (GLbitfield))                                                  \
    X(void, glGenVertexArrays, (GLsizei, GLuint*))                                          \
    X(void, glBindVertexArray, (GLuint))

#define DECLARE_GL(ret, name, args) \
    typedef ret(APIENTRY* PFN_##name) args; \
    static PFN_##name name = nullptr;
GL_FUNCS(DECLARE_GL)

typedef HGLRC(WINAPI* PFN_wglCreateContextAttribsARB)(HDC, HGLRC, const int*);
typedef BOOL(WINAPI* PFN_wglSwapIntervalEXT)(int);

static void fatal(const char* msg) {
    std::fprintf(stderr, "\nERROR: %s\n", msg);
    MessageBoxA(nullptr, msg, "Poincaré series", MB_ICONERROR);
    std::exit(1);
}

static void loadGL() {
#define LOAD_GL(ret, name, args)                                      \
    name = (PFN_##name)wglGetProcAddress(#name);                      \
    if (!name) fatal("OpenGL 4.3 function missing: " #name);
    GL_FUNCS(LOAD_GL)
}

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
// Shaders
// ---------------------------------------------------------------------------
static const char* COMPUTE_SRC = R"GLSL(
layout(local_size_x = 16, local_size_y = 16) in;
#define TILE 256
layout(std430, binding = 0) readonly buffer Mats { VEC4 mats[]; };  // (a, b), (c, d)
layout(std430, binding = 1) buffer Acc { dvec4 acc[]; };            // (P_{-1,-3}, P_{0,-4})
uniform int uStart, uEnd, uW, uH;
uniform dvec2 uCenter;
uniform double uScale;
shared VEC4 tileM[2 * TILE];

VEC2 cm(VEC2 a, VEC2 b) { return VEC2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
#ifdef FP64
// fp64 range is ~1e308, so no scaling is needed; one division per term
// (double division is emulated on NVIDIA GPUs and very expensive).
VEC2 crecip(VEC2 q) { return VEC2(q.x, -q.y) * (1.0LF / dot(q, q)); }
#else
VEC2 crecip(VEC2 q) {  // scaled 1/q, avoids overflow of |q|^2 in fp32
    REAL rs = 1.0 / max(abs(q.x), abs(q.y));
    VEC2 t = q * rs;
    return VEC2(t.x, -t.y) * (rs / dot(t, t));
}
#endif

void main() {
    ivec2 p = ivec2(gl_GlobalInvocationID.xy);
    bool inside = p.x < uW && p.y < uH;
    VEC2 z = VEC2(uCenter + (dvec2(p) + 0.5 - 0.5 * dvec2(uW, uH)) * uScale);
    VEC4 sum = VEC4(0);
    int lid = int(gl_LocalInvocationIndex);
    for (int base = uStart; base < uEnd; base += TILE) {
        int n = min(TILE, uEnd - base);
        if (lid < n) {
            tileM[2 * lid] = mats[2 * (base + lid)];
            tileM[2 * lid + 1] = mats[2 * (base + lid) + 1];
        }
        barrier();
        if (inside) {
            for (int i = 0; i < n; i++) {
                VEC4 ab = tileM[2 * i], cd = tileM[2 * i + 1];
                VEC2 w1 = cm(ab.xy, z) + ab.zw;  // az + b
                VEC2 w2 = cm(cd.xy, z) + cd.zw;  // cz + d
                VEC2 w2s = cm(w2, w2);
                VEC2 r = crecip(cm(w1, cm(w2s, w2s)));  // 1 / ((az+b)(cz+d)^4)
                sum += VEC4(cm(w2, r), cm(w1, r));
            }
        }
        barrier();
    }
    if (inside) acc[p.y * uW + p.x] += dvec4(sum);
}
)GLSL";

static const char* VERT_SRC = R"GLSL(
#version 430
void main() {
    vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
)GLSL";

static const char* FRAG_SRC = R"GLSL(
#version 430
layout(std430, binding = 1) readonly buffer Acc { dvec4 acc[]; };
uniform int uW, uH, uMode, uFunc;
uniform float uDens;
out vec4 o;
const float PI = 3.14159265358979;

dvec2 cdiv(dvec2 a, dvec2 b) {
    double s = max(abs(b.x), abs(b.y));
    a /= s; b /= s;
    return dvec2(a.x * b.x + a.y * b.y, a.y * b.x - a.x * b.y) / dot(b, b);
}
vec3 hsv(float h, float s, float v) {
    vec3 k = clamp(abs(mod(h * 6.0 + vec3(0, 4, 2), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    return v * mix(vec3(1.0), k * k * (3.0 - 2.0 * k), s);
}
float saw(float x) { float w = fwidth(x); return mix(fract(x), 0.5, smoothstep(0.3, 0.8, w)); }
float gridline(float x) {
    float w = max(fwidth(x), 1e-6);
    return 1.0 - smoothstep(0.0, 1.2, abs(fract(x - 0.5) - 0.5) / w);
}

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    dvec4 a = acc[p.y * uW + p.x];
    dvec2 fd = uFunc == 0 ? cdiv(a.xy, a.zw) : (uFunc == 1 ? a.zw : a.xy);
    double s = max(abs(fd.x), abs(fd.y));
    if (!(s > 0.0) || isinf(s) || isnan(s)) { o = vec4(0, 0, 0, 1); return; }
    vec2 dir = vec2(fd / s);                       // safe to convert to float
    float lm = log2(float(s)) + 0.5 * log2(dot(dir, dir));
    if (isinf(float(s))) lm = 128.0;
    float h = fract(atan(dir.y, dir.x) / (2.0 * PI));
    vec3 base = hsv(h, 0.9, 1.0);
    vec3 col;
    if (uMode == 0) {
        col = base * (0.55 + 0.45 * sqrt(saw(lm * uDens) * saw(h * 12.0 * uDens)));
    } else if (uMode == 1) {
        col = base;
    } else if (uMode == 2) {
        col = base * (0.45 + 0.55 * saw(lm * uDens));
    } else {
        // f / 2^floor(log2|f|): grid spacing adapts to the magnitude of f
        vec2 fs = dir * exp2(log2(float(s)) - floor(lm)) * 4.0 * uDens;
        float g = max(gridline(fs.x), gridline(fs.y));
        col = mix(base * 0.95, vec3(0.05), 0.75 * g);
    }
    col *= 1.0 - smoothstep(14.0, 24.0, abs(lm));   // fade at poles / zeros
    o = vec4(pow(col, vec3(1.0 / 1.1)), 1.0);
}
)GLSL";

static GLuint compile(GLenum type, const std::string& src) {
    GLuint s = glCreateShader(type);
    const char* p = src.c_str();
    glShaderSource(s, 1, &p, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[8192];
        glGetShaderInfoLog(s, sizeof log, nullptr, log);
        std::fprintf(stderr, "%s\n", log);
        fatal("shader compilation failed (see console)");
    }
    return s;
}
static GLuint link(std::initializer_list<GLuint> shaders) {
    GLuint p = glCreateProgram();
    for (GLuint s : shaders) glAttachShader(p, s);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[8192];
        glGetProgramInfoLog(p, sizeof log, nullptr, log);
        std::fprintf(stderr, "%s\n", log);
        fatal("program link failed (see console)");
    }
    for (GLuint s : shaders) glDeleteShader(s);
    return p;
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
        z.push_back(n & 0xFF); z.push_back(n >> 8);
        z.push_back(~n & 0xFF); z.push_back((~n >> 8) & 0xFF);
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
    HDC hdc = nullptr;
    int W = 0, H = 0;
    double cx = 0.0, cy = -0.35, span = 5.0;  // span = width of the shorter screen side in z-units
    size_t N = size_t(1) << 17;
    bool fp64 = false;
    int mode = 0, func = 0;
    float dens = 1.0f;

    std::vector<Mat> mats;
    GLuint matF = 0, matD = 0, acc = 0, vao = 0;
    GLuint progF = 0, progD = 0, progDisp = 0;
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
    if (!app.acc) return;
    GLuint zero = 0;
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, app.acc);
    glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_R32UI, GL_RED_INTEGER, GL_UNSIGNED_INT, &zero);
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
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, app.matF);
    glBufferData(GL_SHADER_STORAGE_BUFFER, f.size() * sizeof(float), f.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, app.matD);
    glBufferData(GL_SHADER_STORAGE_BUFFER, d.size() * sizeof(double), d.data(), GL_STATIC_DRAW);
    std::printf("group: %d elements, max ||g|| = %.1f  (enumerated in %.0f ms)\n", app.count,
                std::sqrt(norm2(app.mats.back())), msSince(t0));
    restart();
}

static void resize(int w, int h) {
    if (w <= 0 || h <= 0 || (w == app.W && h == app.H)) return;
    app.W = w;
    app.H = h;
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, app.acc);
    glBufferData(GL_SHADER_STORAGE_BUFFER, GLsizeiptr(w) * h * 4 * sizeof(double), nullptr, GL_DYNAMIC_COPY);
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
    int px = std::clamp(app.mx, 0, app.W - 1), py = std::clamp(app.H - 1 - app.my, 0, app.H - 1);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, app.acc);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, GLintptr(size_t(py) * app.W + px) * 4 * sizeof(double), sizeof a, a);
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

static PFN_wglSwapIntervalEXT swapInterval = nullptr;
static int currentInterval = -1;
static std::string snapshotPath;  // --snapshot <file.png>: render to completion, save, exit

static void savePNG(const char* path = nullptr) {
    std::vector<uint8_t> px(size_t(app.W) * app.H * 4), rgb(size_t(app.W) * app.H * 3);
    glReadPixels(0, 0, app.W, app.H, GL_RGBA, GL_UNSIGNED_BYTE, px.data());
    for (int y = 0; y < app.H; y++)
        for (int x = 0; x < app.W; x++)
            for (int k = 0; k < 3; k++)
                rgb[(size_t(app.H - 1 - y) * app.W + x) * 3 + k] = px[(size_t(y) * app.W + x) * 4 + k];
    SYSTEMTIME t;
    GetLocalTime(&t);
    char name[128];
    std::snprintf(name, sizeof name, "poincare_%04d%02d%02d_%02d%02d%02d.png", t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond);
    if (path) std::snprintf(name, sizeof name, "%s", path);
    if (writePNG(name, app.W, app.H, rgb))
        std::printf("saved %s (%dx%d, %d/%d terms)\n", name, app.W, app.H, app.done, app.count);
}

static void frame() {
    if (app.done < app.count) {
        GLuint prog = app.fp64 ? app.progD : app.progF;
        glUseProgram(prog);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, app.fp64 ? app.matD : app.matF);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, app.acc);
        glUniform1i(glGetUniformLocation(prog, "uW"), app.W);
        glUniform1i(glGetUniformLocation(prog, "uH"), app.H);
        glUniform2d(glGetUniformLocation(prog, "uCenter"), app.cx, app.cy);
        glUniform1d(glGetUniformLocation(prog, "uScale"), app.scale());
        GLint locS = glGetUniformLocation(prog, "uStart"), locE = glGetUniformLocation(prog, "uEnd");
        auto t0 = Clock::now();
        // Keep each dispatch short (Windows resets the GPU if one runs > ~2 s)
        // and the whole frame under ~30 ms so the window stays responsive.
        while (app.done < app.count && msSince(t0) < 30.0) {
            int end = std::min(app.count, app.done + app.chunk);
            glUniform1i(locS, app.done);
            glUniform1i(locE, end);
            auto td = Clock::now();
            glDispatchCompute((app.W + 15) / 16, (app.H + 15) / 16, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            glFinish();
            double dt = msSince(td);
            app.done = end;
            if (dt < 4.0) app.chunk = std::min(app.chunk * 2, 1 << 16);
            else if (dt > 12.0) app.chunk = std::max(app.chunk / 2, 16);
        }
        if (app.done >= app.count) {
            app.renderMs = msSince(app.startT);
            std::printf("render complete: %d terms, %dx%d, %s, %.0f ms\n", app.count, app.W, app.H, app.fp64 ? "fp64" : "fp32", app.renderMs);
        }
    }
    glViewport(0, 0, app.W, app.H);
    glUseProgram(app.progDisp);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, app.acc);
    glUniform1i(glGetUniformLocation(app.progDisp, "uW"), app.W);
    glUniform1i(glGetUniformLocation(app.progDisp, "uH"), app.H);
    glUniform1i(glGetUniformLocation(app.progDisp, "uMode"), app.mode);
    glUniform1i(glGetUniformLocation(app.progDisp, "uFunc"), app.func);
    glUniform1f(glGetUniformLocation(app.progDisp, "uDens"), app.dens);
    glBindVertexArray(app.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    if (app.saveRequested) { savePNG(); app.saveRequested = false; }
    if (!snapshotPath.empty() && app.done >= app.count) {
        invarianceCheck();
        savePNG(snapshotPath.c_str());
        PostQuitMessage(0);
        snapshotPath.clear();
    }
    // No vsync while summing (don't wait 16 ms between compute batches); vsync once idle.
    int want = app.done < app.count ? 0 : 1;
    if (swapInterval && want != currentInterval) { swapInterval(want); currentInterval = want; }
    SwapBuffers(app.hdc);
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
        if (app.acc) resize(LOWORD(lp), HIWORD(lp));
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
    for (int i = 1; i + 1 < argc; i++) {
        if (std::string(argv[i]) == "--snapshot") snapshotPath = argv[i + 1];
        if (std::string(argv[i]) == "--terms") app.N = std::strtoull(argv[i + 1], nullptr, 10);
    }
    if (argc > 1 && std::string(argv[argc - 1]) == "--fp64") app.fp64 = true;
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
    app.hdc = GetDC(app.hwnd);

    PIXELFORMATDESCRIPTOR pfd = {sizeof pfd, 1, PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER, PFD_TYPE_RGBA, 32};
    pfd.cDepthBits = 0;
    SetPixelFormat(app.hdc, ChoosePixelFormat(app.hdc, &pfd), &pfd);
    HGLRC legacy = wglCreateContext(app.hdc);
    wglMakeCurrent(app.hdc, legacy);
    auto createCtx = (PFN_wglCreateContextAttribsARB)wglGetProcAddress("wglCreateContextAttribsARB");
    if (!createCtx) fatal("wglCreateContextAttribsARB not available");
    const int attribs[] = {WGL_CONTEXT_MAJOR_VERSION_ARB, 4, WGL_CONTEXT_MINOR_VERSION_ARB, 3,
                           WGL_CONTEXT_PROFILE_MASK_ARB, WGL_CONTEXT_CORE_PROFILE_BIT_ARB, 0};
    HGLRC ctx = createCtx(app.hdc, nullptr, attribs);
    if (!ctx) fatal("could not create an OpenGL 4.3 core context");
    wglMakeCurrent(app.hdc, ctx);
    wglDeleteContext(legacy);
    loadGL();
    swapInterval = (PFN_wglSwapIntervalEXT)wglGetProcAddress("wglSwapIntervalEXT");

    app.gpu = (const char*)glGetString(GL_RENDERER);
    std::printf("OpenGL %s\nGPU    %s\n", (const char*)glGetString(GL_VERSION), app.gpu.c_str());
    if (app.gpu.find("NVIDIA") == std::string::npos && app.gpu.find("Quadro") == std::string::npos)
        std::printf("WARNING: not running on the NVIDIA GPU. Set this exe to 'High performance' in\n"
                    "         Windows Settings > System > Display > Graphics.\n");
    std::printf("\n");

    std::string csF = std::string("#version 430\n#define REAL float\n#define VEC2 vec2\n#define VEC4 vec4\n") + COMPUTE_SRC;
    std::string csD = std::string("#version 430\n#define FP64\n#define REAL double\n#define VEC2 dvec2\n#define VEC4 dvec4\n") + COMPUTE_SRC;
    app.progF = link({compile(GL_COMPUTE_SHADER, csF)});
    app.progD = link({compile(GL_COMPUTE_SHADER, csD)});
    app.progDisp = link({compile(GL_VERTEX_SHADER, VERT_SRC), compile(GL_FRAGMENT_SHADER, FRAG_SRC)});
    glGenVertexArrays(1, &app.vao);
    glGenBuffers(1, &app.matF);
    glGenBuffers(1, &app.matD);
    glGenBuffers(1, &app.acc);

    RECT cr;
    GetClientRect(app.hwnd, &cr);
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
        frame();
    }
}
