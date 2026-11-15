#pragma once
// Shared Metal Shading Language frame, compiled at startup by the Metal framework.
//
// shaderSource(groupMsl, fp64) returns SHADER_HEAD + groupMsl + SHADER_TAIL,
// preceded by `#define FP64 0` or `1`. Each source is compiled into its own
// library holding the compute kernel and the display shaders, so the D key
// switches the precision of both the sums and the displayed value.
//
// Apple GPUs have no double precision. With FP64 set, the real type R is the
// double-float `df` below: an unevaluated sum hi + lo of two floats, about
// 48 significant bits, with error-free transformations (two-sum, fma). Without
// FP64, R is float. C is the complex type over R (float2 or cdf). Either way
// the pixel coordinate and the accumulated sums are kept as double-floats.
//
// The group's MSL must define, using the helpers in SHADER_HEAD:
//
//   #define FADE_POLES 0 or 1          darken near |value| = 2^±14..24
//   void term(C w1, C w2, uint index, thread C& t0, thread C& t1)
//       Contribution of group element `index` to the two per-pixel sums,
//       with w1 = a z + b and w2 = c z + d from its uploaded matrix.
//   C plotValue(constant View& v, C a0, C a1, C z)
//       Value to colour from the accumulated sums a0, a1 at point z (uses v.func).
//   #define PLOT_OVERLAY 1             optional, with
//   float3 plotOverlay(constant View& v, float3 col, C z)
//       Final colour at z given the domain colouring `col`, e.g. to draw a curve.
//   #define PLOT_CUSTOM_COLOR 1        optional, with
//   bool plotCustomColor(constant View& v, C z, thread float3& col)
//       Return true to bypass phase colouring, e.g. for categorical basins.
//
// Floating-point literals need an `f` suffix: MSL has no double.

#include <string>

static const char* SHADER_HEAD = R"MSL(
#include <metal_stdlib>
using namespace metal;

// Must match ViewConstants in metal.h. Pairs are double-floats (hi, lo).
struct View {
    uint start, end, w, h;
    float cx[2], cy[2], scale[2];
    uint mode, func;
    float dens;
    uint sphere;
    float radius[2];
    float rot[9][2];  // 3x3 rotation, row j = image of view axis j (x, y, z)
    uint palette;     // phase colours: 0 hue wheel, 1 GreenPinkTones
    uint count;       // total number of summed items; when drawing, `end` of them are done
};

// ---------------------------------------------------------------------------
// Double-float arithmetic (Dekker, Knuth; as in QD's dd_real). Correct only if
// the compiler keeps the rounding of each operation, so the library is built
// with fast math off (MTLMathModeSafe), and contraction into fma is disabled.
#pragma clang fp contract(off)

struct df {
    float hi, lo;
    df() {}
    df(float h) : hi(h), lo(0.0f) {}
    df(float h, float l) : hi(h), lo(l) {}
};
inline df df2(constant float* p) { return df(p[0], p[1]); }
inline df quickTwoSum(float a, float b) { float s = a + b; return df(s, b - (s - a)); }
inline df twoSum(float a, float b) {
    float s = a + b, v = s - a;
    return df(s, (a - (s - v)) + (b - v));
}
inline df operator-(df a) { return df(-a.hi, -a.lo); }
inline df operator+(df a, df b) {
    df s = twoSum(a.hi, b.hi), t = twoSum(a.lo, b.lo);
    s.lo += t.hi;
    s = quickTwoSum(s.hi, s.lo);
    s.lo += t.lo;
    return quickTwoSum(s.hi, s.lo);
}
inline df operator-(df a, df b) { return a + (-b); }
inline df operator*(df a, df b) {
    float p = a.hi * b.hi;
    float e = fma(a.hi, b.hi, -p);  // exact error of the product
    e = fma(a.hi, b.lo, e);
    e = fma(a.lo, b.hi, e);
    return quickTwoSum(p, e);
}
inline df operator/(df a, df b) {
    float q1 = a.hi / b.hi;
    df r = a - df(q1) * b;
    float q2 = r.hi / b.hi;
    r = r - df(q2) * b;
    float q3 = r.hi / b.hi;
    return quickTwoSum(q1, q2) + df(q3);
}
inline df dsqrt(df a) {
    if (!(a.hi > 0.0f)) return df(0.0f);
    float s = sqrt(a.hi);
    df x(s), r = a - x * x;  // one Newton step doubles the precision
    return x + df(r.hi * (0.5f / s));
}
inline df ldexpdf(df a, int k) { return df(ldexp(a.hi, k), ldexp(a.lo, k)); }
inline bool operator<(df a, df b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
inline bool operator>(df a, df b) { return b < a; }
inline bool operator>=(df a, df b) { return !(a < b); }
inline float toF(df a) { return a.hi; }
inline float toF(float a) { return a; }

struct cdf {
    df x, y;
    cdf() {}
    cdf(df a, df b) : x(a), y(b) {}
    cdf(float a, float b) : x(a), y(b) {}
};
inline cdf operator+(cdf a, cdf b) { return cdf(a.x + b.x, a.y + b.y); }
inline cdf operator-(cdf a, cdf b) { return cdf(a.x - b.x, a.y - b.y); }
inline cdf operator-(cdf a) { return cdf(-a.x, -a.y); }
inline cdf operator*(df s, cdf a) { return cdf(s * a.x, s * a.y); }
inline cdf operator*(cdf a, df s) { return cdf(a.x * s, a.y * s); }
inline cdf operator/(cdf a, df s) { return cdf(a.x / s, a.y / s); }
inline float2 toF2(cdf a) { return float2(a.x.hi, a.y.hi); }
inline float2 toF2(float2 a) { return a; }
inline cdf toDF(cdf a) { return a; }
inline cdf toDF(float2 a) { return cdf(a.x, a.y); }

// Complex helpers for both float2 and cdf.
template <typename T> inline T cm(T a, T b) { return T(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
template <typename T> inline T conjugate(T a) { return T(a.x, -a.y); }
inline float abs2(float2 a) { return a.x * a.x + a.y * a.y; }
inline df abs2(cdf a) { return a.x * a.x + a.y * a.y; }
inline float2 crecip(float2 q) {
    float rs = 1.0f / max(abs(q.x), abs(q.y));
    float2 t = q * rs;
    return float2(t.x, -t.y) * (rs / dot(t, t));
}
inline cdf crecip(cdf q) {
    // Scale by a power of two (exact) so |q|^2 neither overflows nor underflows.
    int e = -int(clamp(floor(log2(max(abs(q.x.hi), abs(q.y.hi)))), -120.0f, 120.0f));
    cdf t = cdf(ldexpdf(q.x, e), ldexpdf(q.y, e));
    df n = abs2(t);
    return cdf(ldexpdf(t.x / n, e), ldexpdf(-t.y / n, e));
}
template <typename T> inline T cdiv(T a, T b) { return cm(a, crecip(b)); }

#if FP64
typedef df R;
typedef cdf C;
inline C toC(cdf z) { return z; }
#else
typedef float R;
typedef float2 C;
inline C toC(cdf z) { return float2(z.x.hi, z.y.hi); }
#endif

// exp(i p) for complex p, and sqrt, in the precision of R. Shared by the
// Pearcey integrand and the q = exp(2 pi i tau) of the modular programs.
#if FP64
// Double-float exp, and exp(i p) for complex p. ln 2 and pi/2 are given as
// double-floats (hi, lo).
df dexp(df x) {  // x = k ln 2 + r, |r| <= 0.35, Taylor series in r
    if (x.hi < -80.0f) return df(0.0f);
    float k = rint(x.hi * 1.44269504f);
    df r = x - df(k) * df(0.6931472f, -1.9046542e-09f);
    df p = df(1.0f);
    for (int n = 13; n >= 1; --n) p = df(1.0f) + r * p / df(float(n));
    return ldexpdf(p, int(k));
}
cdf expi(cdf p) {
    // p.x = k pi/2 + r, |r| <= pi/4, Taylor series in r
    float k = rint(p.x.hi * 0.636619772f);
    df r = p.x - df(k) * df(1.5707964f, -4.371139e-08f), r2 = r * r;
    df c = df(1.0f), s = df(1.0f);
    for (int n = 16; n >= 2; n -= 2) {
        c = df(1.0f) - r2 * c / df(float(n * (n - 1)));
        s = df(1.0f) - r2 * s / df(float(n * (n + 1)));
    }
    s = s * r;
    int q = int(k) & 3;
    cdf e = q == 0 ? cdf(c, s) : (q == 1 ? cdf(-s, c) : (q == 2 ? cdf(-c, -s) : cdf(s, -c)));
    return e * dexp(-p.y);
}
#else
float dsqrt(float d) { return sqrt(d); }
float2 expi(float2 p) { float c; float s = sincos(p.x, c); return float2(c, s) * exp(-p.y); }
#endif

// Complex coordinate of pixel p (top-left origin); App::toZ in app.h mirrors it.
// Plane view: z = centre + offset * scale.
// Sphere view: orthographic view of a sphere of radius v.radius pixels, turned by v.rot, and
//   z = (X + iY) / (1 + Z), stereographic projection from the south pole, so
//   0 is the north pole and infinity the south pole. Seen from outside, this
//   map preserves orientation, so the sphere and the plane are not mirrored.
//   Returns false off the sphere. `front` is the view-space depth (1 at the
//   centre of the disc, 0 at the rim) and `cover` the antialiased disc coverage.
// Computed in double-float in both modes, to keep deep zooms sharp.
bool pixelZ(constant View& v, uint2 p, thread cdf& z, thread float& front, thread float& cover) {
    // Pixel offsets are exact in float at supported window sizes.
    float2 offset = float2(float(p.x) + 0.5f - 0.5f * float(v.w), 0.5f * float(v.h) - 0.5f - float(p.y));
    front = 1.0f; cover = 1.0f;
    df scale = df2(v.scale);
    if (v.sphere == 0) {
        z = cdf(df2(v.cx) + df(offset.x) * scale, df2(v.cy) + df(offset.y) * scale);
        return true;
    }
    df rad = df2(v.radius), rad2 = rad * rad;
    df o2 = df(offset.x) * df(offset.x) + df(offset.y) * df(offset.y);
    cover = saturate(rad.hi - length(offset) + 0.5f);
    if (o2 >= rad2) { z = cdf(0.0f, 0.0f); front = 0.0f; return false; }
    df S = dsqrt(rad2 - o2);
    front = toF(S / rad);
    df V[3];
    for (int i = 0; i < 3; i++)
        V[i] = df(offset.x) * df2(v.rot[i]) + df(offset.y) * df2(v.rot[3 + i]) + S * df2(v.rot[6 + i]);
    // z = (X + iY) / (1 + Z) = V.xy / (R + V.z). Near the south pole R + V.z
    // cancels, so use the equal form V.xy (R - V.z) / (V.x^2 + V.y^2) there.
    cdf Vxy = cdf(V[0], V[1]);
    if (V[2].hi >= 0.0f) z = Vxy / (rad + V[2]);
    else z = Vxy * ((rad - V[2]) / (V[0] * V[0] + V[1] * V[1]));
    return true;
}
)MSL";

static const char* SHADER_TAIL = R"MSL(
// Matrix k (0..3 = a, b, c, d) of element i. fp32: 8 floats per element;
// FP64: 16 floats, each entry's real and imaginary parts as (hi, lo).
inline C loadEntry(device const float* m, uint i, uint k) {
#if FP64
    device const float* e = m + 16 * i + 4 * k;
    return cdf(df(e[0], e[1]), df(e[2], e[3]));
#else
    return float2(m[8 * i + 2 * k], m[8 * i + 2 * k + 1]);
#endif
}

// Each pixel's two complex sums are kept as double-floats: acc[2 i] holds the
// hi parts (s0.x, s0.y, s1.x, s1.y), acc[2 i + 1] the lo parts.
kernel void computeKernel(uint2 p [[thread_position_in_grid]],
                          constant View& v [[buffer(0)]],
                          device const float* mats [[buffer(1)]],
                          device float4* acc [[buffer(2)]]) {
    if (p.x >= v.w || p.y >= v.h) return;
    cdf zd; float front, cover;
    if (!pixelZ(v, p, zd, front, cover)) return;
    C z = toC(zd);
    C s0 = C(0.0f, 0.0f), s1 = C(0.0f, 0.0f);
    for (uint i = v.start; i < v.end; ++i) {
        C w1 = cm(loadEntry(mats, i, 0), z) + loadEntry(mats, i, 1);
        C w2 = cm(loadEntry(mats, i, 2), z) + loadEntry(mats, i, 3);
        C t0, t1;
        term(w1, w2, i, t0, t1);
        s0 = s0 + t0;
        s1 = s1 + t1;
    }
    cdf a0 = toDF(s0), a1 = toDF(s1);
    uint index = p.y * v.w + p.x;
    if (v.start != 0) {
        float4 hi = acc[2 * index], lo = acc[2 * index + 1];
        a0 = a0 + cdf(df(hi.x, lo.x), df(hi.y, lo.y));
        a1 = a1 + cdf(df(hi.z, lo.z), df(hi.w, lo.w));
    }
    acc[2 * index] = float4(a0.x.hi, a0.y.hi, a1.x.hi, a1.y.hi);
    acc[2 * index + 1] = float4(a0.x.lo, a0.y.lo, a1.x.lo, a1.y.lo);
}

struct VOut { float4 pos [[position]]; };
vertex VOut vertexMain(uint id [[vertex_id]]) {
    float2 p = float2(float((id << 1) & 2), float(id & 2));
    VOut o;
    o.pos = float4(p * 2.0f - 1.0f, 0.0f, 1.0f);
    return o;
}
float3 hsv(float h) {
    float3 k = saturate(abs((fract((h * 6.0f + float3(0, 4, 2)) / 6.0f) * 6.0f) - 3.0f) - 1.0f);
    return mix(float3(1.0f), k * k * (3.0f - 2.0f * k), 0.9f);
}
// Mathematica's ColorData["GreenPinkTones"] at 0, 1/64, ..., 1: dark green,
// light green, near white, pink, dark magenta.
constant float3 GREEN_PINK[65] = {
float3(0.000,0.242,0.023),float3(0.000,0.292,0.028),float3(0.000,0.342,0.033),float3(0.000,0.391,0.038),float3(0.000,0.441,0.042),float3(0.001,0.493,0.048),float3(0.003,0.547,0.054),float3(0.005,0.600,0.061),
float3(0.007,0.653,0.067),float3(0.015,0.702,0.077),float3(0.031,0.745,0.092),float3(0.046,0.788,0.106),float3(0.061,0.831,0.121),float3(0.082,0.870,0.139),float3(0.123,0.889,0.171),float3(0.163,0.907,0.204),
float3(0.204,0.926,0.236),float3(0.245,0.945,0.269),float3(0.296,0.952,0.315),float3(0.347,0.958,0.361),float3(0.398,0.964,0.408),float3(0.450,0.970,0.454),float3(0.502,0.973,0.505),float3(0.555,0.974,0.557),
float3(0.609,0.975,0.609),float3(0.662,0.976,0.661),float3(0.711,0.975,0.710),float3(0.754,0.968,0.753),float3(0.797,0.961,0.796),float3(0.840,0.955,0.839),float3(0.880,0.946,0.879),float3(0.902,0.923,0.901),
float3(0.924,0.901,0.923),float3(0.946,0.878,0.945),float3(0.968,0.855,0.967),float3(0.976,0.816,0.974),float3(0.982,0.775,0.980),float3(0.988,0.735,0.986),float3(0.994,0.694,0.991),float3(0.996,0.649,0.993),
float3(0.995,0.600,0.992),float3(0.995,0.552,0.991),float3(0.995,0.504,0.990),float3(0.992,0.459,0.987),float3(0.986,0.419,0.981),float3(0.979,0.379,0.975),float3(0.973,0.339,0.968),float3(0.965,0.300,0.961),
float3(0.943,0.278,0.940),float3(0.922,0.255,0.919),float3(0.900,0.232,0.898),float3(0.878,0.209,0.877),float3(0.838,0.194,0.837),float3(0.793,0.180,0.793),float3(0.748,0.166,0.748),float3(0.704,0.153,0.704),
float3(0.654,0.140,0.654),float3(0.601,0.129,0.601),float3(0.548,0.117,0.547),float3(0.495,0.105,0.494),float3(0.443,0.094,0.441),float3(0.393,0.083,0.391),float3(0.343,0.073,0.341),float3(0.294,0.062,0.291),
float3(0.244,0.051,0.241)};
// Colour of the phase h = arg/2pi in [0, 1). GreenPinkTones is not cyclic, so
// it runs from arg = -pi (dark green) through 0 (near white) to +pi (dark
// magenta), with its seam on the negative real axis.
float3 phaseColor(constant View& v, float h) {
    if (v.palette == 0) return hsv(h);
    float t = fract(h + 0.5f) * 64.0f;
    uint i = min(uint(t), 63u);
    return mix(GREEN_PINK[i], GREEN_PINK[i + 1], t - float(i));
}
float saw(float x) { return mix(fract(x), 0.5f, smoothstep(0.3f, 0.8f, fwidth(x))); }
float gridline(float x) {
    return 1.0f - smoothstep(0.0f, 1.2f, abs(fract(x - 0.5f) - 0.5f) / max(fwidth(x), 1e-6f));
}
fragment float4 fragmentMain(VOut in [[stage_in]],
                             constant View& v [[buffer(0)]],
                             device const float4* acc [[buffer(1)]]) {
    // Clamp: quads on the bottom/right edge run helper lanes (needed for
    // fwidth) one pixel outside the target, which would read past acc.
    uint2 p = min(uint2(in.pos.xy), uint2(v.w - 1, v.h - 1));
    cdf zd; float front, cover;
    // Off the sphere acc is not written; colour a constant instead (no early
    // return, so neighbouring on-sphere lanes keep valid derivatives).
    bool on = pixelZ(v, p, zd, front, cover);
    C z = toC(zd);
#ifdef PLOT_CUSTOM_COLOR
    float3 custom;
    if (on && plotCustomColor(v, z, custom)) return float4(custom, 1.0f);
#endif
    C fd = C(1.0f, 0.0f);
    if (on) {
        uint index = p.y * v.w + p.x;
        float4 hi = acc[2 * index], lo = acc[2 * index + 1];
#if FP64
        C a0 = cdf(df(hi.x, lo.x), df(hi.y, lo.y)), a1 = cdf(df(hi.z, lo.z), df(hi.w, lo.w));
#else
        C a0 = hi.xy, a1 = hi.zw;
#endif
        fd = plotValue(v, a0, a1, z);
    }
    float2 f = toF2(fd);
    float s = max(abs(f.x), abs(f.y));
    if (!(s > 0.0f)) return float4(0, 0, 0, 1);
    float2 dir = f / s;
    if (any(isnan(dir)) || any(isinf(dir))) return float4(0, 0, 0, 1);
    float lm = log2(s) + 0.5f * log2(dot(dir, dir));
    if (isinf(s)) lm = 128.0f;
    float h = fract(atan2(dir.y, dir.x) / (2.0f * 3.14159265358979f));
    float3 base = phaseColor(v, h), col;
    float dens = v.dens;
    if (v.mode == 0) col = base * (0.55f + 0.45f * sqrt(saw(lm * dens) * saw(h * 12.0f * dens)));
    else if (v.mode == 1) col = base;
    else if (v.mode == 2) col = base * (0.45f + 0.55f * saw(lm * dens));
    else {
        float2 fs = dir * exp2(log2(s) - floor(lm)) * 4.0f * dens;
        col = mix(base * 0.95f, float3(0.05f), 0.75f * max(gridline(fs.x), gridline(fs.y)));
    }
#if FADE_POLES
    col *= 1.0f - smoothstep(14.0f, 24.0f, abs(lm));
#endif
#ifdef PLOT_OVERLAY
    col = plotOverlay(v, col, z);
#endif
    if (v.sphere != 0) col = mix(float3(0.045f, 0.05f, 0.065f), col * (0.55f + 0.45f * front), cover);
    return float4(pow(abs(col), 1.0f / 1.1f), 1.0f);
}
)MSL";

inline std::string shaderSource(const char* groupMsl, bool fp64) {
    return std::string(fp64 ? "#define FP64 1\n" : "#define FP64 0\n") + SHADER_HEAD + groupMsl + SHADER_TAIL;
}
