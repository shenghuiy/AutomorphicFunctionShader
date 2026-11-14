// Native Metal renderer for the Pearcey integral
//     Psi(x, y) = integral over the real line of exp(i (t^4 + y t^2 + x t)) dt
// as a function of z = x + i y. Each pixel sums a midpoint rule on a path bent
// into the complex plane so that the integrand decays (../../dx12/pearcey/pearcey.h).
//
// The window, renderer and controls are in ../common/app.h. There the summed
// items are called group elements; here they are quadrature nodes.

#include "app.h"
#include "pearcey.h"

// Each node is uploaded as [[1, 0], [0, u + i w]] (pearcey::nodes), so w1 = z = (x, y)
// and w2 = (u, w): position in (-1, 1) and weight.
// Sums: t0 = Psi,  t1 = dPsi/dx
static const char* PEARCEY_MSL = R"MSL(
#define FADE_POLES 0
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

// Path t(s) = s + i LAMBDA g / sqrt(1 + g^2), g = phi'(s) / kappa; see pearcey.h.
void term(C z, C node, uint index, thread C& t0, thread C& t1) {
    R x = z.x, y = z.y, one = R(1.0f);
    // pearcey::halfWidth, in float like the CPU
    R S = R(3.0f + max(sqrt(abs(toF(y))), pow(0.5f * abs(toF(x)), 1.0f / 3.0f)));
    R s = S * node.x, s2 = s * s, y2 = y + y;
    R kappa = R(4.0f) * S;  // 16 S LAMBDA^2 with LAMBDA = 1/2
    R g = (R(4.0f) * s2 * s + y2 * s + x) / kappa, dg = (R(12.0f) * s2 + y2) / kappa;
    R q = dsqrt(one + g * g), q3 = q * q * q;
    C t = C(s, g / (q + q)), dt = C(one, dg / (q3 + q3));
    C t2 = cm(t, t);
    C v = cm(expi(cm(t2, t2) + y * t2 + x * t), dt) * (S * node.y);
    t0 = v;
    t1 = cm(C(-t.y, t.x), v);
}
// F key: Psi, Psi with the caustic drawn, dPsi/dx
// While the sum is incomplete, the first `end` of `count` nodes give end/count
// of the integral (pearcey::nodes), so scale it back up.
C plotValue(constant View& v, C a0, C a1, C z) {
    return (v.func == 2 ? a1 : a0) * (R(float(v.count)) / R(float(max(v.end, 1u))));
}
// The caustic 27 x^2 + 8 y^3 = 0: the (x, y) where two real stationary points
// of the phase merge. Drawn where |c| / |grad c|, the distance to the curve to
// first order, is under about 1.5 pixels.
#define PLOT_OVERLAY 1
float3 plotOverlay(constant View& v, float3 col, C zc) {
    if (v.func != 1) return col;
    float2 z = toF2(zc);
    float c = 27.0f * z.x * z.x + 8.0f * z.y * z.y * z.y;
    float d = abs(c) / max(length(float2(54.0f * z.x, 24.0f * z.y * z.y)), 1e-30f);
    return col * smoothstep(1.0f, 2.0f, d / v.scale[0]);
}
)MSL";

static size_t nodeCount = 0;

static std::vector<Mat> build(size_t requested) {
    size_t n = 1;
    while (n < requested) n *= 2;  // pearcey::nodes needs a power of two
    nodeCount = n;
    std::printf("quadrature: %zu midpoint nodes per pixel on the deformed path\n", n);
    return pearcey::nodes(n);
}

// A much finer rule on the same path; the GPU sum has converged if it agrees.
static void reference(cd z, cd out[2]) { pearcey::psi(z.real(), z.imag(), 20000, out); }

static cd value(cd psi, cd dpsi, int func, cd) { return func == 2 ? dpsi : psi; }

// I key. There is no group here; check the quadrature instead.
static void invariance() {
    auto err = [](cd a, cd b) { return std::abs(a - b) / std::abs(b); };
    cd p[2], q[2];
    pearcey::psi(0, 0, nodeCount, p);
    cd exact = 2 * std::tgamma(1.25) * std::polar(1.0, 3.14159265358979323846 / 8);
    std::printf("\nquadrature check, %zu nodes\n", nodeCount);
    std::printf("  Psi(0,0)      = %+.12f %+.12fi   vs 2 Gamma(5/4) e^{i pi/8}: rel.err %.1e\n", p[0].real(), p[0].imag(), err(p[0], exact));
    const double pts[][2] = {{1.5, -2.0}, {-3.0, 1.0}, {2.0, -4.0}, {7.0, -12.0}};
    for (const auto& xy : pts) {
        pearcey::psi(xy[0], xy[1], nodeCount, p);
        std::printf("  Psi(%g,%g) = %+.12f %+.12fi", xy[0], xy[1], p[0].real(), p[0].imag());
        if (xy[1] >= -4) std::printf("   vs rotated path: rel.err %.1e", err(p[0], pearcey::psiRotated(xy[0], xy[1])));
        pearcey::psi(xy[0], xy[1], 20000, q);
        std::printf("   vs 20000 nodes: rel.err %.1e", err(p[0], q[0]));
        // Consistency check only: for -x the path is the mirror image t -> -t and the
        // nodes are symmetric, so this sums the same terms in reverse order and
        // agrees to rounding error by construction.
        pearcey::psi(-xy[0], xy[1], nodeCount, q);
        std::printf("   Psi(-x,y) = Psi(x,y): rel.err %.1e\n", err(q[0], p[0]));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Pearcey integral",
        "Psi(x, y) = int exp(i (t^4 + y t^2 + x t)) dt,  z = x + i y",
        "pearcey",
        "pearcey",
        {"Psi", "Psi + caustic", "dPsi/dx"},
        {"Psi", "dPsi/dx"},
        0.0, -5.0, 20.0,  // |x| < 10, -15 < y < 5 on the shorter side
        512,              // quadrature nodes; --terms and +/-
        PEARCEY_MSL,
        build, reference, value, invariance,
        false,  // no sphere view
        128,    // partial sums over 128, 256, ... nodes are valid coarser rules
    };
    return runApp(group, argc, argv);
}
