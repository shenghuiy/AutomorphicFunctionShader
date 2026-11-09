// Native DirectX 12 renderer for the automorphic function
//     f(z) = P_{-1,-3}(z) / P_{0,-4}(z),   P_{k,l}(z) = sum_{g in G} (az+b)^k (cz+d)^l
// for the Kleinian group G = <g1, g2> (kleinian.h).
//
// The Poincaré series are summed per pixel in a compute shader (fp32 or fp64),
// accumulated progressively in a double-precision buffer, and displayed with
// domain colouring. The window, renderer and controls are in ../common/app.h.

#include "app.h"
#include "kleinian.h"

// Accumulator: .xy = P_{-1,-3}(z) = sum 1/((az+b)(cz+d)^3),  .zw = P_{0,-4}(z) = sum 1/(cz+d)^4
static const char* KLEINIAN_HLSL = R"HLSL(
#define FADE_POLES 1
#ifdef COMPUTE
VEC4 term(VEC2 w1, VEC2 w2, uint index) {
    VEC2 w2s=cm(w2,w2);
    VEC2 r=crecip(cm(w1,cm(w2s,w2s)));
    return VEC4(cm(w2,r),cm(w1,r));
}
#else
// F key: f = P_{-1,-3}/P_{0,-4}, P_{0,-4}, P_{-1,-3}
double2 plotValue(double4 a, double2 z) {
    return uFunc==0 ? cdiv(a.xy,a.zw) : (uFunc==1 ? a.zw : a.xy);
}
#endif
)HLSL";

static std::vector<Mat> elements;  // the N group elements currently summed

static std::vector<Mat> build(size_t N) {
    auto t0 = Clock::now();
    elements = kleinian::enumerate(N);
    std::printf("group: %zu elements, max ||g|| = %.1f  (enumerated in %.0f ms)\n", elements.size(),
                std::sqrt(norm2(elements.back())), msSince(t0));
    return elements;
}

static void reference(cd z, cd out[2]) {
    auto s = kleinian::evalSeries(elements, z);
    out[0] = s.P1;
    out[1] = s.P2;
}

static cd value(cd P1, cd P2, int func, cd) { return func == 0 ? P1 / P2 : (func == 1 ? P2 : P1); }

static void invariance() {
    cd z0(0.3, 0.6);
    cd fz = kleinian::evalSeries(elements, z0).f;
    std::printf("\ninvariance check with %zu elements, z0 = 0.3+0.6i\n  f(z0)        = %+.8f %+.8fi\n", elements.size(), fz.real(), fz.imag());
    for (int k = 0; k < 4; k++) {
        cd fg = kleinian::evalSeries(elements, mobius(kleinian::GENS[k], z0)).f;
        std::printf("  f(%-5s(z0)) = %+.8f %+.8fi   rel.err %.1e\n", kleinian::GEN_NAMES[k], fg.real(), fg.imag(), std::abs(fg - fz) / std::abs(fz));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Poincaré series",
        "f(z) = P_{-1,-3}(z) / P_{0,-4}(z),  G = <g1, g2>",
        "poincare.exe",
        "poincare",
        {"f", "P0,-4", "P-1,-3"},
        {"P-1,-3", "P0,-4"},
        0.0, -0.35, 5.0,
        size_t(1) << 17,
        KLEINIAN_HLSL,
        build, reference, value, invariance,
    };
    return runApp(group, argc, argv);
}
