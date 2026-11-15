// Native Metal renderer for the automorphic function
//     f(z) = P_{-1,-3}(z) / P_{0,-4}(z),   P_{k,l}(z) = sum_{g in G} (az+b)^k (cz+d)^l
// for the Kleinian group G = <g1, g2> (../../dx12/poincare/kleinian.h).
//
// The Poincaré series are summed per pixel in a compute kernel (fp32 or
// double-float), accumulated progressively in a double-float buffer, and
// displayed with domain colouring. The window, renderer and controls are in
// ../common/app.h.

#include "app.h"
#include "kleinian.h"

// Sums: t0 = P_{-1,-3}(z) = sum 1/((az+b)(cz+d)^3),  t1 = P_{0,-4}(z) = sum 1/(cz+d)^4
static const char* KLEINIAN_MSL = R"MSL(
#define FADE_POLES 1
void term(C w1, C w2, uint index, thread C& t0, thread C& t1) {
    C w2s = cm(w2, w2);
    C r = crecip(cm(w1, cm(w2s, w2s)));
    t0 = cm(w2, r);
    t1 = cm(w1, r);
}
// F key: f = P_{-1,-3}/P_{0,-4}, P_{0,-4}, P_{-1,-3}
C plotValue(constant View& v, C a0, C a1, C z) {
    return v.func == 0 ? cdiv(a0, a1) : (v.func == 1 ? a1 : a0);
}
)MSL";

static std::vector<Mat> elements;  // the N group elements currently summed

static std::vector<Mat> build(size_t N) {
    auto t0 = Clock::now();
    elements = kleinian::enumerate(N);
    consoleLog("group: %zu elements, max ||g|| = %.1f  (enumerated in %.0f ms)\n", elements.size(),
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
    consoleLog("\ninvariance check with %zu elements, z0 = 0.3+0.6i\n  f(z0)        = %+.8f %+.8fi\n", elements.size(), fz.real(), fz.imag());
    for (int k = 0; k < 4; k++) {
        cd fg = kleinian::evalSeries(elements, mobius(kleinian::GENS[k], z0)).f;
        consoleLog("  f(%-5s(z0)) = %+.8f %+.8fi   rel.err %.1e\n", kleinian::GEN_NAMES[k], fg.real(), fg.imag(), std::abs(fg - fz) / std::abs(fz));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Poincaré series",
        "f(z) = P_{-1,-3}(z) / P_{0,-4}(z),  G = <g1, g2>",
        "poincare",
        "poincare",
        {"f", "P0,-4", "P-1,-3"},
        {"P-1,-3", "P0,-4"},
        0.0, -0.35, 5.0,
        size_t(1) << 17,
        KLEINIAN_MSL,
        build, reference, value, invariance,
    };
    return runApp(group, argc, argv);
}
