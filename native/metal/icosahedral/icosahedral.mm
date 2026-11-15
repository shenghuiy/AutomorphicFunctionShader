// Native Metal renderer for the icosahedral automorphic function
//     F(z) = sum_{g in G} h(g(z)),   h(z) = (z - 0.05)^-5,
// where G is the group of 60 rotations of the icosahedron acting on the
// Riemann sphere (../../dx12/icosahedral/icosahedral.h). G is finite, so the sum is exact.
//
// The window, renderer and controls are in ../common/app.h. See docs/coxeter.md.

#include "app.h"
#include "icosahedral.h"

// The uploaded matrices are pre-multiplied by [[1, -0.05], [0, 1]]
// (icosahedral::gpuMatrices), so w2/w1 = 1/(g(z) - 0.05) and each term is h(g(z)).
// Sums: t0 = F(z) = sum h(g(z)),  t1 = h(z) (the identity term only)
static const char* ICOSAHEDRAL_MSL = R"MSL(
// F is typically around 2^13, and its poles read clearly without darkening.
#define FADE_POLES 0
void term(C w1, C w2, uint index, thread C& t0, thread C& t1) {
    C r = cm(w2, crecip(w1));
    C r2 = cm(r, r);
    C t = cm(cm(r2, r2), r);
    t0 = t;
    t1 = index == 0 ? t : C(0.0f, 0.0f);
}
// Klein's invariant J = H^3/(1728 f^5). z -> -1/z is in G, so evaluate
// inside the unit disk to keep the powers small.
C kleinJ(C z) {
    if (toF(abs2(z)) > 1.0f) z = -crecip(z);
    C one = C(1.0f, 0.0f);
    C z5 = cm(cm(cm(z, z), cm(z, z)), z), z10 = cm(z5, z5), z15 = cm(z10, z5), z20 = cm(z10, z10);
    C f = cm(z, z10 + R(11.0f) * z5 - one);
    C H = -(z20 + one) + R(228.0f) * (z15 - z5) - R(494.0f) * z10;
    C f2 = cm(f, f);
    return cdiv(cm(cm(H, H), H), R(1728.0f) * cm(cm(f2, f2), f));
}
// F key: F = sum h(g(z)), h, Klein's J
C plotValue(constant View& v, C a0, C a1, C z) {
    return v.func == 0 ? a0 : (v.func == 1 ? a1 : kleinJ(z));
}
)MSL";

static std::vector<Mat> rotations;  // the 60 rotations (unshifted)

static std::vector<Mat> build(size_t) {
    rotations = icosahedral::group();
    if (rotations.size() != icosahedral::ORDER) fatal("Icosahedral group closure did not give 60 elements");
    consoleLog("group: icosahedral rotations, %zu elements (exact sum)\n", rotations.size());
    return icosahedral::gpuMatrices(rotations);
}

static void reference(cd z, cd out[2]) {
    out[0] = icosahedral::F(rotations, z);
    out[1] = icosahedral::h(z);
}

static cd value(cd F, cd h, int func, cd z) { return func == 0 ? F : (func == 1 ? h : icosahedral::kleinJ(z)); }

static void invariance() {
    cd z0(0.3, 0.6);
    cd Fz = icosahedral::F(rotations, z0);
    consoleLog("\ninvariance check, icosahedral group (exact), z0 = 0.3+0.6i\n  F(z0)        = %+.8g %+.8gi\n", Fz.real(), Fz.imag());
    for (int k = 0; k < 4; k++) {
        cd Fg = icosahedral::F(rotations, mobius(icosahedral::GENS[k], z0));
        consoleLog("  F(%-5s(z0)) = %+.8g %+.8gi   rel.err %.1e\n", icosahedral::GEN_NAMES[k], Fg.real(), Fg.imag(), std::abs(Fg - Fz) / std::abs(Fz));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Icosahedral group",
        "F(z) = sum_g h(g(z)),  h(z) = (z - 0.05)^-5,  G = 60 rotations of the icosahedron",
        "icosahedral",
        "icosahedral",
        {"F", "h", "J"},
        {"F", "h"},
        0.0, 0.0, 5.0,
        0,  // finite group: always all 60 elements
        ICOSAHEDRAL_MSL,
        build, reference, value, invariance,
        true,  // G acts on the whole sphere: offer the Riemann sphere view
    };
    return runApp(group, argc, argv);
}
