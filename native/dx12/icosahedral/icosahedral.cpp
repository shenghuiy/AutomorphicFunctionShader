// Native DirectX 12 renderer for the icosahedral automorphic function
//     F(z) = sum_{g in G} h(g(z)),   h(z) = (z - 0.05)^-5,
// where G is the group of 60 rotations of the icosahedron acting on the
// Riemann sphere (icosahedral.h). G is finite, so the sum is exact.
//
// The window, renderer and controls are in ../common/app.h. See docs/coxeter.md.

#include "app.h"
#include "icosahedral.h"

// The uploaded matrices are pre-multiplied by [[1, -0.05], [0, 1]]
// (icosahedral::gpuMatrices), so w2/w1 = 1/(g(z) - 0.05) and each term is h(g(z)).
// Accumulator: .xy = F(z) = sum h(g(z)),  .zw = h(z) (the identity term only)
static const char* ICOSAHEDRAL_HLSL = R"HLSL(
// F is typically around 2^13, and its poles read clearly without darkening.
#define FADE_POLES 0
#ifdef COMPUTE
VEC4 term(VEC2 w1, VEC2 w2, uint index) {
    VEC2 r=cm(w2,crecip(w1));
    VEC2 r2=cm(r,r);
    VEC2 t=cm(cm(r2,r2),r);
    VEC4 s=VEC4(t,(VEC2)0);
    if (index==0) s.zw=t;
    return s;
}
#else
// Klein's invariant J = H^3/(1728 f^5). z -> -1/z is in G, so evaluate
// inside the unit disk to keep the powers small.
double2 kleinJ(double2 z) {
    if (z.x*z.x+z.y*z.y>1.0L) z=-cdiv(double2(1,0),z);
    double2 z5=cmd(cmd(cmd(z,z),cmd(z,z)),z), z10=cmd(z5,z5), z15=cmd(z10,z5), z20=cmd(z10,z10);
    double2 f=cmd(z,z10+11.0L*z5-double2(1,0));
    double2 H=-(z20+double2(1,0))+228.0L*(z15-z5)-494.0L*z10;
    double2 f2=cmd(f,f);
    return cdiv(cmd(cmd(H,H),H),1728.0L*cmd(cmd(f2,f2),f));
}
// F key: F = sum h(g(z)), h, Klein's J
double2 plotValue(double4 a, double2 z) {
    return uFunc==0 ? a.xy : (uFunc==1 ? a.zw : kleinJ(z));
}
#endif
)HLSL";

static std::vector<Mat> rotations;  // the 60 rotations (unshifted)

static std::vector<Mat> build(size_t) {
    rotations = icosahedral::group();
    if (rotations.size() != icosahedral::ORDER) fatal("Icosahedral group closure did not give 60 elements");
    std::printf("group: icosahedral rotations, %zu elements (exact sum)\n", rotations.size());
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
    std::printf("\ninvariance check, icosahedral group (exact), z0 = 0.3+0.6i\n  F(z0)        = %+.8g %+.8gi\n", Fz.real(), Fz.imag());
    for (int k = 0; k < 4; k++) {
        cd Fg = icosahedral::F(rotations, mobius(icosahedral::GENS[k], z0));
        std::printf("  F(%-5s(z0)) = %+.8g %+.8gi   rel.err %.1e\n", icosahedral::GEN_NAMES[k], Fg.real(), Fg.imag(), std::abs(Fg - Fz) / std::abs(Fz));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Icosahedral group",
        "F(z) = sum_g h(g(z)),  h(z) = (z - 0.05)^-5,  G = 60 rotations of the icosahedron",
        "icosahedral.exe",
        "icosahedral",
        {"F", "h", "J"},
        {"F", "h"},
        0.0, 0.0, 5.0,
        0,  // finite group: always all 60 elements
        ICOSAHEDRAL_HLSL,
        build, reference, value, invariance,
        true,  // G acts on the whole sphere: offer the Riemann sphere view
    };
    return runApp(group, argc, argv);
}
