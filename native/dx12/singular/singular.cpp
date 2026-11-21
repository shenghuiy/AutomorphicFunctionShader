#include "app.h"

// The Cantor measure is approximated by equal masses at the midpoints of
// the surviving intervals of a ternary Cantor construction. Each node has
// weight 1/N. As N=2^depth grows, these discrete measures converge weakly
// to the non-atomic Cantor probability measure pushed onto the unit circle.
static unsigned cantorDepth = 10;
static unsigned atomCount = 1u << cantorDepth;
static constexpr double tau = 6.283185307179586;

static const char* SINGULAR_HLSL = R"HLSL(
#define FADE_POLES 1
#ifdef COMPUTE
VEC4 term(VEC2 delta, VEC2 numerator, uint index) {
    if (delta.x == 0 && delta.y == 0) return (VEC4)0;
    VEC2 exponent = cm(numerator, crecip(delta));
    return VEC4(exponent, index == 0 ? exponent : (VEC2)0);
}
#else
double2 plotValue(double4 accumulator, double2 z) {
    if (uFunc == 2) {
        double2 alpha = double2(0.35L, 0.2L);
        return cdiv(z-alpha, double2(1,0)-cmd(double2(alpha.x,-alpha.y),z));
    }

    double2 exponent;
    if (uFunc == 0) exponent = accumulator.xy / (double)uCount;
    else exponent = cdiv(z+double2(1,0),z-double2(1,0));
    float phase = frac((float)(exponent.y / 6.283185307179586L)) * 6.28318530718;
    float magnitude = exp(clamp((float)exponent.x, -70.0, 70.0));
    return magnitude * float2(cos(phase), sin(phase));
}
#endif
)HLSL";

static std::vector<cd> support;

static std::vector<Mat> build(size_t requested) {
    cantorDepth = 0;
    atomCount = 1;
    while (atomCount <= requested/2 && cantorDepth < 20) {
        atomCount *= 2;
        ++cantorDepth;
    }
    std::printf("Cantor quadrature: depth %u, %u equal weights, total mass 1\n",
                cantorDepth, atomCount);
    support.clear();
    std::vector<Mat> matrices;
    for (unsigned index = 0; index < atomCount; ++index) {
        double left = 0;
        double width = 1;
        for (unsigned digit = 0; digit < cantorDepth; ++digit) {
            width /= 3;
            if ((index >> digit) & 1u) left += 2 * width;
        }
        const cd w = std::polar(1.0, tau * (left + width / 2));
        support.push_back(w);
        matrices.push_back({1, -w, 1, w});
    }
    return matrices;
}

static void reference(cd z, cd out[2]) {
    out[0] = out[1] = 0;

    for (cd w : support) out[0] += (z + w) / (z - w);
    out[1] = (z + support.front()) / (z - support.front());
}

static cd value(cd cantor, cd firstNode, int function, cd z) {

    const cd alpha(0.35, 0.2);
    if (function == 2) return (z - alpha) / (1.0 - std::conj(alpha) * z);
    return std::exp(function == 0 ? cantor / double(atomCount) : (z+1.0)/(z-1.0));
}

static void checkCenter() {
    cd sums[2];
    reference(0, sums);
    std::printf("At z=0, Cantor modulus = %.9f, atomic modulus = %.9f (expected exp(-1))\n",
                std::abs(value(sums[0], sums[1], 0, 0)),
                std::abs(value(sums[0], sums[1], 1, 0)));
}

int main(int argc, char** argv) {
    const Group group = {
        "Cantor singular inner function",
        "S_mu(z) = exp(integral (z+w)/(z-w) dmu(w)), mu = Cantor probability measure",
        "singular.exe", "singular",
        {"Cantor", "Atomic", "Blaschke"},
        {"Cantor exponent sum", "First quadrature node"},
        0, 0, 3.4, 1024,
        SINGULAR_HLSL, build, reference, value, checkCenter
    };
    return runApp(group, argc, argv);
}
