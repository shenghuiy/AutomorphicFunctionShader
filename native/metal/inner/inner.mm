// Native Metal port of ../../dx12/inner/inner.cpp; see docs/inner.md.
#include "app.h"

// The Cantor measure is approximated by equal masses at the midpoints of
// the surviving intervals of a ternary Cantor construction. Each node has
// weight 1/N. As N=2^depth grows, these discrete measures converge weakly
// to a non-atomic Cantor measure. Three scaled copies cover the circle,
// followed by a smooth angular distortion to adjust cluster spacing.
static unsigned cantorDepth = 10;
static unsigned atomCount = 1u << cantorDepth;
static constexpr double measureMass = 1.0;
static constexpr unsigned sectors = 3;
static constexpr double rotation = 0.12;
static constexpr double tau = 6.283185307179586;

static const char* MSL = R"MSL(
#define FADE_POLES 0
void term(C delta, C numerator, uint index, thread C& t0, thread C& t1) {
    t0 = t1 = C(0.0f, 0.0f);
    if (toF(delta.x) == 0.0f && toF(delta.y) == 0.0f) return;
    t0 = cm(numerator, crecip(delta));
    if (index == 0) t1 = t0;
}
C plotValue(constant View& v, C a0, C a1, C z) {
    if (v.func == 2) {
        // Split constants retain the CPU alpha in double-float mode.
        C alpha = toC(cdf(df(0.35f, 5.9604646e-09f), df(0.2f, -2.9802323e-09f)));
        return cdiv(z - alpha, C(1.0f, 0.0f) - cm(conjugate(alpha), z));
    }
    C exponent = v.func == 0 ? a0 / R(float(v.count))
                            : cdiv(z + C(1.0f, 0.0f), z - C(1.0f, 0.0f));
    // Match the DX12 display clamp and reduce phase before sin/cos.
    float phase = fract(toF(exponent.y / R(6.283185307179586f))) * 6.283185307179586f;
    float magnitude = exp(clamp(toF(exponent.x), -70.0f, 70.0f));
    return C(R(magnitude * cos(phase)), R(magnitude * sin(phase)));
}
)MSL";

static std::vector<cd> support;

static std::vector<Mat> build(size_t requested) {
    cantorDepth = 0;
    atomCount = 1;
    while (atomCount <= requested/(2*sectors) && cantorDepth < 20) {
        atomCount *= 2;
        ++cantorDepth;
    }
    consoleLog("Cantor quadrature: depth %u, %u equal weights, total mass 1\n",
                cantorDepth, atomCount*sectors);
    support.clear();
    std::vector<Mat> matrices;
    for (unsigned sector = 0; sector < sectors; ++sector) {
        for (unsigned index = 0; index < atomCount; ++index) {
            double left = 0;
            double width = 1;
            for (unsigned digit = 0; digit < cantorDepth; ++digit) {
                width /= 3;
                if ((index >> digit) & 1u) left += 2 * width;
            }
            double angle = rotation + tau * (sector + left + width/2) / sectors;
            angle += 0.10 * std::sin(2*angle + 0.4);
            const cd w = std::polar(1.0, angle);
            support.push_back(w);
            matrices.push_back({1, -w, 1, w});
        }
    }
    return matrices;
}

static void reference(cd z, cd out[2]) {
    out[0] = out[1] = 0;

    for (cd w : support) out[0] += (z + w) / (z - w);
    out[1] = (z + support.front()) / (z - support.front());
}

static cd value(cd cantor, cd, int function, cd z) {

    const cd alpha(0.35, 0.2);
    if (function == 2) return (z - alpha) / (1.0 - std::conj(alpha) * z);
    return std::exp(function == 0 ? measureMass * cantor / double(support.size()) : (z+1.0)/(z-1.0));
}

static void checkCenter() {
    cd sums[2];
    reference(0, sums);
    consoleLog("At z=0, Cantor modulus = %.9f, atomic modulus = %.9f (Cantor expected exp(-1), atomic exp(-1))\n",
                std::abs(value(sums[0], sums[1], 0, 0)),
                std::abs(value(sums[0], sums[1], 1, 0)));
}

int main(int argc, char** argv) {
    Group group = {
        "Inner function",
        "S_mu(z) = exp(integral (z+w)/(z-w) dmu(w)), mu = Cantor probability measure",
        "inner", "inner",
        {"Cantor", "Atomic", "Blaschke"},
        {"Cantor exponent sum", "First quadrature node"},
        0, -0.04, 2.25, 3072,
        MSL, build, reference, value, checkCenter,
        false, 0
    };
    group.defaultPalette = 1;
    return runApp(group, argc, argv);
}
