// Native Metal renderer for Klein's modular j-invariant
//     j(tau) = E4(tau)^3 / Delta(tau) = 1/q + 744 + 196884 q + ...,   q = exp(2 pi i tau),
// on the upper half-plane, the automorphic function of the modular group
// SL(2, Z): j((a tau + b) / (c tau + d)) = j(tau).
//
// Each pixel reduces tau into the standard fundamental domain
// |Re tau| <= 1/2, |tau| >= 1 with tau -> tau - n and tau -> -1/tau. There
// |q| <= exp(-pi sqrt 3) < 0.0044, so 12 terms of the q-series reach double
// precision. Delta = q prod (1 - q^n)^24 is computed as a product, which
// avoids the cancellation in (E4^3 - E6^2) / 1728 high up the half-plane.
// The Eisenstein series E4 and E6 are modular forms of weight 4 and 6, so
// their values at the reduced point are carried back with the factor
// tau^-k of each inversion.
//
// The window, renderer and controls are in ../common/app.h. There is one
// summed item (the identity matrix), as in the Kantorovich program.

#include "app.h"

#include <array>

// Sums: t0 = E4(tau), t1 = E6(tau). j itself is computed in plotValue from z.
static const char* J_MSL = R"MSL(
// j grows like exp(2 pi Im tau) towards the cusp; darkening it would hide that.
#define FADE_POLES 0
// tau -> tau - n and tau -> -1/tau until tau lies in the fundamental domain.
// f is the product of the tau that were inverted, so a modular form of weight
// k satisfies E(tau0) = E(tau) / f^k. False if 64 steps were not enough
// (only within about 1e-30 of the real axis).
bool reduceTau(thread C& tau, thread C& f) {
    f = C(1.0f, 0.0f);
    for (int i = 0; i < 64; i++) {
        tau.x = tau.x - R(rint(toF(tau.x)));
        if (toF(abs2(tau)) >= 1.0f) return true;
        f = cm(f, tau);
        tau = -crecip(tau);
    }
    return false;
}
// q = exp(2 pi i tau); 2 pi as a double-float in the FP64 frame.
C qOf(C tau) {
#if FP64
    return expi(df(6.2831855f, -1.7484555e-07f) * tau);
#else
    return expi(6.28318531f * tau);
#endif
}
// Lambert series E4 = 1 + 240 sum n^3 q^n / (1 - q^n), E6 = 1 - 504 sum n^5 q^n / (1 - q^n),
// and Delta = q prod (1 - q^n)^24, for |q| < 0.0044.
void eisenstein(C q, thread C& e4, thread C& e6, thread C& delta) {
    C one = C(1.0f, 0.0f), qn = q, s3 = C(0.0f, 0.0f), s5 = C(0.0f, 0.0f), p = one;
    for (int n = 1; n <= 12; n++) {
        C t = cdiv(qn, one - qn);
        float n3 = float(n * n * n);
        s3 = s3 + R(n3) * t;
        s5 = s5 + R(n3 * float(n * n)) * t;
        p = cm(p, one - qn);
        qn = cm(qn, q);
    }
    e4 = one + R(240.0f) * s3;
    e6 = one - R(504.0f) * s5;
    C p2 = cm(p, p), p4 = cm(p2, p2), p8 = cm(p4, p4), p16 = cm(p8, p8);
    delta = cm(q, cm(p8, p16));
}
void term(C z, C unused, uint index, thread C& t0, thread C& t1) {
    C tau = z, f, e4, e6, delta;
    if (!(toF(z.y) > 0.0f) || !reduceTau(tau, f)) {
        t0 = t1 = C(0.0f, 0.0f);
        return;
    }
    eisenstein(qOf(tau), e4, e6, delta);
    C f2 = cm(f, f), f4 = cm(f2, f2);
    t0 = cdiv(e4, f4);
    t1 = cdiv(e6, cm(f4, f2));
}
// j is invariant, so it needs only the reduced point.
C jOf(C z) {
    C tau = z, f, e4, e6, delta;
    if (!reduceTau(tau, f)) return C(0.0f, 0.0f);
    // Above Im tau = 12, |j| ~ exp(2 pi Im tau) leaves float range (the discs
    // touching the real axis at the cusps). There j = 1/q to 30 digits: keep
    // its phase, -2 pi Re tau, and cap the modulus.
    if (toF(tau.y) > 12.0f) {
        float phase = -6.28318531f * toF(tau.x);
        return C(R(1e37f * cos(phase)), R(1e37f * sin(phase)));
    }
    C q = qOf(tau);
    eisenstein(q, e4, e6, delta);
    return cdiv(cm(cm(e4, e4), e4), delta);
}
// F key: j, E4, E6
C plotValue(constant View& v, C a0, C a1, C z) {
    return v.func == 0 ? jOf(z) : (v.func == 1 ? a0 : a1);
}
// The lower half-plane is not in the domain.
#define PLOT_CUSTOM_COLOR 1
bool plotCustomColor(constant View& v, C z, thread float3& col) {
    col = float3(0.045f, 0.05f, 0.065f);
    return !(toF(z.y) > 0.0f);
}
)MSL";

// ---------------------------------------------------------------------------
// CPU side, in double: the same reduction and series for the title readout,
// and an independent direct q-series (no reduction) for the checks.

namespace modular {

constexpr double PI = 3.14159265358979323846;

inline cd qOf(cd tau) { return std::exp(cd(0, 2 * PI) * tau); }

// tau -> tau - n, tau -> -1/tau into the fundamental domain; f = product of the inverted tau.
inline bool reduce(cd& tau, cd& f) {
    f = 1;
    for (int i = 0; i < 64; i++) {
        tau -= std::nearbyint(tau.real());
        if (std::norm(tau) >= 1) return true;
        f *= tau;
        tau = -1.0 / tau;
    }
    return false;
}

// E4, E6 and Delta from q with `terms` terms of each series.
inline void series(cd q, int terms, cd& e4, cd& e6, cd& delta) {
    cd qn = q, s3 = 0, s5 = 0, p = 1;
    for (int n = 1; n <= terms; n++) {
        cd t = qn / (1.0 - qn);
        double n3 = double(n) * n * n;
        s3 += n3 * t;
        s5 += n3 * n * n * t;
        p *= 1.0 - qn;
        qn *= q;
    }
    e4 = 1.0 + 240.0 * s3;
    e6 = 1.0 - 504.0 * s5;
    delta = q * std::pow(p, 24);
}

// E4(tau), E6(tau) and j(tau) through the fundamental domain, as on the GPU.
struct Values { cd e4, e6, j; };
inline Values evaluate(cd tau) {
    cd f, e4, e6, delta;
    if (!(tau.imag() > 0) || !reduce(tau, f)) return {NAN, NAN, NAN};
    series(qOf(tau), 12, e4, e6, delta);
    return {e4 / std::pow(f, 4), e6 / std::pow(f, 6), e4 * e4 * e4 / delta};
}

// Directly from the q-series at tau itself, with enough terms for Im tau >= 0.2.
inline Values direct(cd tau) {
    cd e4, e6, delta;
    series(qOf(tau), 600, e4, e6, delta);
    return {e4, e6, e4 * e4 * e4 / delta};
}

}  // namespace modular

static std::vector<Mat> build(size_t) {
    consoleLog("modular group SL(2,Z): each pixel reduced to the fundamental domain, 12-term q-series\n");
    return {{1., 0., 0., 1.}};
}

// Independent of the GPU's method where the direct series converges.
static void reference(cd z, cd out[2]) {
    modular::Values v = z.imag() >= 0.2 ? modular::direct(z) : modular::evaluate(z);
    out[0] = v.e4;
    out[1] = v.e6;
}

static cd value(cd e4, cd e6, int func, cd z) {
    return func == 0 ? modular::evaluate(z).j : (func == 1 ? e4 : e6);
}

// I key: known values of j, invariance under SL(2,Z), and the reduction
// against the direct series.
static void invariance() {
    using namespace modular;
    auto err = [](cd a, cd b) { return std::abs(a - b) / std::max(std::abs(b), 1e-300); };
    const cd rho = std::polar(1.0, 2 * PI / 3), tau0(0.3, 0.6);
    consoleLog("j(i) = %.10g%+.3gi (1728)   j(i sqrt 2) = %.10g%+.3gi (8000)   |j(rho)| = %.2g (0)\n",
               evaluate(cd(0, 1)).j.real(), evaluate(cd(0, 1)).j.imag(),
               evaluate(cd(0, std::sqrt(2.0))).j.real(), evaluate(cd(0, std::sqrt(2.0))).j.imag(),
               std::abs(evaluate(rho).j));
    const cd j0 = evaluate(tau0).j;
    const std::array<std::array<double, 4>, 3> gammas = {{{1, 1, 0, 1}, {0, -1, 1, 0}, {2, 1, 7, 4}}};
    consoleLog("j(0.3+0.6i) = %+.10g %+.10gi   vs direct series: rel.err %.1e\n", j0.real(), j0.imag(),
               err(j0, direct(tau0).j));
    for (const auto& g : gammas) {
        cd t = (g[0] * tau0 + g[1]) / (g[2] * tau0 + g[3]);
        consoleLog("  j(%g,%g;%g,%g . tau) rel.err %.1e   E4 weight 4: rel.err %.1e\n", g[0], g[1], g[2], g[3],
                   err(evaluate(t).j, j0),
                   err(evaluate(t).e4, std::pow(g[2] * tau0 + g[3], 4) * evaluate(tau0).e4));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Modular j-invariant",
        "j(tau) = E4^3 / Delta = 1/q + 744 + 196884 q + ...,  G = SL(2,Z)",
        "invariant_j",
        "invariant_j",
        {"j", "E4", "E6"},
        {"E4", "E6"},
        0.0, 0.9, 2.0,  // -0.1 < Im tau < 1.9 on the shorter side
        0,              // one summed item
        J_MSL,
        build, reference, value, invariance,
    };
    return runApp(group, argc, argv);
}
