// Harmonic Newton example from the supplied Kantorovich theorem illustration.
// Native Metal version of ../../dx12/kantorovich/kantorovich.cpp; see docs/kantorovich.md.
#include "app.h"

static cd evaluate(cd z, cd* derivative = nullptr) {
    cd d = z * z * z - 0.216;
    if (derivative) *derivative = -z * (z * z * z + 0.432) / (d * d);
    return z - std::conj(z * z / d);
}
static void reference(cd z, cd out[2]) { out[0] = evaluate(z, &out[1]); }
static cd value(cd f, cd g, int mode, cd z) {
    if (mode == 0) return f;
    if (mode == 1) return -(f + std::conj(g) * std::conj(f)) / (1.0 - std::norm(g));
    for (int i = 0; i < 80; i++) {
        f = evaluate(z, &g);
        if (std::abs(f) < 1e-6) return z;
        double j = 1.0 - std::norm(g);
        if (!std::isfinite(j) || std::abs(j) < 1e-10) break;
        z -= (f + std::conj(g) * std::conj(f)) / j;
    }
    return {NAN, NAN};
}
static void check() {
    double worst = 0; int positive = 0, negative = 0;
    // Independent real finite differences validate the Wirtinger derivative.
    for (cd z : {cd(.3, .4), cd(-.8, .2), cd(1.2, -.7)}) {
        cd g; evaluate(z, &g); double e = 1e-6;
        cd fx = (evaluate(z + e) - evaluate(z - e)) / (2 * e);
        cd fy = (evaluate(z + cd(0, e)) - evaluate(z - cd(0, e))) / (2 * e);
        worst = std::max(worst, std::abs(fx - (cd(1) - std::conj(g))));
        worst = std::max(worst, std::abs(fy - cd(0, 1) * (cd(1) + std::conj(g))));
    }
    std::vector<cd> roots;
    for (int x = -15; x <= 15; x++) for (int y = -15; y <= 15; y++) {
        cd z(x * .1, y * .1), f, g; f = evaluate(z, &g);
        cd r = value(f, g, 2, z);
        if (!std::isfinite(r.real()) || std::abs(evaluate(r)) > 1e-6) continue;
        bool seen = false; for (cd old : roots) if (std::abs(old - r) < 1e-4) seen = true;
        if (!seen) roots.push_back(r);
    }
    for (cd r : roots) { cd g; evaluate(r, &g); if (std::norm(g) < 1) positive++; else negative++; }
    std::printf("Derivative finite-difference error %.3g; roots %zu (%d preserving, %d reversing)\n", worst, roots.size(), positive, negative);
    if (worst > 1e-7 || roots.size() != 10 || positive != 4 || negative != 6) fatal("Kantorovich mathematical check failed");
}
static std::vector<Mat> build(size_t) { check(); return {{1., 0., 0., 1.}}; }

// The single uploaded matrix is the identity, so w1 = z.
// Sums: t0 = f(z),  t1 = g(z) = -z (z^3 + 0.432) / (z^3 - 0.216)^2, the derivative the Newton step uses
static const char* MSL = R"MSL(
#define FADE_POLES 0
void term(C z, C unused, uint index, thread C& t0, thread C& t1) {
    C z2 = cm(z, z), z3 = cm(z2, z), d = z3 - C(R(27.0f) / R(125.0f), R(0.0f));
    C r = cm(z2, crecip(d));
    t0 = z - conjugate(r);
    t1 = -cm(cm(z, z3 + C(R(54.0f) / R(125.0f), R(0.0f))), crecip(cm(d, d)));
}
bool evaluate(C z, thread C& f, thread C& g, thread R& j) {
    C z2 = cm(z, z), z3 = cm(z2, z), d = z3 - C(R(27.0f) / R(125.0f), R(0.0f));
    f = g = C(0.0f, 0.0f); j = R(0.0f);
    if (toF(abs2(d)) < 1e-20f || toF(abs2(z)) > 1e16f) return false;
    f = z - conjugate(cdiv(z2, d));
    g = -cdiv(cm(z, z3 + C(R(54.0f) / R(125.0f), R(0.0f))), cm(d, d));
    j = R(1.0f) - abs2(g);
    return all(isfinite(toF2(f))) && all(isfinite(toF2(g)));
}
C plotValue(constant View& v, C a0, C a1, C z) {
    if (v.func == 0) return a0;
    R j = R(1.0f) - abs2(a1);
    if (abs(toF(j)) < 1e-10f) return C(0.0f, 0.0f);
    return -(a0 + cm(conjugate(a1), conjugate(a0))) / j;
}
#define PLOT_OVERLAY 1
float3 plotOverlay(constant View& v, float3 col, C z) {
    C f, g; R j;
    bool ok = evaluate(z, f, g, j);
    float jf = toF(j);
    // fwidth outside the branch, so every lane of the quad takes part.
    float contour = 1.0f - smoothstep(0.0f, 1.5f, abs(jf) / max(fwidth(jf), 1e-6f));
    if (!ok) return float3(0.0f);
    // The orientation shading refers to f in both phase portraits.
    return mix(col * (jf > 0.0f ? 1.0f : 0.48f), float3(0.04f), contour * 0.8f);
}
#define PLOT_CUSTOM_COLOR 1
bool plotCustomColor(constant View& v, C z, thread float3& col) {
    col = float3(0.025f, 0.03f, 0.045f);
    if (v.func != 2) return false;
    for (int k = 0; k < 80; k++) {
        C f, g; R j;
        if (!evaluate(z, f, g, j)) return true;
        if (toF(abs2(f)) < 1e-12f) {
            float2 zf = toF2(z);
            float r = length(zf);
            int id = 0;
            if (r > 0.05f) {
                float theta = atan2(zf.y, zf.x);
                int sector = int(floor((theta - (r < 1.0f ? 1.0471975512f : 0.0f)) / 2.0943951024f + 0.5f));
                sector = (sector % 3 + 3) % 3;
                id = 1 + sector + (r < 0.5f ? 0 : (r < 1.0f ? 3 : 6));
            }
            float h = fract(float(id) * 0.61803398875f);
            float3 rgb = saturate(abs(fract(h + float3(0.0f, 2.0f / 3.0f, 1.0f / 3.0f)) * 6.0f - 3.0f) - 1.0f);
            col = mix(float3(1.0f), rgb, 0.78f) * (0.45f + 0.55f * exp(-0.025f * float(k)));
            return true;
        }
        if (abs(toF(j)) < 1e-10f) return true;
        z = z - (f + cm(conjugate(g), conjugate(f))) / j;
    }
    return true;
}
)MSL";

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") { check(); return 0; }
    const Group group = {"Kantorovich: harmonic Newton iteration",
        "f(z) = z - conjugate(z^2 / (z^3 - 0.6^3))",
        "kantorovich", "kantorovich", {"f", "N(z)-z", "Newton basins"}, {"f", "g prime"},
        0, 0, 3, 0, MSL, build, reference, value, check};
    return runApp(group, argc, argv);
}
