#include "app.h"

// a_n = (1 - 0.8/(n+1)^2) exp(i n goldenAngle), n >= 1.
// Sum (1-|a_n|) is finite. The normalized factors converge to 1.
static constexpr double goldenAngle = 2.39996322972865332;
static std::vector<cd> zeros;

static const char* BLASCHKE_HLSL = R"HLSL(
#define FADE_POLES 0
#ifdef COMPUTE
// This renderer evaluates the product directly in the pixel shader.
VEC4 term(VEC2 w1, VEC2 w2, uint index) {
    return (VEC4)0;
}
#else
#define PLOT_CUSTOM_COLOR
bool plotCustomColor(double2 z, out float3 color) {
    color = float3(0.045, 0.05, 0.065);
    return z.x*z.x + z.y*z.y >= 1.0L;
}
double2 plotValue(double4 unused, double2 z) {
    float2 product = float2(1,0);
    float2 coordinate = (float2)z;
    for (uint n = 1; n <= uCount; ++n) {
        float radius = 1.0 - 0.8 / ((float)(n+1)*(float)(n+1));
        float angle = (float)n * 2.39996322972865332;
        float2 direction = normalize(float2(cos(angle), sin(angle)));
        float2 a = radius * direction;
        float2 numerator = a-coordinate;
        float2 denominator = float2(1,0)-float2(a.x*coordinate.x+a.y*coordinate.y, a.x*coordinate.y-a.y*coordinate.x);
        float2 factor = float2(dot(numerator,denominator),
            numerator.y*denominator.x-numerator.x*denominator.y)/dot(denominator,denominator);
        factor = float2(direction.x*factor.x+direction.y*factor.y,
                        direction.x*factor.y-direction.y*factor.x);
        product = float2(product.x*factor.x-product.y*factor.y,
                         product.x*factor.y+product.y*factor.x);
    }
    return product;
}
#endif
)HLSL";

static std::vector<Mat> build(size_t count) {
    zeros.clear();
    for (size_t n = 1; n <= count; ++n) {
        float angle = float(double(n) * goldenAngle);
        cd direction(std::cos(angle), std::sin(angle));
        direction /= std::abs(direction);
        zeros.push_back((1.0 - 0.8 / (double(n+1)*double(n+1))) * direction);
    }
    std::printf("Blaschke product: %zu factors, sum(1-|a_n|) < 0.516\n", count);
    return std::vector<Mat>(count, Mat{1,0,0,1});
}

static void reference(cd, cd out[2]) { out[0] = out[1] = 0; }

static cd value(cd, cd, int, cd z) {
    cd product = 1;
    for (cd a : zeros)
        product *= (std::abs(a)/a) * (a-z)/(1.0-std::conj(a)*z);
    return product;
}

static void check() {
    double worst = 0;
    for (int k = 0; k < 32; ++k) {
        cd z = std::polar(1.0, 6.283185307179586 * (k+0.5)/32);
        worst = std::max(worst, std::abs(std::abs(value(0,0,0,z))-1));
    }
    std::printf("CPU boundary modulus error %.3e; B_N(0) = %.12g\n",
                worst, value(0,0,0,0).real());
    if (worst > 1e-10) fatal("Blaschke boundary check failed");
}

int main(int argc, char** argv) {
    const Group group = {
        "Infinite Blaschke product",
        "B(z)=prod_n (|a_n|/a_n)(a_n-z)/(1-conj(a_n)z), a_n=(1-0.8/(n+1)^2)exp(i*n*goldenAngle)",
        "blaschke.exe", "blaschke",
        {"B_N", "B_N", "B_N"}, {"Unused", "Unused"},
        0,0,2.2,256, BLASCHKE_HLSL, build, reference, value, check
    };
    return runApp(group, argc, argv);
}

