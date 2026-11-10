// Native DirectX 12 renderer for the Pearcey integral
//     Psi(x, y) = integral over the real line of exp(i (t^4 + y t^2 + x t)) dt
// as a function of z = x + i y. Each pixel sums a midpoint rule on a path bent
// into the complex plane so that the integrand decays (pearcey.h).
//
// The window, renderer and controls are in ../common/app.h. There the summed
// items are called group elements; here they are quadrature nodes.

#include "app.h"
#include "pearcey.h"

// Each node is uploaded as [[1, 0], [0, u + i w]] (pearcey::nodes), so w1 = z = (x, y)
// and w2 = (u, w): position in (-1, 1) and weight.
// Accumulator: .xy = Psi,  .zw = dPsi/dx
static const char* PEARCEY_HLSL = R"HLSL(
#define FADE_POLES 0
#ifdef COMPUTE
// Small integer n as a runtime value. On at least one NVIDIA driver double
// immediates read as zero in some instructions (see pixelZ in shaders.h), which
// here made the fp64 sums NaN, so the fp64 code has no double constants at all.
uint num(uint n) { return n*min(uW,1u); }
#ifdef FP64
#define REAL double
// A double constant as a sum of three floats, each scaled by the runtime 1.
double split(float a, float b, float c) {
    float one=(float)num(1);
    return (double)(a*one)+(double)(b*one)+(double)(c*one);
}
// Shader model 5 has no double sqrt, exp, sin or cos.
double dsqrt(double d) {  // float estimate, two Newton steps
    double r=(double)sqrt((float)d);
    r=(r*r+d)/(r+r);
    r=(r*r+d)/(r+r);
    return r;
}
double dexp(double x) {  // x = k ln 2 + r, |r| <= 0.35, Taylor series in r
    float k=round(max((float)x,-80.0)*1.44269504);
    double r=x-(double)k*split(0.6931472,-1.9046542e-09,-1.110223e-16);
    double one=(double)num(1), p=one;
    for (uint n=13; n>=1; --n) p=one+r*p/(double)num(n);
    return (float)x < -80.0 ? one-one : p*(double)exp2(k);
}
// exp(i p) for complex p
double2 expi(double2 p) {
    // p.x = k pi/2 + r, |r| <= pi/4, Taylor series in r
    float k=round((float)p.x*0.636619772);
    double r=p.x-(double)k*split(1.5707964,-4.371139e-08,-1.7763568e-15), r2=r*r;
    double one=(double)num(1), c=one, s=one;
    for (uint n=16; n>=2; n-=2) {
        c=one-r2*c/(double)num(n*(n-1));
        s=one-r2*s/(double)num(n*(n+1));
    }
    s*=r;
    int q=((int)k)&3;
    double2 e=q==0 ? double2(c,s) : (q==1 ? double2(-s,c) : (q==2 ? double2(-c,-s) : double2(s,-c)));
    return e*dexp(-p.y);
}
#else
#define REAL float
float dsqrt(float d) { return sqrt(d); }
float2 expi(float2 p) { float s, c; sincos(p.x,s,c); return float2(c,s)*exp(-p.y); }
#endif

// Path t(s) = s + i LAMBDA g / sqrt(1 + g^2), g = phi'(s) / kappa; see pearcey.h.
VEC4 term(VEC2 z, VEC2 node, uint index) {
    REAL x=z.x, y=z.y, one=(REAL)num(1);
    REAL S=(REAL)(3.0+max(sqrt(abs((float)y)),pow(0.5*abs((float)x),1.0/3.0)));  // pearcey::halfWidth
    REAL s=S*node.x, s2=s*s, y2=y+y;
    REAL kappa=(REAL)num(4)*S;  // 16 S LAMBDA^2 with LAMBDA = 1/2
    REAL g=((REAL)num(4)*s2*s+y2*s+x)/kappa, dg=((REAL)num(12)*s2+y2)/kappa;
    REAL q=dsqrt(one+g*g), q3=q*q*q;
    VEC2 t=VEC2(s,g/(q+q)), dt=VEC2(one,dg/(q3+q3));
    VEC2 t2=cm(t,t);
    VEC2 v=cm(expi(cm(t2,t2)+y*t2+x*t),dt)*(S*node.y);
    return VEC4(v,cm(VEC2(-t.y,t.x),v));
}
#else
// F key: Psi, Psi with the caustic drawn, dPsi/dx
double2 plotValue(double4 a, double2 z) {
    return uFunc==2 ? a.zw : a.xy;
}
// The caustic 27 x^2 + 8 y^3 = 0: the (x, y) where two real stationary points
// of the phase merge. Drawn where |c| / |grad c|, the distance to the curve to
// first order, is under about 1.5 pixels.
#define PLOT_OVERLAY 1
float3 plotOverlay(float3 col, double2 zd) {
    if (uFunc!=1) return col;
    float2 z=(float2)zd;
    float c=27.0*z.x*z.x+8.0*z.y*z.y*z.y;
    float d=abs(c)/max(length(float2(54.0*z.x,24.0*z.y*z.y)),1e-30);
    return col*smoothstep(1.0,2.0,d/(float)uScale);
}
#endif
)HLSL";

static size_t nodeCount = 0;

static std::vector<Mat> build(size_t n) {
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
        pearcey::psi(-xy[0], xy[1], nodeCount, q);
        std::printf("   Psi(-x,y) = Psi(x,y): rel.err %.1e\n", err(q[0], p[0]));
    }
}

int main(int argc, char** argv) {
    const Group group = {
        "Pearcey integral",
        "Psi(x, y) = int exp(i (t^4 + y t^2 + x t)) dt,  z = x + i y",
        "pearcey.exe",
        "pearcey",
        {"Psi", "Psi + caustic", "dPsi/dx"},
        {"Psi", "dPsi/dx"},
        0.0, -5.0, 20.0,  // |x| < 10, -15 < y < 5 on the shorter side
        512,              // quadrature nodes; --terms and +/-
        PEARCEY_HLSL,
        build, reference, value, invariance,
        false,  // no sphere view
        1,      // start with the GreenPinkTones phase colours
    };
    return runApp(group, argc, argv);
}
