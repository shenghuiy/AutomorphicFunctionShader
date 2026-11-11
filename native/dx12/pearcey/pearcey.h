#pragma once
// Pearcey integral
//     Psi(x, y) = integral over the real line of exp(i (t^4 + y t^2 + x t)) dt,
// plotted as a function of z = x + i y.
//
// The integrand does not decay on the real line, so the path is moved into the
// complex plane. With phi(t) = t^4 + y t^2 + x t, the path is
//     t(s) = s + i delta(s),   delta(s) = LAMBDA g / sqrt(1 + g^2),   g = phi'(s) / kappa,
// for real s. To first order Im phi(t(s)) = delta(s) phi'(s) >= 0, so
// |exp(i phi)| <= 1 and it decays like exp(-LAMBDA |phi'(s)|) away from the
// stationary points of phi. Exactly, Im phi = delta (phi' - 4 s delta^2)
// = (LAMBDA g^2 / q^3) (kappa (1 + g^2) - 4 LAMBDA^2 s g) with q = sqrt(1 + g^2),
// which is non-negative for all g iff kappa >= 2 LAMBDA^2 |s|, so on |s| <= S
// once kappa >= 2 S LAMBDA^2. kappa is 8 times that: a gentler bend keeps the
// square root's branch points away from the real axis, and the midpoint rule
// then needs about half as many nodes as with kappa = 4 S LAMBDA^2. Nothing
// large is cancelled, so single precision is already usable.
//
// The integrand is analytic in s, so the midpoint rule on [-S, S] converges
// exponentially in the number of nodes. S depends on (x, y): the real roots of
// phi' lie in |s| <= max(sqrt|y|, cbrt(|x|/2)), and MARGIN beyond them the
// integrand is below double precision.
//
// This header holds the CPU side. The GPU side (pearcey.cpp) uses the same path.

#include <cmath>
#include <complex>
#include <vector>

#include "mobius.h"

namespace pearcey {

constexpr double LAMBDA = 0.5;  // height of the path above / below the real axis
constexpr double MARGIN = 3.0;

inline double halfWidth(double x, double y) {
    // float, like the shader; any value near this one gives the same integral
    return double(float(MARGIN) + std::max(std::sqrt(std::fabs(float(y))), std::pow(0.5f * std::fabs(float(x)), 1.0f / 3.0f)));
}

// Midpoint nodes u_j in (-1, 1) and weights, uploaded as the d entry of a
// matrix [[1, 0], [0, u + i w]], so the shader sees w1 = z and w2 = (u, w).
//
// n must be a power of two. The nodes are in bit-reversed order, so the first
// m = 2^k of them are every (n/m)-th node: an evenly spaced rule of its own.
// A partial sum over m nodes is then m/n times a coarser but valid estimate of
// the integral, and a frame drawn before the sum is complete shows the right picture.
inline std::vector<Mat> nodes(size_t n) {
    std::vector<Mat> out;
    out.reserve(n);
    for (size_t i = 0; i < n; i++) {
        size_t j = 0;
        for (size_t bit = 1, rev = n >> 1; bit < n; bit <<= 1, rev >>= 1)
            if (i & bit) j |= rev;
        out.push_back({1, 0, 0, cd(-1 + (2.0 * j + 1) / n, 2.0 / n)});
    }
    return out;
}

// out[0] = Psi(x, y), out[1] = dPsi/dx = integral of i t exp(i phi) dt, with n nodes.
inline void psi(double x, double y, size_t n, cd out[2]) {
    const cd I(0, 1);
    double S = halfWidth(x, y), kappa = 16 * S * LAMBDA * LAMBDA;
    out[0] = out[1] = 0;
    for (size_t j = 0; j < n; j++) {
        double s = S * (-1 + (2.0 * j + 1) / n);
        double g = (4 * s * s * s + 2 * y * s + x) / kappa, dg = (12 * s * s + 2 * y) / kappa;
        double q = std::sqrt(1 + g * g);
        cd t(s, LAMBDA * g / q), dt(1, LAMBDA * dg / (q * q * q));
        cd t2 = t * t;
        cd v = std::exp(I * (t2 * t2 + y * t2 + x * t)) * dt * (2 * S / n);
        out[0] += v;
        out[1] += I * t * v;
    }
}

// Independent check on the straight path t = e^{i pi/8} s, where exp(i t^4) = exp(-s^4).
// Only for small |x|, |y|: for y << 0 the integrand grows to exp(y^2 / 8) before it decays.
inline cd psiRotated(double x, double y) {
    const cd I(0, 1), w = std::polar(1.0, 3.14159265358979323846 / 8);
    const int n = 4000;
    const double S = 8;
    cd sum = 0;
    for (int j = 0; j < n; j++) {
        cd t = w * (S * (-1 + (2.0 * j + 1) / n));
        cd t2 = t * t;
        sum += std::exp(I * (t2 * t2 + y * t2 + x * t));
    }
    return sum * w * (2 * S / n);
}

}  // namespace pearcey
