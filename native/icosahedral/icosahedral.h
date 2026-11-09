#pragma once
// Icosahedral automorphic function on the Riemann sphere
//     F(z) = sum_{g in G} h(g(z)),   h(z) = (z - 0.05)^-5,
// where G is the group of 60 rotations of the icosahedron: the rotation
// subgroup of the Coxeter group H3 = [5,3], i.e. the von Dyck group D(2,3,5).
// Because G is finite, the sum is exact. See docs/coxeter.md.
//
// This header holds the CPU side: group, F, h and Klein's invariant J.
// The GPU side is in icosahedral.cpp.

#include <algorithm>
#include <cmath>
#include <vector>

#include "mobius.h"

namespace icosahedral {

constexpr double POLE = 0.05;  // h has a pole of order 5 here
constexpr size_t ORDER = 60;

inline cd eps(int k) { return std::polar(1.0, 2 * 3.14159265358979323846 * k / 5); }

// Klein's generators in SU(2): S(z) = eps*z is the 72-degree turn about the
// axis through the vertices 0 and infinity; T is a half-turn (T^2 = -I).
inline const Mat S = {eps(3), 0, 0, eps(2)};
inline const Mat T = [] {
    double r = 1 / std::sqrt(5.0);
    cd p = (eps(1) - eps(4)) * r, q = (eps(2) - eps(3)) * r;
    return Mat{-p, q, q, p};
}();
inline const Mat GENS[4] = {S, inv(S), T, inv(T)};
inline const char* const GEN_NAMES[4] = {"S", "S^-1", "T", "T^-1"};

// Closure of {S, T}; matrices differing by -1 are the same rotation.
// Returns ORDER elements; the caller checks the count.
inline std::vector<Mat> group() {
    auto same = [](const Mat& p, const Mat& q) {
        double plus = norm2({p.a - q.a, p.b - q.b, p.c - q.c, p.d - q.d});
        double minus = norm2({p.a + q.a, p.b + q.b, p.c + q.c, p.d + q.d});
        return std::min(plus, minus) < 1e-12;
    };
    std::vector<Mat> out = {{1, 0, 0, 1}};
    for (size_t i = 0; i < out.size(); i++) {
        for (const Mat& g : {S, T}) {
            Mat m = mul(g, out[i]);
            if (std::none_of(out.begin(), out.end(), [&](const Mat& o) { return same(o, m); })) out.push_back(m);
        }
    }
    return out;
}

// Matrices for the GPU: each g pre-multiplied by [[1, -POLE], [0, 1]], so with
// w1 = a z + b, w2 = c z + d of the shifted matrix, h(g(z)) = (w2 / w1)^5.
inline std::vector<Mat> gpuMatrices(const std::vector<Mat>& g) {
    const Mat shift = {1, -POLE, 0, 1};
    std::vector<Mat> out;
    out.reserve(g.size());
    for (const Mat& m : g) out.push_back(mul(shift, m));
    return out;
}

inline cd h(cd w) { cd r = 1.0 / (w - POLE); cd r2 = r * r; return r2 * r2 * r; }

inline cd F(const std::vector<Mat>& g, cd z) {
    cd sum = 0;
    for (const Mat& m : g) sum += h(mobius(m, z));
    return sum;
}

// Klein's invariant J = H^3 / (1728 f^5) with
//   f = z (z^10 + 11 z^5 - 1),   H = -(z^20 + 1) + 228 (z^15 - z^5) - 494 z^10.
// z -> -1/z is in G, so evaluate inside |z| <= 1 to keep the powers small.
inline cd kleinJ(cd z) {
    if (std::norm(z) > 1) z = -1.0 / z;
    cd z5 = z * z * z * z * z, z10 = z5 * z5, z15 = z10 * z5, z20 = z10 * z10;
    cd f = z * (z10 + 11.0 * z5 - 1.0);
    cd H = -(z20 + 1.0) + 228.0 * (z15 - z5) - 494.0 * z10;
    return H * H * H / (1728.0 * f * f * f * f * f);
}

}  // namespace icosahedral
