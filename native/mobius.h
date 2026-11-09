#pragma once
// SL(2,C) matrices acting on the Riemann sphere as Möbius transformations.
// Shared by the Kleinian group (poincare.cpp) and the icosahedral group
// (icosahedral.h).

#include <complex>

using cd = std::complex<double>;
struct Mat { cd a, b, c, d; };

inline Mat mul(const Mat& p, const Mat& q) {
    return {p.a * q.a + p.b * q.c, p.a * q.b + p.b * q.d, p.c * q.a + p.d * q.c, p.c * q.b + p.d * q.d};
}
inline Mat inv(const Mat& m) { return {m.d, -m.b, -m.c, m.a}; }  // det = 1
inline double norm2(const Mat& m) { return std::norm(m.a) + std::norm(m.b) + std::norm(m.c) + std::norm(m.d); }
inline cd mobius(const Mat& m, cd z) { return (m.a * z + m.b) / (m.c * z + m.d); }
