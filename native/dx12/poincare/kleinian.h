#pragma once
// Kleinian group G = <g1, g2>, treated as free on its two generators, and its
// truncated Poincaré series
//     P_{k,l}(z) = sum_{g in G} (az+b)^k (cz+d)^l.
// See docs/background.md.

#include <queue>
#include <vector>

#include "mobius.h"

namespace kleinian {

inline const cd B(-0.10997580175045825, -0.0011537600196144);
inline const cd Cc(0.07184238006772357, -1.8193355782014326);
inline const Mat G1 = {cd(1, 0.1), B, Cc, cd(1, 0.1)};
inline const Mat G2 = {cd(1.1, -1), cd(1.1, 0), cd(1.1, 0), cd(1.1, 1)};
inline const Mat GENS[4] = {G1, inv(G1), G2, inv(G2)};  // g1, g1^-1, g2, g2^-1
inline const char* const GEN_NAMES[4] = {"g1", "g1^-1", "g2", "g2^-1"};
inline const int INVERSE[4] = {1, 0, 3, 2};

// The N reduced words of smallest Frobenius norm, smallest first. G is free
// on g1, g2, so each reduced word is a distinct group element.
inline std::vector<Mat> enumerate(size_t N) {
    struct Node {
        double n2;
        Mat m;
        int last;
        bool operator<(const Node& o) const { return n2 > o.n2; }  // min-heap
    };
    std::priority_queue<Node> pq;
    pq.push({2.0, {1, 0, 0, 1}, -1});
    std::vector<Mat> out;
    out.reserve(N);
    while (out.size() < N && !pq.empty()) {
        Node e = pq.top();
        pq.pop();
        out.push_back(e.m);
        for (int k = 0; k < 4; k++) {
            if (e.last >= 0 && k == INVERSE[e.last]) continue;
            Mat m = mul(GENS[k], e.m);
            pq.push({norm2(m), m, k});
        }
    }
    return out;
}

// P1 = P_{-1,-3}, P2 = P_{0,-4}, f = P1 / P2 (the automorphic function).
struct Series { cd P1, P2, f; };
inline Series evalSeries(const std::vector<Mat>& mats, cd z) {
    cd P1 = 0, P2 = 0;
    for (const Mat& m : mats) {
        cd w1 = m.a * z + m.b, w2 = m.c * z + m.d;
        cd r = 1.0 / (w1 * w2 * w2 * w2 * w2);
        P1 += w2 * r;  // 1 / ((az+b)(cz+d)^3)
        P2 += w1 * r;  // 1 / (cz+d)^4
    }
    return {P1, P2, P1 / P2};
}

}  // namespace kleinian
