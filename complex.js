"use strict";
// Complex arithmetic (double precision, CPU side) and SL(2,C) matrices acting
// on the Riemann sphere as Möbius transformations. Shared by the Kleinian group
// (index.html) and the icosahedral group (icosahedral.js). Mirrors native/mobius.h.

const C = (re, im = 0) => [re, im];
const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const div = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
const abs2 = a => a[0] * a[0] + a[1] * a[1];
const neg = a => [-a[0], -a[1]];

// Matrix [a,b,c,d] of complex numbers
const mmul = (P, Q) => [
  add(mul(P[0], Q[0]), mul(P[1], Q[2])), add(mul(P[0], Q[1]), mul(P[1], Q[3])),
  add(mul(P[2], Q[0]), mul(P[3], Q[2])), add(mul(P[2], Q[1]), mul(P[3], Q[3])),
];
const minv = M => [M[3], neg(M[1]), neg(M[2]), M[0]]; // det = 1
const mnorm2 = M => abs2(M[0]) + abs2(M[1]) + abs2(M[2]) + abs2(M[3]);
const mobius = (M, z) => div(add(mul(M[0], z), M[1]), add(mul(M[2], z), M[3]));
