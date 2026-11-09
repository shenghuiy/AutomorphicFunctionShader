"use strict";
// Icosahedral automorphic function on the Riemann sphere
//     F(z) = Σ_{g∈G} h(g(z)),   h(z) = (z − 0.05)^−5,
// where G is the group of 60 rotations of the icosahedron: the rotation
// subgroup of the Coxeter group H3 = [5,3], i.e. the von Dyck group D(2,3,5).
// Because G is finite, the sum is exact. See docs/coxeter.md.
//
// CPU side (group, F, h, Klein's invariant J), the GLSL snippets index.html
// splices into its shaders, and the panel text. Mirrors native/icosahedral.h.
// Requires complex.js.

const Icosahedral = (() => {
  const POLE = 0.05; // h has a pole of order 5 here
  const ORDER = 60;
  const eps = k => [Math.cos(2 * Math.PI * k / 5), Math.sin(2 * Math.PI * k / 5)];

  // Klein's generators in SU(2): S(z) = εz is the 72° turn about the axis
  // through the vertices 0 and ∞; T is a half-turn (T² = −I).
  const S = [eps(3), C(0), C(0), eps(2)];
  const T = (() => {
    const r5 = 1 / Math.sqrt(5);
    const p = add(eps(1), neg(eps(4))), q = add(eps(2), neg(eps(3)));
    const sc = a => [a[0] * r5, a[1] * r5];
    return [sc(neg(p)), sc(q), sc(q), sc(p)];
  })();
  const GENS = [S, minv(S), T, minv(T)];
  const GEN_NAMES = ['S', 'S⁻¹', 'T', 'T⁻¹'];

  // Close {S, T} under multiplication; matrices that differ by −1 are the same
  // rotation. Returns ORDER elements; the caller checks the count.
  function group() {
    const key = M => {
      const lead = M.find(x => abs2(x) > 1e-18), s = Math.sqrt(abs2(lead));
      const u = [lead[0] / s, -lead[1] / s]; // conj(lead)/|lead| removes the overall sign
      return M.map(x => mul(x, u)).flat().map(x => Math.round(x * 1e7) || 0).join(',');
    };
    const I = [C(1), C(0), C(0), C(1)];
    const seen = new Set([key(I)]), out = [I];
    for (let i = 0; i < out.length; i++) {
      for (const G of [S, T]) {
        const M = mmul(G, out[i]), k = key(M);
        if (!seen.has(k)) { seen.add(k); out.push(M); }
      }
    }
    return out;
  }

  // Matrices for the GPU: each g pre-multiplied by [[1, −POLE], [0, 1]], so with
  // w1 = az+b, w2 = cz+d of the shifted matrix, h(g(z)) = (w2/w1)⁵.
  const gpuMatrices = g => { const shift = [C(1), C(-POLE), C(0), C(1)]; return g.map(M => mmul(shift, M)); };

  const h = w => { const r = div(C(1), add(w, C(-POLE))); const r2 = mul(r, r); return mul(mul(r2, r2), r); };
  function F(g, z) { let s = [0, 0]; for (const M of g) s = add(s, h(mobius(M, z))); return s; }

  // Klein's invariant J = H³ / (1728 f⁵), f = z(z¹⁰+11z⁵−1), H = −(z²⁰+1)+228(z¹⁵−z⁵)−494z¹⁰.
  // z ↦ −1/z is in the group, so evaluate inside the unit disk to keep powers small.
  function kleinJ(z) {
    if (abs2(z) > 1) z = neg(div(C(1), z));
    const z5 = mul(mul(mul(z, z), mul(z, z)), z), z10 = mul(z5, z5), z15 = mul(z10, z5), z20 = mul(z10, z10);
    const f = mul(z, add(add(z10, [11 * z5[0], 11 * z5[1]]), C(-1)));
    const H = add(add(neg(add(z20, C(1))), mul(C(228), add(z15, neg(z5)))), mul(C(-494), z10));
    const f2 = mul(f, f);
    return div(mul(mul(H, H), H), mul(C(1728), mul(mul(f2, f2), f)));
  }

  // GLSL for the accumulation pass: one term h(g(z)) = (w2/w1)⁵ from the shifted
  // matrix. Uses cm and cinv from ACC_FS.
  const GLSL_ACC = `
vec2 icosaTerm(vec2 w1, vec2 w2){
  vec2 r = cm(w2, cinv(w1));
  vec2 r2 = cm(r, r);
  return cm(cm(r2, r2), r);
}`;

  // GLSL for the display pass: Klein's invariant J. Uses cm and cdiv from DISP_FS.
  const GLSL_DISP = `
vec2 kleinJ(vec2 z){
  if (dot(z, z) > 1.0) z = -cdiv(vec2(1.0, 0.0), z);
  vec2 z5 = cm(cm(cm(z, z), cm(z, z)), z), z10 = cm(z5, z5), z15 = cm(z10, z5), z20 = cm(z10, z10);
  vec2 f = cm(z, z10 + 11.0 * z5 - vec2(1.0, 0.0));
  vec2 H = -(z20 + vec2(1.0, 0.0)) + 228.0 * (z15 - z5) - 494.0 * z10;
  vec2 f2 = cm(f, f);
  return cdiv(cm(cm(H, H), H), 1728.0 * cm(cm(f2, f2), f));
}`;

  // Panel text for index.html
  const TEXT = {
    title: 'F(z) = Σ<sub>g∈G</sub> h(g(z))',
    formula: 'h(z) = (z − 0.05)<sup>−5</sup>, G = icosahedral rotations (|G| = 60)',
    funcs: ['F = Σ h∘g (invariant)', 'h (single term)', 'Klein’s J (invariant)'],
  };

  return { POLE, ORDER, GENS, GEN_NAMES, group, gpuMatrices, h, F, kleinJ, GLSL_ACC, GLSL_DISP, TEXT };
})();
