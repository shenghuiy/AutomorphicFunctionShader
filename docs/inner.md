# Inner function phase portrait

The portrait illustrates the singular factor in the canonical factorization
`I = exp(iγ) B Sμ` of an inner function, where `B` is a Blaschke product.
The Axler reference image describes a non-atomic singular measure on a Cantor
subset of the unit circle. The renderer uses one example from this family,
with parameters chosen visually to match the reference.

## Cantor measure

The measure has total mass 1, divided equally among three scaled copies of
the middle-thirds Cantor probability measure. For a Cantor-distributed
point `t` in [0, 1] and a sector `k` in {0, 1, 2}, set
`θ = 0.12 + 2π(k + t)/3`. Place the point on the unit circle at angle
`θ + 0.10 sin(2θ + 0.4)`. This smooth change of angle preserves the
singular, non-atomic character of the measure.

The shader evaluates `Sμ(z) = exp(∫ (z + w)/(z − w) dμ(w))`. The limiting
function is zero-free in the open disk and continues analytically across gaps
in the Cantor support. The same formula gives the exterior portrait and
satisfies `Sμ(1/conj(z)) = 1/conj(Sμ(z))`.

## Rendering and resolution

Build with `native\dx12\build.bat inner` and run `native\build\inner.exe`.
The default view is centered at `(0, −0.04)` with span 1.95. Press **F** to
cycle through the Cantor portrait, a single atom, and a Blaschke factor.
Press **P** to save a PNG.

At depth 10, each copy contributes 1,024 surviving interval midpoints,
giving 3,072 quadrature nodes of mass 1/3,072 each. Use `--terms N` or the
`+` and `-` keys to change resolution. The renderer rounds the count down
to three times a power of two, with a minimum of three nodes and a maximum
depth of 20 per copy.

The discrete measures converge to the Cantor measure as depth increases.
Features near its boundary support are especially sensitive to resolution.
The shader clamps the real part of the displayed exponent to [−70, 70] to
avoid overflow and underflow. It disables extreme-magnitude fading to keep
boundary phase detail bright.

## Numerical checks

GPU exponent sums passed comparisons with the CPU evaluator in FP32 and FP64.
Increasing the quadrature depth from 10 to 11 changed the function value by
at most 1.89 × 10⁻⁸ at the four sampled points `0`, `0.3 + 0.2i`,
`0.8 + 0.1i`, and `0.97 + 0.01i`. This check measures convergence only at
those points. A uniform error bound near the boundary remains unknown.
