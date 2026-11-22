# Inner function phase portrait

The portrait illustrates the singular factor in the canonical factorization
$`I = e^{i\gamma} B S_\mu`$ of an inner function, where $`B`$ is a Blaschke
product. The factorization is given in Garcia and Ross, §2.4, Theorem 2.12. [1](#references)

The Axler reference image describes a non-atomic singular measure on a Cantor
subset of the unit circle. The renderer uses one example from this family,
with parameters chosen visually to match the reference.

## Inner functions and Blaschke products

An inner function is a bounded analytic function on the open unit disk whose
radial boundary values have modulus one almost everywhere on the unit circle.
Its canonical factorization separates its zeros into a Blaschke product and
its singular measure into a zero-free singular inner factor. [1](#references)

A Blaschke product is a finite or convergent infinite product of disk factors.
The renderer's Blaschke example uses a single factor
$`b_\alpha(z) = \frac{z-\alpha}{1-\overline{\alpha}z}`$, with
$`\alpha = 0.35 + 0.2i`$. For finite products and their properties, see the
survey by Garcia, Mashreghi, and Ross. [2](#references)

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

On Windows, build with `native\dx12\build.bat inner` and run
`native\build\inner.exe`. On macOS, build with `native/metal/build.sh inner`
and run `native/build/inner`.
The default view is centered at `(0, −0.04)` with span 2.25. The default palette is the shared GreenPinkTones scheme.
Press **C** to switch to the hue wheel. Press **F** to
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

## References

1. Stephan Ramon Garcia and William T. Ross. *Model spaces: a survey* (2013),
   §2.4, especially Theorem 2.12 on canonical factorization.
   [arXiv:1312.5018](https://arxiv.org/html/1312.5018#S2.SS4).
2. Stephan Ramon Garcia, Javad Mashreghi, and William T. Ross.
   *Finite Blaschke products: a survey* (2015, revised 2016).
   [arXiv:1512.05444](https://arxiv.org/abs/1512.05444).
