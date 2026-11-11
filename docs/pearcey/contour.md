# The integration path

[← The integral](definition.md) · [overview](README.md) · next: [Quadrature](quadrature.md)

## Why not the real line

On the real line the integrand $`e^{i\varphi(t)}`$ has modulus 1 for every $`t`$. The
integral converges only because the oscillations cancel, and truncating it at
$`\pm T`$ leaves an error of order $`1/T^3`$. Any rule that sums the integrand on the
real line has to sum a huge number of nodes, and those nodes cancel each other
almost completely.

The integrand is entire in $`t`$. By Cauchy's theorem the path can therefore be moved
to any curve with the same ends at infinity, as long as the integrand decays at
those ends. For $`t^4`$ the decaying sectors at infinity are centred on $`\arg t = \pi/8`$
and $`\arg t = 7\pi/8`$. A path that stays near the real axis in the middle and
tends into those sectors gives the same integral.

## The simplest choice and its problem

The simplest deformation is the straight line $`t = e^{i\pi/8}s`$. On it,
$`i t^4 = -s^4`$, and the integrand is a Gaussian-like bump. The program keeps this as
an independent check (`psiRotated` in `pearcey.h`). It fails for large negative
$`y`$, however. On the same line,

```math
\big\lvert e^{i y t^2}\big\rvert = e^{-y \sin(\pi/4)\, s^2} = e^{\lvert y\rvert s^2/\sqrt2}
\quad (y < 0),
```

and the integrand $`e^{-s^4 + \lvert y\rvert s^2/\sqrt2}`$ peaks at
$`e^{y^2/8}`$ before it decays. For $`y = -15`$, the bottom of the default view,
that peak is $`e^{28} \approx 10^{12}`$. The result is of order 1, so the sum loses
about twelve digits to cancellation. The rotated line only works near the origin:

```cpp
// Independent check on the straight path t = e^{i pi/8} s, where exp(i t^4) = exp(-s^4).
// Only for small |x|, |y|: for y << 0 the integrand grows to exp(y^2 / 8) before it decays.
```

## A path that follows the phase

The renderer instead uses a path adapted to each $`(x, y)`$. For real $`s`$, set

```math
t(s) = s + i\,\delta(s),
\qquad
\delta(s) = \lambda\,\frac{g}{\sqrt{1 + g^2}},
\qquad
g = \frac{\varphi'(s)}{\kappa},
```

with constants $`\lambda = \tfrac12`$ (`LAMBDA`) and $`\kappa = 16 S\lambda^2`$, where
$`S`$ is the half-width of the interval ([quadrature.md](quadrature.md#the-interval)).

The path lifts the integration off the real axis in the direction in which
$`e^{i\varphi}`$ decreases:

- Near a real stationary point, $`\varphi'(s) \approx 0`$, so $`\delta \approx 0`$. The path
  passes through it on the real axis, where the main contribution comes from.
- Where $`\lvert\varphi'(s)\rvert \gg \kappa`$, $`g/\sqrt{1+g^2} \to \pm 1`$, and the path
  runs at height $`\pm\lambda`$, with the sign of $`\varphi'`$.
- $`\delta`$ is smooth (analytic in $`s`$ near the real axis), and $`\lvert\delta\rvert < \lambda`$.

### Why the integrand decays on it

To first order in $`\delta`$, $`\varphi(s + i\delta) \approx \varphi(s) + i\delta\,\varphi'(s)`$, so

```math
\operatorname{Im}\varphi(t(s)) \approx \delta(s)\,\varphi'(s)
= \frac{\lambda\kappa\,g^2}{\sqrt{1+g^2}} \;\ge\; 0 ,
```

because $`\delta`$ and $`\varphi'`$ always have the same sign. The modulus of the
integrand is $`e^{-\operatorname{Im}\varphi}`$, so it is at most 1, and once
$`\lvert g\rvert \gg 1`$ it falls off like $`e^{-\lambda\lvert\varphi'(s)\rvert}`$.

The exact imaginary part follows from expanding $`(s + i\delta)^4`$,
$`(s + i\delta)^2`$ and $`x(s + i\delta)`$:

```math
\operatorname{Im}\varphi = (4s^3\delta - 4s\delta^3) + y\,(2s\delta) + x\,\delta
= \delta\,\big(\varphi'(s) - 4s\,\delta^2\big).
```

With $`q = \sqrt{1+g^2}`$ this is

```math
\operatorname{Im}\varphi = \frac{\lambda g^2}{q^3}\Big(\kappa\,(1 + g^2) - 4\lambda^2\, s\, g\Big).
```

The bracket, a quadratic in $`g`$, is non-negative for all $`g`$ if and only if
$`\kappa \ge 2\lambda^2\lvert s\rvert`$. On the integration interval $`\lvert s\rvert \le S`$
it is therefore enough that $`\kappa \ge 2S\lambda^2`$. The comment in `pearcey.h`
uses the simpler sufficient bound $`\kappa \ge 4S\lambda^2`$, and the code takes
$`\kappa`$ four times that:

```cpp
// Exactly, Im phi = delta (phi' - 4 s delta^2), which
// is non-negative for |s| <= S once kappa >= 4 S LAMBDA^2. kappa is 4 times
// that: a gentler bend keeps the square root's branch points away from the real
// axis, and the midpoint rule then needs about half as many nodes.
double S = halfWidth(x, y), kappa = 16 * S * LAMBDA * LAMBDA;
```

So $`\lvert e^{i\varphi(t(s))}\rvert \le 1`$ holds exactly on the whole interval.
Nothing in the sum is larger than the answer by more than a modest factor. That is
why fp32 summation already gives a usable picture ([gpu.md](gpu.md#fp32-and-fp64)).

### Why κ is larger than it needs to be

The path has a cost: it brings in the square root $`q = \sqrt{1 + g(s)^2}`$. As a
function of complex $`s`$, $`q`$ has branch points where $`g(s) = \pm i`$, that is, where
$`\varphi'(s) = \pm i\kappa`$. Those singularities limit how fast the quadrature
converges ([quadrature.md](quadrature.md#why-it-converges-exponentially)): the farther
they are from the real $`s`$-axis, the faster. A larger $`\kappa`$ moves them
farther out. The bend becomes gentler, so the decay starts a little later, but
the integrand becomes smoother. With $`\kappa = 16 S\lambda^2`$, the midpoint rule
needs about half as many nodes as at the minimum.

### The Jacobian

With $`dt = t'(s)\,ds`$, the integral becomes a real-line integral over $`s`$:

```math
\Psi(x, y) = \int e^{i\varphi(t(s))}\, t'(s)\, ds,
\qquad
t'(s) = 1 + i\lambda\,\frac{g'(s)}{(1 + g^2)^{3/2}},
\qquad
g'(s) = \frac{\varphi''(s)}{\kappa} = \frac{12s^2 + 2y}{\kappa}.
```

The derivative uses $`\frac{d}{dg}\frac{g}{\sqrt{1+g^2}} = (1+g^2)^{-3/2}`$. In code
(`pearcey::psi`):

```cpp
double g = (4 * s * s * s + 2 * y * s + x) / kappa, dg = (12 * s * s + 2 * y) / kappa;
double q = std::sqrt(1 + g * g);
cd t(s, LAMBDA * g / q), dt(1, LAMBDA * dg / (q * q * q));
```

The ends of the path are not on the rays $`\arg t = \pi/8, 7\pi/8`$. They sit at
height $`\pm\lambda`$ above and below the real axis. That does not matter: the
integrand is already below double precision where the interval is cut off,
because $`\operatorname{Im}\varphi`$ grows like $`\lambda\lvert\varphi'(s)\rvert`$
there. The next note picks the cut-off.
