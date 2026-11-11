# The integral

[← Pearcey overview](README.md) · next: [The integration path](contour.md)

## Definition and notation

The program plots

```math
\Psi(x, y) = \int_{-\infty}^{\infty} e^{i\varphi(t)}\,dt,
\qquad
\varphi(t) = t^4 + y\,t^2 + x\,t,
```

with $`(x, y)`$ read from the pixel as $`z = x + iy`$. This is the cusp case $`K = 2`$
of the canonical integrals in [DLMF §36.2](https://dlmf.nist.gov/36.2):

```math
\Psi_K(\mathbf{x}) = \int_{-\infty}^{\infty} \exp\big(i\,\Phi_K(t;\mathbf{x})\big)\,dt,
\qquad
\Phi_K(t;\mathbf{x}) = t^{K+2} + \sum_{m=1}^{K} x_m t^m
\qquad\text{(DLMF 36.2.1, 36.2.4)}.
```

So $`\Psi(x, y) = \Psi_2(x_1 = x,\; x_2 = y)`$: the program's $`x`$ multiplies $`t`$ and
its $`y`$ multiplies $`t^2`$. The Pearcey function $`P`$ of DLMF 36.2.14 takes its
arguments the other way round, $`\Psi_2(x_1, x_2) = P(x_2, x_1)`$, so
$`\Psi(x, y) = P(y, x)`$.

The integral converges only conditionally: $`\lvert e^{i\varphi}\rvert = 1`$ on the real
line, and it converges because the oscillations speed up like $`t^3`$. The DLMF notes
that convergence can be shown by rotating the contour. The renderer does the same
thing in a more careful way ([contour.md](contour.md)).

## What the code plots

The `F` key cycles through three displays (`funcNames` in `main`):

| `F` | Display | Accumulator |
|---|---|---|
| 0 | $`\Psi`$ | `acc.xy` |
| 1 | $`\Psi`$ with the caustic drawn in black | `acc.xy` |
| 2 | $`\partial\Psi/\partial x`$ | `acc.zw` |

Differentiating under the integral sign gives

```math
\frac{\partial\Psi}{\partial x} = \int_{-\infty}^{\infty} i\,t\,e^{i\varphi(t)}\,dt,
```

which is what the second accumulator sums (see [gpu.md](gpu.md)).

$`\Psi`$ is a smooth function of the two real variables $`x`$ and $`y`$, but it is not
an analytic function of $`z = x + iy`$. A phase portrait of it therefore does not
behave like one of a holomorphic function: the zeros inside the cusp do not all
turn the colours the same way round.

## Symmetry and special values

**Mirror symmetry in x.** Substituting $`t \to -t`$ maps $`\varphi(t; x, y)`$ to
$`\varphi(t; -x, y)`$, so

```math
\Psi(-x, y) = \Psi(x, y),
```

which is DLMF 36.2.22 for $`K = 2`$. The picture is symmetric about the imaginary
axis, and $`\partial\Psi/\partial x`$ is odd in $`x`$. The `I` key checks the first
relation numerically ([verification.md](verification.md)).

**Value at the origin.** On the ray $`t = e^{i\pi/8}u`$, $`u \ge 0`$, we have
$`i t^4 = i \cdot i u^4 = -u^4`$, and the arcs at infinity contribute nothing, so

```math
\Psi(0, 0) = 2\int_0^\infty e^{i t^4}\,dt
= 2 e^{i\pi/8}\int_0^\infty e^{-u^4}\,du
= 2\,\Gamma\!\big(\tfrac54\big)\,e^{i\pi/8}
\approx 1.67481 + 0.69373\,i,
```

the value listed in DLMF 36.2.16. The `I` key compares the quadrature against
this value:

```cpp
cd exact = 2 * std::tgamma(1.25) * std::polar(1.0, 3.14159265358979323846 / 8);
```

**The imaginary axis.** On $`x = 0`$ the integral reduces to Bessel functions
(DLMF 36.2.19). In the program's variables:

```math
\Psi(0, y) = \frac{\pi}{2}\sqrt{\frac{\lvert y\rvert}{2}}\;e^{-iy^2/8}
\Big[e^{i\pi/8} J_{-1/4}\!\big(\tfrac{y^2}{8}\big)
- \operatorname{sgn}(y)\, e^{-i\pi/8} J_{1/4}\!\big(\tfrac{y^2}{8}\big)\Big].
```

This form was checked numerically against the integral to about $`10^{-17}`$. The
program does not use it, but it is a useful independent test for large $`\lvert y\rvert`$,
where the program's own rotated-path check cannot be used
([verification.md](verification.md)).

## The cusp caustic

By stationary phase, the main contributions to $`\Psi`$ come from the real
**stationary points** of the phase, the real roots of

```math
\varphi'(t) = 4t^3 + 2y\,t + x = 0 .
```

A cubic has either one or three real roots. The two cases are separated by the
curve where two roots merge, that is, where $`\varphi'(t) = \varphi''(t) = 0`$:

```math
\varphi''(t) = 12t^2 + 2y = 0 \;\Rightarrow\; t^2 = -\frac{y}{6},
\qquad
x = -t\,(4t^2 + 2y) = -\frac{4y}{3}\,t .
```

Squaring the second equation and substituting the first:

```math
x^2 = \frac{16y^2}{9}\cdot\Big(-\frac{y}{6}\Big) = -\frac{8y^3}{27}
\quad\Longleftrightarrow\quad
27x^2 + 8y^3 = 0 .
```

This is the **cusp caustic**, the discriminant of the cubic $`\varphi'`$. It has its
point at the origin and opens towards $`y < 0`$:

- Inside the cusp, where $`27x^2 + 8y^3 < 0`$, there are three real stationary
  points. Their contributions interfere, which gives the lattice of zeros seen in
  the plot.
- Outside it there is one real stationary point, and $`\Psi`$ has no interference
  pattern. Away from the caustic its modulus varies slowly.

The default view, centre $`(0, -5)`$ with span 20, shows the cusp point and a long
stretch of both branches. `F = 1` draws the curve; [gpu.md](gpu.md#the-caustic-overlay)
explains how.

The stationary points also set how far along the real axis the integrand matters.
That is the subject of the next note.
