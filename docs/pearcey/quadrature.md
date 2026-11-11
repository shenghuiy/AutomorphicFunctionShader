# Quadrature

[← The integration path](contour.md) · [overview](README.md) · next: [GPU evaluation](gpu.md)

After the change of variables in [contour.md](contour.md), each pixel computes a
real-line integral of a smooth, rapidly decaying function:

```math
\Psi(x, y) = \int_{-\infty}^{\infty} f(s)\,ds,
\qquad
f(s) = e^{i\varphi(t(s))}\,t'(s).
```

The program cuts it off at $`\lvert s\rvert \le S`$ and applies the midpoint rule.

## The interval

`halfWidth` picks $`S`$ from $`(x, y)`$:

```math
S(x, y) = \underbrace{3}_{\texttt{MARGIN}} + \max\!\Big(\sqrt{\lvert y\rvert},\; \sqrt[3]{\lvert x\rvert/2}\Big).
```

```cpp
inline double halfWidth(double x, double y) {
    // float, like the shader; any value near this one gives the same integral
    return double(float(MARGIN) + std::max(std::sqrt(std::fabs(float(y))), std::pow(0.5f * std::fabs(float(x)), 1.0f / 3.0f)));
}
```

**The max term bounds the stationary points.** Let
$`R = \max(\sqrt{\lvert y\rvert}, \sqrt[3]{\lvert x\rvert/2})`$. If $`\lvert s\rvert > R`$, then
$`2\lvert s\rvert^3 > \lvert x\rvert`$ and $`2\lvert s\rvert^3 = 2\lvert s\rvert\, s^2 > 2\lvert y\rvert\lvert s\rvert`$, so

```math
\lvert 4s^3\rvert = 2\lvert s\rvert^3 + 2\lvert s\rvert^3 > \lvert x\rvert + 2\lvert y\rvert\lvert s\rvert \ge \lvert 2ys + x\rvert,
```

and $`\varphi'(s) = 4s^3 + 2ys + x \ne 0`$. All real stationary points therefore lie in
$`\lvert s\rvert \le R`$.

**The margin makes the tail negligible.** Past the stationary points
$`\lvert\varphi'\rvert`$ grows like $`4\lvert s\rvert^3`$, and the integrand's modulus is
about $`e^{-\lambda\lvert\varphi'(s)\rvert}`$ ([contour.md](contour.md#why-the-integrand-decays-on-it)).
At the origin, for example, $`S = 3`$, $`\varphi'(3) = 108`$, and the integrand at the
ends is about $`e^{-54} \approx 4\times10^{-24}`$, far below double precision. Larger
$`\lvert x\rvert`$ or $`\lvert y\rvert`$ only make $`\varphi'(S)`$ larger.

**Why $`S`$ is computed in float.** The shader computes $`S`$ in single precision
(`float` `sqrt` and `pow`). The CPU copy rounds the same way, so both sides use
nearly the same path. The exact value does not matter much: by Cauchy's theorem
any $`S`$ in this range gives the same integral, and only the rounding and
discretisation errors change.

## The midpoint rule

With $`n`$ nodes on $`[-S, S]`$:

```math
\Psi \approx h\sum_{j=0}^{n-1} f(s_j),
\qquad
s_j = S\,u_j,
\qquad
u_j = -1 + \frac{2j+1}{n},
\qquad
h = \frac{2S}{n}.
```

The CPU reference does exactly this, and also accumulates the integrand of
$`\partial\Psi/\partial x`$, which is $`i\,t\,f`$:

```cpp
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
```

$`\varphi(t)`$ is evaluated as `t2*t2 + y*t2 + x*t` with `t2 = t*t`. This costs two
complex multiplies instead of three.

### Why it converges exponentially

For a general smooth function, the midpoint rule's error is $`O(h^2)`$. Here the
error is much smaller, because of two properties of $`f`$:

1. $`f`$ and all its derivatives are negligible at $`\pm S`$. The Euler–Maclaurin
   correction terms all involve derivatives at the endpoints, so they vanish to
   machine precision. The rule then behaves like the trapezoidal rule on the
   whole real line.
2. $`f`$ is analytic in a strip $`\lvert\operatorname{Im} s\rvert < d`$ around the real axis.

For such integrands the error of an equispaced rule decays like

```math
\lvert\text{error}\rvert = O\!\big(e^{-2\pi d/h}\big) = O\!\big(e^{-\pi d\, n/S}\big),
```

see Trefethen and Weideman, *The exponentially convergent trapezoidal rule*,
SIAM Review 56 (2014). The strip width $`d`$ is set by the nearest singularity of
$`f`$, and $`f`$ is built from entire functions except for $`q = \sqrt{1 + g^2}`$. Its
branch points, where $`\varphi'(s) = \pm i\kappa`$, set $`d`$. This is why the code
takes $`\kappa`$ larger than the decay condition requires ([contour.md](contour.md#why-κ-is-larger-than-it-needs-to-be)).

The node count needed grows with $`S`$, and so with the distance from the origin.
Here is the measured relative error of the rule against a 40 000-node result on
the same path:

| $`(x, y)`$ | $`S`$ | 64 nodes | 128 nodes | 256 nodes | 512 nodes |
|---|---|---|---|---|---|
| $`(0, 0)`$ | 3.00 | $`2\times10^{-15}`$ | $`1\times10^{-16}`$ | $`4\times10^{-16}`$ | $`1\times10^{-16}`$ |
| $`(1.5, -2)`$ | 4.41 | $`3\times10^{-8}`$ | $`1\times10^{-16}`$ | $`2\times10^{-16}`$ | $`2\times10^{-16}`$ |
| $`(2, -4)`$ | 5.00 | $`2\times10^{-6}`$ | $`3\times10^{-14}`$ | $`0`$ | $`1\times10^{-16}`$ |
| $`(7, -12)`$ | 6.46 | $`3\times10^{-2}`$ | $`5\times10^{-5}`$ | $`2\times10^{-12}`$ | $`5\times10^{-16}`$ |
| $`(10, -15)`$ | 6.87 | $`2\times10^{-1}`$ | $`2\times10^{-3}`$ | $`8\times10^{-9}`$ | $`2\times10^{-15}`$ |

The error falls by orders of magnitude with each doubling, and the default of 512
nodes reaches double precision across the default view. Its worst corner is
$`(\pm10, -15)`$. Further out it does not: at $`(0, -30)`$, 512 nodes give
$`4\times10^{-9}`$ and 1024 nodes $`2\times10^{-14}`$. Press `+` (or pass `--terms`)
when zoomed out.

## Bit-reversed node order and progressive rendering

The app sums the nodes in chunks across several frames, so that the window stays
responsive while dragging or resizing. A partial midpoint sum over the first $`m`$
nodes in the natural order covers only the left part of $`[-S, S]`$. Drawn as is,
it would show a meaningless picture.

`pearcey::nodes` fixes this with the order of the nodes. $`n`$ is a power of two,
and node $`i`$ of the upload holds midpoint $`j = \operatorname{rev}(i)`$, where
$`\operatorname{rev}`$ reverses the $`\log_2 n`$ bits of $`i`$:

```cpp
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
```

If $`m = 2^k`$, the indices $`i < m`$ use only the low $`k`$ bits. After reversal these
become the high $`k`$ bits, so $`j = (n/m)\,r`$ for $`r = 0, \dots, m-1`$. The first $`m`$
nodes are every $`(n/m)`$-th midpoint:

```math
u = -1 + \frac1n + \frac{2r}{m},
\qquad r = 0, \dots, m-1 .
```

This is an equispaced rule with spacing $`2/m`$. It is shifted from the $`m`$-point
midpoint rule by a fraction of a step, but that does not matter for the
exponential convergence above. The partial sum carries weight $`2/n`$ instead of
$`2/m`$, so it equals $`m/n`$ times an $`m`$-point estimate of $`\Psi`$:

```math
\sum_{i<m} h\, f(s_{\operatorname{rev}(i)}) \approx \frac{m}{n}\,\Psi .
```

The app therefore stops only at $`m = 128, 256, 512, \dots`$ nodes (`minPartial = 128`
in `main`, with each step doubling `done` in `frame()`). The pixel shader scales
the partial sum back up by $`n/m`$:

```hlsl
double2 plotValue(double4 a, double2 z) {
    return (uFunc==2 ? a.zw : a.xy)*((double)uCount/(double)max(uEnd,1u));
}
```

A frame drawn after 128 nodes is a complete, if coarser, quadrature. From the
table above it is already accurate to the eye over most of the default view.

## Packing nodes into the group-sum frame

The shared compute shader was written to sum $`\mathrm{term}(w_1, w_2)`$ over
matrices $`\begin{pmatrix} a & b \\ c & d\end{pmatrix}`$, with $`w_1 = az + b`$ and
$`w_2 = cz + d`$ ([background.md](../background.md)). Each node is uploaded as

```math
\begin{pmatrix} 1 & 0 \\ 0 & u_j + i\,w \end{pmatrix},
\qquad w = \frac{2}{n},
```

so the frame hands the Pearcey `term` $`w_1 = z = (x, y)`$ and
$`w_2 = (u_j, w)`$, the node position on $`(-1, 1)`$ and its weight. The shader
multiplies both by $`S`$ to get $`s_j = S u_j`$ and $`h = S w`$. The same uploaded
nodes serve every pixel, whatever its $`S`$.
