# Mathematical background

This note explains what the renderers plot and why the plot looks the way it
does. The formulas match the code in `native/dx12/poincare/` (`kleinian.h` and
`poincare.cpp`), `native/dx12/common/shaders.h` and `index.html`.

## The function

Both renderers plot

```math
f(z) = \frac{P_{-1,-3}(z)}{P_{0,-4}(z)},
\qquad
P_{k,l}(z) = \sum_{g \in G} (a_g z + b_g)^{k}\,(c_g z + d_g)^{l},
```

where the sum runs over every element

```math
g = \begin{pmatrix} a_g & b_g \\ c_g & d_g \end{pmatrix}
```

of a discrete group $`G \subset \mathrm{SL}(2,\mathbb{C})`$. Written out, the two series are

```math
P_{0,-4}(z) = \sum_{g \in G} \frac{1}{(c_g z + d_g)^4},
\qquad
P_{-1,-3}(z) = \sum_{g \in G} \frac{1}{(a_g z + b_g)(c_g z + d_g)^3}.
```

Each $`g`$ acts on the Riemann sphere as the Möbius transformation
$`g(z) = (a_g z + b_g)/(c_g z + d_g)`$. Because $`\det g = 1`$, its derivative is
$`g'(z) = (c_g z + d_g)^{-2}`$. The series can therefore be written without
matrix entries:

```math
P_{0,-4}(z) = \sum_{g \in G} g'(z)^2,
\qquad
P_{-1,-3}(z) = \sum_{g \in G} \frac{g'(z)^2}{g(z)}.
```

$`P_{0,-4}`$ is the classical Poincaré theta series of weight 4. It is the
series Poincaré used in the 1880s to construct automorphic functions for
Fuchsian and Kleinian groups.

## Why $`f`$ is automorphic

Fix a group element

```math
h = \begin{pmatrix} \alpha & \beta \\ \gamma & \delta \end{pmatrix}.
```

For any $`g \in G`$, put $`g' = gh`$. Then

```math
a_g\,h(z) + b_g = \frac{a_{g'} z + b_{g'}}{\gamma z + \delta},
\qquad
c_g\,h(z) + d_g = \frac{c_{g'} z + d_{g'}}{\gamma z + \delta}.
```

The map $`g \mapsto gh`$ is a bijection of $`G`$, so reindexing the sum gives

```math
P_{k,l}\big(h(z)\big) = (\gamma z + \delta)^{-(k+l)}\,P_{k,l}(z).
```

Both series have $`k + l = -4`$, so both pick up the same factor
$`(\gamma z + \delta)^{4}`$. In the ratio the factors cancel:

```math
f\big(h(z)\big) = f(z) \qquad \text{for every } h \in G.
```

So $`f`$ is a meromorphic function on the domain of discontinuity that is
invariant under $`G`$: an *automorphic function*. It descends to a meromorphic
function on the quotient surface $`\Omega(G)/G`$. The plot is the same
picture tiled over and over, shrinking toward the limit set.

You can check the invariance numerically. Press `I` in the native renderer, or
use the invariance check in the browser renderer. Either one evaluates $`f`$ at
$`z_0 = 0.3 + 0.6i`$ and at $`g_1^{\pm 1}(z_0)`$, $`g_2^{\pm 1}(z_0)`$, and prints
the relative errors. The errors get smaller as more terms are added.

## The group

$`G = \langle g_1, g_2 \rangle`$ has two generators:

```math
g_1 = \begin{pmatrix} 1 + 0.1i & B \\ C & 1 + 0.1i \end{pmatrix},
\qquad
g_2 = \begin{pmatrix} 1.1 - i & 1.1 \\ 1.1 & 1.1 + i \end{pmatrix},
```

with

```math
B \approx -0.1099758 - 0.0011538\,i,
\qquad
C \approx 0.0718424 - 1.8193356\,i.
```

$`B`$ and $`C`$ are chosen so that $`BC = (1 + 0.1i)^2 - 1`$, which makes
$`\det g_1 = 1`$. For $`g_2`$, $`\det g_2 = (1.1 - i)(1.1 + i) - 1.1^2 = 1`$.

- $`\mathrm{tr}\, g_1 = 2 + 0.2i`$ is not real, so $`g_1`$ is **loxodromic**:
  it spirals points toward its attracting fixed point.
- $`\mathrm{tr}\, g_2 = 2.2`$ is real and greater than 2, so $`g_2`$ is
  **hyperbolic**: it moves points along circular arcs between its two fixed
  points.

The code treats $`G`$ as a free group on $`g_1, g_2`$, as in a Schottky group. In
a free group every *reduced word* (a product of generators with no
$`g\,g^{-1}`$ next to each other) is a different group element. That makes it
possible to list the group without checking for duplicates.

## Convergence

Poincaré's classical result says that for a Kleinian group the series
$`\sum_g |g'(z)|^{q}`$ converges for $`q \ge 2`$, uniformly on compact subsets of
the domain of discontinuity $`\Omega(G)`$. Here $`|g'(z)|^2 = |c_g z + d_g|^{-4}`$,
which is exactly the borderline case $`q = 2`$. So $`P_{0,-4}`$ converges
absolutely on $`\Omega(G)`$. $`P_{-1,-3}`$ has the extra factor $`1/g(z)`$. As long
as $`0`$ is not in the limit set, $`g(z)`$ stays away from $`0`$ for all but
finitely many $`g`$, and that series converges as well.

Near the limit set $`\Lambda(G)`$ the terms decay more and more slowly, and the
series diverges on $`\Lambda(G)`$ itself. In the plot this shows up as fine
detail that keeps changing as more terms are added close to the limit set.

## Truncating the sum

The renderers cannot sum over an infinite group, so they keep the $`N`$
elements with the smallest Frobenius norm

```math
\|g\|^2 = |a_g|^2 + |b_g|^2 + |c_g|^2 + |d_g|^2.
```

For a fixed $`z`$ in $`\Omega(G)`$, $`|c_g z + d_g|`$ grows like $`\|g\|`$, so terms
with small norm contribute the most. The enumeration (`kleinian::enumerate` in
`native/dx12/poincare/kleinian.h`) is a best-first search over reduced words:

1. Start from the identity.
2. Repeatedly take the element with the smallest norm from a min-heap and
   output it.
3. Multiply it on the left by each of the four generators
   $`g_1^{\pm 1}, g_2^{\pm 1}`$, skipping the one that would cancel the last letter, and push
   the results onto the heap.

Doubling `N` with `+` adds the next shell of larger-norm elements. Away from
the limit set the picture barely changes. Near the limit set more detail
appears.

## Numerical precision

Large-norm group elements produce terms that are tiny compared with the
leading ones. Summing many of them in single precision loses their
contribution and can also cancel badly near zeros of $`P_{0,-4}`$. That is why
both summation modes accumulate into a double-precision buffer, and why the
native renderer lets you switch the per-term arithmetic between fp32 and fp64
with `D`. The browser renderer uses fp32 throughout.

## Reading the domain coloring

Each pixel shows the complex value $`w = f(z)`$ (or $`P_{0,-4}`$ or
$`P_{-1,-3}`$ if you switch with `F`):

- **Hue** is $`\arg w`$. Red, yellow, green, cyan, blue and magenta go once
  around the color wheel.
- **Brightness bands** follow $`\log_2 |w|`$. Each band is a doubling of the
  modulus, and `[` / `]` change how dense the bands are.

Some features are easy to spot:

- A **zero** of $`f`$ is a point where all hues meet, going around the wheel
  counterclockwise. A **pole** has the hues in the opposite order. Because $`f`$
  is automorphic, every zero and pole is repeated in every copy of the
  fundamental domain.
- Zeros of $`P_{0,-4}`$ are poles of $`f`$, and zeros of $`P_{-1,-3}`$ are zeros of
  $`f`$. Switching to the individual series with `F` shows which is which.
  Unlike $`f`$, the individual series are *not* invariant. They are automorphic
  forms, so their pictures are distorted by the factor $`(\gamma z + \delta)^4`$
  between copies.
- The shrinking copies of the pattern pile up along the **limit set**, a
  fractal curve or dust where the series stop converging.

## Further reading

- [coxeter.md](coxeter.md): the finite counterpart in this repo. It covers the
  icosahedral group, Coxeter and triangle groups, and how averaging over a finite
  group relates to these Poincaré series.

- H. Poincaré, "Mémoire sur les groupes kleinéens", *Acta Mathematica* 3
  (1883).
- D. Mumford, C. Series and D. Wright, *Indra's Pearls: The Vision of Felix
  Klein*, Cambridge University Press (2002). Schottky groups, limit sets, and
  enumerating group words.
- I. Kra, *Automorphic Forms and Kleinian Groups*, W. A. Benjamin (1972).
  Convergence of Poincaré series and automorphic forms.
- E. Wegert, *Visual Complex Functions: An Introduction with Phase Portraits*,
  Birkhäuser (2012). Domain coloring.
