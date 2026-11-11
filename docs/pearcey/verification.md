# Verification

[← GPU evaluation](gpu.md) · [overview](README.md)

The Pearcey program has three ways to check its numbers: the `I` key, the
`--verify-gpu` option, and the independent checks below, which were used to
write these notes.

## The `I` key: checking the quadrature on the CPU

Other programs in this repository use `I` to test invariance under their group.
The Pearcey program has no group, so `invariance()` in `pearcey.cpp` tests the
quadrature instead. It runs on the CPU with the current node count $`n`$ and
prints, for $`\Psi(0, 0)`$ and four sample points:

| Comparison | What it tests |
|---|---|
| $`\Psi(0,0)`$ vs $`2\,\Gamma(\tfrac54)\,e^{i\pi/8}`$ | Exact value ([definition.md](definition.md#symmetry-and-special-values)): path, Jacobian and weights. |
| vs `psiRotated` (only for $`y \ge -4`$) | Agreement with a different contour, $`t = e^{i\pi/8}s`$, 4000 nodes on $`\lvert s\rvert \le 8`$. By Cauchy's theorem both give the same integral. |
| vs 20 000 nodes | Convergence: the rule with $`n`$ nodes against a much finer rule on the same path. |
| $`\Psi(-x, y)`$ vs $`\Psi(x, y)`$ | The mirror symmetry (DLMF 36.2.22). Changing $`x \to -x`$ flips the sign of $`\varphi'(-s)`$, so the path becomes its mirror image $`t \to -t`$, and the midpoint nodes are symmetric. The two sums have the same terms in reverse order, so this is a consistency check of the code: it should agree to rounding error. |

```cpp
const double pts[][2] = {{1.5, -2.0}, {-3.0, 1.0}, {2.0, -4.0}, {7.0, -12.0}};
for (const auto& xy : pts) {
    pearcey::psi(xy[0], xy[1], nodeCount, p);
    ...
    if (xy[1] >= -4) std::printf("   vs rotated path: rel.err %.1e", err(p[0], pearcey::psiRotated(xy[0], xy[1])));
    pearcey::psi(xy[0], xy[1], 20000, q);
    std::printf("   vs 20000 nodes: rel.err %.1e", err(p[0], q[0]));
    pearcey::psi(-xy[0], xy[1], nodeCount, q);
    std::printf("   Psi(-x,y) = Psi(x,y): rel.err %.1e\n", err(q[0], p[0]));
}
```

The rotated path is skipped for $`y < -4`$ because its integrand grows to $`e^{y^2/8}`$
there ([contour.md](contour.md#the-simplest-choice-and-its-problem)). At
$`y = -12`$ that is $`e^{18}`$, so the reference itself would lose about eight
digits.

The check also runs automatically before `--snapshot` saves its PNG.

## `--verify-gpu`: GPU against CPU

```powershell
.\native\build\pearcey.exe --fp64 --verify-gpu --snapshot check-pearcey.png
```

When the render is complete, `verifyGpuValues` in `app.h` reads back the raw
accumulator (`acc`, both sums) at three pixels and compares it with
`reference()`:

```cpp
// A much finer rule on the same path; the GPU sum has converged if it agrees.
static void reference(cd z, cd out[2]) { pearcey::psi(z.real(), z.imag(), 20000, out); }
```

This checks the whole GPU chain at once: node upload, the shared summation loop,
the HLSL integrand, the home-made `dsqrt`, `dexp` and `expi`, and convergence of
the GPU's $`n`$-node rule. It fails with an error if the relative error of either
$`\Psi`$ or $`\partial\Psi/\partial x`$ exceeds $`10^{-9}`$ in fp64 or $`5\times10^{-3}`$
in fp32.

## Independent checks

These are not part of the program. They were run in Python
(NumPy and mpmath) while writing these notes, against a NumPy port of
`pearcey::psi`:

- **Against DLMF 36.2.16.** $`\Psi(0,0) = 1.674813393538173 + 0.693730422047619\,i`$,
  which agrees with $`2\Gamma(\tfrac54)e^{i\pi/8}`$ to $`10^{-16}`$.
- **Against arbitrary precision.** mpmath integration on the rotated contour at
  30 digits, at $`(0,0)`$, $`(1.5,-2)`$, $`(-3,1)`$ and $`(2,-4)`$: relative differences
  of $`1`$ to $`2\times10^{-16}`$.
- **Against the Bessel form on $`x = 0`$** ([definition.md](definition.md#symmetry-and-special-values)),
  which works for every $`y`$ and so tests the region the rotated path cannot reach:

  | $`y`$ | 256 nodes | 512 nodes | 1024 nodes |
  |---|---|---|---|
  | 5 | $`1\times10^{-16}`$ | $`0`$ | $`1\times10^{-16}`$ |
  | −12 | $`7\times10^{-14}`$ | $`1\times10^{-14}`$ | $`6\times10^{-16}`$ |
  | −15 | $`7\times10^{-11}`$ | $`3\times10^{-15}`$ | $`2\times10^{-15}`$ |
  | −30 | $`5\times10^{-3}`$ | $`4\times10^{-9}`$ | $`2\times10^{-14}`$ |

- **Convergence in n.** See the table in
  [quadrature.md](quadrature.md#why-it-converges-exponentially).

The $`y = -30`$ row shows the default 512 nodes running out of accuracy outside the
default view. Doubling the node count with `+` brings it back to near double
precision, as the exponential convergence predicts.
