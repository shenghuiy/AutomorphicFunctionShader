# Metal renderer (macOS)

Native Metal versions of the four programs of the DirectX 12 renderer in
[`../dx12/`](../dx12/), with the same plots, controls and command-line options,
and a fifth program, `invariant_j`, that so far exists only for Metal.

| Program | Plot | Default center |
|---|---|---|
| `native/build/poincare` | Kleinian ⟨g₁, g₂⟩, truncated Poincaré series | `0 -0.35` |
| `native/build/icosahedral` | 60 icosahedral rotations, exact sum | `0 0` |
| `native/build/pearcey` | Pearcey integral, quadrature per pixel | `0 -5` |
| `native/build/kantorovich` | Harmonic function and Newton basins | `0 0` |
| `native/build/invariant_j` | Modular j-invariant, SL(2,ℤ) | `0 0.9` |

## Build

Requires a Mac with a Metal GPU and Apple clang from Xcode or the Command Line
Tools. Tested on an Apple M1 Max with macOS 27 and Apple clang 21. The shaders are compiled at
startup by the Metal framework, so the offline `metal` compiler is not needed.

```sh
native/metal/build.sh               # all five programs
native/metal/build.sh pearcey       # one program
native/metal/build.sh poincare p2   # native/build/p2, e.g. while poincare is running
```

Outputs go to `native/build/`, next to the Windows executables.

## Run

```sh
native/build/poincare --center 0.2 -0.6 --span 1.5 --terms 20000 --fp64 --snapshot zoom.png
native/build/icosahedral --sphere --center 1.5 0.8
native/build/pearcey --function 1 --colors greenpink
native/build/kantorovich --fp64 --function 2
native/build/invariant_j --fp64 --center 0.3 0.01 --span 0.02
native/build/invariant_j --disk --function 1
```

The options and keys are those described in the [top-level README](../../README.md)
(`--help` lists them). The window's drawable is 1400 × 1000 pixels, so snapshots
match the Windows ones. Pinch on a trackpad also zooms. Run from a terminal: the
console output (render times, invariance checks, saved files) goes to stdout.
In an interactive terminal it shows only the last 5 lines below the start-up
banner, redrawn in place as new lines arrive; piped output and `--snapshot` runs
print every line.

## Precision: fp64 is double-float

Apple GPUs have no double-precision arithmetic, so `--fp64` and the `D` key select
**double-float** arithmetic: each real number is an unevaluated sum `hi + lo` of
two floats, about 48 significant bits, computed with error-free two-sum and
`fma`-based two-product. `--verify-gpu` checks against the CPU reference with the
same thresholds as on Windows (1e-9 for fp64, 5e-3 for fp32); observed relative
errors in fp64 are 1e-15 to 3e-13 for all four programs, including a zoom to a
span of 1e-9.

Differences from the DirectX 12 version:

- The pixel coordinate and the per-pixel accumulator are double-floats in both
  modes, like the doubles on Windows. The accumulator's `hi` part is a float, so
  sums beyond about 3·10³⁸ overflow, where a double would not.
- `D` switches the precision of the displayed value too (Klein's `J`, the Newton
  basins, `N(z)-z`): fp32 display is fast, fp64 display matches the Windows
  double-precision display. On Windows the display is always double.
- The library is compiled with fast math off (`MTLMathModeSafe`) and floating-point
  contraction disabled, which the double-float arithmetic depends on.
- The window is redrawn only when something changes.

## Modular j-invariant

`invariant_j` plots Klein's j-invariant on the upper half-plane,

`j(τ) = E₄(τ)³ / Δ(τ) = 1/q + 744 + 196884 q + …`,  `q = e^{2πiτ}`,

the automorphic function of the modular group SL(2,ℤ):
`j((aτ + b)/(cτ + d)) = j(τ)`. The default view is `−0.1 < Im τ < 1.9` on the
shorter window side; the lower half-plane is drawn dark.

Each pixel moves τ into the standard fundamental domain `|Re τ| ≤ ½, |τ| ≥ 1`
with `τ → τ − n` and `τ → −1/τ`. There `|q| ≤ e^{−π√3} < 0.0044`, so 12 terms of
each q-series reach double precision:

- `E₄ = 1 + 240 Σ n³qⁿ/(1 − qⁿ)` and `E₆ = 1 − 504 Σ n⁵qⁿ/(1 − qⁿ)` (Lambert series);
- `Δ = q Π (1 − qⁿ)²⁴`, a product, which avoids the cancellation in
  `(E₄³ − E₆²)/1728` high up the half-plane.

`F` cycles j, `E₄` and `E₆`. The Eisenstein series are modular forms of weight
4 and 6, so their values at the reduced point are carried back with the factor
`τ⁻ᵏ` of each inversion. `E₄` vanishes at `ρ = e^{2πi/3}` and its images, `E₆` at
`i` and its images, and j has a triple zero at `ρ`.

Close to a cusp (a rational point on the real axis), the reduced point has
`Im τ > 12` and `|j| ≈ e^{2π Im τ}` exceeds float range. There j equals `1/q` to
30 digits, so its phase `−2π Re τ` is kept and only the modulus is capped.

`H` (or `--disk`) maps the upper half-plane conformally onto the unit disk,
`τ = i(1 + w)/(1 − w)`, the Poincaré disk model: `w = 0` is `τ = i`, the boundary
circle is the real axis and `w = 1` is `τ = ∞`. The fundamental-domain tiling then
shows as ideal triangles accumulating at the rim. The wheel zooms about the
centre; the window title reports `τ` under the cursor.

`I` prints j(i) = 1728, j(i√2) = 8000 and j(ρ) = 0, j at `0.3 + 0.6i` against
the direct q-series (no reduction), and the invariance of j and the weight-4
transformation of `E₄` under `τ + 1`, `−1/τ` and `(2τ + 1)/(7τ + 4)`.
`--verify-gpu` compares the GPU's `E₄` and `E₆` with the direct series.

## Source files

- `common/app.h`: Cocoa window, input, command line, progressive rendering, GPU checks.
- `common/metal.h`: Metal device, pipelines, buffers, and readback.
- `common/shaders.h`: shared MSL frame: double-float and complex arithmetic,
  pixel coordinates, the summation kernel, and domain colouring.
- `<program>/<program>.mm`: each program's MSL term and display value, and `main`.
  `invariant_j/invariant_j.mm` also holds its CPU reference.
- `build.sh`: builds one or all programs.

The CPU-side math is shared with the Windows programs and included from
`../dx12`: `common/mobius.h`, `common/png.h`, `poincare/kleinian.h`,
`icosahedral/icosahedral.h` and `pearcey/pearcey.h`. Kantorovich's CPU code is
small and is repeated in `kantorovich/kantorovich.mm`; `invariant_j` has no
DirectX 12 counterpart.
