# Metal renderer (macOS)

Native Metal versions of the four programs of the DirectX 12 renderer in
[`../dx12/`](../dx12/): the same plots, controls and command-line options.

| Program | Plot | Default center |
|---|---|---|
| `native/build/poincare` | Kleinian ⟨g₁, g₂⟩, truncated Poincaré series | `0 -0.35` |
| `native/build/icosahedral` | 60 icosahedral rotations, exact sum | `0 0` |
| `native/build/pearcey` | Pearcey integral, quadrature per pixel | `0 -5` |
| `native/build/kantorovich` | Harmonic function and Newton basins | `0 0` |

## Build

Requires a Mac with a Metal GPU and Apple clang from Xcode or the Command Line
Tools. Tested on an Apple M1 Max with macOS 27 and Apple clang 21. The shaders are compiled at
startup by the Metal framework, so the offline `metal` compiler is not needed.

```sh
native/metal/build.sh               # all four programs
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

## Source files

- `common/app.h`: Cocoa window, input, command line, progressive rendering, GPU checks.
- `common/metal.h`: Metal device, pipelines, buffers, and readback.
- `common/shaders.h`: shared MSL frame: double-float and complex arithmetic,
  pixel coordinates, the summation kernel, and domain colouring.
- `<program>/<program>.mm`: each program's MSL term and display value, and `main`.
- `build.sh`: builds one or all programs.

The CPU-side math is shared with the Windows programs and included from
`../dx12`: `common/mobius.h`, `common/png.h`, `poincare/kleinian.h`,
`icosahedral/icosahedral.h` and `pearcey/pearcey.h`. Kantorovich's CPU code is
small and is repeated in `kantorovich/kantorovich.mm`.
