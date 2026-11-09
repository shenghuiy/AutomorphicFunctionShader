# Complex Plot Shader

Interactive domain coloring of the automorphic function
`f(z) = P_{-1,-3}(z) / P_{0,-4}(z)` for a two-generator Kleinian group.
Both renderers approximate the Poincaré series with a finite number of group elements.

![Domain coloring of the automorphic function rendered with 20,000 group elements](docs/poincare.png)

DirectX 12 renderer preview using 20,000 terms at 1400 × 1000 pixels.

See [docs/background.md](docs/background.md) for the math: why `f` is invariant
under the group, how the series is truncated, and how to read the colors.

## Icosahedral group

Both renderers can also plot an automorphic function for the 60 rotations of the
icosahedron acting on the Riemann sphere:

`F(z) = Σ_{g∈G} h(g(z))`, with `h(z) = (z − 0.05)⁻⁵`.

The group is finite, so the sum is exact. The default view is `|Re z|, |Im z| < 2.5`.
Choose **Icosahedral** in the browser's Group menu, or run the separate native program
`native\icosahedral\icosahedral.exe`.
Its Function choices are `F`, the single term `h`, and Klein's invariant `J`.

The icosahedral plot can also be shown on the Riemann sphere. Choose
**View → Riemann sphere** in the browser, or press `S` in either the browser or
`icosahedral.exe` (or start it with `--sphere`). The sphere uses stereographic
projection from the south pole, `z = (X + iY)/(1 + Z)`, so 0 is the north pole and
∞ the south pole. Seen from outside, this keeps the plane's orientation rather than
mirroring it. Switching views keeps the point at the centre of the window. On the
sphere, drag rotates it and the scroll wheel zooms.

![Phase portrait of the icosahedral automorphic function](docs/icosahedral.png)

See [docs/coxeter.md](docs/coxeter.md) for Coxeter and triangle groups, Klein's
icosahedral invariants, and how this finite case relates to automorphic functions
on other Riemann surfaces.

## Native Windows renderer

Requires Visual Studio or Build Tools with the MSVC x64 C++ toolchain, the Windows SDK, and a DirectX 12 GPU with double-precision shader support (Windows 10/11).

There are two programs, one per group:

| Program | Group | Default center |
|---|---|---|
| `native\poincare\poincare.exe` | Kleinian ⟨g₁, g₂⟩, truncated Poincaré series | `0 -0.35` |
| `native\icosahedral\icosahedral.exe` | 60 icosahedral rotations, exact sum | `0 0` |

Run `native\build.bat` to build both, or `native\build.bat poincare` /
`native\build.bat icosahedral` to build one.
Group elements are enumerated on the CPU; HLSL compute shaders sum the series per pixel.
DirectX 12 handles GPU buffers, domain coloring, and presentation. Both summation modes
accumulate into a double-precision buffer.

Close the running executable before rebuilding, or pass an alternate output name,
for example `native\build.bat poincare poincare-dx12.exe`.

Optionally set the initial view from the command line:

```powershell
.\native\poincare\poincare.exe --center 0 -0.35 --span 5
.\native\poincare\poincare.exe --center 0.2 -0.6 --span 1.5 --terms 20000 --fp64 --snapshot zoom.png
.\native\icosahedral\icosahedral.exe --snapshot icosahedral.png
.\native\icosahedral\icosahedral.exe --sphere --center 1.5 0.8 --snapshot sphere.png
```

`--center <real> <imag>` sets the complex center.
`--span <value>` sets the extent of the shorter window dimension in complex-plane
units (default: `5`); it must be finite and positive. The other dimension follows
the window aspect ratio. Both options are optional. `--terms` exists only for
`poincare.exe`, and `--sphere` (start in the Riemann sphere view, with `--center`
facing the viewer) only for `icosahedral.exe`. Run either program with `--help` for all options.

Controls (both programs):

- Drag to pan; scroll to zoom.
- `+` / `-`: double / halve the number of terms (`poincare.exe` only).
- `D`: toggle fp32 / fp64 summation.
- `M`: cycle plot styles; `F`: cycle displayed functions.
- `[` / `]`: adjust contour density; `R`: reset the view.
- `S`: toggle complex plane / Riemann sphere (`icosahedral.exe` only). On the
  sphere, drag rotates and scroll zooms.
- `I`: print an invariance check; `P`: save a PNG.
- `Esc`: quit.

To validate GPU sums against the CPU evaluator at three sample pixels:

```powershell
.\native\poincare\poincare.exe --terms 20000 --verify-gpu --snapshot check.png
.\native\poincare\poincare.exe --terms 1024 --fp64 --verify-gpu --snapshot check-fp64.png
.\native\icosahedral\icosahedral.exe --fp64 --verify-gpu --snapshot check-icosahedral.png
```

## Browser renderer

Open `index.html` in a browser supporting WebGL2 and the `EXT_color_buffer_float` extension.
Keep `complex.js` and `icosahedral.js` in the same folder; the page loads them.
Use the on-screen controls to adjust the plot.

## Source files

- `native/poincare/`: the Kleinian program.
  - `kleinian.h`: generators, reduced-word enumeration, and the CPU Poincaré series.
  - `poincare.cpp`: the group's HLSL term and display value, and `main`.
- `native/icosahedral/`: the icosahedral program.
  - `icosahedral.h`: the 60 rotations, `F`, `h`, and Klein's `J` on the CPU.
  - `icosahedral.cpp`: the group's HLSL term and display value, and `main`.
- `native/common/`: code shared by both programs.
  - `app.h`: Win32 window, controls, command line, progressive rendering, GPU checks.
  - `shaders.h`: shared HLSL frame (summation loop and domain coloring) around each group's functions.
  - `dx12.h`: DirectX 12 device, pipelines, resource management, and readback.
  - `mobius.h`: SL(2,C) matrices and Möbius transformations.
  - `png.h`: dependency-free PNG writer.
- `native/build.bat`: builds one or both programs.
- `index.html`: browser implementation (WebGL2 shaders, UI, Kleinian group).
- `complex.js`: complex numbers and Möbius matrices for the browser, shared by both groups.
- `icosahedral.js`: icosahedral group, `F`, `h`, Klein's `J`, and their GLSL snippets for the browser.
- `docs/background.md`: mathematical background.
- `docs/coxeter.md`: Coxeter groups and the icosahedral automorphic function.
