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
Choose **Icosahedral** in the browser's Group menu, press `G` in the native renderer,
or start it with `--group icosahedral`.
In this mode the Function choices are `F`, the single term `h`, and Klein's invariant `J`.

![Phase portrait of the icosahedral automorphic function](docs/icosahedral.png)

See [docs/coxeter.md](docs/coxeter.md) for Coxeter and triangle groups, Klein's
icosahedral invariants, and how this finite case relates to automorphic functions
on other Riemann surfaces.

## Native Windows renderer

Requires Visual Studio or Build Tools with the MSVC x64 C++ toolchain, the Windows SDK, and a DirectX 12 GPU with double-precision shader support (Windows 10/11).

Run `native\build.bat`, then launch `native\poincare.exe`.
Group elements are enumerated on the CPU; HLSL compute shaders sum the series per pixel.
DirectX 12 handles GPU buffers, domain coloring, and presentation. Both summation modes
accumulate into a double-precision buffer.

Close the running executable before rebuilding, or use `native\build.bat poincare-dx12.exe`
to build an alternate executable.

Optionally set the initial view from the command line:

```powershell
.\native\poincare.exe --center 0 -0.35 --span 5
.\native\poincare.exe --center 0.2 -0.6 --span 1.5 --terms 20000 --fp64 --snapshot zoom.png
.\native\poincare.exe --group icosahedral --snapshot icosahedral.png
```

`--group kleinian|icosahedral` selects the group (default: `kleinian`). The
icosahedral group defaults to center `0 0`.

`--center <real> <imag>` sets the complex center (default: `0 -0.35`).
`--span <value>` sets the extent of the shorter window dimension in complex-plane
units (default: `5`); it must be finite and positive. The other dimension follows
the window aspect ratio. Both options are optional. Run with `--help` for all options.

Controls:

- Drag to pan; scroll to zoom.
- `+` / `-`: double / halve the number of terms (Kleinian group only).
- `G`: toggle between the Kleinian and icosahedral groups.
- `D`: toggle fp32 / fp64 summation.
- `M`: cycle plot styles; `F`: cycle displayed functions.
- `[` / `]`: adjust contour density; `R`: reset the view.
- `I`: print an invariance check; `P`: save a PNG.
- `Esc`: quit.

To validate GPU sums against the CPU evaluator at three sample pixels:

```powershell
.\native\poincare.exe --terms 20000 --verify-gpu --snapshot check.png
.\native\poincare.exe --terms 1024 --fp64 --verify-gpu --snapshot check-fp64.png
```

## Browser renderer

Open `index.html` in a browser supporting WebGL2 and the `EXT_color_buffer_float` extension.
Keep `complex.js` and `icosahedral.js` in the same folder; the page loads them.
Use the on-screen controls to adjust the plot.

## Source files

- `native/poincare.cpp`: Win32 application, Kleinian group enumeration, and PNG export.
- `native/mobius.h`: SL(2,C) matrices and Möbius transformations, shared by both groups.
- `native/icosahedral.h`: icosahedral group, `F`, `h`, Klein's `J`, and their HLSL snippets.
- `native/dx12.h`: DirectX 12 device, pipelines, resource management, and readback.
- `native/shaders.h`: embedded HLSL compute and display shaders.
- `native/build.bat`: native build script.
- `index.html`: browser implementation (WebGL2 shaders, UI, Kleinian group).
- `complex.js`: complex numbers and Möbius matrices for the browser, shared by both groups.
- `icosahedral.js`: icosahedral group, `F`, `h`, Klein's `J`, and their GLSL snippets for the browser.
- `docs/background.md`: mathematical background.
- `docs/coxeter.md`: Coxeter groups and the icosahedral automorphic function.
