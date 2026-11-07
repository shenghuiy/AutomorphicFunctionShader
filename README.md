# Complex Plot Shader

Interactive domain coloring of the automorphic function
`f(z) = P_{-1,-3}(z) / P_{0,-4}(z)` for a two-generator Kleinian group.
Both renderers approximate the Poincaré series with a finite number of group elements.

## Native Windows renderer

Requires Visual Studio or Build Tools with the MSVC x64 C++ toolchain, and a GPU supporting OpenGL 4.3.

Run `native\build.bat`, then launch `native\poincare.exe`.
Group elements are enumerated on the CPU; a compute shader sums the series per pixel.

Controls:

- Drag to pan; scroll to zoom.
- `+` / `-`: double / halve the number of terms.
- `D`: toggle fp32 / fp64 summation.
- `M`: cycle plot styles; `F`: cycle displayed functions.
- `[` / `]`: adjust contour density; `R`: reset the view.
- `I`: print an invariance check; `P`: save a PNG.
- `Esc`: quit.

## Browser renderer

Open `index.html` in a browser supporting WebGL2 and the `EXT_color_buffer_float` extension.
Use the on-screen controls to adjust the plot.

## Source files

- `native/poincare.cpp`: Win32/OpenGL renderer and group enumeration.
- `native/build.bat`: native build script.
- `index.html`: standalone browser implementation.
