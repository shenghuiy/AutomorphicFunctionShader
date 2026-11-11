# GPU evaluation

[← Quadrature](quadrature.md) · [overview](README.md) · next: [Verification](verification.md)

The HLSL in `PEARCEY_HLSL` (`pearcey.cpp`) is spliced into the shared frame of
`common/shaders.h`. It supplies two functions. `term` runs in the compute shader
and returns one node's contribution. `plotValue` and `plotOverlay` run in the
pixel shader and turn the accumulated sums into a colour.

## The summation loop

The shared compute shader runs one thread per pixel in 16×16 groups. It loads
the nodes into group-shared memory 256 at a time, and every thread then loops
over the tile:

```hlsl
for (uint base = uStart; base < uEnd; base += 256) {
    ...                                   // load tileM from mats[base ...]
    if (inside) {
        for (uint i=0; i<n; ++i) {
            VEC4 ab=tileM[2*i], cd=tileM[2*i+1];
            VEC2 w1=cm(ab.xy,z)+ab.zw, w2=cm(cd.xy,z)+cd.zw;
            sum += term(w1, w2, base+i);
        }
    }
}
// acc[index] = sum on the first dispatch, acc[index] += sum after that
```

`VEC2`/`VEC4`/`REAL` are `float` or `double` types depending on the `D` key
(`FP64` defined or not). The accumulator buffer `acc` is always `double4`. Each
dispatch covers the nodes `[uStart, uEnd)`. Its fp32 or fp64 sum is converted to
double and added to the pixel's running total.

## The integrand at one node

`term` evaluates the same formulas as `pearcey::psi` ([quadrature.md](quadrature.md#the-midpoint-rule)).
With $`\lambda = \tfrac12`$ fixed, it folds the constant into the code:
$`\kappa = 16S\lambda^2 = 4S`$, $`\lambda g/q = g/(2q)`$, and
$`\lambda g'/q^3 = g'/(2q^3)`$.

```hlsl
// Path t(s) = s + i LAMBDA g / sqrt(1 + g^2), g = phi'(s) / kappa; see pearcey.h.
VEC4 term(VEC2 z, VEC2 node, uint index) {
    REAL x=z.x, y=z.y, one=(REAL)num(1);
    REAL S=(REAL)(3.0+max(sqrt(abs((float)y)),pow(0.5*abs((float)x),1.0/3.0)));  // pearcey::halfWidth
    REAL s=S*node.x, s2=s*s, y2=y+y;
    REAL kappa=(REAL)num(4)*S;  // 16 S LAMBDA^2 with LAMBDA = 1/2
    REAL g=((REAL)num(4)*s2*s+y2*s+x)/kappa, dg=((REAL)num(12)*s2+y2)/kappa;
    REAL q=dsqrt(one+g*g), q3=q*q*q;
    VEC2 t=VEC2(s,g/(q+q)), dt=VEC2(one,dg/(q3+q3));
    VEC2 t2=cm(t,t);
    VEC2 v=cm(expi(cm(t2,t2)+y*t2+x*t),dt)*(S*node.y);
    return VEC4(v,cm(VEC2(-t.y,t.x),v));
}
```

Line by line:

| Code | Math |
|---|---|
| `s=S*node.x` | $`s_j = S u_j`$ |
| `g`, `dg` | $`g = \varphi'(s)/\kappa`$, $`g' = \varphi''(s)/\kappa`$ |
| `t=VEC2(s,g/(q+q))` | $`t = s + i\lambda g/q`$ |
| `dt=VEC2(one,dg/(q3+q3))` | $`t'(s) = 1 + i\lambda g'/q^3`$ |
| `cm(t2,t2)+y*t2+x*t` | $`\varphi(t) = t^4 + yt^2 + xt`$ |
| `expi(p)` | $`e^{ip}`$ for complex $`p`$ |
| `*(S*node.y)` | weight $`h = S\cdot 2/n`$ |
| `cm(VEC2(-t.y,t.x),v)` | $`i\,t\,v`$, the $`\partial\Psi/\partial x`$ integrand |

The returned `VEC4` is $`(\Psi\text{ term}, \partial_x\Psi\text{ term})`$ as two complex
numbers, `acc.xy` and `acc.zw`. `cm` is the complex multiply from `shaders.h`.

## fp32 and fp64

In fp32 mode `expi` is direct:

```hlsl
float2 expi(float2 p) { float s, c; sincos(p.x,s,c); return float2(c,s)*exp(-p.y); }
```

since $`e^{i(p_x + i p_y)} = e^{-p_y}(\cos p_x + i\sin p_x)`$. On the path,
$`p_y = \operatorname{Im}\varphi \ge 0`$, so every term has modulus at most $`\lvert t'\rvert h`$,
and the sum has no large cancellation ([contour.md](contour.md#why-the-integrand-decays-on-it)).
The fp32 error is then set mainly by the phase. $`\operatorname{Re}\varphi`$ reaches
a few thousand near the ends of the interval, and a float carries it with an
absolute error of about $`10^{-4}`$. The picture is right, and `--verify-gpu`
accepts fp32 up to a relative error of $`5\times10^{-3}`$.

fp64 mode needs `exp`, `sin`, `cos` and `sqrt` in double. Shader model 5 has
none of them, so `pearcey.cpp` builds them from additions, multiplications and
divisions.

### Square root

The float result as a starting point, then two Newton steps for $`r^2 = d`$. Each
step roughly doubles the number of correct digits, $`24 \to 48 \to 96`$ bits:

```hlsl
double dsqrt(double d) {  // float estimate, two Newton steps
    double r=(double)sqrt((float)d);
    r=(r*r+d)/(r+r);
    r=(r*r+d)/(r+r);
    return r;
}
```

### Exponential

Argument reduction, $`x = k\ln 2 + r`$ with integer $`k`$ and $`\lvert r\rvert \le \tfrac12\ln 2 \approx 0.35`$,
then a Taylor polynomial in Horner form:

```math
e^x = 2^k e^r,
\qquad
e^r \approx 1 + r\Big(1 + \frac r2\Big(1 + \frac r3\Big(\cdots\Big(1 + \frac r{13}\Big)\Big)\Big)\Big).
```

The truncation error is about $`0.35^{14}/14! \approx 4\times10^{-18}`$, below double
rounding.

```hlsl
double dexp(double x) {  // x = k ln 2 + r, |r| <= 0.35, Taylor series in r
    float k=round(max((float)x,-80.0)*1.44269504);
    double r=x-(double)k*split(0.6931472,-1.9046542e-09,-1.110223e-16);
    double one=(double)num(1), p=one;
    for (uint n=13; n>=1; --n) p=one+r*p/(double)num(n);
    return (float)x < -80.0 ? one-one : p*(double)exp2(k);
}
```

- $`\ln 2`$ is given as the sum of three floats (`split`). Their sum is $`\ln 2`$ to
  well beyond double precision, so $`r = x - k\ln 2`$ stays accurate. This is a
  Cody–Waite style reduction.
- $`2^k`$ comes from the float `exp2`. That is exact for an integer $`k`$ in the
  normal float range.
- Arguments below $`-80`$ return 0. $`e^{-80} \approx 2\times10^{-35}`$ is negligible
  next to the terms near the stationary points, which have size of order $`h`$.
  The cut also keeps $`2^k`$ inside the float range.

### Sine and cosine

`expi` reduces $`p_x`$ modulo $`\pi/2`$, $`p_x = k\pi/2 + r`$ with $`\lvert r\rvert \le \pi/4`$.
It then evaluates both Taylor series in Horner form up to $`r^{16}`$ and $`r^{17}`$:

```math
\cos r \approx 1 - \frac{r^2}{2\cdot1}\Big(1 - \frac{r^2}{4\cdot3}\Big(\cdots\Big(1 - \frac{r^2}{16\cdot15}\Big)\Big)\Big),
\qquad
\sin r \approx r\Big(1 - \frac{r^2}{3\cdot2}\Big(\cdots\Big(1 - \frac{r^2}{17\cdot16}\Big)\Big)\Big).
```

The truncation error is about $`(\pi/4)^{18}/18! \approx 2\times10^{-18}`$. The quadrant
$`k \bmod 4`$ then rotates $`(\cos r, \sin r)`$ into place:
$`e^{i(k\pi/2 + r)} = i^k e^{ir}`$.

```hlsl
double2 expi(double2 p) {
    // p.x = k pi/2 + r, |r| <= pi/4, Taylor series in r
    float k=round((float)p.x*0.636619772);
    double r=p.x-(double)k*split(1.5707964,-4.371139e-08,-1.7763568e-15), r2=r*r;
    double one=(double)num(1), c=one, s=one;
    for (uint n=16; n>=2; n-=2) {
        c=one-r2*c/(double)num(n*(n-1));
        s=one-r2*s/(double)num(n*(n+1));
    }
    s*=r;
    int q=((int)k)&3;
    double2 e=q==0 ? double2(c,s) : (q==1 ? double2(-s,c) : (q==2 ? double2(-c,-s) : double2(s,-c)));
    return e*dexp(-p.y);
}
```

$`k`$ is found in float. If rounding picks a neighbouring $`k`$, $`\lvert r\rvert`$ is a little
larger than $`\pi/4`$, which the series still handles. The phase is at most a few
thousand, so $`k`$ fits exactly in a float.

### No double constants

Each constant in the fp64 code is built at run time from floats:

```hlsl
// Small integer n as a runtime value. On at least one NVIDIA driver double
// immediates read as zero in some instructions (see pixelZ in shaders.h), which
// here made the fp64 sums NaN, so the fp64 code has no double constants at all.
uint num(uint n) { return n*min(uW,1u); }
...
// A double constant as a sum of three floats, each scaled by the runtime 1.
double split(float a, float b, float c) {
    float one=(float)num(1);
    return (double)(a*one)+(double)(b*one)+(double)(c*one);
}
```

`uW` is the window width from the constant buffer, so `min(uW,1u)` is 1, but the
compiler cannot fold it. This is why the shader writes `(REAL)num(4)` instead of
`4.0`, and `one-one` instead of `0.0`.

## Display

The pixel shader reads `acc` and calls `plotValue`. That returns $`\Psi`$ or
$`\partial\Psi/\partial x`$, scaled up for a partial sum
([quadrature.md](quadrature.md#bit-reversed-node-order-and-progressive-rendering)).
The shared domain colouring in `shaders.h` then maps the phase to a colour and
draws the modulus contours. `pearcey.exe` starts with
Mathematica's `GreenPinkTones` palette (`palette = 1` in `main`). It runs from dark
green at $`\arg = -\pi`$ through near white at 0 to dark magenta at $`+\pi`$. It is not
cyclic, so it has a seam on the negative real axis. `C` switches to the hue wheel.

## The caustic overlay

With `F = 1` the pixel shader darkens the pixels near the caustic
$`c(x, y) = 27x^2 + 8y^3 = 0`$ ([definition.md](definition.md#the-cusp-caustic)). The
distance from $`(x, y)`$ to the zero set of $`c`$ is, to first order,

```math
d \approx \frac{\lvert c\rvert}{\lvert\nabla c\rvert},
\qquad
\nabla c = (54x,\; 24y^2).
```

Dividing by the pixel size `uScale` gives the distance in pixels. A `smoothstep`
from 1 to 2 pixels draws an antialiased line about 1.5 pixels wide:

```hlsl
float3 plotOverlay(float3 col, double2 zd) {
    if (uFunc!=1) return col;
    float2 z=(float2)zd;
    float c=27.0*z.x*z.x+8.0*z.y*z.y*z.y;
    float d=abs(c)/max(length(float2(54.0*z.x,24.0*z.y*z.y)),1e-30);
    return col*smoothstep(1.0,2.0,d/(float)uScale);
}
```

At the cusp point both $`c`$ and $`\nabla c`$ vanish. The `max(..., 1e-30)` keeps the
division finite, and the pixel there is drawn dark, as it should be.
