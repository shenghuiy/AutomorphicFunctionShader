#pragma once
// Shared HLSL frame, compiled at startup with the Windows SDK's HLSL compiler.
//
// shaderSource(groupHlsl) returns SHADER_HEAD + groupHlsl + SHADER_TAIL. The
// group's HLSL must define, using the helpers in SHADER_HEAD:
//
//   #define FADE_POLES 0 or 1          darken near |value| = 2^±14..24
//   #ifdef COMPUTE
//   VEC4 term(VEC2 w1, VEC2 w2, uint index)
//       Contribution of group element `index` to the per-pixel accumulator,
//       with w1 = a z + b and w2 = c z + d from its uploaded matrix.
//   #else
//   double2 plotValue(double4 a, double2 z)
//       Value to colour from the accumulated sums `a` at point z (uses uFunc).
//   #endif

#include <string>

static const char* SHADER_HEAD = R"HLSL(
cbuffer View : register(b0) {
    uint uStart, uEnd, uW, uH;
    double2 uCenter;
    double uScale;
    uint uMode, uFunc;
    float uDens;
};

#ifdef COMPUTE
StructuredBuffer<VEC4> mats : register(t0);
RWStructuredBuffer<double4> acc : register(u0);
groupshared VEC4 tileM[512];

VEC2 cm(VEC2 a, VEC2 b) {
    return VEC2(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);
}
VEC2 crecip(VEC2 q) {
#ifdef FP64
    // Divide directly: on at least one NVIDIA driver, `* (1.0L / |q|^2)` read
    // back as zero in the icosahedral shader, though fxc's output looked correct.
    return VEC2(q.x, -q.y) / (q.x*q.x + q.y*q.y);
#else
    float rs = 1.0 / max(abs(q.x), abs(q.y));
    float2 t = q * rs;
    return float2(t.x, -t.y) * (rs / dot(t,t));
#endif
}
#else
StructuredBuffer<double4> acc : register(t1);
double2 cmd(double2 a,double2 b) { return double2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x); }
double2 cdiv(double2 a,double2 b) {
    double s=max(abs(b.x),abs(b.y));
    a/=s; b/=s;
    return double2(a.x*b.x+a.y*b.y,a.y*b.x-a.x*b.y)/(b.x*b.x+b.y*b.y);
}
#endif
)HLSL";

static const char* SHADER_TAIL = R"HLSL(
#ifdef COMPUTE
[numthreads(16,16,1)]
void compute(uint3 p : SV_DispatchThreadID, uint lid : SV_GroupIndex) {
    bool inside = p.x < uW && p.y < uH;
    // Pixel offsets are exact in float at supported window sizes; transform
    // to the complex plane in double precision to retain deep zoom accuracy.
    float2 offset = float2((float)p.x+0.5-0.5*uW, 0.5*uH-0.5-(float)p.y);
    VEC2 z = (VEC2)(uCenter + (double2)offset*uScale);
    VEC4 sum = (VEC4)0;
    for (uint base = uStart; base < uEnd; base += 256) {
        uint n = min(256u, uEnd-base);
        if (lid < n) {
            tileM[2*lid] = mats[2*(base+lid)];
            tileM[2*lid+1] = mats[2*(base+lid)+1];
        }
        GroupMemoryBarrierWithGroupSync();
        if (inside) {
            for (uint i=0; i<n; ++i) {
                VEC4 ab=tileM[2*i], cd=tileM[2*i+1];
                VEC2 w1=cm(ab.xy,z)+ab.zw, w2=cm(cd.xy,z)+cd.zw;
                sum += term(w1, w2, base+i);
            }
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (inside) {
        uint index=p.y*uW+p.x;
        if (uStart==0) acc[index]=(double4)sum;
        else acc[index]+=(double4)sum;
    }
}
#else
float4 vertex(uint id : SV_VertexID) : SV_Position {
    float2 p=float2((id<<1)&2,id&2);
    return float4(p*2.0-1.0,0,1);
}
float3 hsv(float h) {
    float3 k=saturate(abs((frac((h*6.0+float3(0,4,2))/6.0)*6.0)-3.0)-1.0);
    return lerp(float3(1,1,1),k*k*(3.0-2.0*k),0.9);
}
float saw(float x) { return lerp(frac(x),0.5,smoothstep(0.3,0.8,fwidth(x))); }
float gridline(float x) {
    return 1.0-smoothstep(0.0,1.2,abs(frac(x-0.5)-0.5)/max(fwidth(x),1e-6));
}
float4 pixel(float4 position : SV_Position) : SV_Target {
    // Clamp: quads on the bottom/right edge run helper lanes (needed for
    // fwidth) one pixel outside the target, which would read past acc.
    uint2 p=min((uint2)position.xy,uint2(uW-1,uH-1));
    float2 offset=float2((float)p.x+0.5-0.5*uW,0.5*uH-0.5-(float)p.y);
    double2 fd=plotValue(acc[p.y*uW+p.x],uCenter+(double2)offset*uScale);
    double s=max(abs(fd.x),abs(fd.y));
    if (!(s>0.0L)) return float4(0,0,0,1);
    float2 dir=(float2)(fd/s);
    if (any(isnan(dir)) || any(isinf(dir))) return float4(0,0,0,1);
    float lm=log2((float)s)+0.5*log2(dot(dir,dir));
    if (isinf((float)s)) lm=128.0;
    float h=frac(atan2(dir.y,dir.x)/(2.0*3.14159265358979));
    float3 base=hsv(h), col;
    if (uMode==0) col=base*(0.55+0.45*sqrt(saw(lm*uDens)*saw(h*12.0*uDens)));
    else if (uMode==1) col=base;
    else if (uMode==2) col=base*(0.45+0.55*saw(lm*uDens));
    else {
        float2 fs=dir*exp2(log2((float)s)-floor(lm))*4.0*uDens;
        col=lerp(base*0.95,float3(0.05,0.05,0.05),0.75*max(gridline(fs.x),gridline(fs.y)));
    }
#if FADE_POLES
    col*=1.0-smoothstep(14.0,24.0,abs(lm));
#endif
    return float4(pow(abs(col),1.0/1.1),1);
}
#endif
)HLSL";

inline std::string shaderSource(const char* groupHlsl) {
    return std::string(SHADER_HEAD) + groupHlsl + SHADER_TAIL;
}
