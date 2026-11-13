// Harmonic Newton example from the supplied Kantorovich theorem illustration.
#include "app.h"

static cd evaluate(cd z, cd* derivative = nullptr) {
    cd d=z*z*z-0.216;
    if (derivative) *derivative=-z*(z*z*z+0.432)/(d*d);
    return z-std::conj(z*z/d);
}
static void reference(cd z, cd out[2]) { out[0]=evaluate(z,&out[1]); }
static cd value(cd f, cd g, int mode, cd z) {
    if(mode==0) return f;
    if(mode==1) return -(f+std::conj(g)*std::conj(f))/(1.0-std::norm(g));
    for(int i=0;i<80;i++) {
        f=evaluate(z,&g);
        if(std::abs(f)<1e-6) return z;
        double j=1.0-std::norm(g);
        if(!std::isfinite(j)||std::abs(j)<1e-10) break;
        z-=(f+std::conj(g)*std::conj(f))/j;
    }
    return {NAN,NAN};
}
static void check() {
    double worst=0; int positive=0, negative=0;
    // Independent real finite differences validate the Wirtinger derivative.
    for(cd z : {cd(.3,.4),cd(-.8,.2),cd(1.2,-.7)}) {
        cd g; evaluate(z,&g); double e=1e-6;
        cd fx=(evaluate(z+e)-evaluate(z-e))/(2*e);
        cd fy=(evaluate(z+cd(0,e))-evaluate(z-cd(0,e)))/(2*e);
        worst=std::max(worst,std::abs(fx-(cd(1)-std::conj(g))));
        worst=std::max(worst,std::abs(fy-cd(0,1)*(cd(1)+std::conj(g))));
    }
    std::vector<cd> roots;
    for(int x=-15;x<=15;x++) for(int y=-15;y<=15;y++) {
        cd z(x*.1,y*.1), f,g; f=evaluate(z,&g);
        cd r=value(f,g,2,z);
        if(!std::isfinite(r.real()) || std::abs(evaluate(r))>1e-6) continue;
        bool seen=false; for(cd old:roots) if(std::abs(old-r)<1e-4) seen=true;
        if(!seen) roots.push_back(r);
    }
    for(cd r:roots) {cd g; evaluate(r,&g); if(std::norm(g)<1) positive++; else negative++;}
    std::printf("Derivative finite-difference error %.3g; roots %zu (%d preserving, %d reversing)\n",worst,roots.size(),positive,negative);
    if(worst>1e-7 || roots.size()!=10 || positive!=4 || negative!=6) fatal("Kantorovich mathematical check failed");
}
static std::vector<Mat> build(size_t) { check(); return {{1.,0.,0.,1.}}; }

static const char* HLSL=R"HLSL(
#define FADE_POLES 0
#ifdef COMPUTE
#ifdef FP64
#define VEC1 double
#else
#define VEC1 float
#endif
VEC4 term(VEC2 z,VEC2 unused,uint index) {
    VEC2 z2=cm(z,z), z3=cm(z2,z), d=z3-VEC2((VEC1)27/(VEC1)125,0);
    VEC2 r=cm(z2,crecip(d));
    VEC2 g=-cm(cm(z,z3+VEC2((VEC1)54/(VEC1)125,0)),crecip(cm(d,d)));
    return VEC4(z-VEC2(r.x,-r.y),g);
}
#else
double2 conjugate(double2 z) {return double2(z.x,-z.y);}
bool evaluate(double2 z,out double2 f,out double2 g,out double j) {
    double2 z2=cmd(z,z), z3=cmd(z2,z), d=z3-double2((double)27/(double)125,0);
    f=g=(double2)0; j=0;
    if (d.x*d.x+d.y*d.y<1e-20 || z.x*z.x+z.y*z.y>1e16) return false;
    f=z-conjugate(cdiv(z2,d));
    g=-cdiv(cmd(z,z3+double2((double)54/(double)125,0)),cmd(d,d));
    j=1.0-(g.x*g.x+g.y*g.y);
    return all(isfinite((float2)f)) && all(isfinite((float2)g));
}
double2 plotValue(double4 a,double2 z) {
    if(uFunc==0) return a.xy;
    double j=1.0-(a.z*a.z+a.w*a.w);
    if(abs(j)<1e-10) return (double2)0;
    return -(a.xy+cmd(conjugate(a.zw),conjugate(a.xy)))/j;
}
#define PLOT_OVERLAY 1
float3 plotOverlay(float3 col,double2 z) {
    double2 f,g; double j;
    if(!evaluate(z,f,g,j)) return (float3)0;
    // The orientation shading refers to f in both phase portraits.
    float contour=1-smoothstep(0.0,1.5,abs((float)j)/max(fwidth((float)j),1e-6));
    return lerp(col*(j>0 ? 1.0 : 0.48),float3(0.04,0.04,0.04),contour*.8);
}
#define PLOT_CUSTOM_COLOR 1
bool plotCustomColor(double2 z,out float3 col) {
    col=float3(.025,.03,.045);
    if(uFunc!=2) return false;
    [loop] for(int k=0;k<80;k++) {
        double2 f,g; double j;
        if(!evaluate(z,f,g,j)) return true;
        if(f.x*f.x+f.y*f.y<1e-12) {
            float r=length((float2)z);
            int id=0;
            if(r>.05) {
                float theta=atan2((float)z.y,(float)z.x);
                int sector=(int)floor((theta-(r<1 ? 1.0471975512 : 0))/2.0943951024+.5);
                sector=(sector%3+3)%3;
                id=1+sector+(r<.5 ? 0 : (r<1 ? 3 : 6));
            }
            float h=frac(id*.61803398875);
            float3 rgb=saturate(abs(frac(h+float3(0,2.0/3,1.0/3))*6-3)-1);
            col=lerp(float3(1,1,1),rgb,.78)*(.45+.55*exp(-.025*k));
            return true;
        }
        if(abs(j)<1e-10) return true;
        z-=(f+cmd(conjugate(g),conjugate(f)))/j;
    }
    return true;
}
#endif
)HLSL";
int main(int argc,char** argv) {
    if(argc==2 && std::string(argv[1])=="--self-test") {check();return 0;}
    const Group group={"Kantorovich: harmonic Newton iteration",
        "f(z) = z - conjugate(z^2 / (z^3 - 0.6^3))",
        "kantorovich.exe","kantorovich",{"f","N(z)-z","Newton basins"},{"f","g prime"},
        0,0,3,0,HLSL,build,reference,value,check};
    return runApp(group,argc,argv);
}

