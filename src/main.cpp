#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <glm/glm.hpp>
#include <vector>
#include <thread>
#include <initializer_list>
#include <chrono>
#include <iostream>

using namespace glm;

struct Ray
{
    vec3 o;
    vec3 d;
    Ray(vec3 o_, vec3 d_)
    : o(o_)
    , d(d_)
    {}
};

enum Refl_t
{
    DIFF,
    SPEC,
    REFR
};  // material types, used in radiance()

struct Sphere {
    float rad;       // radius
    vec3 p, e, c;      // position, emission, color
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    
    Sphere(float rad_, vec3 p_, vec3 e_, vec3 c_, Refl_t refl_)
    :rad(rad_)
    , p(p_)
    , e(e_)
    , c(c_)
    , refl(refl_)
    {}
    
    float intersect(const Ray &r) const
    { // returns distance, 0 if nohit
        vec3 op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        float t, eps=1e-4, b=dot(op, r.d), det=b*b-dot(op, op)+rad*rad;
        if (det < 0)
            return 0;
        else
            det = ::sqrtf(det);
        return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
    }
};
Sphere spheres[] = {//Scene: radius, position, emission, color, material
    Sphere(1e3, vec3( 1e3+1,40.8,81.6), vec3(0),vec3(.75,.25,.25),DIFF),//Left
    Sphere(1e3, vec3(-1e3+99,40.8,81.6),vec3(0),vec3(.25,.25,.75),DIFF),//Rght
    Sphere(1e3, vec3(50,40.8, 1e3),     vec3(0),vec3(.75,.75,.75),DIFF),//Back
    Sphere(1e3, vec3(50,40.8,-1e3+170), vec3(0),vec3(),           DIFF),//Frnt
    Sphere(1e3, vec3(50, 1e3, 81.6),    vec3(0),vec3(.75,.75,.75),DIFF),//Botm
    Sphere(1e3, vec3(50,-1e3+81.6,81.6),vec3(0),vec3(.75,.75,.75),DIFF),//Top
    Sphere(16.5,vec3(27,16.5,47),       vec3(0),vec3(1,1,1)*0.999f, SPEC),//Mirr
    Sphere(16.5,vec3(73,16.5,78),       vec3(0),vec3(1,1,1)*0.999f, REFR),//Glas
    Sphere(2.5, vec3(50, 81.6 - 6.0, 81.6),vec3(4,4,4)*50.0f,  vec3(), DIFF),//Lite
};

int numSpheres = sizeof(spheres) / sizeof(Sphere);

inline double clamp(double x)
{
    return x<0 ? 0 : x>1 ? 1 : x;
}

inline int toInt(double x)
{
    return int(pow(clamp(x),1/2.2)*255+.5);
}

inline bool intersect(const Ray &r, float &t, int &id)
{
    float inf = std::numeric_limits<float>::max();
    t = inf;
    for(int i = 0; i < numSpheres; ++i)
    {
        float d = spheres[i].intersect(r);
        if(d != 0.0 && d < t)
        {
            t = d;
            id = i;
        }
    }
    return t < inf;
}

vec3 radiance(const Ray &r, int depth, unsigned short *Xi,int E=1)
{
    float t;                               // distance to intersection
    int id = 0;                               // id of intersected object
    if (depth > 10)
        return vec3();
    if (!intersect(r, t, id))
        return vec3(); // if miss, return black
    const Sphere &obj = spheres[id];        // the hit object
    vec3 rayIntersectionPoint = r.o + r.d * t;
    vec3 sphereNormal = normalize(rayIntersectionPoint - obj.p);
    vec3 orientedNormal = dot(sphereNormal, r.d) < 0 ? sphereNormal : -sphereNormal;
    vec3 f = obj.c;
    float p = (f.x > f.y && f.x > f.z) ? f.x :(f.y > f.z ? f.y : f.z); // max refl
    if (++depth>5||!p)
    {
        if (erand48(Xi)<p)
            f=f*(1.0f/p);
        else
            return obj.e * (float)E;
    }
    
    if (obj.refl == DIFF)
    {
        float angleAround = 2 * M_PI * erand48(Xi);
        float distFromCenter = erand48(Xi);
        float distFromCentersq = ::sqrtf(distFromCenter);
        vec3 w = orientedNormal;
        vec3 u = normalize(cross(fabsf(w.x) > 0.1 ? vec3(0, 1, 0) : vec3(1, 0, 0), w));
        vec3 v = cross(w, u);
        vec3 randomReflection = normalize(u * cosf(angleAround) * distFromCentersq + v * sinf(angleAround) * distFromCentersq + w * sqrtf(1.0f - distFromCenter));
        
        vec3 e;
        for (int i = numSpheres - 1; i < numSpheres; ++i)
        {
            const Sphere& s = spheres[i];
            if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0)
                continue; // skip no light
            vec3 sw = s.p - rayIntersectionPoint;
            vec3 su = normalize(cross(fabs(sw.x) > 0.1f ? vec3(0.0f, 1.0f, 0.0f) : vec3(1.0f, 0.0f, 0.0f), sw));
            vec3 sv = cross(sw, su);
            float cos_a_max = sqrtf(1.0f - s.rad * s.rad / dot(rayIntersectionPoint - s.p, rayIntersectionPoint - s.p));
            float eps1 = erand48(Xi), eps2 = erand48(Xi);
            float cos_a = 1.0f - eps1 + eps1 * cos_a_max;
            float sin_a = sqrtf(1.0f - cos_a * cos_a);
            float phi = 2 * M_PI * eps2;
            vec3 l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
            l = normalize(l);
            if (intersect(Ray(rayIntersectionPoint, l), t, id) && id==i)
            {  // shadow ray
                float omega = 2.0f * M_PI * (1.0f - cos_a_max);
                e = e + f * (s.e * dot(l, orientedNormal) * omega) * (float) M_1_PI;  // 1/pi for brdf
            }
        }
//        return f;
        return obj.e * (float) E + e + f * (radiance(Ray(rayIntersectionPoint, randomReflection), depth, Xi, 0));
    }
    else if (obj.refl == SPEC)              // Ideal SPECULAR reflection
    {
        return obj.e + f * (radiance(Ray(rayIntersectionPoint, r.d - sphereNormal * 2.0f * dot(sphereNormal, r.d)),depth, Xi));
    }
    
    Ray reflRay(rayIntersectionPoint, r.d- sphereNormal * 2.0f * dot(sphereNormal, r.d));     // Ideal dielectric REFRACTION
    bool into = dot(sphereNormal, orientedNormal) > 0;                // Ray from outside going in?
    float nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=dot(r.d, orientedNormal), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj.e + f * (radiance(reflRay,depth,Xi));
    vec3 tdir = normalize(r.d*nnt - sphereNormal *((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
    float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1.0f-(into?-ddn:dot(tdir, sphereNormal));
    float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj.e + f * (depth>2 ? (erand48(Xi)<P ?   // Russian roulette
                                     radiance(reflRay,depth,Xi)*RP:radiance(Ray(rayIntersectionPoint,tdir),depth,Xi)*TP) :
                          radiance(reflRay,depth,Xi)*Re+radiance(Ray(rayIntersectionPoint,tdir),depth,Xi)*Tr);
}

void render(int fromRow, int toRow, int w, int h, int samps, const Ray& cam, const vec3& cx, const vec3& cy, vec3* c)
{
    for (int y = fromRow; y < toRow; ++y)
    {                       // Loop over image rows
//        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        unsigned short Xi[3] = {0, 0, static_cast<unsigned short>(y*y*y)};
        for (unsigned short x = 0; x < w; ++x)   // Loop cols
            for (int sy=0, i = (h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
                for (int sx=0; sx<2; sx++){        // 2x2 subpixel cols
                    vec3 r = vec3();
                    for (int s=0; s<samps; s++){
                        float r1=2*erand48(Xi), dx=r1<1 ? ::sqrtf(r1)-1: 1-::sqrtf(2-r1);
                        float r2=2*erand48(Xi), dy=r2<1 ? ::sqrtf(r2)-1: 1-::sqrtf(2-r2);
                        vec3 d = cx*( ( (sx+.5f + dx)/2.0f + x)/w - .5f) + cy*( ( (sy+.5f + dy)/2.0f + y)/h - .5f) + cam.d;
                        r = r + radiance(Ray(cam.o+d * 140.0f, normalize(d)),0,Xi)*(1.0f/samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + vec3(clamp(r.x),clamp(r.y),clamp(r.z))*.25f;
                }
    }
}

int main(int argc, char *argv[])
{
    int w=1024, h=768;
    int samps = argc==2 ? atoi(argv[1])/4 : 100; // # samples
    Ray cam(vec3(50.0f, 35.0f, 270.0f), normalize(vec3(0.0f, 0.0f, -1.0f))); // cam pos, dir
    vec3 cx=vec3(w * 0.5135f / h, 0.0f, 0.0f);
    vec3 cy=normalize(cross(cx, cam.d))*0.5135f;
    vec3 r, *c=new vec3[w*h];
    int parts = 4;
    std::thread *tt = new std::thread[parts];
    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < parts; ++i)
    {
        int from = i * (h / parts);
        int to = (i + 1) * (h / parts);
        tt[i] = std::thread(render, from, to, w, h, samps, cam, cx, cy, c);
    }
    for (int i = 0; i < parts; ++i)
    {
        tt[i].join();
    }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i=0; i<w*h; i++)
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
