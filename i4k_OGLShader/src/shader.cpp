#include "shader.h"

static const char *vsh = R"SHADER(
#version 430

layout (location=0) in vec2 inVer;
out vec2 p;

out gl_PerVertex 
{
  vec4 gl_Position;
};


void main()
{
  gl_Position=vec4(inVer,0.0,1.0);
  p=inVer;
}
)SHADER";

static const char * fsh = R"SHADER(
#version 430

layout (location=0) uniform vec4 fpar[];
layout (location=0) out vec4 co;
in vec2 p;

vec2 iResolution = vec2(1.0);
float iTime = 1.0;

// -----------------------------------------------------------------------------


#define TOLERANCE       0.00001
#define MAX_RAY_LENGTH  64.0
#define MAX_BOUNCES     10
#define MAX_RAY_MARCHES 80

#define PI              3.141592654
#define TAU             (2.0*PI)

const float fadeInTime  = 2.0;
const float fadeOutTime = 28.0;

float maxComp(in vec3 p)
{
  return max(p.x,max(p.y,p.z));
}

float maxcomp(in vec2 p)
{
  return max(p.x, p.y);
}

float smin(float a, float b, float k)
{
  float res = exp( -k*a ) + exp( -k*b );
  return -log( res )/k;
}

float sdSphere(in vec3 p, in float r)
{
  return length(p) - r;
}

float sdBox(vec3 p, vec3 b)
{
  vec3  di = abs(p) - b;
  float mc = maxComp(di);
  return min(mc,length(max(di,0.0)));
}

float sdCross(in vec3 p, float s)
{
  float da = maxcomp(abs(p.xy));
  float db = maxcomp(abs(p.yz));
  float dc = maxcomp(abs(p.zx));
  return min(da,min(db,dc))-s;
}

float lengthN(in vec3 v, in float n)
{
  vec3 vv = pow(v, vec3(n));
  return pow(vv.x + vv.y + vv.z, 1.0/n);
}

float sdRoundCube(in vec3 p, float r)
{
  return lengthN(abs(p), 8.0) - r;
}

float pModMirror1(inout float p, float size) 
{
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize,size) - halfsize;
  p *= mod(c, 2.0)*2.0 - 1.0;
  return c;
}

float sdRepeatedCross(in vec3 p, in float r)
{
  p -= r;
  p = mod(p, 2.0*r) - r;
  return sdCross(p, r/3.0);
}

float pMod1(inout float p, float size) 
{
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize, size) - halfsize;
  return c;
}

vec3 pMod3(inout vec3 p, vec3 size) 
{
  vec3 c = floor((p + size*0.5)/size);
  p = mod(p + size*0.5, size) - size*0.5;
  return c;
}

float sgn(float x) 
{
  return (x<0.0)?-1.0:1.0;
}

float pReflect(inout vec3 p, vec3 planeNormal, float offset) 
{
  float t = dot(p, planeNormal)+offset;
  if (t < 0.0) 
  {
    p = p - (2.0*t)*planeNormal;
  }
  return sgn(t);
}

float impulse1(in vec3 p, out vec3 col, out float ref)
{
  col = vec3(1.0);  
  ref = 0.9;

  float sb = sdBox(p, vec3(0.4));
  
  vec3 pp = p - vec3(0.1) - vec3(0.0, iTime*0.15 + 10.0, 0.0);
  pReflect(pp, normalize(vec3(1.0, 0.5, 0.2)), 0.3);
  pReflect(pp, normalize(vec3(0.2, 0.5, 1.0)), 0.2);
  //pReflect(pp, normalize(vec3(0.5, 0.2, 1.0)), 0.07);
  pMod3(pp, vec3(0.5, 0.3, 0.4));
  
  vec3 ppp = p - vec3(0.2) - vec3(0.0, iTime*0.05 + 10.0, 0.0);
  pReflect(ppp, normalize(vec3(0.7, 0.5, 0.4)), 0.3);
  pReflect(ppp, normalize(vec3(0.5, 0.4, 0.7)), 0.1);
  pMod3(ppp, vec3(0.7, 0.6, 0.4));

  float ss = sdSphere(pp, 0.05);
  float sss = sdSphere(ppp, 0.1);

  float st = smin(ss, sss, 20.0);
  
  return max(sb, -st);
}

float impulse2(in vec3 p, out vec3 col, out float ref)
{
  col = vec3(1.0);  
  ref = 0.9;
    
  float s1 = sdBox(p, vec3(0.4));
  pMod1(p.x, 0.4);
  pMod1(p.y, 0.4);
  pMod1(p.z, 0.4);
  float s2 = sdSphere(p, 0.18);  
  float s = max(s1, -s2);
 
  return s;
}

void pR(inout vec2 p, float a) 
{
  p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}


float impulse3(in vec3 p, out vec3 col, out float ref)
{
  col = vec3(1.0);  
  ref = 0.9;
  float s1 = sdRoundCube(p, 0.33);
  float s2 = sdSphere(p, 0.40);
  vec3 pp = p;
  pp -= vec3(0.6, -0.31, 0.0);
  
  float s3 = sdSphere(pp, 0.1);
 
  return min(min(s1, s2), s3);
}

void pR45(inout vec2 p) 
{
  p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

float distanceEstimator(in vec3 p, out vec3 col, out float ref)
{
  float b = sdRoundCube(p - vec3(0.0, -1.41, 0.0), 1.0); 
  float i = impulse2(p, col, ref);
  float f = min(b, i);
  if (f == b)
  {
    col = vec3(1.0);
    ref = 0.0;
  }
  return f;
}

vec3 saturate(in vec3 a) { return clamp(a, 0.0, 1.0); }
vec2 saturate(in vec2 a) { return clamp(a, 0.0, 1.0); }
float saturate(in float a) { return clamp(a, 0.0, 1.0); }

const vec3 lightPos1 = 100.0*vec3(-0.3, 0.15, 1.0);
const vec3 lightPos2 = 100.0*vec3(-0.33,  -0.2, -1.0);

vec3 getSkyColor(vec3 rayDir)
{
  vec3 lightDir1 = normalize(lightPos1);
  vec3 lightDir2 = normalize(lightPos2);
  vec3 lightCol1 = vec3(8.0/8.0,7.0/8.0,6.0/8.0);
  vec3 lightCol2 = vec3(8.0/8.0,6.0/8.0,7.0/8.0);
  float ld1      = max(dot(lightDir1, rayDir), 0.0);
  float ld2      = max(dot(lightDir2, rayDir), 0.0);
  vec3 final     = vec3(0.125);

  if ((rayDir.y > abs(rayDir.x)*1.0) && (rayDir.y > abs(rayDir.z*0.25))) final = vec3(2.0)*rayDir.y;
  float roundBox = length(max(abs(rayDir.xz/max(0.0,rayDir.y))-vec2(0.9, 4.0),0.0))-0.1;
  final += vec3(0.8)* pow(saturate(1.0 - roundBox*0.5), 6.0);
  
  final += 1.0*pow(lightCol1, vec3(2.0, 1.5, 1.5)) * pow(ld1, 8.0);
  final += 1.0*lightCol1 * pow(ld1, 200.0);
  final += 1.0*pow(lightCol2, vec3(2.0, 1.5, 1.5)) * pow(ld2, 8.0);
  final += 1.0*lightCol2 * pow(ld2, 200.0);
  return final;
}

vec3 normal(in vec3 pos)
{
  vec3 col;
  float ref;
  vec3  eps = vec3(.001,0.0,0.0);
  vec3 nor;
  nor.x = distanceEstimator(pos+eps.xyy, col, ref) - distanceEstimator(pos-eps.xyy, col, ref);
  nor.y = distanceEstimator(pos+eps.yxy, col, ref) - distanceEstimator(pos-eps.yxy, col, ref);
  nor.z = distanceEstimator(pos+eps.yyx, col, ref) - distanceEstimator(pos-eps.yyx, col, ref);
  return normalize(nor);
}

float rayMarch(in float dmod, in vec3 ro, inout vec3 rd, float mint, float maxt, out int rep, out vec3 col, out float ref)
{
  float t = mint;
  int i = 0;
  for (i = 0; i < MAX_RAY_MARCHES; i++)
  {
    float distance_ = distanceEstimator(ro + rd*t, col, ref);
    float distance = dmod*distance_;
    if (distance < TOLERANCE || t > maxt) break;
    t += max(distance, 0.001);
  }
  rep = i;
  return t;
}

float ambientOcclusion(vec3 p, vec3 n)
{
  vec3 col;
  float ref;
  float stepSize = 0.01;
  float t = stepSize;

  float oc = 0.0;

  for(int i = 0; i < 10; i++)
  {
    float d = distanceEstimator(p + n * t, col, ref);
    oc += t - d; // Actual distance to surface - distance field value
    t += stepSize;
  }

  return clamp(oc, 0.0, 1.0);
}

float softShadow(in vec3 pos, in vec3 ld, in float ll, float mint, float k)
{
  vec3 col;
  float ref;
  const float minShadow = 0.25;
  float res = 1.0;
  float t = mint;
  for (int i=0; i<24; i++)
  {
    float distance = distanceEstimator(pos + ld*t, col, ref);
    res = min(res, k*distance/t);
    if (ll <= t) break;
    if(res <= minShadow) break;
    t += max(mint*0.2, distance);
  }
  return clamp(res,minShadow,1.0);
}

float specular(in vec3 nor, in vec3 ld, in vec3 rd)
{
  return pow(max(dot(reflect(-ld, nor), -rd), 0.), 75.);
}

vec3 postProcess(in vec3 col, in vec2 q) 
{
  col=pow(clamp(col,0.0,1.0),vec3(0.45)); 
  col=col*0.6+0.4*col*col*(3.0-2.0*col);  // contrast
  col=mix(col, vec3(dot(col, vec3(0.33))), -0.4);  // satuation
  col*=0.5+0.5*pow(19.0*q.x*q.y*(1.0-q.x)*(1.0-q.y),0.7);  // vigneting
  return col;
}


vec3 render(in vec3 ro, in vec3 rd)
{
  vec3 lightPos = 2.0*vec3(1.5, 3.0, 1.0);

  vec3 col    = vec3(0.0);

  vec3 ragg2 = vec3(1.0);
  
  float tdist = 0.0;

  float refraction = 0.9;
  
  bool inside = false;
    
  for (int i = 0; i < MAX_BOUNCES && maxComp(ragg2) > 0.01; ++i)
  {
    vec3 mat    = vec3(0.0);
    float dmod  = inside ? -1.0 : 1.0;
    float rscale= 0.0;
    int rep     = 0;
    float t     = rayMarch(dmod, ro, rd, 0.01, MAX_RAY_LENGTH, rep, mat, rscale);
    tdist       += t;
  
    vec3 pos    = ro + t*rd;

    vec3 nor = vec3(0.0, 1.0, 0.0);
    
    if (t < MAX_RAY_LENGTH)
    {
      // Ray intersected object
      nor = normal(pos);
    }
    else
    {
      // Ray intersected sky
      col += ragg2*getSkyColor(rd);
      break;
    }

    vec3 refr;
    vec3 refl;
    
    if (!inside)
    {
      refl = reflect(rd, nor);
      refr = refract(rd, nor, refraction);
    }
    else
    {
      refl = reflect(rd, -nor);
      refr = refract(rd, -nor, 1.0/refraction);
    }
    

    vec3 lv   = lightPos - pos;
    //float ll2 = dot(lv, lv);
    vec3  ld  = normalize(lv);
    float ll  = length(lv);
    // TODO: Rework shadow to "work" with transparent objects
    float sha = 1.0;
    if (!inside)
    {
      sha = softShadow(pos, ld, ll, 0.01, 64.0);
    }
    //float sha = 1.0;

    float dif = max(dot(nor,ld),0.0);
    //float occ = 1.0 - ambientOcclusion(pos, nor);
    //float spe = specular(nor, ld, rd);
    float occ = 1.0 - float(rep)/float(MAX_RAY_MARCHES);
    float l   = dif*occ*sha;

    float lin = mix(0.2, 1.0, l);
    
    vec3 mcol = 0.8*lin*mat + 0.2*getSkyColor(refl);

    vec3 beer = vec3(1.0);
    
    if (inside)
    {
      const vec3 color = vec3(2.0, 2.0, 1.0);
      beer = exp(-color*t);
    }
    col        += (1.0 - rscale)*ragg2*beer*mcol;
    ragg2      *= rscale*beer;

    if (abs(dot(rd, nor)) < 0.1)
    {
      // Ray bounce will follow the surface if flat. Assume
      // sky intersect
      col += ragg2*getSkyColor(reflect(rd, nor));
      break;
    }
      
    ro        = pos;           
    
    if (refr == vec3(0.0))
    {
       rd = refl;
    }
    else
    {
      rd = refr;
    }
      
    if (dot(refl, rd) < 0.9)
    {
      inside = !inside;
    }
       //rd = refl;
  }
    
 
  return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
  vec2 q=fragCoord.xy/iResolution.xy; 
  vec2 p = -1.0 + 2.0*q;
  p.x *= iResolution.x/iResolution.y;

  float ctime = iTime/2.0;
  //ctime = 0.0;

  // camera
  float z  = 2.0;
  vec3 ro = 0.4*vec3(z*2.0, 1.0, 0.0);
  //pR(ro.xy, -PI/2.0*(0.5 + 0.5*sin(ctime)) + 0.1);
  pR(ro.xz, ctime);
  vec3 ww = normalize(vec3(0.0, 0.0, 0.0) - ro);
  vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww ));
  vec3 vv = normalize(cross(ww,uu));
  vec3 rd = normalize( p.x*uu + p.y*vv + 2.5*ww );

  vec3 col = render(ro, rd);

  float fadeIn = smoothstep(0.0, fadeInTime, iTime);
  float fadeOut = 1.0 - smoothstep(fadeOutTime, fadeOutTime + 2.0, iTime);
  fadeIn = 1.0;
  fadeOut = 1.0;

  fragColor = vec4(postProcess(col, q)*fadeIn*fadeOut,1.0);
}

// -----------------------------------------------------------------------------

void main()
{
  iTime = fpar[0].x;
  iResolution.x = fpar[0].y;
  iResolution.y = fpar[0].z;
  vec2 pp = (p + 1.0)*0.5*iResolution.xy;

  mainImage(co, pp);
}
)SHADER";


extern "C"
{
  char * getVertexShader() { return const_cast<char*>(vsh);  }
  char * getFragmentShader() { return const_cast<char*>(fsh); }
}