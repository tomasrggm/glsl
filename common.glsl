/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);

    //Calculate eye_offset and ray direction
    vec3 ps;  //pixel samole or ray in camera coordinates
    ps.x = cam.width * (pixel_sample.x / iResolution.x - 0.5);
    ps.y = cam.height * (pixel_sample.y / iResolution.y - 0.5);
    ps.z = -cam.planeDist;

    //para o dof
    vec3 p;
    p.x = ps.x * cam.focusDist;
    p.y = ps.y * cam.focusDist;
    p.z = ps.z * cam.focusDist;

    //calculates ray in world coordinates
    // vec3 vX = cam.u * ps.x;
    // vec3 vY = cam.v * ps.y;
    // vec3 vZ = cam.n * ps.z;
    vec3 vX = cam.u * (p.x - ls.x);
    vec3 vY = cam.v * (p.y - ls.y);
    vec3 vZ = cam.n * p.z;
    vec3 ray_dir = normalize((vX + vY + vZ));

    //vec3 ray_direction = vec3(  p.x, cam.focusDist * cam.height * (pixel_sample.y / iResolution.y - 0.5), -cam.focusDist * cam.planeDist);
    vec3 eye_offset = cam.eye + cam.u * ls.x + cam.v * ls.y;

    //return createRay(eye_offset, normalize(ray_direction.x * cam.u + ray_direction.y * cam.v + ray_direction.z * cam.n), time);
    return createRay(eye_offset, ray_dir, time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};

//baseeie-me no  manual ray tracing one weekend e funciona
bool refract (vec3 v, vec3 n, float niOverNt, out vec3 refracted){
    vec3 uv = normalize(v);
    float dt = dot(uv,n);
    float discriminant = 1.0 - niOverNt * niOverNt * (1.0 - dt * dt);
    if(discriminant > 0.0){
        refracted = niOverNt * (uv - n * dt) - n * sqrt(discriminant);
        return true;
    }
    return false;
}


//vimos na lab que esta bem
float schlick(float cosine, float refIdx)
{

    float x = (1.0 - refIdx) / (1.0 + refIdx);
    x = x*x;
    float temp = clamp(1.0 - cosine, 0.0, 1.0);
    
    return x + (1.0 - x) * pow(clamp(1.0 - cosine,0.0,1.0),5.0);
    //INSERT YOUR CODE HERE
}




bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if(rec.material.type == MT_DIFFUSE)
    {
        //INSERT CODE HERE,
        vec3 target = rec.pos + rec.normal + randomUnitVector(gSeed);
        rScattered = createRay(rec.pos + epsilon * rec.normal, normalize(target-rec.pos),rIn.t);
        //rScattered = createRay(rec.pos + epsilon * rec.normal, normalize(refl + rec.material.roughness * randomInUnitSphere(gSeed)),rIn.t);
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
       //INSERT CODE HERE, consider fuzzy reflections
       vec3 refl = reflect(rIn.d,rec.normal);
       float cosi = -dot(rIn.d,rec.normal);
       float white = 1.0;
       rScattered = createRay(rec.pos + epsilon * rec.normal, normalize(refl + rec.material.roughness * randomInUnitSphere(gSeed)),rIn.t);
       
       //link que usei para resolver https://github.com/RayTracing/raytracing.github.io/issues/844 aqui faz sentido que seja a especular visto que e uma carateristica inerente a metais enquanto que albedo e mais para componentes difusas.
       // alias bastava olhar para as funcoes de create material e percebia-se logo. albedo e 0 para metais.
       atten = rec.material.specColor + ( 1.0 - rec.material.specColor ) * pow(clamp(1.0-cosi,0.0,1.0),5.0); 
       return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        atten = vec3(1.0);
        vec3 outwardNormal;
        float niOverNt;
        float cosine;
        cosine = -dot(rIn.d, rec.normal);
        if(cosine < 0.0) //hit inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
             
            atten = exp(-rec.material.refractColor * rec.t);
        }
        else  //hit from outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray

        float reflectProb;
        vec3 refracted; //PRECISO DE VER ESTE VALOR AQUI

        //if no total reflection  reflectProb = schlick(cosine, rec.material.refIdx);  
        //else reflectProb = 1.0;
        if(refract(rIn.d,outwardNormal,niOverNt,refracted)){
            if(niOverNt > 1.0){
                cosine =sqrt(1.0-niOverNt*niOverNt*(1.0-cosine*cosine));
            }
            reflectProb = schlick(cosine,rec.material.refIdx);
        }
        else{
            reflectProb = 1.0;
        }
        

        if( hash1(gSeed) < reflectProb){  //Reflection
            vec3 reflected = reflect(rIn.d, rec.normal);
            rScattered = createRay(rec.pos + epsilon * rec.normal, normalize(reflected),rIn.t);

        // rScattered = calculate reflected ray
          // atten *= vec3(reflectProb); not necessary since we are only scattering reflectProb rays and not all reflected rays
        
        //else  //Refraction
        // rScattered = calculate refracted ray
           // atten *= vec3(1.0 - reflectProb); not necessary since we are only scattering 1-reflectProb rays and not all refracted rays
        }else{
            rScattered = createRay(rec.pos - epsilon * outwardNormal, normalize(refracted),rIn.t);
        }

        return true;
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    //calculate a valid t and normal
    //calculate a valid t and normal
    vec3 top = t.b - t.a;
    vec3 left = t.c - t.a;
    vec3 cro = cross(r.d,left);
    float det = dot(top,cro);
    rec.normal = normalize(cross(top,left)); //formar normal como no proj1
    if(abs(det) < 0.000001) return false;

    vec3 pointToOrigin = r.o - t.a;
    float invDet = 1.0/ det;

    float u = dot(pointToOrigin,cro) * invDet;
    if(u < 0.0 || u > 1.0){
        return false;
    } 

    vec3 cro2 = cross(pointToOrigin , top);
    float v = dot(r.d,cro2) * invDet;
    if(v < 0.0 || u+v > 1.0) {
        return false;
    }

    float temp = (dot(top,cross(left,pointToOrigin))) * invDet;
    if(temp < tmax && temp > tmin)
    {
        rec.t = temp;
        //rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    //center captures the sphere moving linearly from one point to the other
    return mvsphere.center0 + ((time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0)) * (mvsphere.center1 - mvsphere.center0);
}


bool solveQuadratic( float a, float b,  float c, float x0, float x1)
{
	float discr = b * b - 4.0 * a * c;
	if (discr < 0.0) return false;
	else if (discr == 0.0) {
		x0 = x1 = -0.5 * b / a;
	}
	else {
		float q = (b > 0.0) ?
			-0.5 * (b + sqrt(discr)) :
			-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}

	return true;
}

/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    //calculate a valid t and normal

    vec3 oc = r.o - s.center;
    float a = dot(r.d,r.d);
    float b = dot(oc,r.d);
    float c = dot(oc,oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;
    if(discriminant > 0.0){
        float sqrtDiscriminant = sqrt(discriminant);
        float temp = (-b - sqrtDiscriminant) / a;
        if(temp < tmax && temp > tmin){
            rec.t = temp;
            rec.pos = pointOnRay(r,rec.t);
            rec.normal = (rec.pos - s.center) / s.radius;
            return true;
        }
        temp = (-b + sqrtDiscriminant) / a;
        if(temp < tmax && temp > tmin){
            rec.t = temp;
            rec.pos = pointOnRay(r,rec.t);
            rec.normal = (rec.pos - s.center) / s.radius;
            return true;
        }
    }
    return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 sphereCenter = center(s,r.t);
    vec3 oc = r.o - sphereCenter;
    float a = dot(r.d,r.d);
    float b = dot(oc,r.d);
    float c = dot(oc,oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;
    if(discriminant > 0.0){
        float sqrtDiscriminant = sqrt(discriminant);
        float temp = (-b - sqrtDiscriminant) / a;
        if(temp < tmax && temp > tmin){
            rec.t = temp;
            rec.pos = pointOnRay(r,rec.t);
            rec.normal = (rec.pos - sphereCenter) / s.radius;
            return true;
        }
        temp = (-b + sqrtDiscriminant) / a;
        if(temp < tmax && temp > tmin){
            rec.t = temp;
            rec.pos = pointOnRay(r,rec.t);
            rec.normal = (rec.pos - sphereCenter) / s.radius;
            return true;
        }
    }
    return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}