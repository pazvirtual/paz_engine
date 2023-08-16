uniform mat4 invProjection;
uniform vec4 sunDir;
uniform vec4 sunIll;
uniform uint numLights;
uniform vec4 light0[16];
uniform vec4 light1[16];
in vec3 _pos;
in vec3 _nor;
in vec3 _loc;
layout(location = 0) out vec4 accum;
layout(location = 1) out vec4 r;
const vec3 diffCol = vec3(0.7, 0.8, 1.); //TEMP
const int specPower = 32;
const float specFac = 0.3;
const int edgePower = 2;
const float minAlpha = 0.01;
const float glareFac = 0.1;
float luminance(in vec3 v)
{
    return dot(v, vec3(0.2126, 0.7152, 0.0722));
}
void main()
{
    vec2 uv = 0.5*_loc.xy/_loc.z + 0.5;
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 color = 0.01*diffCol;
    vec3 nor = normalize(_nor);
    float cosViewAngle = abs(dot(nor, dir));
    float cosSunAngle = abs(dot(nor, sunDir));
    float diff = cosSunAngle;
    vec3 halfwayDir = normalize(sunDir - dir);
    float spec = cosSunAngle*pow(abs(dot(nor, halfwayDir)), specPower);
    color += sunIll*mix(diff, spec, specFac)*diffCol;
    for(uint i = 0u; i < numLights; ++i)
    {
        vec3 lightPos = light0[i].xyz;
        vec3 intens = vec3(light0[i].w, light1[i].xy);
        float falloff = light1[i].z;
        vec3 lightDir = lightPos - _pos;
        float lightDist = length(lightDir);
        lightDir /= lightDist;
        float cosLightAngle = abs(dot(nor, lightDir));
        float diff = cosLightAngle;
        vec3 halfwayDir = normalize(lightDir - dir);
        float spec = cosLightAngle*pow(abs(dot(nor, halfwayDir)), specPower);
        vec3 ill = intens*exp(-falloff*lightDist); //TEMP
        color += ill*mix(diff, spec, specFac)*diffCol;
    }
    float alpha = max(minAlpha, mix(pow(1. - cosViewAngle, edgePower),
        luminance(color), glareFac));
    float w = pow(abs(_pos.z), -0.1);
    accum = vec4(color*alpha, alpha)*w;
    r = vec4(0., 0., 0., alpha);
}
