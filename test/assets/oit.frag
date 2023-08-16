uniform mat4 invProjection;
in vec3 _pos;
in vec3 _nor;
in vec3 _loc;
in vec3 _lightPos;
in vec3 _intens;
in float _falloff;
layout(location = 0) out vec4 accum;
layout(location = 1) out vec4 r;
const vec3 diffCol = vec3(0.1, 0.8, 1.); //TEMP
void main()
{
    vec2 uv = 0.5*_loc.xy/_loc.z + 0.5;
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 lightDir = _lightPos - _pos;
    float lightDist = length(lightDir);
    lightDir /= lightDist;
    vec3 nor = normalize(_nor);
    float cosLightAngle = abs(dot(nor, lightDir));
    float cosViewAngle = abs(dot(nor, dir));
    float diff = cosLightAngle;
    vec3 halfwayDir = normalize(lightDir - dir);
    float spec = cosLightAngle*pow(abs(dot(nor, halfwayDir)), 32);
    vec3 ill = _intens*exp(-_falloff*lightDist); //TEMP
    vec3 color = (0.01 + ill*mix(diff, spec, 0.3))*diffCol;
    float alpha = max(0.01, mix((1. - cosViewAngle)*(1. - cosViewAngle), spec,
        0.3));
    float w = pow(abs(_pos.z), -0.1);
    accum = vec4(color*alpha, alpha)*w;
    r = vec4(0., 0., 0., alpha);
}
