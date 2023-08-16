in vec3 lightPos;
in vec3 intens;
in float falloff;
in vec3 loc;
uniform sampler2D diffuseMap;
// ...
uniform sampler2D emissMap;
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
layout(location = 0) out vec4 color;
const float zNear = 0.1; //TEMP
const float zFar = 1e3; //TEMP
void main()
{
    vec2 uv = 0.5*loc.xy/loc.z + 0.5;
    vec3 diffCol = texture(diffuseMap, uv).rgb;
    vec3 emissCol = texture(emissMap, uv).rgb;
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    float depth = texture(depthMap, uv).r;
    float dist = zFar*zNear/(zFar + depth*(zNear - zFar))/-dir.z;
    vec3 pos = dist*dir;
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    vec3 lightDir = lightPos - pos;
    float lightDist = length(lightDir);
    lightDir /= lightDist;
    float cosAngle = max(0., dot(nor, lightDir));
    float diff = cosAngle;
    vec3 halfwayDir = normalize(lightDir - dir);
    float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
    vec3 ill = intens*exp(-falloff*lightDist); //TEMP
    color = vec4(ill*mix(diff, spec, 0.3)*diffCol, 1.);
}
