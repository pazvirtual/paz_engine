in vec2 uv;
uniform sampler2D diffuseMap;
// ...
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
uniform vec4 sun;
layout(location = 0) out vec4 color;
void main()
{
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    float cosAngle = max(0., dot(nor, sun.xyz));
    float diff = cosAngle;
    vec3 halfwayDir = normalize(sun.xyz - dir);
    float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
    vec3 diffCol = texture(diffuseMap, uv).rgb;
    color = vec4((0.1 + mix(diff, spec, 0.3))*diffCol, 1.);
}
