in vec3 loc;
uniform sampler2D diffuseMap;
// ...
uniform sampler2D emissMap;
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
uniform vec4 lightDir;
uniform vec4 ill;
layout(location = 0) out vec4 color;
void main()
{
    vec2 uv = loc.xy/loc.z;
    vec3 diffCol = texture(diffuseMap, uv).rgb;
    vec3 emissCol = texture(emissMap, uv).rgb;
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    float cosAngle = max(0., dot(nor, lightDir.xyz));
    float diff = cosAngle;
    vec3 halfwayDir = normalize(lightDir.xyz - dir);
    float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
    color = vec4((0.01 + ill.rgb*mix(diff, spec, 0.3))*diffCol + emissCol, 1.);
}
