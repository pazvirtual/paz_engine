in vec3 loc;
uniform sampler2D diffuseMap;
// ...
uniform sampler2D emissMap;
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
layout(location = 0) out vec4 color;
void main()
{
    vec2 uv = loc.xy/loc.z;
    vec3 diffCol = texture(diffuseMap, uv).rgb;
    vec3 emissCol = texture(emissMap, uv).rgb;
    color = vec4(0.01*diffCol + emissCol, 1.);
}
