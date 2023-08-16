in vec2 uv;
//uniform usampler2D materialMap;
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 projection;
uniform vec4 sun;
layout(location = 0) out vec4 color;
void main()
{
    vec3 dir = normalize((inverse(projection)*vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    float diff = max(0., dot(nor, sun.xyz));
    vec3 halfwayDir = normalize(sun.xyz - dir);
    float spec = pow(max(0., dot(nor, halfwayDir)), 32);
    color = vec4(mix(diff, spec, 0.3) + 5e-3);
}
