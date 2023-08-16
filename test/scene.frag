in vec2 uv;
uniform usampler2D materialMap;
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
uniform vec4 sun;
layout(location = 0) out vec4 color;
void main()
{
    vec3 dir = normalize((invProjection*vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 nor = texture(normalMap, uv).rgb;
    float cosAngle = max(0., dot(nor, sun.xyz));
    float diff = cosAngle;
    vec3 halfwayDir = normalize(sun.xyz - dir);
    float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
    color.rgb = mix(diff, spec, 0.3)*vec3(0.8, 1, 0.5) + vec3(0.1, 0, 0);
    uint mtl = texture(materialMap, uv).r;
    color.rgb = mix(vec3(0, 1, 1), color.rgb, float(mtl > 0u));
}
