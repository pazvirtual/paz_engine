in vec2 uv;
uniform usampler2D materialMap;
uniform sampler2D normalMap;
uniform sampler2D coordMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
uniform vec4 sun;
layout(location = 0) out vec4 color;
void main()
{
    vec3 dir = normalize((invProjection*vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    vec2 coord = texture(coordMap, uv).rg;
    float lat = (coord.y - 0.5)*3.14159;
    float lon = (coord.x - 0.5)*2.*3.14159;
    float cosAngle = max(0., dot(nor, sun.xyz));
    float diff = cosAngle;
    vec3 halfwayDir = normalize(sun.xyz - dir);
    float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
    color.rgb = mix(diff, spec, 0.3)*vec3(0.8, 1, 0.5) + vec3(0.1, 0, 0);
    if(sin(100.*lat) > 0.95) color.r = 1.;
    if(sin(100.*lon) > 0.95 || coord.y > 0.9 || coord.y < 0.1) color.g = 1.;
    uint mtl = texture(materialMap, uv).r;
    color.rgb = mix(vec3(0, 1, 1), color.rgb, float(mtl > 0u));
}
