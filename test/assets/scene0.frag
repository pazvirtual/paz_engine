in vec3 loc;
uniform sampler2D diffuseMap;
// ...
uniform sampler2D emissMap;
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
uniform vec4 lightDir;
uniform vec4 ill;
uniform depthSampler2D shadowMap;
uniform mat4 invView;
uniform mat4 lightProjection;
uniform mat4 lightView;
uniform float zFar;
uniform float zNear;
uniform vec4 relLightPos;
layout(location = 0) out vec4 color;
const uint res = 8u;
void main()
{
    ivec2 texSize = textureSize(shadowMap, 0);
    vec2 texOffset = vec2(1./texSize.x, 1./texSize.y);
    vec2 uv = loc.xy/loc.z;
    vec3 diffCol = texture(diffuseMap, uv).rgb;
    vec3 emissCol = texture(emissMap, uv).rgb;
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    float cosAngle = max(0., dot(nor, lightDir.xyz));
    float diff = cosAngle;
    vec3 halfwayDir = normalize(lightDir.xyz - dir);
    float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
    float depth = texture(depthMap, uv).r;
    float dist = zFar*zNear/(zFar + depth*(zNear - zFar))/-dir.z;
    vec3 pos = dist*dir;
    vec3 temp0 = mul(invView, vec4(pos, 1.)).xyz - relLightPos.xyz;
    vec3 posLs = mul(lightView, vec4(temp0, 1.)).xyz;
    vec4 temp1 = mul(lightProjection, vec4(posLs, 1.));
    temp1 = 0.5*temp1/temp1.w + 0.5;
    vec2 shadowUv = temp1.xy;
    float depthLs = temp1.z;
    float bias = min(3e-3*sqrt(1. - cosAngle*cosAngle)/cosAngle, 3e-2); // should compute in shadow map
    uint c1 = 0;
    for(uint i = 0u; i < res; ++i) // should be gaussian kernel
    {
        for(uint j = 0u; j < res; ++j)
        {
            vec2 offset = (2./float(res - 1u)*vec2(i, j) - 1.)*texOffset;
            if(depthLs > texture(shadowMap, shadowUv + offset).r + bias)
            {
                ++c1;
            }
        }
    }
    float c = 1. - float(c1)/(float(res)*float(res));
    color = vec4((0.01 + c*ill.rgb*mix(diff, spec, 0.3))*diffCol + emissCol,
        1.);
}
