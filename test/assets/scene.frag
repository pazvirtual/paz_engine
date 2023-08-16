in vec2 uv;
uniform sampler2D diffuseMap;
// ...
uniform sampler2D normalMap;
uniform depthSampler2D depthMap;
uniform mat4 invProjection;
uniform sampler2D lights; //TEMP - should be 1D
layout(location = 0) out vec4 color;
const float zNear = 0.1; //TEMP
const float zFar = 1e3; //TEMP
void main()
{
    vec3 diffCol = texture(diffuseMap, uv).rgb;
    vec3 dir = normalize(mul(invProjection, vec4(2.*uv - 1., 1., 1.)).xyz);
    float depth = texture(depthMap, uv).r;
    float dist = zFar*zNear/(zFar + depth*(zNear - zFar))/-dir.z;
    vec3 pos = dist*dir;
    vec3 nor = normalize(texture(normalMap, uv).rgb);
    color = vec4(0.01*diffCol, 1.);
    int numLights = textureSize(lights, 0).y;
    for(int i = 0; i < numLights; ++i)
    {
        vec2 c = vec2(0.5, (i + 0.5)/numLights); //TEMP should use `texelFetch`
        vec4 curLight = texture(lights, c);
        vec3 lightPos = curLight.xyz;
        vec3 lightDir = lightPos - pos;
        float lightDist = length(lightDir);
        lightDir /= lightDist;
        float cosAngle = max(0., dot(nor, lightDir));
        float diff = cosAngle;
        vec3 halfwayDir = normalize(lightDir - dir);
        float spec = cosAngle*pow(max(0., dot(nor, halfwayDir)), 32);
        float ill = curLight.w*2.*exp(-0.1*lightDist); //TEMP
        color.rgb += ill*mix(diff, spec, 0.3)*diffCol;
    }
}
