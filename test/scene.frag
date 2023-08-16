in vec2 uv;
//uniform usampler2D materialMap;
uniform sampler2D normalMap;
//uniform sampler2D directionMap;
//uniform depthSampler2D depthMap;
layout(location = 0) out vec4 color;
void main()
{
    color = vec4(abs(texture(normalMap, uv)));
}
