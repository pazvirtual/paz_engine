uniform sampler2D hdrRender;
uniform depthSampler2D depthMap;
uniform float minDepth;
uniform float maxDepth;
in vec2 uv;
layout(location = 0) out vec4 color;
void main()
{
    float depth = texture(depthMap, uv).r;
    color = texture(hdrRender, uv);
    if(depth < minDepth || depth > maxDepth)
    {
        color += texture(hdrRender, uv + 1e-3);
        color += texture(hdrRender, uv - 1e-3);
        color += texture(hdrRender, uv + vec2(1e-3, -1e-3));
        color += texture(hdrRender, uv - vec2(1e-3, -1e-3));
        color /= 5.;
    }
}
