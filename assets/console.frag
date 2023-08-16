uniform float u;
uniform float v;
uniform float width;
uniform float height;
in vec2 uv;
layout(location = 0) out vec4 color;
void main()
{
    color = vec4(0., 0., 0., 0.);
    if(uv.x > u && uv.x < width + u && uv.y > v && uv.y < height + v)
    {
        color.a = 0.95;
    }
}
