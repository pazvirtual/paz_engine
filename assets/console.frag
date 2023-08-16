uniform float width;
uniform float height;
in vec2 uv;
layout(location = 0) out vec4 color;
void main()
{
    color = vec4(0., 0., 0., 0.);
    if(uv.y < height && uv.x < width)
    {
        color.a = 0.5;
    }
}
