layout(location = 0) in vec2 pos;
uniform float x;
uniform float y;
uniform int scale;
uniform int width;
uniform int height;
out vec2 uv;
void main()
{
    vec2 origin = 2.*vec2(x, y) - 1.;
    gl_Position = vec4(float(scale)*pos/vec2(width, height) + origin, 0., 1.);
    uv = 0.5*pos + 0.5;
}
