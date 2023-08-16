layout(location = 0) in vec2 pos;
uniform float x;
uniform float y;
uniform int scale;
uniform int width;
uniform int height;
uniform int idx;
uniform float aspect;
out vec2 uv;
void main()
{
    vec2 origin = 2.*vec2(x, y) - 1.;
    float a = floor(0.5*(scale*pos.x + origin.x*float(width)))*2./float(width);
    float b = floor(0.5*(scale*pos.y + origin.y*float(height)))*2./float(
        height);
    gl_Position = vec4(a, b, 0., 1.);
    uv = 0.5*pos + 0.5;
    uv.x = (uv.x + idx)*aspect;
}
