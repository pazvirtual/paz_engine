uniform float u;
uniform float v;
uniform int charWidth;
uniform int baseWidth;
uniform int baseHeight;
uniform int scale;
uniform int height;
uniform int width;
layout(location = 0) in vec2 pos;
layout(location = 1) in float highlight [[instance]];
layout(location = 2) in int character [[instance]];
layout(location = 3) in int col [[instance]];
layout(location = 4) in int row [[instance]];
out vec2 uv;
flat out float h;
void main()
{
    float a = 1./float(charWidth);
    uv = 0.5*pos + 0.5;
    float x = 2.*u - 1.;
    float y = 2.*v - 1.;
    float xx = x + (uv.x + a + float(col)*(1. + a))*float(scale)*float(
        charWidth)*2./float(width);
    float yy = y + (uv.y + float(row))*float(scale)*float(baseHeight)*2./float(
        height);
    gl_Position = vec4(floor(0.5*float(width)*xx)*2./float(width), floor(0.5*
        float(height)*yy)*2./float(height), 0., 1.);
    uv.x = (uv.x + float(character))*float(charWidth)/float(baseWidth);
    h = highlight;
}
