uniform sampler2D img;
in vec2 uv;
layout(location = 0) out float lum;
void main()
{
    lum = dot(texture(img, uv).rgb, vec3(0.2126, 0.7152, 0.0722));
}
