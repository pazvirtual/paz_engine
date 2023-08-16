uniform sampler2D font;
in vec2 uv;
flat in float h;
layout(location = 0) out vec4 color;
void main()
{
    float t = texture(font, uv).x;
    color = vec4(mix(vec3(0.8, 0.8, 0.8), vec3(1., 1., 0.), h), t);
}
