in vec2 uv;
uniform sampler2D tex;
uniform float h;
layout(location = 0) out vec4 color;
void main()
{
    float t = texture(tex, uv).r;
    color = vec4(mix(vec3(0.8, 0.8, 0.8), vec3(1., 1., 0.), h), t);
}
