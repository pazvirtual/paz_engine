uniform sampler2D tex;
uniform float h;
in vec2 uv;
layout(location = 0) out vec4 color;
void main()
{
    float t = texture(tex, uv).r;
    vec3 col = mix(vec3(0.8, 0.8, 0.8), vec3(1., 1., 0.), max(0., h));
    col = mix(col, vec3(0.1, 0.1, 0.1), -min(0., h));
    color = vec4(col, t);
}
