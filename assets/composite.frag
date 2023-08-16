uniform sampler2D accumTex;
uniform sampler2D revealTex;
in vec2 uv;
layout(location = 0) out vec4 color;
void main()
{
    vec4 accum = texture(accumTex, uv);
    float r = texture(revealTex, uv).a;
    color = vec4(accum.rgb/max(accum.a, 1e-5), r);
}
