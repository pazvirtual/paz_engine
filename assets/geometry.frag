in vec4 posCs;
in vec4 norCs;
in vec2 uv;
uniform sampler2D diffTex;
uniform vec3 emiss;
layout(location = 0) out vec4 diffCol;
// ...
layout(location = 1) out vec4 emissCol;
layout(location = 2) out vec4 normal;
void main()
{
    diffCol = texture(diffTex, uv);
    // ...
    emissCol = vec4(emiss, 1.);
    normal = vec4(normalize(norCs.xyz), 0.);
}
