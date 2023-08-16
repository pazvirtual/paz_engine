in vec4 norCs;
layout(location = 0) out vec4 normal;
void main()
{
    normal = vec4(normalize(norCs.xyz), 0.);
}
