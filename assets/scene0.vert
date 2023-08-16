layout(location = 0) in vec2 pos;
out vec3 loc;
void main()
{
    gl_Position = vec4(pos, 0, 1);
    loc = vec3(0.5*pos + 0.5, 1.);
}
