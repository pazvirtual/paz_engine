uniform vec4 col;
in float _z;
layout(location = 0) out vec4 accum;
layout(location = 1) out vec4 r;
void main()
{
    vec3 ci = col.a*col.rgb;
    float ai = col.a;
    float w = 1./(_z*_z);
    accum = vec4(ci, ai)*w;
    r = vec4(0., 0., 0., ai);
}
