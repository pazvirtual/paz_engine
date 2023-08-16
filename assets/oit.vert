layout(location = 0) in vec4 position;
layout(location = 1) in vec4 normal;
uniform vec4 model0;
uniform vec2 model1;
uniform mat4 projection;
uniform mat4 view;
uniform vec4 light0;
uniform vec4 light1;
out vec3 _pos;
out vec3 _nor;
out vec3 _loc;
out vec3 _lightPos;
out vec3 _intens;
out float _falloff;
void main()
{
    vec4 att = vec4(model0.xyz, sqrt(1. - dot(model0.xyz, model0.xyz)));
    float xx = att.x*att.x;
    float yy = att.y*att.y;
    float zz = att.z*att.z;
    float xy = att.x*att.y;
    float zw = att.z*att.w;
    float xz = att.x*att.z;
    float yw = att.y*att.w;
    float yz = att.y*att.z;
    float xw = att.x*att.w;
    mat4 model = mat4(1. - 2.*(yy + zz), 2.*(xy + zw), 2.*(xz - yw), 0.,
                      2.*(xy - zw), 1. - 2.*(xx + zz), 2.*(yz + xw), 0.,
                      2.*(xz + yw), 2.*(yz - xw), 1. - 2.*(xx + yy), 0.,
                      model0.w, model1.x, model1.y, 1.);
    mat4 mv = mul(view, model);
    _pos = mul(mv, position).xyz;
    gl_Position = mul(projection, vec4(_pos, 1.));
    _nor = mul(mv, normal).xyz;
    _loc = gl_Position.xyw;
    _lightPos = light0.xyz;
    _intens = vec3(light0.w, light1.xy);
    _falloff = light1.z;
}
