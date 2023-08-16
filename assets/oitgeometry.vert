layout(location = 0) in vec4 position;
layout(location = 1) in vec4 normal;
uniform vec4 model0;
uniform vec2 model1;
uniform mat4 projection;
uniform mat4 view;
out vec4 norCs;
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
    vec4 posCs = mul(mv, position);
    norCs = mul(mv, normal);
    gl_Position = mul(projection, posCs);
}
