layout(location = 0) in vec4 position;
layout(location = 1) in vec4 light0 [[instance]];
layout(location = 2) in vec4 light1 [[instance]];
uniform mat4 projection;
out vec3 lightPos;
out vec3 intens;
out float falloff;
const float minIll = 1e-3;
void main()
{
    lightPos = light0.xyz;
    intens = vec3(light0.w, light1.xy);
    falloff = light1.z;
    float totalIntens = dot(intens, vec3(0.2126, 0.7152, 0.0722));
    float r = log(totalIntens/minIll)/falloff;
    gl_Position = mul(projection, vec4(r*position.xyz + lightPos, 1.));
}
