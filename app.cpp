#include "object.hpp"
#include "PAZ_Engine"
#include <cmath>

#define _cameraGrounded _cameraObject->grounded()
#define _cameraX _cameraObject->x()
#define _cameraY _cameraObject->y()
#define _cameraZ (_cameraObject->z() + _cameraObject->height())
#define _cameraXVel _cameraObject->xVel()
#define _cameraYVel _cameraObject->yVel()
//#define _cameraZVel _cameraObject->zVel()

static constexpr double Radius = 0.2;
static constexpr double CosMaxAngle = 0.6; // 53.13 deg

////
static paz::Framebuffer _geometryBuffer;
static paz::Framebuffer _renderBuffer;

//static paz::RenderTarget _materialMap(1, paz::TextureFormat::R8UInt, paz::MinMagFilter::Nearest, paz::MinMagFilter::Nearest);
static paz::RenderTarget _normalMap(1, paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Nearest, paz::MinMagFilter::Nearest);
static paz::RenderTarget _directionMap(1, paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Nearest, paz::MinMagFilter::Nearest);
static paz::RenderTarget _depthMap(1, paz::TextureFormat::Depth32Float, paz::MinMagFilter::Nearest, paz::MinMagFilter::Nearest);
static paz::RenderTarget _hdrRender(1, paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Nearest, paz::MinMagFilter::Nearest);

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass;
static paz::RenderPass _postPass;

static paz::VertexBuffer _quadVertices;

static double _cameraPitch = 0.5*3.14159;
static double _cameraYaw = 0.;

static /*const*/ paz::Object* _cameraObject = nullptr;
////

static std::array<float, 16> mat_mult(const std::array<float, 16>& a, const
    std::array<float, 16>& b)
{
    std::array<float, 16> m = {};
    for(std::size_t i = 0; i < 4; ++i)
    {
        for(std::size_t j = 0; j < 4; ++j)
        {
            for(std::size_t k = 0; k < 4; ++k)
            {
                m[4*j + i] += a[4*k + i]*b[4*j + k];
            }
        }
    }
    return m;
}

inline double fract(const double n)
{
    return n - std::floor(n);
}

inline double normalize_angle(const double n)
{
    return fract(n/(2.*3.14159))*2.*3.14159;
}

static const std::string QuadVertSrc = 1 + R"===(
layout(location = 0) in vec2 pos;
out vec2 uv;
void main()
{
    gl_Position = vec4(pos, 0, 1);
    uv = 0.5*pos + 0.5;
}
)===";

static const std::string GeometryVertSrc = 1 + R"====(
layout(location = 0) in vec4 position;
layout(location = 1) in vec4 normal;
//layout(location = 2) in uint material;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
//flat out uint mtl;
out vec4 posCs;
out vec4 norCs;
void main()
{
    mat4 mv = view*model;
//    mtl = material;
    posCs = mv*position;
    norCs = mv*normal;
    gl_Position = projection*posCs;
}
)====";

static const std::string GeometryFragSrc = 1 + R"====(
//flat in uint mtl;
in vec4 posCs;
in vec4 norCs;
//layout(location = 0) out uint material;
layout(location = 0/*1*/) out vec4 normal;
layout(location = 1/*2*/) out vec4 direction;
void main()
{
//    material = mtl;
    normal = norCs;
    direction = vec4(normalize(posCs.xyz), 0);
}
)====";

static const std::string PostFragSrc = 1 + R"===(
const float eps = 1e-6;
uniform sampler2D hdrRender;
uniform float whitePoint;
in vec2 uv;
layout(location = 0) out vec4 color;
float luminance(in vec3 v)
{
    return dot(v, vec3(0.2126, 0.7152, 0.0722));
}
vec3 reinhard(in vec3 col, in float w)
{
    float lOld = max(eps, luminance(col));
    float lNew = lOld*(1. + lOld/max(eps, w*w))/(1. + lOld);
    return col*lNew/lOld;
}
void main()
{
    vec3 hdrCol = texture(hdrRender, uv).rgb;
    color = vec4(pow(reinhard(hdrCol, whitePoint), vec3(0.4545)), 1.);
}
)===";

static constexpr std::array<float, 8> QuadPos = {1, -1, 1, 1, -1, -1, -1, 1};

void paz::App::Init(const std::string& sceneShaderPath)
{
//    _geometryBuffer.attach(_materialMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_directionMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

    const VertexFunction geometryVert(GeometryVertSrc);
    const VertexFunction quadVert(QuadVertSrc);
    const FragmentFunction geometryFrag(GeometryFragSrc);
    const FragmentFunction sceneFrag(load_file(sceneShaderPath).str());
    const FragmentFunction postFrag(PostFragSrc);

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass = RenderPass(_renderBuffer, quadVert, sceneFrag);
    _postPass = RenderPass(quadVert, postFrag);

    _quadVertices.attribute(2, QuadPos);

    Window::SetCursorMode(CursorMode::Disable);
}

void paz::App::Run()
{
    while(!Window::Done())
    {
        physics();

        _cameraYaw = normalize_angle(_cameraYaw - 0.1*Window::MousePos().first*
            Window::FrameTime());
        _cameraPitch = normalize_angle(_cameraPitch + 0.1*Window::MousePos().
            second*Window::FrameTime());
        _cameraPitch = std::max(0.05*3.14159, std::min(0.95*3.14159,
            _cameraPitch));

        const double cosYaw = std::cos(_cameraYaw);
        const double sinYaw = std::sin(_cameraYaw);

        if(_cameraGrounded)
        {
            _cameraXVel = 0.;
            _cameraYVel = 0.;
            if(Window::KeyDown(Key::A))
            {
                _cameraXVel += 3.*-cosYaw;
                _cameraYVel += 3.*-sinYaw;
            }
            if(Window::KeyDown(Key::D))
            {
                _cameraXVel -= 3.*-cosYaw;
                _cameraYVel -= 3.*-sinYaw;
            }
            if(Window::KeyDown(Key::W))
            {
                _cameraXVel += 3.*-sinYaw;
                _cameraYVel += 3.*cosYaw;
            }
            if(Window::KeyDown(Key::S))
            {
                _cameraXVel -= 3.*-sinYaw;
                _cameraYVel -= 3.*cosYaw;
            }
        }

        if(Window::KeyPressed(Key::E)) //TEMP
        {
            for(const auto& n : objects())
            {
                reinterpret_cast<Object*>(n.first)->onInteract(*_cameraObject);
            }
        }

        const auto projection = perspective(1., Window::AspectRatio(), 0.1,
            100.);
        const std::array<float, 16> transMat =
        {
                    1,         0,         0, 0,
                    0,         1,         0, 0,
                    0,         0,         1, 0,
            static_cast<float>(-_cameraX), static_cast<float>(-_cameraY), static_cast<float>(-_cameraZ), 1
        };
        const double cosPitch = std::cos(_cameraPitch);
        const double sinPitch = std::sin(_cameraPitch);
        const std::array<float, 16> pitchMat =
        {
            1,               0,                0, 0,
            0, static_cast<float>(cosPitch), -static_cast<float>(sinPitch), 0,
            0, static_cast<float>(sinPitch),  static_cast<float>(cosPitch), 0,
            0,               0,                0, 1
        };
        const std::array<float, 16> yawMat =
        {
            static_cast<float>(cosYaw), -static_cast<float>(sinYaw), 0, 0,
            static_cast<float>(sinYaw),  static_cast<float>(cosYaw), 0, 0,
                        0,              0, 1, 0,
                        0,              0, 0, 1
        };
        const auto view = mat_mult(pitchMat, mat_mult(yawMat, transMat));

        for(const auto& n : objects())
        {
            Object* o = reinterpret_cast<Object*>(n.first);
            if(o->collisionType() == CollisionType::Default)
            {
                o->zVel() -= 9.81*Window::FrameTime();
            }
        }

        for(const auto& n : objects())
        {
            reinterpret_cast<Object*>(n.first)->update();
        }

////////
for(auto& a0 : objects())
{
    Object* a = reinterpret_cast<Object*>(a0.first);
    if(a->collisionType() != CollisionType::Default)
    {
        continue;
    }

    bool collided = false;

    double minDist = std::numeric_limits<double>::infinity();
    double gx = 0.;
    double gy = 0.;
    double gz = 1.;

    Object* collidedWith = nullptr;
    for(auto& b0 : objects())
    {
        Object* b = reinterpret_cast<Object*>(b0.first);
        if(b->collisionType() != CollisionType::World)
        {
            continue;
        }

        double x = a->x() - b->x();
        double y = a->y() - b->y();
        double z = a->z() - b->z() + Radius;

        double hx, hy, hz;
        const double c = b->model().collide(x, y, z, hx, hy, hz, Radius);
        if(c < minDist)
        {
            collided = true;
            minDist = c;
            gx = hx;
            gy = hy;
            gz = hz;
            collidedWith = b;
        }
    }

    if(collided)
    {
        a->onCollide(*collidedWith);
        collidedWith->onCollide(*a);

        const double normNormal = std::sqrt(gx*gx + gy*gy + gz*gz);
        a->localNorX() = gx/normNormal;
        a->localNorY() = gy/normNormal;
        a->localNorZ() = gz/normNormal;
#if 1
        a->x() += gx;
        a->y() += gy;
        a->z() += gz;
#else
        const double gDotV = gx*a->xVel() + gy*a->yVel() + gz*a->zVel();
        if(gDotV < 0.)
        {
            const double normVel = std::sqrt(a->xVel()*a->xVel() + a->yVel()*a->
                yVel() + a->zVel()*a->zVel());
            x -= gDotV*a->xVel()/normVel/normVel;
            y -= gDotV*a->yVel()/normVel/normVel;
            z -= gDotV*a->zVel()/normVel/normVel;
        }
#endif
        const double norVel = a->xVel()*a->localNorX() + a->yVel()*a->
            localNorY() + a->zVel()*a->localNorZ();
        if(norVel < 0.)
        {
#if 1
            a->xVel() -= norVel*a->localNorX();
            a->yVel() -= norVel*a->localNorY();
            a->zVel() -= norVel*a->localNorZ();
            // friction here ...
#else
            a->xVel() = 0.;
            a->yVel() = 0.;
            a->zVel() = 0.;
#endif
        }
    }

    a->grounded() = collided && a->localNorZ() > CosMaxAngle;
}
////////

        _geometryPass.begin(std::vector<LoadAction>(4, LoadAction::Clear),
            LoadAction::Clear);
        _geometryPass.cull(CullMode::Back);
        _geometryPass.depth(DepthTestMode::Less);
        _geometryPass.uniform("projection", projection);
        _geometryPass.uniform("view", view);
        _geometryPass.uniform("model", std::array<float, 16>{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});
        for(const auto& n : objects())
        {
            const Object* o = reinterpret_cast<const Object*>(n.first);
            if(o->model()._vc)
            {
                _geometryPass.uniform("model", std::array<float, 16>{
1, 0, 0, 0,
0, 1, 0, 0,
0, 0, 1, 0,
static_cast<float>(o->x()), static_cast<float>(o->y()), static_cast<float>(o->z() + o->height()), 1});
                _geometryPass.indexed(PrimitiveType::Triangles, o->model()._vg, o->model()._ig);
            }
        }
        _geometryPass.end();

        _renderPass.begin();
//        _renderPass.read("materialMap", _materialMap);
        _renderPass.read("normalMap", _normalMap);
//        _renderPass.read("directionMap", _directionMap);
//        _renderPass.read("depthMap", _depthMap);
        _renderPass.primitives(PrimitiveType::TriangleStrip, _quadVertices);
        _renderPass.end();

        _postPass.begin();
        _postPass.read("hdrRender", _hdrRender);
        _postPass.uniform("whitePoint", 1.f);
        _postPass.primitives(PrimitiveType::TriangleStrip, _quadVertices);
        _postPass.end();

        Window::EndFrame();
    }
}

void paz::App::AttachCamera(/*const*/ Object& o)
{
    _cameraObject = &o;
}

double paz::App::CameraYaw()
{
    return _cameraYaw;
}
