#include "object.hpp"
#include "PAZ_Engine"
#include <cmath>
#include <limits>

#define _cameraGrounded _cameraObject->grounded()
#define _cameraX _cameraObject->x()
#define _cameraY _cameraObject->y()
#define _cameraZ (_cameraObject->z() + _cameraObject->height())
#define _cameraXVel _cameraObject->xVel()
#define _cameraYVel _cameraObject->yVel()
//#define _cameraZVel _cameraObject->zVel()

static constexpr double Radius = 0.2;
static constexpr double CosMaxAngle = 0.6; // 53.13 deg
static constexpr std::array<float, 4> SunVec = {0.57735, 0.57735, 0.57735, 0.};
static constexpr double InteractRangeBehindSq = 4.;
static constexpr double InteractRangeInFrontSq = 9.;

////
static paz::Framebuffer _geometryBuffer;
static paz::Framebuffer _renderBuffer;
static paz::Framebuffer _postBuffer;
static paz::Framebuffer _lumBuffer;

//static paz::RenderTarget _materialMap(1, paz::TextureFormat::R8UInt);
static paz::RenderTarget _normalMap(1, paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _depthMap(1, paz::TextureFormat::Depth32Float);
static paz::RenderTarget _hdrRender(1, paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _finalRender(1, paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _finalLumMap(1, paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass;
static paz::RenderPass _postPass;
static paz::RenderPass _lumPass;
static paz::RenderPass _fxaaPass;

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

static std::array<float, 4> mat_mult(const std::array<float, 16>& a, const std::
    array<float, 4>& b)
{
    std::array<float, 4> m = {};
    for(std::size_t i = 0; i < 4; ++i)
    {
        for(std::size_t j = 0; j < 4; ++j)
        {
            m[i] += a[4*j + i]*b[j];
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
void main()
{
//    material = mtl;
    normal = vec4(normalize(norCs.xyz), 0.);
}
)====";

static const std::string LumFragSrc = 1 + R"===(
uniform sampler2D img;
in vec2 uv;
layout(location = 0) out float lum;
void main()
{
    lum = dot(texture(img, uv).rgb, vec3(0.2126, 0.7152, 0.0722));
}
)===";

static const std::string FxaaFragSrc = 1 + R"===(
// Linear LDR -> Antialiased gamma-corrected LDR
const float edgeSteps[] = float[](1., 1.5, 2., 2., 2., 2., 2., 2., 2., 4.);
const float edgeGuess = 8.;
uniform sampler2D img;
uniform sampler2D lum;
in vec2 uv;
layout(location = 0) out vec4 color;
struct LumData
{
    float m;
    float n;
    float e;
    float s;
    float w;
    float ne;
    float nw;
    float se;
    float sw;
    float minLum;
    float maxLum;
    float contrast;
};
LumData lum_neighborhood(in sampler2D lum, in vec2 uv)
{
    vec2 texOffset = 1./textureSize(lum, 0);

    LumData l;
    l.m = texture(lum, uv).r;
    l.n = texture(lum, uv + texOffset*vec2(0, 1)).r;
    l.e = texture(lum, uv + texOffset*vec2(1, 0)).r;
    l.s = texture(lum, uv + texOffset*vec2(0, -1)).r;
    l.w = texture(lum, uv + texOffset*vec2(-1, 0)).r;
    l.ne = texture(lum, uv + texOffset*vec2(1, 1)).r;
    l.nw = texture(lum, uv + texOffset*vec2(-1, 1)).r;
    l.se = texture(lum, uv + texOffset*vec2(1, -1)).r;
    l.sw = texture(lum, uv + texOffset*vec2(-1, -1)).r;
    l.minLum = min(min(min(min(l.n, l.e), l.s), l.w), l.m);
    l.maxLum = max(max(max(max(l.n, l.e), l.s), l.w), l.m);
    l.contrast = l.maxLum - l.minLum;
    return l;
}
float blend_fac(in LumData l)
{
    float fac = (2.*(l.n + l.e + l.s + l.w) + l.ne + l.nw + l.se + l.sw)/12.;
    fac = abs(fac - l.m);
    fac = clamp(fac/l.contrast, 0., 1.);
    fac = smoothstep(0., 1., fac);
    return fac*fac;
}
struct EdgeData
{
    bool isHorizontal;
    float pixelStep;
    float oppositeLum;
    float grad;
};
EdgeData determine_edge(in sampler2D lum, in LumData l)
{
    vec2 texOffset = 1./textureSize(lum, 0);

    EdgeData e;

    float horizontal = 2.*abs(l.n + l.s - 2.*l.m) + abs(l.ne + l.se - 2.*l.e) +
        abs(l.nw + l.sw - 2.*l.w);
    float vertical = 2.*abs(l.e + l.w - 2.*l.m) + abs(l.ne + l.nw - 2.*l.n) +
        abs(l.se + l.sw - 2.*l.s);
    e.isHorizontal = horizontal >= vertical;

    float pLum = mix(l.e, l.n, float(e.isHorizontal));
    float nLum = mix(l.w, l.s, float(e.isHorizontal));
    float pGrad = abs(pLum - l.m);
    float nGrad = abs(nLum - l.m);
    e.pixelStep = mix(texOffset.x, texOffset.y, float(e.isHorizontal));

    e.pixelStep *= sign(pGrad - nGrad);
    e.oppositeLum = mix(pLum, nLum, float(pGrad < nGrad));
    e.grad = mix(pGrad, nGrad, float(pGrad < nGrad));

    return e;
}
float edge_blend_fac(in sampler2D lum, in LumData l, in EdgeData e, in vec2 uv)
{
    vec2 texOffset = 1./textureSize(lum, 0);

    vec2 uvEdge = uv + mix(vec2(0.5*e.pixelStep, 0), vec2(0, 0.5*e.pixelStep),
        float(e.isHorizontal));
    vec2 edgeStep = mix(vec2(0, texOffset.y), vec2(texOffset.x, 0), float(e.
        isHorizontal));

    float edgeLum = 0.5*(l.m + e.oppositeLum);
    float gradThresh = 0.25*e.grad;

    vec2 puv = uvEdge + edgeStep;
    float pLumDelta = texture(lum, puv).r - edgeLum;
    bool pAtEnd = abs(pLumDelta) >= gradThresh;
    for(int i = 0; i + 1 < edgeSteps.length() && !pAtEnd; ++i)
    {
        puv += edgeStep*edgeSteps[i];
        pLumDelta = texture(lum, puv).r - edgeLum;
        pAtEnd = abs(pLumDelta) >= gradThresh;
    }
    puv += (1. - float(pAtEnd))*edgeStep*edgeGuess;
    float pDist = mix(puv.y - uv.y, puv.x - uv.x, float(e.isHorizontal));

    vec2 nuv = uvEdge - edgeStep;
    float nLumDelta = texture(lum, nuv).r - edgeLum;
    bool nAtEnd = abs(nLumDelta) >= gradThresh;
    for(int i = 0; i + 1 < edgeSteps.length() && !nAtEnd; ++i)
    {
        nuv -= edgeStep*edgeSteps[i];
        nLumDelta = texture(lum, nuv).r - edgeLum;
        nAtEnd = abs(nLumDelta) >= gradThresh;
    }
    nuv -= (1. - float(nAtEnd))*edgeStep*edgeGuess;
    float nDist = mix(uv.y - nuv.y, uv.x - nuv.x, float(e.isHorizontal));

    float shortestDist = mix(nDist, pDist, float(pDist < nDist));
    bool deltaSign = (pDist < nDist && pLumDelta >= 0.) || (pDist >= nDist &&
        nLumDelta >= 0.);

    return float(deltaSign == (l.m - edgeLum >= 0.))*edgeSteps.length()*
        shortestDist;
}
void main()
{
    LumData l = lum_neighborhood(lum, uv);
    float pixelFac = blend_fac(l);
    EdgeData e = determine_edge(lum, l);
    float edgeFac = edge_blend_fac(lum, l, e, uv);
    float fac = max(pixelFac, edgeFac);
    vec2 deltaUv = mix(vec2(e.pixelStep*fac, 0.), vec2(0., e.pixelStep*fac),
        float(e.isHorizontal));
    color = texture(img, uv + deltaUv);
    color.rgb = pow(color.rgb, vec3(0.4545));
}
)===";

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
    color = vec4(reinhard(hdrCol, whitePoint), 1.);
}
)===";

static constexpr std::array<float, 8> QuadPos = {1, -1, 1, 1, -1, -1, -1, 1};

void paz::App::Init(const std::string& sceneShaderPath)
{
//    _geometryBuffer.attach(_materialMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

    _postBuffer.attach(_finalRender);

    _lumBuffer.attach(_finalLumMap);

    const VertexFunction geometryVert(GeometryVertSrc);
    const VertexFunction quadVert(QuadVertSrc);
    const FragmentFunction geometryFrag(GeometryFragSrc);
    const FragmentFunction sceneFrag(load_file(sceneShaderPath).str());
    const FragmentFunction lumFrag(LumFragSrc);
    const FragmentFunction fxaaFrag(FxaaFragSrc);
    const FragmentFunction postFrag(PostFragSrc);

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass = RenderPass(_renderBuffer, quadVert, sceneFrag);
    _postPass = RenderPass(_postBuffer, quadVert, postFrag);
    _lumPass = RenderPass(_lumBuffer, quadVert, lumFrag);
    _fxaaPass = RenderPass(quadVert, fxaaFrag);

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

        if(Window::KeyPressed(Key::E))
        {
            for(const auto& n : objects())
            {
                Object* o = reinterpret_cast<Object*>(n.first);
                const double relX = o->x() - _cameraX;
                const double relY = o->y() - _cameraY;
                const double relZ = o->z() - _cameraZ;
                const double distSq = relX*relX + relY*relY + relZ*relZ;
                const double dotProd = -sinYaw*relX + cosYaw*relY;
                if((distSq < InteractRangeBehindSq || (distSq <
                    InteractRangeInFrontSq && dotProd > 0.)))
                {

                    o->onInteract(*_cameraObject);
                }
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

        // Get geometry map.
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
            if(!o->model()._i.empty())
            {
                _geometryPass.uniform("model", std::array<float, 16>{1, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 1, 0, static_cast<float>(o->x()),
                    static_cast<float>(o->y()), static_cast<float>(o->z() + o->
                    height()), 1});
                _geometryPass.indexed(PrimitiveType::Triangles, o->model()._v,
                    o->model()._i);
            }
        }
        _geometryPass.end();

        // Render in HDR.
        _renderPass.begin();
//        _renderPass.read("materialMap", _materialMap);
        _renderPass.read("normalMap", _normalMap);
        _renderPass.read("depthMap", _depthMap);
        _renderPass.uniform("projection", projection);
        _renderPass.uniform("sun", mat_mult(view, SunVec));
        _renderPass.primitives(PrimitiveType::TriangleStrip, _quadVertices);
        _renderPass.end();

        // Tonemap to linear LDR.
        _postPass.begin();
        _postPass.read("hdrRender", _hdrRender);
        _postPass.uniform("whitePoint", 1.f);
        _postPass.primitives(PrimitiveType::TriangleStrip, _quadVertices);
        _postPass.end();

        // Get luminance map.
        _lumPass.begin();
        _lumPass.read("img", _finalRender);
        _lumPass.primitives(PrimitiveType::TriangleStrip, _quadVertices);
        _lumPass.end();

        // Anti-alias and correct gamma.
        _fxaaPass.begin();
        _fxaaPass.read("img", _finalRender);
        _fxaaPass.read("lum", _finalLumMap);
        _fxaaPass.primitives(PrimitiveType::TriangleStrip, _quadVertices);
        _fxaaPass.end();

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
