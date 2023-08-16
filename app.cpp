#include "object.hpp"
#include "shared.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>

//#define DO_FXAA

#define _cameraBasePos Vec{{_cameraObject->x(), _cameraObject->y(), _cameraObject->z()}}
#define _cameraAtt Vec{{_cameraObject->xAtt(), _cameraObject->yAtt(), _cameraObject->zAtt(), std::sqrt(1. - _cameraObject->xAtt()*_cameraObject->xAtt() - _cameraObject->yAtt()*_cameraObject->yAtt() - _cameraObject->zAtt()*_cameraObject->zAtt())}}

static constexpr double CosMaxAngle = 0.6; // 53.13 deg
static constexpr std::size_t NumSteps = 100;
static const paz::Vec SunVec{{0.57735, 0.57735, 0.57735, 0.}};
static constexpr double InteractRangeBehindSq = 4.;
static constexpr double InteractRangeInFrontSq = 9.;

static paz::Framebuffer _geometryBuffer;
static paz::Framebuffer _renderBuffer;
#ifdef DO_FXAA
static paz::Framebuffer _postBuffer;
static paz::Framebuffer _lumBuffer;
#endif

static paz::RenderTarget _materialMap(1, paz::TextureFormat::R8UInt);
static paz::RenderTarget _normalMap(1, paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _depthMap(1, paz::TextureFormat::Depth32Float);
static paz::RenderTarget _hdrRender(1, paz::TextureFormat::RGBA16Float);
#ifdef DO_FXAA
static paz::RenderTarget _finalRender(1, paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _finalLumMap(1, paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
#endif

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass;
static paz::RenderPass _postPass;
#ifdef DO_FXAA
static paz::RenderPass _lumPass;
static paz::RenderPass _fxaaPass;
#endif

static paz::VertexBuffer _quadVertices;

static double _cameraPitch = 0.;

static /*const*/ paz::Object* _cameraObject = nullptr;

static paz::Vec _mousePos = paz::Vec::Zero(2);

static paz::Mat convert_mat(const std::array<float, 16>& m)
{
    paz::Mat res(4, 4);
    std::copy(m.begin(), m.end(), res.begin());
    return res;
}

static std::array<float, 16> convert_mat(const paz::Mat& m)
{
    if(m.rows() != 4 || m.cols() != 4)
    {
        throw std::runtime_error("Must be a 4x4 matrix.");
    }
    std::array<float, 16> res;
    std::copy(m.begin(), m.end(), res.begin());
    return res;
}

static std::array<float, 4> convert_vec(const paz::Vec& v)
{
    if(v.size() != 4)
    {
        throw std::runtime_error("Must be a four-vector.");
    }
    std::array<float, 4> res;
    std::copy(v.begin(), v.end(), res.begin());
    return res;
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
layout(location = 2) in uint material;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
flat out uint mtl;
out vec4 posCs;
out vec4 norCs;
void main()
{
    mat4 mv = view*model;
    mtl = material;
    posCs = mv*position;
    norCs = mv*normal;
    gl_Position = projection*posCs;
}
)====";

static const std::string GeometryFragSrc = 1 + R"====(
flat in uint mtl;
in vec4 posCs;
in vec4 norCs;
layout(location = 0) out uint material;
layout(location = 1) out vec4 normal;
void main()
{
    material = mtl;
    normal = vec4(normalize(norCs.xyz), 0.);
}
)====";

#ifdef DO_FXAA
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
const int numEdgeSteps = 10;
const float edgeSteps[] = {1., 1.5, 2., 2., 2., 2., 2., 2., 2., 4.};
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
LumData lum_neighborhood(in sampler2D lum, in vec2 texOffset, in vec2 uv)
{
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
EdgeData determine_edge(in sampler2D lum, in vec2 texOffset, in LumData l)
{
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
float edge_blend_fac(in sampler2D lum, in vec2 texOffset, in LumData l, in
    EdgeData e, in vec2 uv)
{
    vec2 uvEdge = uv + mix(vec2(0.5*e.pixelStep, 0), vec2(0, 0.5*e.pixelStep),
        float(e.isHorizontal));
    vec2 edgeStep = mix(vec2(0, texOffset.y), vec2(texOffset.x, 0), float(e.
        isHorizontal));

    float edgeLum = 0.5*(l.m + e.oppositeLum);
    float gradThresh = 0.25*e.grad;

    vec2 puv = uvEdge + edgeStep;
    float pLumDelta = texture(lum, puv).r - edgeLum;
    bool pAtEnd = abs(pLumDelta) >= gradThresh;
    for(int i = 0; i + 1 < numEdgeSteps && !pAtEnd; ++i)
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
    for(int i = 0; i + 1 < numEdgeSteps && !nAtEnd; ++i)
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

    return float(deltaSign == (l.m - edgeLum >= 0.))*numEdgeSteps*shortestDist;
}
void main()
{
    ivec2 texSize = textureSize(lum, 0);
    vec2 texOffset = vec2(1./texSize.x, 1./texSize.y);
    LumData l = lum_neighborhood(lum, texOffset, uv);
    float pixelFac = blend_fac(l);
    EdgeData e = determine_edge(lum, texOffset, l);
    float edgeFac = edge_blend_fac(lum, texOffset, l, e, uv);
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
#else
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
    color.rgb = pow(color.rgb, vec3(0.4545));
}
)===";
#endif

static constexpr std::array<float, 8> QuadPos = {1, -1, 1, 1, -1, -1, -1, 1};

void paz::App::Init(const std::string& sceneShaderPath)
{
    _geometryBuffer.attach(_materialMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

#ifdef DO_FXAA
    _postBuffer.attach(_finalRender);

    _lumBuffer.attach(_finalLumMap);
#endif

    const VertexFunction geometryVert(GeometryVertSrc);
    const VertexFunction quadVert(QuadVertSrc);
    const FragmentFunction geometryFrag(GeometryFragSrc);
    const FragmentFunction sceneFrag(getAsset(sceneShaderPath).str());
#ifdef DO_FXAA
    const FragmentFunction lumFrag(LumFragSrc);
    const FragmentFunction fxaaFrag(FxaaFragSrc);
#endif
    const FragmentFunction postFrag(PostFragSrc);

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass = RenderPass(_renderBuffer, quadVert, sceneFrag);
#ifdef DO_FXAA
    _postPass = RenderPass(_postBuffer, quadVert, postFrag);
    _lumPass = RenderPass(_lumBuffer, quadVert, lumFrag);
    _fxaaPass = RenderPass(quadVert, fxaaFrag);
#else
    _postPass = RenderPass(quadVert, postFrag);
#endif

    _quadVertices.attribute(2, QuadPos);

    Window::SetCursorMode(CursorMode::Disable);
}

void paz::App::Run()
{
    while(!Window::Done())
    {
        physics();
        gravity();

        // Identify all objects that can collide and precompute as much as possible.
        std::vector<Object*> a;
        std::vector<Object*> b;
        for(auto& n : objects())
        {
            Object* o = reinterpret_cast<Object*>(n.first);
            if(o->collisionType() == CollisionType::Default)
            {
                a.push_back(o);
            }
            if(o->collisionType() == CollisionType::World)
            {
                b.push_back(o);
            }
        }

        std::vector<double> offsetX(a.size());
        std::vector<double> offsetY(a.size());
        std::vector<double> offsetZ(a.size());
        for(std::size_t i = 0; i < a.size(); ++i)
        {
            const double q0 = a[i]->xAtt();
            const double q1 = a[i]->yAtt();
            const double q2 = a[i]->zAtt();
            const auto q3 = -std::sqrt(1. - q0*q0 - q1*q1 - q2*q2);
            const auto xx = q0*q0;
            const auto yy = q1*q1;
            const auto xz = q0*q2;
            const auto yw = q1*q3;
            const auto yz = q1*q2;
            const auto xw = q0*q3;
            offsetX[i] = a[i]->collisionRadius()*2.*(xz + yw);
            offsetY[i] = a[i]->collisionRadius()*2.*(yz - xw);
            offsetZ[i] = a[i]->collisionRadius()*(1. - 2.*(xx + yy));
        }

        std::vector<std::vector<double>> bX(b.size(), std::vector<double>(NumSteps));
        std::vector<std::vector<double>> bY(b.size(), std::vector<double>(NumSteps));
        std::vector<std::vector<double>> bZ(b.size(), std::vector<double>(NumSteps));
        for(std::size_t i = 0; i < b.size(); ++i)
        {
            for(std::size_t j = 0; j < NumSteps; ++j)
            {
                bX[i][j] = b[i]->xPrev() + (j + 1)*(b[i]->x() - b[i]->xPrev())/NumSteps;
                bY[i][j] = b[i]->yPrev() + (j + 1)*(b[i]->y() - b[i]->yPrev())/NumSteps;
                bZ[i][j] = b[i]->zPrev() + (j + 1)*(b[i]->z() - b[i]->zPrev())/NumSteps;
            }
        }

        // Find and handle collisions. (Time is the outer loop to ensure collision repsonses occur in the correct order.)
std::vector<bool> tempDone(a.size(), false);
        for(std::size_t i = 0; i < NumSteps; ++i)
        {
            for(std::size_t j = 0; j < a.size(); ++j)
            {
if(tempDone[j]){ continue; }
                for(std::size_t k = 0; k < b.size(); ++k)
                {
                    const double x = offsetX[i] + a[j]->xPrev() + (i + 1)*(a[j]->x() - a[j]->xPrev())/NumSteps - bX[k][i];
                    const double y = offsetY[i] + a[j]->yPrev() + (i + 1)*(a[j]->y() - a[j]->yPrev())/NumSteps - bY[k][i];
                    const double z = offsetZ[i] + a[j]->zPrev() + (i + 1)*(a[j]->z() - a[j]->zPrev())/NumSteps - bZ[k][i];

                    double xNew, yNew, zNew, xNor, yNor, zNor;
                    const double c = b[k]->model().collide(x, y, z, a[j]->collisionRadius(), xNew, yNew, zNew, xNor, yNor, zNor);
                    if(c < a[j]->collisionRadius())
                    {
                        const double xVel = a[j]->xVel() - b[k]->xVel();
                        const double yVel = a[j]->yVel() - b[k]->yVel();
                        const double zVel = a[j]->zVel() - b[k]->zVel();
                        const double norVel = xVel*xNor + yVel*yNor + zVel*zNor;
                        if(norVel < 0.)
                        {
                            a[j]->xVel() -= norVel*xNor;
                            a[j]->yVel() -= norVel*yNor;
                            a[j]->zVel() -= norVel*zNor;

                            // Apply friction.
                            a[j]->xVel() = 0.5*a[j]->xVel() + 0.5*b[k]->xVel();
                            a[j]->yVel() = 0.5*a[j]->yVel() + 0.5*b[k]->yVel();
                            a[j]->zVel() = 0.5*a[j]->zVel() + 0.5*b[k]->zVel();
                        }
                        a[j]->x() = xNew - offsetX[j] + (NumSteps - i - 1)*a[j]->xVel()/NumSteps*paz::Window::FrameTime();
                        a[j]->y() = yNew - offsetY[j] + (NumSteps - i - 1)*a[j]->yVel()/NumSteps*paz::Window::FrameTime();
                        a[j]->z() = zNew - offsetZ[j] + (NumSteps - i - 1)*a[j]->zVel()/NumSteps*paz::Window::FrameTime();
                        a[j]->localNorX() = xNor;
                        a[j]->localNorY() = yNor;
                        a[j]->localNorZ() = zNor;
                        a[j]->grounded() = true;// > CosMaxAngle;
                        a[j]->onCollide(*b[k]);
                        b[k]->onCollide(*a[j]);
tempDone[j] = true;
                    }
                }
            }
        }
//std::cout << static_cast<bool>(_cameraObject->grounded()) << " " << std::sqrt(_cameraObject->xVel()*_cameraObject->xVel() +_cameraObject->yVel()*_cameraObject->yVel() +_cameraObject->zVel()*_cameraObject->zVel()) << std::endl;
{
const double r = std::sqrt(_cameraObject->x()*_cameraObject->x() + _cameraObject->y()*_cameraObject->y() + _cameraObject->z()*_cameraObject->z());
const double lat = std::asin(_cameraObject->z()/r);
const double lon = std::atan2(_cameraObject->y(), _cameraObject->x());
std::cout << lat*180./M_PI << " " << lon*180./M_PI << std::endl;
}


        update();

        const Vec gravDir = -_cameraBasePos.normalized();

        Vec cameraAtt = _cameraAtt;

        // Get basis vectors (rows of rotation matrix).
        Mat cameraRot = to_mat(cameraAtt);
        Vec cameraForward = cameraRot.row(1).trans();
        Vec cameraRight = gravDir.cross(cameraForward).normalized();

        _cameraObject->yAngRate() = 0.;
        if(_cameraObject->grounded())
        {
            _mousePos = Vec::Zero(2);
            _cameraObject->xAngRate() = 0.;
            _cameraObject->zAngRate() = -0.1*Window::MousePos().first;
            const Vec baseForward = cameraRight.cross(gravDir).normalized();
            const double deltaPitch = std::acos(std::max(0., std::min(1.,
                baseForward.dot(cameraForward))))*(cameraForward.dot(gravDir) >
                0. ? -1. : 1.) + 0.1*Window::MousePos().second*Window::
                FrameTime();
            if(std::abs(deltaPitch) > 1e-6)
            {
                _cameraPitch = std::max(-0.45*M_PI, std::min(0.45*
                    M_PI, _cameraPitch + deltaPitch));
            }
            cameraRot.setRow(0, cameraRight.trans());
            cameraRot.setRow(1, baseForward.trans());
            cameraRot.setRow(2, -gravDir.trans());
            Vec cameraBaseAtt = to_quat(cameraRot);
            if(cameraBaseAtt(3) < 0.)
            {
                cameraBaseAtt = -cameraBaseAtt;
            }
            _cameraObject->xAtt() = cameraBaseAtt(0);
            _cameraObject->yAtt() = cameraBaseAtt(1);
            _cameraObject->zAtt() = cameraBaseAtt(2);
            cameraAtt = qmult(axis_angle(Vec{{1, 0, 0}}, _cameraPitch + 0.5*
                M_PI), _cameraAtt);
        }
        else
        {
            if(_cameraPitch)
            {
                cameraAtt = qmult(axis_angle(Vec{{1, 0, 0}}, _cameraPitch),
                    cameraAtt);
                _cameraObject->xAtt() = cameraAtt(0);
                _cameraObject->yAtt() = cameraAtt(1);
                _cameraObject->zAtt() = cameraAtt(2);
                cameraRot = to_mat(_cameraAtt);
                _cameraPitch = 0.;
            }
            cameraAtt = qmult(axis_angle(Vec{{1, 0, 0}}, 0.5*M_PI), cameraAtt);
            _mousePos(0) += Window::MousePos().first;
            _mousePos(1) += Window::MousePos().second;
            const double norm = _mousePos.norm();
            if(norm > 100.)
            {
                _mousePos *= 100./norm;
            }
            else if(norm > 0.1)
            {
                _mousePos -= 50.0/norm*Window::FrameTime()*_mousePos;
            }
            else
            {
                _mousePos = Vec::Zero(2);
            }
            _cameraObject->xAngRate() = 0.006*_mousePos(1);
            _cameraObject->zAngRate() = -0.006*_mousePos(0);
        }

        cameraRight = cameraRot.row(0).trans();
        cameraForward = cameraRot.row(1).trans();
        const Vec cameraUp = cameraRot.row(2).trans();

        if(_cameraObject->grounded())
        {
            if(Window::KeyDown(Key::A))
            {
                _cameraObject->xVel() -= 3.*cameraRight(0);
                _cameraObject->yVel() -= 3.*cameraRight(1);
                _cameraObject->zVel() -= 3.*cameraRight(2);
            }
            if(Window::KeyDown(Key::D))
            {
                _cameraObject->xVel() += 3.*cameraRight(0);
                _cameraObject->yVel() += 3.*cameraRight(1);
                _cameraObject->zVel() += 3.*cameraRight(2);
            }
            if(Window::KeyDown(Key::W))
            {
                _cameraObject->xVel() += 3.*cameraForward(0);
                _cameraObject->yVel() += 3.*cameraForward(1);
                _cameraObject->zVel() += 3.*cameraForward(2);
            }
            if(Window::KeyDown(Key::S))
            {
                _cameraObject->xVel() -= 3.*cameraForward(0);
                _cameraObject->yVel() -= 3.*cameraForward(1);
                _cameraObject->zVel() -= 3.*cameraForward(2);
            }
        }
        else
        {
            if(Window::KeyDown(Key::A))
            {
                _cameraObject->xVel() -= 12.*cameraRight(0)*Window::FrameTime();
                _cameraObject->yVel() -= 12.*cameraRight(1)*Window::FrameTime();
                _cameraObject->zVel() -= 12.*cameraRight(2)*Window::FrameTime();
            }
            if(Window::KeyDown(Key::D))
            {
                _cameraObject->xVel() += 12.*cameraRight(0)*Window::FrameTime();
                _cameraObject->yVel() += 12.*cameraRight(1)*Window::FrameTime();
                _cameraObject->zVel() += 12.*cameraRight(2)*Window::FrameTime();
            }
            if(Window::KeyDown(Key::W))
            {
                _cameraObject->xVel() += 12.*cameraForward(0)*Window::FrameTime();
                _cameraObject->yVel() += 12.*cameraForward(1)*Window::FrameTime();
                _cameraObject->zVel() += 12.*cameraForward(2)*Window::FrameTime();
            }
            if(Window::KeyDown(Key::S))
            {
                _cameraObject->xVel() -= 12.*cameraForward(0)*Window::FrameTime();
                _cameraObject->yVel() -= 12.*cameraForward(1)*Window::FrameTime();
                _cameraObject->zVel() -= 12.*cameraForward(2)*Window::FrameTime();
            }
        }
        if(Window::KeyDown(Key::LeftShift))
        {
            _cameraObject->xVel() += 12.*cameraUp(0)*Window::FrameTime();
            _cameraObject->yVel() += 12.*cameraUp(1)*Window::FrameTime();
            _cameraObject->zVel() += 12.*cameraUp(2)*Window::FrameTime();
        }
        if(Window::KeyDown(Key::LeftControl))
        {
            _cameraObject->xVel() -= 12.*cameraUp(0)*Window::FrameTime();
            _cameraObject->yVel() -= 12.*cameraUp(1)*Window::FrameTime();
            _cameraObject->zVel() -= 12.*cameraUp(2)*Window::FrameTime();
        }

        /*if(Window::KeyPressed(Key::E))
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
        }*/

        const auto projection = perspective(1., Window::AspectRatio(), 0.1,
            1e3);
        Mat view = Mat::Zero(4);
        view(3, 3) = 1.;
        view.setBlock(0, 0, 3, 3, to_mat(cameraAtt));

        // Get geometry map.
        _geometryPass.begin(std::vector<LoadAction>(4, LoadAction::Clear),
            LoadAction::Clear);
        _geometryPass.cull(CullMode::Back);
        _geometryPass.depth(DepthTestMode::Less);
        _geometryPass.uniform("projection", projection);
        _geometryPass.uniform("view", convert_mat(view));
        _geometryPass.uniform("model", convert_mat(Mat::Identity(4)));
        for(const auto& n : objects())
        {
            const Object* o = reinterpret_cast<const Object*>(n.first);
            if(!o->model()._i.empty())
            {
                const double q0 = o->xAtt();
                const double q1 = o->yAtt();
                const double q2 = o->zAtt();
                const auto q3 = -std::sqrt(1. - q0*q0 - q1*q1 - q2*q2);
                const auto xx = q0*q0;
                const auto yy = q1*q1;
                const auto xz = q0*q2;
                const auto yw = q1*q3;
                const auto yz = q1*q2;
                const auto xw = q0*q3;
                const auto upX = 2.*(xz + yw);
                const auto upY = 2.*(yz - xw);
                const auto upZ = 1. - 2.*(xx + yy);
                Mat model = Mat::Identity(4);
                model(0, 3) = o->x() + o->height()*upX - _cameraObject->x() - _cameraObject->height()*cameraUp(0);
                model(1, 3) = o->y() + o->height()*upY - _cameraObject->y() - _cameraObject->height()*cameraUp(1);
                model(2, 3) = o->z() + o->height()*upZ - _cameraObject->z() - _cameraObject->height()*cameraUp(2);
                model.setBlock(0, 0, 3, 3, to_mat(Vec{{q0, q1, q2, q3}}));
                _geometryPass.uniform("model", convert_mat(model));
                _geometryPass.draw(PrimitiveType::Triangles, o->model()._v, o->
                    model()._i);
            }
        }
        _geometryPass.end();

        // Render in HDR.
        _renderPass.begin();
        _renderPass.read("materialMap", _materialMap);
        _renderPass.read("normalMap", _normalMap);
        _renderPass.read("depthMap", _depthMap);
        _renderPass.uniform("invProjection", convert_mat(convert_mat(projection).inv()));
        _renderPass.uniform("sun", convert_vec(view*SunVec));
        _renderPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _renderPass.end();

        // Tonemap to linear LDR.
        _postPass.begin();
        _postPass.read("hdrRender", _hdrRender);
        _postPass.uniform("whitePoint", 1.f);
        _postPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _postPass.end();

#ifdef DO_FXAA
        // Get luminance map.
        _lumPass.begin();
        _lumPass.read("img", _finalRender);
        _lumPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _lumPass.end();

        // Anti-alias and correct gamma.
        _fxaaPass.begin();
        _fxaaPass.read("img", _finalRender);
        _fxaaPass.read("lum", _finalLumMap);
        _fxaaPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _fxaaPass.end();
#endif

        Window::EndFrame();
    }
}

void paz::App::AttachCamera(/*const*/ Object& o)
{
    _cameraObject = &o;
}
