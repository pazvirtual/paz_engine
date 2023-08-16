#include "object.hpp"
#include "shared.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include "timer.hpp"
#include <limits>
#include <iomanip>
#include <deque>
#include <regex>

//#define DO_FXAA
#define NO_FRICTION

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

static paz::RenderTarget _materialMap(1., paz::TextureFormat::R8UInt);
static paz::RenderTarget _normalMap(1., paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _coordMap(1., paz::TextureFormat::RG16Float);
static paz::RenderTarget _depthMap(1., paz::TextureFormat::Depth32Float);
static paz::RenderTarget _hdrRender(1., paz::TextureFormat::RGBA16Float);
#ifdef DO_FXAA
static paz::RenderTarget _finalRender(1., paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _finalLumMap(1., paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
#endif

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass;
static paz::RenderPass _postPass;
#ifdef DO_FXAA
static paz::RenderPass _lumPass;
static paz::RenderPass _fxaaPass;
#endif
static paz::RenderPass _textPass;
static paz::RenderPass _plotPass;

static paz::VertexBuffer _quadVertices;

static paz::Texture _font;
static std::stringstream _msgStream;

static const paz::Object* _cameraObject = nullptr;

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
layout(location = 3) in vec2 coord;
layout(location = 4) in vec4 model0;
layout(location = 5) in vec4 model1;
layout(location = 6) in vec4 model2;
layout(location = 7) in vec4 model3;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
flat out uint mtl;
out vec4 posCs;
out vec4 norCs;
out vec2 uv;
void main()
{
    mat4 model = mat4(model0, model1, model2, model3);
    mat4 mv = view*model;
    mtl = material;
    posCs = mv*position;
    norCs = mv*normal;
    gl_Position = projection*posCs;
    uv = coord;
}
)====";

static const std::string GeometryFragSrc = 1 + R"====(
flat in uint mtl;
in vec4 posCs;
in vec4 norCs;
in vec2 uv;
layout(location = 0) out uint material;
layout(location = 1) out vec4 normal;
layout(location = 2) out vec2 coord;
void main()
{
    material = mtl;
    normal = vec4(normalize(norCs.xyz), 0.);
    coord = uv;
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

static const std::string TextVertSrc = 1 + R"===(
const int charWidth = 5;
uniform int baseWidth;
uniform int baseHeight;
uniform float scale;
uniform int height;
uniform int width;
uniform int character;
uniform int col;
uniform int row;
layout(location = 0) in vec2 position;
out vec2 uv;
void main()
{
    uv = 0.5*position.xy + 0.5;
    gl_Position = vec4((uv.x +
        0.2 + //TEMP
        float(col)
        *1.2 //TEMP
        )*scale*float(charWidth)*2./float(width) - 1., (uv.y - 1.
        - 0.1 //TEMP
        - float(row))*scale*float(baseHeight)*2./float(height) + 1., 0, 1);
    uv.x = (uv.x + float(character))*float(charWidth)/float(baseWidth);
}
)===";

static const std::string TextFragSrc = 1 + R"===(
uniform sampler2D font;
uniform float highlight;
in vec2 uv;
layout(location = 0) out vec4 color;
void main()
{
    float t = texture(font, uv).x;
    float c = 1. - 0.8*highlight;
    color = vec4(1., 1., c, t);
}
)===";

static const std::string PlotVertSrc = 1 + R"===(
layout(location = 0) in vec2 pos;
void main()
{
    gl_Position = vec4(pos, 0, 1);
}
)===";

static const std::string PlotFragSrc = 1 + R"===(
layout(location = 0) out vec4 color;
void main()
{
    color = vec4(1.);
}
)===";

static constexpr std::array<float, 8> QuadPos = {1, -1, 1, 1, -1, -1, -1, 1};

static paz::Texture get_asset_image(const std::string& path)
{
    static const std::regex bmp("bmp", std::regex_constants::icase);
    static const std::regex pbm("pbm", std::regex_constants::icase);
    //static const std::regex jpg("jpe?g", std::regex_constants::icase);
    //static const std::regex png("png", std::regex_constants::icase);
    const std::string ext = paz::split_path(path)[2];
    if(std::regex_match(ext, bmp))
    {
        return paz::Texture(paz::parse_bmp(paz::get_asset(path)));
    }
    if(std::regex_match(ext, pbm))
    {
        return paz::Texture(paz::parse_pbm(paz::get_asset(path)));
    }
    //if(std::regex_match(ext, jpg))
    //{
    //    return paz::Texture(paz::parse_jpg(paz::get_asset(path)));
    //}
    //if(std::regex_match(ext, png))
    //{
    //    return paz::Texture(paz::parse_png(paz::get_asset(path)));
    //}
    throw std::runtime_error("Unrecognized image extension \"" + ext + "\".");
}


void paz::App::Init(const std::string& sceneShaderPath, const std::string&
    fontPath)
{
    _geometryBuffer.attach(_materialMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_coordMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

#ifdef DO_FXAA
    _postBuffer.attach(_finalRender);

    _lumBuffer.attach(_finalLumMap);
#endif

    const VertexFunction geometryVert(GeometryVertSrc);
    const VertexFunction quadVert(QuadVertSrc);
    const VertexFunction textVert(TextVertSrc);
    const VertexFunction plotVert(PlotVertSrc);
    const FragmentFunction geometryFrag(GeometryFragSrc);
    const FragmentFunction sceneFrag(get_asset(sceneShaderPath).str());
#ifdef DO_FXAA
    const FragmentFunction lumFrag(LumFragSrc);
    const FragmentFunction fxaaFrag(FxaaFragSrc);
#endif
    const FragmentFunction postFrag(PostFragSrc);
    const FragmentFunction textFrag(TextFragSrc);
    const FragmentFunction plotFrag(PlotFragSrc);

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass = RenderPass(_renderBuffer, quadVert, sceneFrag);
#ifdef DO_FXAA
    _postPass = RenderPass(_postBuffer, quadVert, postFrag);
    _lumPass = RenderPass(_lumBuffer, quadVert, lumFrag);
    _fxaaPass = RenderPass(quadVert, fxaaFrag);
#else
    _postPass = RenderPass(quadVert, postFrag);
#endif
    _textPass = RenderPass(textVert, textFrag, BlendMode::Blend);
    _plotPass = RenderPass(plotVert, plotFrag, BlendMode::Blend);

    _quadVertices.attribute(2, QuadPos);

    _font = Texture(get_asset_image(fontPath)); //TEMP - note that only red channel is used

    Window::SetCursorMode(CursorMode::Disable);
}

void paz::App::Run()
{
    Window::MakeFullscreen();
    Window::DisableHidpi();

    while(!Window::Done())
    {
        physics();

paz::Timer timer;
        // Identify all objects that can collide and precompute as much as
        // possible.
        std::vector<Object*> a;
        std::vector<Object*> b;
        a.reserve(objects().size());
        b.reserve(objects().size());
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

        // `a[i]` may collide with any `b[c[i][j]]`.
        std::vector<std::vector<std::size_t>> c(a.size());
        for(std::size_t i = 0; i < a.size(); ++i)
        {
            c[i].reserve(b.size());
            for(std::size_t j = 0; j < b.size(); ++j)
            {
                if(b[j]->model().sweepVol(a[i]->xPrev(), a[i]->yPrev(), a[i]->
                    zPrev(), a[i]->x(), a[i]->y(), a[i]->z(), b[j]->xPrev(), b[
                    j]->yPrev(), b[j]->zPrev(), b[j]->x(), b[j]->y(), b[j]->z(),
                    a[i]->collisionRadius()))
                {
                    c[i].push_back(j);
                }
            }
        }
const double colTime0 = timer.getAndRestart();

        std::vector<Mat> bRot(b.size());
        std::vector<std::vector<double>> bX(b.size(), std::vector<double>(
            NumSteps));
        std::vector<std::vector<double>> bY(b.size(), std::vector<double>(
            NumSteps));
        std::vector<std::vector<double>> bZ(b.size(), std::vector<double>(
            NumSteps));
        for(std::size_t i = 0; i < b.size(); ++i)
        {
            const double wAtt = std::sqrt(1. - b[i]->xAtt()*b[i]->xAtt() - b[
                i]->yAtt()*b[i]->yAtt() - b[i]->zAtt()*b[i]->zAtt());
            const Vec att{{b[i]->xAtt(), b[i]->yAtt(), b[i]->zAtt(), wAtt}};
            bRot[i] = to_mat(att);
            for(std::size_t j = 0; j < NumSteps; ++j)
            {
                static constexpr double t0 = 1./NumSteps;
                bX[i][j] = b[i]->xPrev() + (j + 1)*(b[i]->x() - b[i]->xPrev())*
                    t0;
                bY[i][j] = b[i]->yPrev() + (j + 1)*(b[i]->y() - b[i]->yPrev())*
                    t0;
                bZ[i][j] = b[i]->zPrev() + (j + 1)*(b[i]->z() - b[i]->zPrev())*
                    t0;
            }
        }

        // Find and handle collisions. (Time is the outer loop to ensure
        // collision repsonses occur in the correct order.)
std::vector<bool> tempDone(a.size(), false);
        for(std::size_t i = 0; i < NumSteps; ++i)
        {
            for(std::size_t j = 0; j < a.size(); ++j)
            {
if(tempDone[j]){ continue; }
                for(auto n : c[j])
                {
                    double x = a[j]->xPrev() + (i + 1)*(a[j]->x() - a[j]->
                        xPrev())/NumSteps - bX[n][i];
                    double y = a[j]->yPrev() + (i + 1)*(a[j]->y() - a[j]->
                        yPrev())/NumSteps - bY[n][i];
                    double z = a[j]->zPrev() + (i + 1)*(a[j]->z() - a[j]->
                        zPrev())/NumSteps - bZ[n][i];
                    const Vec relPos = bRot[n]*Vec{{x, y, z}};
                    x = relPos(0);
                    y = relPos(1);
                    z = relPos(2);

                    double xNew, yNew, zNew, xNor, yNor, zNor;
                    const double dist = b[n]->model().collide(x, y, z, a[j]->
                        collisionRadius(), xNew, yNew, zNew, xNor, yNor, zNor);
                    if(dist < a[j]->collisionRadius())
                    {
                        const Vec nor = bRot[n].trans()*Vec{{xNor, yNor, zNor}};
                        xNor = nor(0);
                        yNor = nor(1);
                        zNor = nor(2);
                        const Vec newPos = bRot[n].trans()*Vec{{xNew, yNew,
                            zNew}};
                        xNew = newPos(0);
                        yNew = newPos(1);
                        zNew = newPos(2);
                        const double xVel = a[j]->xVel() - b[n]->xVel();
                        const double yVel = a[j]->yVel() - b[n]->yVel();
                        const double zVel = a[j]->zVel() - b[n]->zVel();
                        const double norVel = xVel*xNor + yVel*yNor + zVel*zNor;
                        if(norVel < 0.)
                        {
                            a[j]->xVel() -= norVel*xNor;
                            a[j]->yVel() -= norVel*yNor;
                            a[j]->zVel() -= norVel*zNor;

                            // Apply friction.
#ifndef NO_FRICTION
                            a[j]->xVel() = b[n]->xVel();
                            a[j]->yVel() = b[n]->yVel();
                            a[j]->zVel() = b[n]->zVel();
#endif
                        }
                        a[j]->x() = xNew + bX[n][i] + (NumSteps - i - 1)*a[j]->
                            xVel()/NumSteps*Window::FrameTime();
                        a[j]->y() = yNew + bY[n][i] + (NumSteps - i - 1)*a[j]->
                            yVel()/NumSteps*Window::FrameTime();
                        a[j]->z() = zNew + bZ[n][i] + (NumSteps - i - 1)*a[j]->
                            zVel()/NumSteps*Window::FrameTime();
                        a[j]->onCollide(*b[n]);
                        b[n]->onCollide(*a[j]);
tempDone[j] = true;
                    }
                }
            }
        }
const double colTime1 = timer.getAndRestart();
static double avgColTime0Sq = 0.;
avgColTime0Sq = 0.95*avgColTime0Sq + 0.05*colTime0*colTime0;
static double avgColTime1Sq = 0.;
avgColTime1Sq = 0.95*avgColTime1Sq + 0.05*colTime1*colTime1;
_msgStream << "Avg. collision times: " << std::fixed << std::setprecision(6) << std::setw(8) << std::sqrt(avgColTime0Sq) << " " << std::setw(8) << std::sqrt(avgColTime1Sq) << std::endl;
static double avgFrameTimeSq = 1./(60.*60.);
avgFrameTimeSq = 0.95*avgFrameTimeSq + 0.05*Window::FrameTime()*Window::FrameTime();
_msgStream << "Avg. FPS: " << std::fixed << std::setprecision(2) << std::setw(6) << 1./std::sqrt(avgFrameTimeSq) << std::endl;
static std::deque<double> rHist(60*30, 0.);
static std::deque<double> latHist(60*30, 0.);
{
const double r = std::sqrt(_cameraObject->x()*_cameraObject->x() + _cameraObject->y()*_cameraObject->y() + _cameraObject->z()*_cameraObject->z());
const double lat = std::asin(_cameraObject->z()/r);
rHist.pop_front();
rHist.push_back(r);
latHist.pop_front();
latHist.push_back(lat);
}

        const auto tempObjects = objects(); //TEMP - this prevents missed or multiple updates when `objects()` changes, but is not ideal
        for(const auto& n : tempObjects)
        {
            if(objects().count(n.first))
            {
                reinterpret_cast<Object*>(n.first)->update();
            }
        }

        const double cameraWAtt = std::sqrt(1. - _cameraObject->xAtt()*
            _cameraObject->xAtt() - _cameraObject->yAtt()*_cameraObject->yAtt()
            - _cameraObject->zAtt()*_cameraObject->zAtt());
        const Vec cameraAtt{{_cameraObject->xAtt(), _cameraObject->yAtt(),
            _cameraObject->zAtt(), cameraWAtt}};

        const auto projection = perspective(1., Window::AspectRatio(), 0.1,
            1e3);
        Mat view = Mat::Zero(4);
        view(3, 3) = 1.;
        view.setBlock(0, 0, 3, 3, Mat{{1., 0., 0.}, {0., 0., 1.}, {0., -1.,
            0.}}*to_mat(cameraAtt));

        // Get geometry map.
timer.start();
        _geometryPass.begin(std::vector<LoadAction>(4, LoadAction::Clear),
            LoadAction::Clear);
        _geometryPass.cull(CullMode::Back);
        _geometryPass.depth(DepthTestMode::Less);
        _geometryPass.uniform("projection", projection);
        _geometryPass.uniform("view", convert_mat(view));
        std::unordered_map<void*, std::vector<const Object*>> objectsByModel;
        for(const auto& n : objects())
        {
            const Object* o = reinterpret_cast<const Object*>(n.first);
            if(!o->model()._i.empty())
            {
                // Address held to by `paz::Model::_t` identifies all copies of
                // the same model.
                objectsByModel[o->model()._t.get()].push_back(o);
            }
        }
        for(const auto& m : objectsByModel)
        {
            std::array<std::vector<float>, 4> modelMats;
            for(auto& n : modelMats)
            {
                n.resize(4*m.second.size());
            }
            for(std::size_t i = 0; i < m.second.size(); ++i)
            {
                const double q0 = m.second[i]->xAtt();
                const double q1 = m.second[i]->yAtt();
                const double q2 = m.second[i]->zAtt();
                const auto q3 = -std::sqrt(1. - q0*q0 - q1*q1 - q2*q2);
                Mat model = Mat::Identity(4);
                model(0, 3) = m.second[i]->x() - _cameraObject->x();
                model(1, 3) = m.second[i]->y() - _cameraObject->y();
                model(2, 3) = m.second[i]->z() - _cameraObject->z();
                model.setBlock(0, 0, 3, 3, to_mat(Vec{{q0, q1, q2, q3}}));
                for(int j = 0; j < 4; ++j)
                {
                    std::copy(model.begin() + 4*j, model.begin() + 4*j + 4,
                        modelMats[j].begin() + 4*i);
                }
            }
            const VertexBuffer verts = m.second.back()->model()._v;
            const IndexBuffer inds = m.second.back()->model()._i;
            paz::VertexBuffer buf;
            for(const auto& n : modelMats)
            {
                buf.attribute(4, n);
            }
            _geometryPass.draw(PrimitiveType::Triangles, verts, buf, inds);
        }
        _geometryPass.end();
const double gTime = timer.getAndRestart();
static double avgGTimeSq = 0.;
avgGTimeSq = 0.95*avgGTimeSq + 0.05*gTime*gTime;
_msgStream << "Avg. G-pass time: " << std::fixed << std::setprecision(6) << std::setw(8) << std::sqrt(avgGTimeSq) << std::endl;

        // Render in HDR.
        _renderPass.begin();
        _renderPass.read("materialMap", _materialMap);
        _renderPass.read("normalMap", _normalMap);
        _renderPass.read("coordMap", _coordMap);
        _renderPass.read("depthMap", _depthMap);
        _renderPass.uniform("invProjection", convert_mat(convert_mat(
            projection).inv()));
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

        const std::string msg = _msgStream.str();
        if(!msg.empty())
        {
            std::stringstream{}.swap(_msgStream);

            _textPass.begin();
            _textPass.read("font", _font);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("baseWidth", _font.width());
            _textPass.uniform("baseHeight", _font.height());
            bool highlight = false;
            _textPass.uniform("scale", std::round(2.f*Window::UiScale()));
            int row = 0;
            int col = 0;
            for(auto n : msg)
            {
                if(n == '`')
                {
                    highlight = !highlight;
                    continue;
                }
                if(n == '\n')
                {
                    ++row;
                    col = 0;
                    continue;
                }
                if(n >= '!' && n <= '~')
                {
                    _textPass.uniform("highlight", static_cast<float>(
                        highlight));
                    _textPass.uniform("row", row);
                    _textPass.uniform("col", col);
                    _textPass.uniform("character", static_cast<int>(n - '!'));
                    _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
                }
                ++col;
            }
            _textPass.end();
        }

        /*_plotPass.begin();
        {
            VertexBuffer plotData;
            {
                std::vector<float> v(2*60*30);
                for(std::size_t i = 0; i < 60*30; ++i)
                {
                    v[2*i] = i*2./(60*30 - 1) - 1.;
                    v[2*i + 1] = 0.3*std::log(rHist[i] - 50.) - 1.;
                }
                plotData.attribute(2, v);
            }
            _plotPass.draw(PrimitiveType::LineStrip, plotData);
        }
        {
            VertexBuffer plotData;
            {
                std::vector<float> v(2*60*30);
                for(std::size_t i = 0; i < 60*30; ++i)
                {
                    v[2*i] = i*2./(60*30 - 1) - 1.;
                    v[2*i + 1] = 2.*latHist[i]/M_PI;
                }
                plotData.attribute(2, v);
            }
            _plotPass.draw(PrimitiveType::LineStrip, plotData);
        }
        _plotPass.end();*/

        Window::EndFrame();
    }
}

void paz::App::AttachCamera(const Object& o)
{
    _cameraObject = &o;
}

std::stringstream& paz::App::MsgStream()
{
    return _msgStream;
}

paz::Bytes paz::App::GetAsset(const std::string& path)
{
    return get_asset(path);
}
