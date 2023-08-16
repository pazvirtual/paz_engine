#include "object.hpp"
#include "ui.hpp"
#include "shared.hpp"
#include "io.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>
#include <iomanip>
#include <deque>

#define NO_FRICTION

static constexpr std::size_t NumSteps = 100;
static constexpr double InteractRangeBehindSq = 4.;
static constexpr double InteractRangeInFrontSq = 9.;
static constexpr std::size_t MaxConsoleLines = 1000;
static constexpr float FontScale = 1.5f;
static constexpr int CharWidth = 5;

static paz::Framebuffer _geometryBuffer;
static paz::Framebuffer _oitAccumBuffer;
static paz::Framebuffer _renderBuffer;
static paz::Framebuffer _postBuffer;
static paz::Framebuffer _lumBuffer;

static paz::RenderTarget _diffuseMap(paz::TextureFormat::RGBA16Float);
// ...
static paz::RenderTarget _emissMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _normalMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _depthMap(paz::TextureFormat::Depth32Float);
static paz::RenderTarget _hdrRender(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _finalRender(paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _finalLumMap(paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _oitAccumTex(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _oitRevealTex(paz::TextureFormat::RGBA8UNorm);

static std::unordered_map<void*, paz::InstanceBuffer> _instances;

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass0;
static paz::RenderPass _renderPass1;
static paz::RenderPass _oitAccumPass;
static paz::RenderPass _oitCompositePass;
static paz::RenderPass _postPass0;
static paz::RenderPass _lumPass;
static paz::RenderPass _fxaaPass;
static paz::RenderPass _postPass1;
static paz::RenderPass _consolePass;
static paz::RenderPass _textPass;
static paz::RenderPass _cursorPass;

static paz::VertexBuffer _quadVertices;
static paz::VertexBuffer _sphereVertices;
static paz::IndexBuffer _sphereIndices;

static paz::Texture _font;
static std::stringstream _msgStream;
static std::deque<std::string> _console;
static paz::ConsoleMode _consoleMode = paz::ConsoleMode::Disable;

static std::string _title;

static paz::Texture _defaultDiffTex;

static paz::Texture _cursor;

static const paz::Object* _cameraObject; //TEMP - vector realloc breaks ptrs
static const paz::Object* _micObject; //TEMP - vector realloc breaks ptrs

static const paz::Object* _soundSrc; //TEMP - vector realloc breaks ptrs

static bool _paused;

static double _gravity;

static bool _fxaaEnabled = true;

static paz::Vec _sunDir = paz::Vec::Zero(4);
static std::array<float, 4> _sunIll;

static const std::vector<paz::Button> OptionsButtons =
{
    {
        [](){ return paz::Window::IsFullscreen() ? "Fullscreen: ON" :
            "Fullscreen: OFF"; },
        [](paz::Menu&)
        {
            if(paz::Window::IsFullscreen())
            {
                paz::Window::MakeWindowed();
                paz::save_setting("fullscreen", "0");
            }
            else
            {
                paz::Window::MakeFullscreen();
                paz::save_setting("fullscreen", "1");
            }
        }
    },
    {
        []()
        {
            if(paz::Window::HidpiSupported())
            {
                return paz::Window::HidpiEnabled() ? "HiDPI:      ON" :
                    "HiDPI:      OFF";
            }
            else
            {
                return "HiDPI:      N/A";
            }
        },
        [](paz::Menu&)
        {
            if(paz::Window::HidpiSupported())
            {
                if(paz::Window::HidpiEnabled())
                {
                    paz::Window::DisableHidpi();
                    paz::save_setting("hidpi", "0");
                }
                else
                {
                    paz::Window::EnableHidpi();
                    paz::save_setting("hidpi", "1");
                }
            }
        }
    },
    {
        [](){ return _fxaaEnabled ? "FXAA:       ON" : "FXAA:       OFF"; },
        [](paz::Menu&)
        {
            _fxaaEnabled = !_fxaaEnabled;
            paz::save_setting("fxaa", _fxaaEnabled ? "1" : "0");
        }
    },
    {"Back", [](paz::Menu& m){ m.setState(0, 1); }}
};


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

static constexpr std::array<float, 8> QuadPos = {1, -1, 1, 1, -1, -1, -1, 1};

void paz::App::Init(const std::string& title)
{
    if(load_setting("fullscreen") == "1")
    {
        Window::MakeFullscreen();
    }
    if(load_setting("hidpi") == "0")
    {
        Window::DisableHidpi();
    }
    if(load_setting("fxaa") == "0")
    {
        _fxaaEnabled = false;
    }

    _title = title;

    _geometryBuffer.attach(_diffuseMap);
    // ...
    _geometryBuffer.attach(_emissMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

    _oitAccumBuffer.attach(_oitAccumTex);
    _oitAccumBuffer.attach(_oitRevealTex);
    _oitAccumBuffer.attach(_depthMap);

    _postBuffer.attach(_finalRender);

    _lumBuffer.attach(_finalLumMap);

    const VertexFunction geometryVert(get_builtin("geometry.vert").str());
    const VertexFunction quadVert(get_builtin("quad.vert").str());
    const VertexFunction sceneVert0(get_builtin("scene0.vert").str());
    const VertexFunction sceneVert1(get_builtin("scene1.vert").str());
    const VertexFunction textVert(get_builtin("text.vert").str());
    const VertexFunction cursorVert(get_builtin("cursor.vert").str());
    const VertexFunction oitVert(get_builtin("oit.vert").str());
    const FragmentFunction geometryFrag(get_builtin("geometry.frag").str());
    const FragmentFunction sceneFrag0(get_asset("scene0.frag").str());
    const FragmentFunction sceneFrag1(get_asset("scene1.frag").str());
    const FragmentFunction lumFrag(get_builtin("lum.frag").str());
    const FragmentFunction fxaaFrag(get_builtin("fxaa.frag").str());
    const FragmentFunction postFrag(get_builtin("post.frag").str());
    const FragmentFunction consoleFrag(get_builtin("console.frag").str());
    const FragmentFunction textFrag(get_builtin("text.frag").str());
    const FragmentFunction cursorFrag(get_builtin("cursor.frag").str());
    const FragmentFunction oitFrag(get_asset("oit.frag").str());
    const FragmentFunction compositeFrag(get_builtin("composite.frag").str());

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass0 = RenderPass(_renderBuffer, sceneVert0, sceneFrag0);
    _renderPass1 = RenderPass(_renderBuffer, sceneVert1, sceneFrag1,
        {BlendMode::One_One});
    _oitAccumPass = RenderPass(_oitAccumBuffer, oitVert, oitFrag, {BlendMode::
        One_One, BlendMode::Zero_InvSrcAlpha});
    _oitCompositePass = RenderPass(_renderBuffer, quadVert, compositeFrag,
        {BlendMode::InvSrcAlpha_SrcAlpha});
    _postPass0 = RenderPass(_postBuffer, quadVert, postFrag);
    _lumPass = RenderPass(_lumBuffer, quadVert, lumFrag);
    _fxaaPass = RenderPass(quadVert, fxaaFrag);
    _postPass1 = RenderPass(quadVert, postFrag);
    _consolePass = RenderPass(quadVert, consoleFrag, {BlendMode::
        SrcAlpha_InvSrcAlpha});
    _textPass = RenderPass(textVert, textFrag, {BlendMode::
        SrcAlpha_InvSrcAlpha});
    _cursorPass = RenderPass(cursorVert, cursorFrag, {BlendMode::
        SrcAlpha_InvSrcAlpha});

    _quadVertices.addAttribute(2, QuadPos);

    {
        std::vector<std::string> names;
        std::vector<std::vector<float>> positions;
        std::vector<std::vector<float>> uvs;
        std::vector<std::vector<float>> normals;
        std::vector<std::vector<unsigned int>> materials;
        std::vector<std::string> materialNames;
        std::vector<std::string> materialLibs;
        std::vector<std::vector<unsigned int>> indices;
        parse_obj(get_builtin("icosphere3.obj"), names, positions, uvs, normals,
            materials, materialNames, materialLibs, indices);
        _sphereVertices.addAttribute(4, positions[0]);
        _sphereIndices = IndexBuffer(indices[0]);
    }

    _font = Texture(get_asset_image("font.pbm")); //TEMP - note that only red channel is used

    std::vector<unsigned char> temp(4*512*512);
    for(int i = 0; i < 512; ++i)
    {
        for(int j = 0; j < 512; ++j)
        {
            if((i*8/512)%2 == (j*8/512)%2)
            {
                temp[4*(512*i + j) + 0] = 0;
                temp[4*(512*i + j) + 1] = 255;
                temp[4*(512*i + j) + 2] = 255;
            }
            else
            {
                temp[4*(512*i + j) + 0] = 255;
                temp[4*(512*i + j) + 1] = 0;
                temp[4*(512*i + j) + 2] = 0;
            }
            temp[4*(512*i + j) + 3] = 255;
        }
    }
    _defaultDiffTex = Texture(TextureFormat::RGBA8UNorm_sRGB, 512, 512, temp.
        data(), MinMagFilter::Linear, MinMagFilter::Linear, MipmapFilter::
        Linear, WrapMode::Repeat, WrapMode::Repeat);

    _cursor = Texture(get_asset_image("cursor.pbm")); //TEMP - note that only red channel is used
}

void paz::App::Run()
{
    Font menuFont(_font, 2*FontScale, CharWidth);

    {
        bool done = false;
        Menu startMenu(menuFont, _title,
        {
            {
                {"Start", [&](Menu&){ done = true; }},
                {"Options", [](Menu& m){ m.setState(1, 0); }},
                {"Quit", [](Menu&){ Window::Quit(); }}
            },
            OptionsButtons
        }
        );
        while(!Window::Done() && !done)
        {
            if(startMenu.curPage() == 1 && (Window::KeyPressed(Key::Escape) ||
                Window::GamepadPressed(GamepadButton::Start)))
            {
                startMenu.setState(0, 1);
            }

            startMenu.update();

            _textPass.begin({LoadAction::Clear});
            _textPass.read("font", startMenu.font().tex());
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", startMenu.font().charWidth());
            _textPass.uniform("baseWidth", startMenu.font().tex().width());
            _textPass.uniform("baseHeight", startMenu.font().tex().height());
            _textPass.uniform("scale", startMenu.font().curScale());
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices,
                startMenu.chars());
            _textPass.end();

            if(Window::MouseActive())
            {
                _cursorPass.begin({LoadAction::Load});
                _cursorPass.uniform("x", static_cast<float>(Window::MousePos().
                    first)/(Window::Width() - 1));
                _cursorPass.uniform("y", static_cast<float>(Window::MousePos().
                    second)/(Window::Height() - 1));
                _cursorPass.uniform("h", startMenu.curButton() < 0 ? 0.f : 1.f);
                _cursorPass.uniform("scale", static_cast<int>(std::round(20.*
                    Window::UiScale()))); //TEMP
                _cursorPass.uniform("width", Window::ViewportWidth());
                _cursorPass.uniform("height", Window::ViewportHeight());
                _cursorPass.read("tex", _cursor);
                _cursorPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
                _cursorPass.end();
            }

            Window::EndFrame();
        }
        if(Window::Done())
        {
            return;
        }
    }

    Window::SetCursorMode(CursorMode::Disable);

    Menu pauseMenu(menuFont, "Paused",
    {
        {
            {"Resume", [&](Menu&)
                {
                    _paused = false;
                    Window::SetCursorMode(CursorMode::Disable);
                }
            },
            {"Options", [](Menu& m){ m.setState(1, 0); }},
            {"Quit", [](Menu&){ Window::Quit(); }}
        },
        OptionsButtons
    }
    );

    while(!Window::Done())
    {
        bool justPaused = false;
        if(!_paused && (Window::KeyPressed(Key::Escape) || Window::
            GamepadPressed(GamepadButton::Start)))
        {
            justPaused = true;
            _paused = true;
            pauseMenu.setState(0, 0);
        }

        if(_micObject && !_paused && objects().count(reinterpret_cast<std::
            uintptr_t>(_soundSrc)))
        {
            const Vec relPos{{_soundSrc->x() - _micObject->x(), _soundSrc->y() -
                _micObject->y(), _soundSrc->z() - _micObject->z()}};
            const Vec relVel{{_soundSrc->xVel() - _micObject->xVel(),
                _soundSrc->yVel() - _micObject->yVel(), _soundSrc->zVel() -
                _micObject->zVel()}};
            const double dist = relPos.norm();
            const Vec dir = relPos/dist;
            const Vec micAtt{{_micObject->xAtt(), _micObject->yAtt(),
                _micObject->zAtt(), std::sqrt(1. - _micObject->xAtt()*
                _micObject->xAtt() - _micObject->yAtt()*_micObject->yAtt() -
                _micObject->zAtt()*_micObject->zAtt())}};
            const Mat micRot = to_mat(micAtt);
            const Vec micX = micRot.row(1).trans();
            const Vec micY = -micRot.row(0).trans();
            const Vec lEar = std::sin(1.)*micX + std::cos(1.)*micY;
            const Vec rEar = std::sin(1.)*micX - std::cos(1.)*micY;
            const double vlos = dir.dot(relVel);
            const double txPwr = 400.;
            static constexpr double vs = 343.;
            static constexpr double maxRxPwr = 0.3;
            const double t0 = (vs + vlos)/(800.*Pi*dist); // assuming f = 200
            const double rxPwr = std::min(maxRxPwr, txPwr*t0*t0);
            double lPwr = rxPwr*(0.6 + 0.4*(dir.dot(lEar))); //TEMP
            double rPwr = rxPwr*(0.6 + 0.4*(dir.dot(rEar))); //TEMP
            lPwr /= SqrtTwo*(0.6 + 0.4*std::cos(0.5*Pi - 1.));
            rPwr /= SqrtTwo*(0.6 + 0.4*std::cos(0.5*Pi - 1.));
            const double lVol = std::sqrt(std::sqrt(std::max(0., lPwr)));
            const double rVol = std::sqrt(std::sqrt(std::max(0., rPwr)));
            const double lFreqScale = 1./(1. + vlos/vs);
            const double rFreqScale = 1./(1. + vlos/vs);
            AudioEngine::SetVolume(lVol, 0);
            AudioEngine::SetVolume(rVol, 1);
            AudioEngine::SetFreqScale(lFreqScale, 0);
            AudioEngine::SetFreqScale(rFreqScale, 1);
        }
        else
        {
            AudioEngine::SetVolume(0.);
        }

if(!_paused)
{
        physics(_gravity);

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

        std::vector<double> times(NumSteps);
        for(std::size_t i = 0; i < NumSteps; ++i)
        {
            times[i] = static_cast<double>(i)/(NumSteps - 1);
        }

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
                bX[i][j] = b[i]->xPrev() + times[j]*(b[i]->x() - b[i]->xPrev());
                bY[i][j] = b[i]->yPrev() + times[j]*(b[i]->y() - b[i]->yPrev());
                bZ[i][j] = b[i]->zPrev() + times[j]*(b[i]->z() - b[i]->zPrev());
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
                    double x = a[j]->xPrev() + times[i]*(a[j]->x() - a[j]->
                        xPrev()) - bX[n][i];
                    double y = a[j]->yPrev() + times[i]*(a[j]->y() - a[j]->
                        yPrev()) - bY[n][i];
                    double z = a[j]->zPrev() + times[i]*(a[j]->z() - a[j]->
                        zPrev()) - bZ[n][i];
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
                        a[j]->x() = xNew + bX[n][i];
                        a[j]->y() = yNew + bY[n][i];
                        a[j]->z() = zNew + bZ[n][i];
                        // Rewind world object to time of collison.
                        b[n]->x() = bX[n][i];
                        b[n]->y() = bY[n][i];
                        b[n]->z() = bZ[n][i];
                        a[j]->onCollide(*b[n]);
                        b[n]->onCollide(*a[j]);
                        // Fast forward both objects.
                        b[n]->x() = bX[n].back();
                        b[n]->y() = bY[n].back();
                        b[n]->z() = bZ[n].back();
                        const double extraTime = PhysTime()*(times.back() -
                            times[i]);
                        a[j]->x() += extraTime*a[j]->xVel();
                        a[j]->y() += extraTime*a[j]->yVel();
                        a[j]->z() += extraTime*a[j]->zVel();
tempDone[j] = true;
                    }
                }
            }
        }

        const auto tempObjects = objects(); //TEMP - this prevents missed or multiple updates when `objects()` changes, but is not ideal
        for(const auto& n : tempObjects)
        {
            if(objects().count(n.first))
            {
                reinterpret_cast<Object*>(n.first)->update();
            }
        }
}

        const double cameraWAtt = std::sqrt(1. - _cameraObject->xAtt()*
            _cameraObject->xAtt() - _cameraObject->yAtt()*_cameraObject->yAtt()
            - _cameraObject->zAtt()*_cameraObject->zAtt());
        const Vec cameraAtt{{_cameraObject->xAtt(), _cameraObject->yAtt(),
            _cameraObject->zAtt(), cameraWAtt}};

        const Vec cameraPos{{_cameraObject->x(), _cameraObject->y(),
            _cameraObject->z()}};

        const auto projection = perspective(1., Window::AspectRatio(), 0.1,
            1e3);
        Mat view = Mat::Zero(4);
        view(3, 3) = 1.;
        view.setBlock(0, 0, 3, 3, Mat{{1., 0., 0.}, {0., 0., 1.}, {0., -1.,
            0.}}*to_mat(cameraAtt));

        // Prepare for rendering.
        std::unordered_map<void*, std::vector<const Object*>> objectsByModel;
        std::vector<const Object*> invisibleObjects;
        for(const auto& n : objects())
        {
            const Object* o = reinterpret_cast<const Object*>(n.first);
            if(o->model()._i.empty() && o->model()._transp.empty())
            {
                invisibleObjects.push_back(o);
            }
            else
            {
                // Address held to by `paz::Model::_t` identifies all copies of
                // the same model.
                objectsByModel[o->model()._t.get()].push_back(o);
            }
        }
        for(const auto& n : objectsByModel)
        {
            if(!_instances.count(n.first) || _instances.at(n.first).size() !=
                n.second.size()) //TEMP - should only need to change _numInstances if larger, not reinit whole buffer
            {
                _instances[n.first] = InstanceBuffer(n.second.size());
                _instances.at(n.first).addAttribute(4, DataType::Float);
                _instances.at(n.first).addAttribute(2, DataType::Float);
            }
        }
        std::vector<float> lightsData0;
        std::vector<float> lightsData1;
        for(const auto& n : objectsByModel)
        {
            std::array<std::vector<float>, 2> modelMatData;
            modelMatData[0].resize(4*n.second.size());
            modelMatData[1].resize(2*n.second.size());
            for(std::size_t i = 0; i < n.second.size(); ++i)
            {
                modelMatData[0][4*i + 0] = n.second[i]->xAtt();
                modelMatData[0][4*i + 1] = n.second[i]->yAtt();
                modelMatData[0][4*i + 2] = n.second[i]->zAtt();
                modelMatData[0][4*i + 3] = n.second[i]->x() - cameraPos(0);
                modelMatData[1][2*i + 0] = n.second[i]->y() - cameraPos(1);
                modelMatData[1][2*i + 1] = n.second[i]->z() - cameraPos(2);
                if(!n.second[i]->lights().empty())
                {
                    const Vec model0{{modelMatData[0][4*i + 0], modelMatData[0][4*i + 1], modelMatData[0][4*i + 2], modelMatData[0][4*i + 3]}};
                    const Vec model1{{modelMatData[1][2*i + 0], modelMatData[1][2*i + 1]}};
                    const Vec att{{model0(0), model0(1), model0(2), std::sqrt(1. - model0.head(3).normSq())}};
                    const double xx = att(0)*att(0);
                    const double yy = att(1)*att(1);
                    const double zz = att(2)*att(2);
                    const double xy = att(0)*att(1);
                    const double zw = att(2)*att(3);
                    const double xz = att(0)*att(2);
                    const double yw = att(1)*att(3);
                    const double yz = att(1)*att(2);
                    const double xw = att(0)*att(3);
                    const Mat mv = view*Mat{{1. - 2.*(yy + zz), 2.*(xy - zw), 2.*(xz + yw), model0(3)},
                                            {2.*(xy + zw), 1. - 2.*(xx + zz), 2.*(yz - xw), model1(0)},
                                            {2.*(xz - yw), 2.*(yz + xw), 1. - 2.*(xx + yy), model1(1)},
                                            {0., 0., 0., 1.}};
                    for(const auto& m : n.second[i]->lights())
                    {
                        const Vec p = mv*Vec{{m[0], m[1], m[2], 1.}};
                        lightsData0.push_back(p(0));
                        lightsData0.push_back(p(1));
                        lightsData0.push_back(p(2));
                        lightsData0.push_back(m[3]);
                        lightsData1.push_back(m[4]);
                        lightsData1.push_back(m[5]);
                        lightsData1.push_back(m[6]);
                        lightsData1.push_back(0.f);
                    }
                }
            }
            _instances.at(n.first).subAttribute(0, modelMatData[0]);
            _instances.at(n.first).subAttribute(1, modelMatData[1]);
        }
        for(const auto& n : invisibleObjects)
        {
            if(!n->lights().empty())
            {
                const double xAtt = n->xAtt();
                const double yAtt = n->yAtt();
                const double zAtt = n->zAtt();
                const double wAtt = std::sqrt(1. - xAtt*xAtt - yAtt*yAtt - zAtt*zAtt);
                const double xx = xAtt*xAtt;
                const double yy = yAtt*yAtt;
                const double zz = zAtt*zAtt;
                const double xy = xAtt*yAtt;
                const double zw = zAtt*wAtt;
                const double xz = xAtt*zAtt;
                const double yw = yAtt*wAtt;
                const double yz = yAtt*zAtt;
                const double xw = xAtt*wAtt;
                const Mat mv = view*Mat{{1. - 2.*(yy + zz), 2.*(xy - zw), 2.*(xz + yw), n->x() - cameraPos(0)},
                                        {2.*(xy + zw), 1. - 2.*(xx + zz), 2.*(yz - xw), n->y() - cameraPos(1)},
                                        {2.*(xz - yw), 2.*(yz + xw), 1. - 2.*(xx + yy), n->z() - cameraPos(2)},
                                        {0., 0., 0., 1.}};
                for(const auto& m : n->lights())
                {
                    const Vec p = mv*Vec{{m[0], m[1], m[2], 1.}};
                    lightsData0.push_back(p(0));
                    lightsData0.push_back(p(1));
                    lightsData0.push_back(p(2));
                    lightsData0.push_back(m[3]);
                    lightsData1.push_back(m[4]);
                    lightsData1.push_back(m[5]);
                    lightsData1.push_back(m[6]);
                    lightsData1.push_back(0.f);
                }
            }
        }
        InstanceBuffer lights; //TEMP - should just sub data if same number
        lights.addAttribute(4, lightsData0);
        lights.addAttribute(4, lightsData1);

        // Get geometry map.
        _geometryPass.begin(std::vector<LoadAction>(4, LoadAction::Clear),
            LoadAction::Clear);
        _geometryPass.cull(CullMode::Back);
        _geometryPass.depth(DepthTestMode::Less);
        _geometryPass.uniform("projection", projection);
        _geometryPass.uniform("view", convert_mat(view));
        for(const auto& n : objectsByModel)
        {
            if(n.second.back()->model()._i.empty())
            {
                continue;
            }
            if(n.second.back()->model()._diffTex.width())
            {
                _geometryPass.read("diffTex", n.second.back()->model().
                    _diffTex);
            }
            else
            {
                _geometryPass.read("diffTex", _defaultDiffTex);
            }
            std::array<float, 4> emiss;
            std::copy(n.second.back()->model()._emiss.begin(), n.second.back()->
                model()._emiss.end(), emiss.begin());
            emiss[3] = 1.f;
            _geometryPass.uniform("emiss", emiss);
            _geometryPass.draw(PrimitiveType::Triangles, n.second.back()->
                model()._v, _instances.at(n.first), n.second.back()->model().
                _i);
        }
        _geometryPass.end();

        // Render in HDR.
        _renderPass0.begin({LoadAction::Clear});//, LoadAction::Load);
        _renderPass0.read("diffuseMap", _diffuseMap);
        // ...
        _renderPass0.read("emissMap", _emissMap);
        _renderPass0.read("normalMap", _normalMap);
        _renderPass0.read("depthMap", _depthMap);
        _renderPass0.uniform("invProjection", convert_mat(convert_mat(
            projection).inv()));
        _renderPass0.uniform("lightDir", convert_vec(view*_sunDir));
        _renderPass0.uniform("ill", _sunIll);
        _renderPass0.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _renderPass0.end();

        _renderPass1.begin({LoadAction::Load});//, LoadAction::Load);
        _renderPass1.cull(CullMode::Front);
//        _renderPass1.depth(DepthTestMode::GreaterNoMask);
        _renderPass1.read("diffuseMap", _diffuseMap);
        // ...
        _renderPass1.read("emissMap", _emissMap);
        _renderPass1.read("normalMap", _normalMap);
        _renderPass1.read("depthMap", _depthMap);
        _renderPass1.uniform("projection", projection);
        _renderPass1.uniform("invProjection", convert_mat(convert_mat(
            projection).inv()));
        _renderPass1.draw(PrimitiveType::Triangles, _sphereVertices, lights,
            _sphereIndices);
        _renderPass1.end();

        // Render transparent objects. //TEMP - skip if possible
        unsigned int numLights = lightsData0.size()/4;
        if(numLights > 16)
        {
            throw std::runtime_error("Too many lights (" + std::to_string(numLights) + " > 16).");
        }
        _oitAccumPass.begin({LoadAction::Clear, LoadAction::Clear}, LoadAction::
            Load);
        _oitAccumPass.depth(DepthTestMode::LessNoMask);
        _oitAccumPass.uniform("numLights", numLights);
        _oitAccumPass.uniform("light0", lightsData0);
        _oitAccumPass.uniform("light1", lightsData1);
        _oitAccumPass.uniform("projection", projection);
        _oitAccumPass.uniform("invProjection", convert_mat(convert_mat(
            projection).inv())); //TEMP - precompute - used elsewhere
        _oitAccumPass.uniform("sunDir", convert_vec(view*_sunDir)); //TEMP - precompute
        _oitAccumPass.uniform("sunIll", _sunIll);
        _oitAccumPass.uniform("view", convert_mat(view));
        for(const auto& n : objectsByModel)
        {
            if(!n.second.back()->model()._transp.empty())
            {
                _oitAccumPass.draw(PrimitiveType::Triangles, n.second.back()->
                    model()._transp, _instances.at(n.first));
            }
        }
        _oitAccumPass.end();

        _oitCompositePass.begin({LoadAction::Load});
        _oitCompositePass.read("accumTex", _oitAccumTex);
        _oitCompositePass.read("revealTex", _oitRevealTex);
        _oitCompositePass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _oitCompositePass.end();

        if(_fxaaEnabled)
        {
            // Tonemap to linear LDR.
            _postPass0.begin();
            _postPass0.read("hdrRender", _hdrRender);
            _postPass0.uniform("whitePoint", 1.f);
            _postPass0.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _postPass0.end();

            // Get luminance map.
            _lumPass.begin();
            _lumPass.read("img", _finalRender);
            _lumPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _lumPass.end();

            // Anti-alias.
            _fxaaPass.begin();
            _fxaaPass.read("img", _finalRender);
            _fxaaPass.read("lum", _finalLumMap);
            _fxaaPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _fxaaPass.end();
        }
        else
        {
            // Tonemap to linear LDR.
            _postPass1.begin();
            _postPass1.read("hdrRender", _hdrRender);
            _postPass1.uniform("whitePoint", 1.f);
            _postPass1.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _postPass1.end();
        }

        std::string line;
        std::size_t numNewRows = 0;
        while(std::getline(_msgStream, line))
        {
            _console.push_back(line);
            ++numNewRows;
            if(_console.size() > MaxConsoleLines)
            {
                _console.pop_front();
            }
        }
        std::stringstream{}.swap(_msgStream);
        InstanceBuffer consoleChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(_consoleMode != ConsoleMode::Disable && !_console.empty())
        {
            const float scale = std::round(FontScale*Window::UiScale());
            int maxVisRows = std::ceil(Window::ViewportHeight()/(scale*_font.
                height()));
            maxVisRows = std::min(static_cast<std::size_t>(maxVisRows),
                _console.size());
            if(_consoleMode == ConsoleMode::CurrentFrame)
            {
                maxVisRows = std::min(static_cast<std::size_t>(maxVisRows),
                    numNewRows);
            }

            const float consoleHeight = std::min(1.f, (maxVisRows + 1.f/_font.
                height())*scale*_font.height()/Window::ViewportHeight());
            std::size_t maxCols = 0;
            for(int i = 0; i < maxVisRows; ++i)
            {
                maxCols = std::max(maxCols, _console.rbegin()[i].size());
            }
            const float consoleWidth = std::min(1.f, (1 + maxCols*(CharWidth +
                1))*scale/Window::ViewportWidth());

            _consolePass.begin({LoadAction::Load});
            _consolePass.uniform("width", consoleWidth);
            _consolePass.uniform("height", consoleHeight);
            _consolePass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _consolePass.end();

            std::vector<float> highlightAttr;
            std::vector<int> characterAttr;
            std::vector<int> colAttr;
            std::vector<int> rowAttr;
            bool highlight = false;
            for(int row = 0; row < maxVisRows; ++row)
            {
                int col = 0;
                for(auto n : _console.rbegin()[row])
                {
                    if(n == '`')
                    {
                        highlight = !highlight;
                    }
                    else if(n == ' ')
                    {
                        ++col;
                    }
                    else if(n == '\t')
                    {
                        col = (col/4 + 1)*4;
                    }
                    else if(n >= '!' && n <= '~')
                    {
                        highlightAttr.push_back(highlight);
                        characterAttr.push_back(n - '!');
                        colAttr.push_back(col);
                        rowAttr.push_back(row);
                        ++col;
                    }
                }
            }

            consoleChars.addAttribute(1, highlightAttr);
            consoleChars.addAttribute(1, characterAttr);
            consoleChars.addAttribute(1, colAttr);
            consoleChars.addAttribute(1, rowAttr);

            _textPass.begin({LoadAction::Load});
            _textPass.read("font", _font);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", CharWidth);
            _textPass.uniform("baseWidth", _font.width());
            _textPass.uniform("baseHeight", _font.height());
            _textPass.uniform("scale", static_cast<int>(scale));
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices,
                consoleChars);
            _textPass.end();
        }

        InstanceBuffer menuChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(_paused)
        {
            _consolePass.begin({LoadAction::Load});
            _consolePass.uniform("width", 1.f);
            _consolePass.uniform("height", 1.f);
            _consolePass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _consolePass.end();

            if(!justPaused && (Window::KeyPressed(Key::Escape) || Window::
                GamepadPressed(GamepadButton::Start)))
            {
                if(pauseMenu.curPage() == 0)
                {
                    _paused = false;
                    Window::SetCursorMode(CursorMode::Disable);
                }
                else if(pauseMenu.curPage() == 1)
                {
                    pauseMenu.setState(0, 1);
                }
            }

            pauseMenu.update();

            _textPass.begin({LoadAction::Load});
            _textPass.read("font", pauseMenu.font().tex());
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", pauseMenu.font().charWidth());
            _textPass.uniform("baseWidth", pauseMenu.font().tex().width());
            _textPass.uniform("baseHeight", pauseMenu.font().tex().height());
            _textPass.uniform("scale", pauseMenu.font().curScale());
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices,
                pauseMenu.chars());
            _textPass.end();

            if(Window::MouseActive())
            {
                _cursorPass.begin({LoadAction::Load});
                _cursorPass.uniform("x", static_cast<float>(Window::MousePos().
                    first)/(Window::Width() - 1));
                _cursorPass.uniform("y", static_cast<float>(Window::MousePos().
                    second)/(Window::Height() - 1));
                _cursorPass.uniform("h", pauseMenu.curButton() < 0 ? 0.f : 1.f);
                _cursorPass.uniform("scale", static_cast<int>(std::round(20.*
                    Window::UiScale()))); //TEMP
                _cursorPass.uniform("width", Window::ViewportWidth());
                _cursorPass.uniform("height", Window::ViewportHeight());
                _cursorPass.read("tex", _cursor);
                _cursorPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
                _cursorPass.end();
            }
        }

        Window::EndFrame();
    }
}

void paz::App::AttachCamera(const Object& o)
{
    _cameraObject = &o;
}

void paz::App::AttachMic(const Object& o)
{
    _micObject = &o;
}

std::stringstream& paz::App::MsgStream()
{
    return _msgStream;
}

paz::Bytes paz::App::GetAsset(const std::string& path)
{
    return get_asset(path);
}

double paz::App::PhysTime()
{
    return std::min(0.1, Window::FrameTime());
}

void paz::App::SetConsole(ConsoleMode mode)
{
    _consoleMode = mode;
}

void paz::App::SetGravity(double acc)
{
    _gravity = acc;
}

void paz::App::SetSound(const Object& o, const AudioTrack& sound, bool loop)
{
    _soundSrc = &o;
    AudioEngine::Play(sound, loop);
}

void paz::App::SetSun(const Vec& dir, const Vec& ill)
{
    std::copy(dir.begin(), dir.end(), _sunDir.begin());
    std::copy(ill.begin(), ill.end(), _sunIll.begin());
}
