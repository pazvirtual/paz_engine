#include "object.hpp"
#include "ui.hpp"
#include "shared.hpp"
#include "io.hpp"
#include "threads.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>
#include <deque>

static constexpr double InteractRangeBehindSq = 4.;
static constexpr double InteractRangeInFrontSq = 9.;
static constexpr std::size_t MaxConsoleLines = 1000;
static constexpr float FontScale = 1.5f;
static constexpr int CharWidth = 5;
static constexpr double Timestep = 1./60.;

static paz::Framebuffer _geometryBuffer;
static paz::Framebuffer _oitAccumBuffer;
static paz::Framebuffer _oitDepthBuffer;
static paz::Framebuffer _renderBuffer;
static paz::Framebuffer _dofBuffer;
static paz::Framebuffer _postBuffer;
static paz::Framebuffer _lumBuffer;

static paz::RenderTarget _diffuseMap(paz::TextureFormat::RGBA16Float);
// ...
static paz::RenderTarget _emissMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _normalMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _depthMap(paz::TextureFormat::Depth32Float);
static paz::RenderTarget _hdrRender(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _dofRender(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _finalRender(paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _finalLumMap(paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _oitAccumTex(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _oitRevealTex(paz::TextureFormat::RGBA8UNorm);

static std::unordered_map<void*, paz::InstanceBuffer> _instances;

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass0;
static paz::RenderPass _renderPass1;
static paz::RenderPass _dofPass;
static paz::RenderPass _oitAccumPass;
static paz::RenderPass _oitCompositePass;
static paz::RenderPass _oitDepthPass;
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

static paz::Texture _fontTex;
static paz::Texture _consoleFontTex;
static std::stringstream _msgStream;
static std::deque<std::string> _console;
static paz::ConsoleMode _consoleMode = paz::ConsoleMode::Disable;
static std::queue<std::pair<std::string, double>> _dialogs;
static double _dialogTime;
static bool _dialogSkipFrame;

static std::string _title;

static paz::Texture _defaultDiffTex;

static paz::Texture _cursor;
static paz::Texture _reticule;
static int _reticuleIdx;
static bool _reticuleHighlight;

static paz::ObjectPtr _cameraObject;
static paz::ObjectPtr _micObject;

static paz::ObjectPtr _soundSrc;

static bool _paused;
static double _accumTime;
static int _lookSensitivity = 5;

static double _gravAcc;
static paz::Threadpool _threads;

static bool _fxaaEnabled = true;

static float _dofMinDepth;
static float _dofMaxDepth;

static paz::Vec _sunDir = paz::Vec::Zero(4);
static std::array<float, 4> _sunIll;

static const std::vector<paz::Button> OptionsButtons =
{
    {
        [](){ return paz::Window::IsFullscreen() ? "Fullscreen:  ON" :
            "Fullscreen:  OFF"; },
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
                return paz::Window::HidpiEnabled() ? "HiDPI:       ON" :
                    "HiDPI:       OFF";
            }
            else
            {
                return "HiDPI:       N/A";
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
        },
        []()
        {
            return paz::Window::HidpiSupported();
        }
    },
    {
        [](){ return _fxaaEnabled ? "FXAA:        ON" : "FXAA:        OFF"; },
        [](paz::Menu&)
        {
            _fxaaEnabled = !_fxaaEnabled;
            paz::save_setting("fxaa", _fxaaEnabled ? "1" : "0");
        }
    },
    {
        [](){ return "Sensitivity: " + std::to_string(_lookSensitivity); },
        [](paz::Menu&)
        {
            ++_lookSensitivity;
            if(_lookSensitivity > 10)
            {
                _lookSensitivity = 1;
            }
            paz::save_setting("sensitivity", std::to_string(_lookSensitivity));
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
    {
        const auto str = load_setting("sensitivity");
        if(!str.empty())
        {
            _lookSensitivity = std::stoi(str);
        }
    }

    Window::EnableDithering();

    _title = title;
    Window::SetTitle(_title);

    _geometryBuffer.attach(_diffuseMap);
    // ...
    _geometryBuffer.attach(_emissMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

    _dofBuffer.attach(_dofRender);

    _oitAccumBuffer.attach(_oitAccumTex);
    _oitAccumBuffer.attach(_oitRevealTex);
    _oitAccumBuffer.attach(_depthMap);

    _oitDepthBuffer.attach(_depthMap);

    _postBuffer.attach(_finalRender);

    _lumBuffer.attach(_finalLumMap);

    const VertexFunction geometryVert(get_builtin("geometry.vert").str());
    const VertexFunction quadVert(get_builtin("quad.vert").str());
    const VertexFunction sceneVert0(get_builtin("scene0.vert").str());
    const VertexFunction sceneVert1(get_builtin("scene1.vert").str());
    const VertexFunction textVert(get_builtin("text.vert").str());
    const VertexFunction cursorVert(get_builtin("cursor.vert").str());
    const VertexFunction oitVert(get_builtin("oit.vert").str());
    const VertexFunction oitDepthVert(get_builtin("oitdepth.vert").str());
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
    const FragmentFunction dofFrag(get_builtin("dof.frag").str());
    const FragmentFunction emptyFrag(get_builtin("empty.frag").str());

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass0 = RenderPass(_renderBuffer, sceneVert0, sceneFrag0);
    _renderPass1 = RenderPass(_renderBuffer, sceneVert1, sceneFrag1, {BlendMode::
        One_One});
    _dofPass = RenderPass(_dofBuffer, quadVert, dofFrag);
    _oitAccumPass = RenderPass(_oitAccumBuffer, oitVert, oitFrag, {BlendMode::
        One_One, BlendMode::Zero_InvSrcAlpha});
    _oitCompositePass = RenderPass(_renderBuffer, quadVert, compositeFrag,
        {BlendMode::InvSrcAlpha_SrcAlpha});
    _oitDepthPass = RenderPass(_oitDepthBuffer, oitDepthVert, emptyFrag);
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
        std::vector<float> positions;
        std::vector<float> uvs;
        std::vector<float> normals;
        std::vector<unsigned int> materials;
        std::vector<std::string> materialNames;
        std::vector<std::string> materialLibs;
        std::vector<unsigned int> indices;
        parse_model(get_builtin("icosphere3.pazmodel"), positions, uvs, normals,
            materials, materialNames, materialLibs, indices);
        _sphereVertices.addAttribute(4, positions);
        _sphereIndices = IndexBuffer(indices);
    }

    _consoleFontTex = Texture(get_builtin_image("consolefont.pbm")); //TEMP - note that only red channel is used
    try
    {
        _fontTex = Texture(get_asset_image("font.pbm")); //TEMP - note that only red channel is used
    }
    catch(...)
    {
        _fontTex = _consoleFontTex;
    }

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
    _reticule = Texture(get_asset_image("reticule.pbm")); //TEMP - note that only red channel is used
}

void paz::App::Run()
{
    Font menuFont(_fontTex, 2.*FontScale, CharWidth);

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
            Window::PollEvents();

            if(startMenu.curPage() == 1 && (Window::KeyPressed(Key::Escape) ||
                Window::MousePressed(paz::MouseButton::Back) || Window::
                GamepadPressed(GamepadButton::Start)))
            {
                startMenu.setState(0, 1);
            }

            startMenu.update();

            _textPass.begin({LoadAction::Clear});
            _textPass.read("font", startMenu.font().tex());
            _textPass.uniform("u", 0.f);
            _textPass.uniform("v", 0.f);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", startMenu.font().charWidth());
            _textPass.uniform("baseWidth", startMenu.font().tex().width());
            _textPass.uniform("baseHeight", startMenu.font().tex().height());
            _textPass.uniform("scale", startMenu.font().curScale());
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices, startMenu.
                chars());
            _textPass.end();

            if(Window::MouseActive())
            {
                _cursorPass.begin({LoadAction::Load});
                _cursorPass.uniform("x", static_cast<float>(Window::MousePos().
                    first)/(Window::Width() - 1));
                _cursorPass.uniform("y", static_cast<float>(Window::MousePos().
                    second)/(Window::Height() - 1));
                _cursorPass.uniform("idx", 0);
                _cursorPass.uniform("aspect", static_cast<float>(_cursor.
                    height())/_cursor.width());
                _cursorPass.uniform("h", static_cast<float>(startMenu.
                    curButtonEnabled()));
                _cursorPass.uniform("scale", static_cast<int>(std::max(1.f, std::
                    round(20.f*Window::UiScale())))); //TEMP
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

    InputData input;

    while(!Window::Done())
    {
        bool justPaused = false;
        if(_paused)
        {
            _accumTime = 0.;
            Window::PollEvents();
            input.resetEvents();
        }
        else
        {
            _accumTime += Window::FrameTime();

            Window::PollEvents();
            if(Window::KeyPressed(Key::Escape) || Window::GamepadPressed(
                GamepadButton::Start))
            {
                justPaused = true;
                _paused = true;
                pauseMenu.setState(0, 0);
            }
            input.copyEvents(Timestep, 0.2*_lookSensitivity);

            while(_accumTime > 0.)
            {
                if(_consoleMode == ConsoleMode::LatestStep)
                {
                    std::stringstream{}.swap(_msgStream);
                }

                do_physics(_gravAcc, Timestep);
                do_collisions(_threads, Timestep);
                const auto tempObjects = objects(); //TEMP - this prevents missed or multiple updates when `objects()` changes, but is not ideal
                for(const auto& n : tempObjects)
                {
                    if(objects().count(n.first))
                    {
                        reinterpret_cast<Object*>(n.first)->update(input);
                    }
                }

                input.resetEvents();
                _accumTime -= Timestep;
            }
        }

        const double fac = 1. + _accumTime/Timestep;

        if(!_paused && _micObject && _soundSrc)
        {
            const Vec relPos{{mix(_soundSrc->xPrev() - _micObject->xPrev(),
                _soundSrc->x() - _micObject->x(), fac), mix(_soundSrc->yPrev() -
                _micObject->yPrev(), _soundSrc->y() - _micObject->y(), fac), mix(
                _soundSrc->zPrev() - _micObject->zPrev(), _soundSrc->z() -
                _micObject->z(), fac)}};
            const Vec relVel{{_soundSrc->xVel() - _micObject->xVel(), _soundSrc->
                yVel() - _micObject->yVel(), _soundSrc->zVel() - _micObject->
                zVel()}};
            const double dist = relPos.norm();
            const Vec dir = relPos/dist;
            const double wAttPrev = std::sqrt(1. - _micObject->xAttPrev()*
                _micObject->xAttPrev() - _micObject->yAttPrev()*_micObject->
                yAttPrev() - _micObject->zAttPrev()*_micObject->zAttPrev());
            const double wAtt = std::sqrt(1. - _micObject->xAtt()*_micObject->
                xAtt() - _micObject->yAtt()*_micObject->yAtt() - _micObject->
                zAtt()*_micObject->zAtt());
            const Vec micAtt = nlerp({{_micObject->xAttPrev(), _micObject->
                yAttPrev(), _micObject->zAttPrev(), wAttPrev}}, {{_micObject->
                xAtt(), _micObject->yAtt(), _micObject->zAtt(), wAtt}}, fac);
            const Mat micRot = to_mat(micAtt);
            const Vec micX = micRot.row(0).trans();
            const Vec micY = micRot.row(1).trans();
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

        const double cameraWAttPrev = std::sqrt(1. - _cameraObject->xAttPrev()*
            _cameraObject->xAttPrev() - _cameraObject->yAttPrev()*_cameraObject->
            yAttPrev() - _cameraObject->zAttPrev()*_cameraObject->zAttPrev());
        const double cameraWAtt = std::sqrt(1. - _cameraObject->xAtt()*
            _cameraObject->xAtt() - _cameraObject->yAtt()*_cameraObject->yAtt() -
            _cameraObject->zAtt()*_cameraObject->zAtt());
        const Vec cameraAtt = nlerp(Vec{{_cameraObject->xAttPrev(),
            _cameraObject->yAttPrev(), _cameraObject->zAttPrev(),
            cameraWAttPrev}}, Vec{{_cameraObject->xAtt(), _cameraObject->yAtt(),
            _cameraObject->zAtt(), cameraWAtt}}, fac);

        const Vec cameraPos{{mix(_cameraObject->xPrev(), _cameraObject->x(), fac),
            mix(_cameraObject->yPrev(), _cameraObject->y(), fac), mix(
            _cameraObject->zPrev(), _cameraObject->z(), fac)}};

        const auto projection = perspective(1., Window::AspectRatio(), 0.1,
            1e3);
        Mat view = Mat::Zero(4);
        view(3, 3) = 1.;
        view.setBlock(0, 0, 3, 3, Mat{{0., -1., 0.}, {0., 0., 1.}, {-1., 0.,
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
                // Address held by `paz::Model::_t` identifies all copies of the
                // same model.
                objectsByModel[o->model()._t.get()].push_back(o);
            }
        }
        for(const auto& n : objectsByModel)
        {
            if(!_instances.count(n.first) || _instances.at(n.first).size() != n.
                second.size()) //TEMP - should only need to change _numInstances if larger, not reinit whole buffer
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
                const double wAttPrev = std::sqrt(1. - n.second[i]->xAttPrev()*
                    n.second[i]->xAttPrev() - n.second[i]->yAttPrev()*n.second[
                    i]->yAttPrev() - n.second[i]->zAttPrev()*n.second[i]->
                    zAttPrev());
                const double wAtt = std::sqrt(1. - n.second[i]->xAtt()*n.second[
                    i]->xAtt() - n.second[i]->yAtt()*n.second[i]->yAtt() - n.
                    second[i]->zAtt()*n.second[i]->zAtt());
                const Vec att = nlerp(Vec{{n.second[i]->xAttPrev(), n.second[
                    i]->yAttPrev(), n.second[i]->zAttPrev(), wAttPrev}}, Vec{{
                    n.second[i]->xAtt(), n.second[i]->yAtt(), n.second[i]->
                    zAtt(), wAtt}}, fac);
                modelMatData[0][4*i + 0] = att(0);
                modelMatData[0][4*i + 1] = att(1);
                modelMatData[0][4*i + 2] = att(2);
                modelMatData[0][4*i + 3] = mix(n.second[i]->xPrev(), n.second[
                    i]->x(), fac) - cameraPos(0);
                modelMatData[1][2*i + 0] = mix(n.second[i]->yPrev(), n.second[
                    i]->y(), fac) - cameraPos(1);
                modelMatData[1][2*i + 1] = mix(n.second[i]->zPrev(), n.second[
                    i]->z(), fac) - cameraPos(2);
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
                const double wAttPrev = std::sqrt(1. - n->xAttPrev()*n->
                    xAttPrev() - n->yAttPrev()*n->yAttPrev() - n->zAttPrev()*n->
                    zAttPrev());
                const double wAtt = std::sqrt(1. - n->xAtt()*n->xAtt() - n->
                    yAtt()*n->yAtt() - n->zAtt()*n->zAtt());
                const Vec att = nlerp(Vec{{n->xAttPrev(), n->yAttPrev(), n->
                    zAttPrev(), wAttPrev}}, Vec{{n->xAtt(), n->yAtt(), n->
                    zAtt(), wAtt}}, fac);
                const double xx = att(0)*att(0);
                const double yy = att(1)*att(1);
                const double zz = att(2)*att(2);
                const double xy = att(0)*att(1);
                const double zw = att(2)*att(3);
                const double xz = att(0)*att(2);
                const double yw = att(1)*att(3);
                const double yz = att(1)*att(2);
                const double xw = att(0)*att(3);
                const Mat mv = view*Mat{{1. - 2.*(yy + zz), 2.*(xy - zw), 2.*(xz + yw), mix(n->xPrev(), n->x(), fac) - cameraPos(0)},
                                        {2.*(xy + zw), 1. - 2.*(xx + zz), 2.*(yz - xw), mix(n->yPrev(), n->y(), fac) - cameraPos(1)},
                                        {2.*(xz - yw), 2.*(yz + xw), 1. - 2.*(xx + yy), mix(n->zPrev(), n->z(), fac) - cameraPos(2)},
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
                _geometryPass.read("diffTex", n.second.back()->model()._diffTex);
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
                model()._v, _instances.at(n.first), n.second.back()->model()._i);
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

        if(!lights.empty())
        {
            _renderPass1.begin({LoadAction::Load});//, LoadAction::Load);
            _renderPass1.cull(CullMode::Front);
//            _renderPass1.depth(DepthTestMode::GreaterNoMask);
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
        }

        // Render transparent objects.
        bool oitEnabled = false;
        for(const auto& n : objectsByModel)
        {
            if(!n.second.back()->model()._transp.empty())
            {
                oitEnabled = true;
                break;
            }
        }
        if(oitEnabled)
        {
            unsigned int numLights = lightsData0.size()/4;
            if(numLights > 100)
            {
                throw std::runtime_error("Too many lights (" + std::to_string(numLights) + " > 100).");
            }
            _oitAccumPass.begin({LoadAction::Clear, LoadAction::Clear}, LoadAction::Load);
            _oitAccumPass.depth(DepthTestMode::LessNoMask);
            _oitAccumPass.uniform("numLights", numLights);
            _oitAccumPass.uniform("light0", lightsData0);
            _oitAccumPass.uniform("light1", lightsData1);
            _oitAccumPass.uniform("projection", projection);
            _oitAccumPass.uniform("invProjection", convert_mat(convert_mat(projection).inv())); //TEMP - precompute - used elsewhere
            _oitAccumPass.uniform("sunDir", convert_vec(view*_sunDir)); //TEMP - precompute
            _oitAccumPass.uniform("sunIll", _sunIll);
            _oitAccumPass.uniform("view", convert_mat(view)); //TEMP - precompute
            for(const auto& n : objectsByModel)
            {
                if(!n.second.back()->model()._transp.empty())
                {
                    _oitAccumPass.draw(PrimitiveType::Triangles, n.second.back()->model()._transp, _instances.at(n.first));
                }
            }
            _oitAccumPass.end();

            _oitCompositePass.begin({LoadAction::Load});
            _oitCompositePass.read("accumTex", _oitAccumTex);
            _oitCompositePass.read("revealTex", _oitRevealTex);
            _oitCompositePass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _oitCompositePass.end();
        }

        const bool dofEnabled = _dofMinDepth < _dofMaxDepth;
        if(dofEnabled)
        {
            if(oitEnabled) //TEMP - can we get nearest transparent depth while performing OIT ?
            {
                _oitDepthPass.begin();
                _oitDepthPass.depth(DepthTestMode::Less);
                _oitDepthPass.uniform("project", projection);
                _oitDepthPass.uniform("view", convert_mat(view));
                for(const auto& n : objectsByModel)
                {
                    if(!n.second.back()->model()._transp.empty())
                    {
                        _oitDepthPass.draw(PrimitiveType::Triangles, n.second.back()->model()._transp, _instances.at(n.first));
                    }
                }
                _oitDepthPass.end();
            }

            _dofPass.begin();
            _dofPass.read("hdrRender", _hdrRender);
            _dofPass.read("depthMap", _depthMap);
            _dofPass.uniform("minDepth", _dofMinDepth);
            _dofPass.uniform("maxDepth", _dofMaxDepth);
            _dofPass.uniform("aspectRatio", Window::AspectRatio());
            _dofPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _dofPass.end();
        }

        if(_fxaaEnabled)
        {
            // Tonemap to linear LDR.
            _postPass0.begin();
            _postPass0.read("hdrRender", dofEnabled ? _dofRender : _hdrRender);
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
            _postPass1.read("hdrRender", dofEnabled ? _dofRender : _hdrRender);
            _postPass1.uniform("whitePoint", 1.f);
            _postPass1.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _postPass1.end();
        }

        InstanceBuffer _dialogChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(!_dialogs.empty())
        {
            if(!_dialogSkipFrame)
            {
                const float scale = std::max(1.f, std::round(2.f*FontScale*
                    Window::UiScale()));

                int maxCols = 0;
                std::vector<float> highlightAttr;
                std::vector<int> characterAttr;
                std::vector<int> colAttr;
                std::vector<int> rowAttr;
                bool highlight = false;
                int row = 0;
                int col = 0;
                for(auto n : _dialogs.front().first)
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
                    else if(n == '\n')
                    {
                        --row;
                        col = 0;
                    }
                    maxCols = std::max(maxCols, col);;
                }
                const int maxRows = 1 - row;
                const float _dialogWidth = (1 + maxCols*(CharWidth + 1))*scale/
                    Window::ViewportWidth();
                const float _dialogHeight = (maxRows + 1.f/_fontTex.height())*
                    scale*_fontTex.height()/Window::ViewportHeight();
                const float u = 0.5f - 0.5f*_dialogWidth;
                static constexpr float v = 0.1f;

                _consolePass.begin({LoadAction::Load});
                _consolePass.uniform("u", u);
                _consolePass.uniform("v", v);
                _consolePass.uniform("width", _dialogWidth);
                _consolePass.uniform("height", _dialogHeight);
                _consolePass.draw(PrimitiveType::TriangleStrip, _quadVertices);
                _consolePass.end();

                _dialogChars.addAttribute(1, highlightAttr);
                _dialogChars.addAttribute(1, characterAttr);
                _dialogChars.addAttribute(1, colAttr);
                _dialogChars.addAttribute(1, rowAttr);

                _textPass.begin({LoadAction::Load});
                _textPass.read("font", _fontTex);
                _textPass.uniform("u", u);
                const float v1 = v + (maxRows - 1)*scale*
                    _fontTex.height()/Window::ViewportHeight(); //TEMP - fix text/UI coords
                _textPass.uniform("v", v1);
                _textPass.uniform("width", Window::ViewportWidth());
                _textPass.uniform("height", Window::ViewportHeight());
                _textPass.uniform("charWidth", CharWidth);
                _textPass.uniform("baseWidth", _fontTex.width());
                _textPass.uniform("baseHeight", _fontTex.height());
                _textPass.uniform("scale", static_cast<int>(scale));
                _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices,
                    _dialogChars);
                _textPass.end();
            }

            if(!_paused)
            {
                _dialogTime += Window::FrameTime();
                if(_dialogs.front().second < _dialogTime)
                {
                    if(_dialogSkipFrame)
                    {
                        _dialogSkipFrame = false;
                        _dialogs.pop();
                        _dialogTime = 0.;
                    }
                    else
                    {
                        _dialogSkipFrame = true;
                    }
                }
            }
        }

        _cursorPass.begin({LoadAction::Load});
        _cursorPass.uniform("x", 0.5f);
        _cursorPass.uniform("y", 0.5f);
        _cursorPass.uniform("idx", _reticuleIdx);
        _cursorPass.uniform("aspect", static_cast<float>(_reticule.height())/
            _reticule.width());
        _cursorPass.uniform("h", _reticuleHighlight ? 1.f : -1.f);
        _cursorPass.uniform("scale", static_cast<int>(std::max(1.f, std::round(
            20.f*Window::UiScale())))); //TEMP
        _cursorPass.uniform("width", Window::ViewportWidth());
        _cursorPass.uniform("height", Window::ViewportHeight());
        _cursorPass.read("tex", _reticule);
        _cursorPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _cursorPass.end();

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
        if(_consoleMode == ConsoleMode::LatestStep)
        {
            _msgStream.clear();
            _msgStream.seekg(0);
        }
        else
        {
            std::stringstream{}.swap(_msgStream);
        }
        InstanceBuffer consoleChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(_consoleMode != ConsoleMode::Disable && !_console.empty())
        {
            const float scale = std::max(1.f, std::round(Window::UiScale()));

            int maxVisRows = std::ceil(Window::ViewportHeight()/(scale*
                _consoleFontTex.height()));
            maxVisRows = std::min(static_cast<std::size_t>(maxVisRows),
                _console.size());
            if(_consoleMode == ConsoleMode::CurrentFrame || _consoleMode ==
                ConsoleMode::LatestStep)
            {
                maxVisRows = std::min(static_cast<std::size_t>(maxVisRows),
                    numNewRows);
            }

            const float consoleHeight = std::min(1.f, (maxVisRows + 1.f/
                _consoleFontTex.height())*scale*_consoleFontTex.height()/Window::
                ViewportHeight());
            std::size_t maxCols = 0;
            for(int i = 0; i < maxVisRows; ++i)
            {
                maxCols = std::max(maxCols, _console.rbegin()[i].size());
            }
            const float consoleWidth = std::min(1.f, (1 + maxCols*(CharWidth +
                1))*scale/Window::ViewportWidth());

            _consolePass.begin({LoadAction::Load});
            _consolePass.uniform("u", 0.f);
            _consolePass.uniform("v", 0.f);
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
                    if(n == ' ')
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
            _textPass.read("font", _consoleFontTex);
            _textPass.uniform("u", 0.f);
            _textPass.uniform("v", 0.f);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", CharWidth);
            _textPass.uniform("baseWidth", _consoleFontTex.width());
            _textPass.uniform("baseHeight", _consoleFontTex.height());
            _textPass.uniform("scale", static_cast<int>(scale));
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices,
                consoleChars);
            _textPass.end();
        }

        InstanceBuffer menuChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(_paused)
        {
            _consolePass.begin({LoadAction::Load});
            _consolePass.uniform("u", 0.f);
            _consolePass.uniform("v", 0.f);
            _consolePass.uniform("width", 1.f);
            _consolePass.uniform("height", 1.f);
            _consolePass.draw(PrimitiveType::TriangleStrip, _quadVertices);
            _consolePass.end();

            if(!justPaused)
            {
                if(pauseMenu.curPage() == 0 && (Window::KeyPressed(Key::Escape)
                    || Window::GamepadPressed(GamepadButton::Start)))
                {
                    _paused = false;
                    Window::SetCursorMode(CursorMode::Disable);
                }
                else if(pauseMenu.curPage() == 1 && (Window::KeyPressed(Key::
                    Escape) || Window::MousePressed(paz::MouseButton::Back) ||
                    Window::GamepadPressed(GamepadButton::Start)))
                {
                    pauseMenu.setState(0, 1);
                }
            }

            pauseMenu.update();

            _textPass.begin({LoadAction::Load});
            _textPass.read("font", pauseMenu.font().tex());
            _textPass.uniform("u", 0.f);
            _textPass.uniform("v", 0.f);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", pauseMenu.font().charWidth());
            _textPass.uniform("baseWidth", pauseMenu.font().tex().width());
            _textPass.uniform("baseHeight", pauseMenu.font().tex().height());
            _textPass.uniform("scale", pauseMenu.font().curScale());
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices, pauseMenu.
                chars());
            _textPass.end();

            if(Window::MouseActive())
            {
                _cursorPass.begin({LoadAction::Load});
                _cursorPass.uniform("x", static_cast<float>(Window::MousePos().
                    first)/(Window::Width() - 1));
                _cursorPass.uniform("y", static_cast<float>(Window::MousePos().
                    second)/(Window::Height() - 1));
                _cursorPass.uniform("idx", 0);
                _cursorPass.uniform("aspect", static_cast<float>(_cursor.
                    height())/_cursor.width());
                _cursorPass.uniform("h", static_cast<float>(pauseMenu.
                    curButtonEnabled()));
                _cursorPass.uniform("scale", static_cast<int>(std::max(1.f, std::
                    round(20.f*Window::UiScale())))); //TEMP
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
    _cameraObject.reset(o);
}

void paz::App::AttachMic(const Object& o)
{
    _micObject.reset(o);
}

std::stringstream& paz::App::MsgStream()
{
    return _msgStream;
}

paz::Bytes paz::App::GetAsset(const std::string& path)
{
    return get_asset(path);
}

void paz::App::SetConsole(ConsoleMode mode)
{
    _consoleMode = mode;
}

void paz::App::SetGravity(double acc)
{
    _gravAcc = acc;
}

void paz::App::SetSound(const Object& o, const AudioTrack& sound, bool loop)
{
    _soundSrc.reset(o);
    AudioEngine::Play(sound, loop);
}

void paz::App::SetSun(const Vec& dir, const Vec& ill)
{
    std::copy(dir.begin(), dir.end(), _sunDir.begin());
    std::copy(ill.begin(), ill.end(), _sunIll.begin());
}

void paz::App::PushDialog(const std::string& msg, double time)
{
    _dialogs.emplace(msg, time);
}

void paz::App::SetReticule(int n, bool h)
{
    _reticuleIdx = n;
    _reticuleHighlight = h;
}
