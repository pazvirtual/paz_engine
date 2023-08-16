#include "object.hpp"
#include "shared.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>
#include <iomanip>
#include <deque>

//#define DO_FXAA
#define NO_FRICTION

static constexpr std::size_t NumSteps = 100;
static constexpr double InteractRangeBehindSq = 4.;
static constexpr double InteractRangeInFrontSq = 9.;
static constexpr std::size_t MaxConsoleLines = 1000;
static constexpr float FontScale = 1.5f;
static constexpr int CharWidth = 5;

static paz::Framebuffer _geometryBuffer;
static paz::Framebuffer _renderBuffer;
#ifdef DO_FXAA
static paz::Framebuffer _postBuffer;
static paz::Framebuffer _lumBuffer;
#endif

static paz::RenderTarget _diffuseMap(paz::TextureFormat::RGBA16Float);
// ...
static paz::RenderTarget _emissMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _normalMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget _depthMap(paz::TextureFormat::Depth32Float);
static paz::RenderTarget _hdrRender(paz::TextureFormat::RGBA16Float);
#ifdef DO_FXAA
static paz::RenderTarget _finalRender(paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget _finalLumMap(paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
#endif

static std::unordered_map<void*, paz::InstanceBuffer> _instances;

static paz::RenderPass _geometryPass;
static paz::RenderPass _renderPass;
static paz::RenderPass _postPass;
#ifdef DO_FXAA
static paz::RenderPass _lumPass;
static paz::RenderPass _fxaaPass;
#endif
static paz::RenderPass _consolePass;
static paz::RenderPass _textPass;

static paz::VertexBuffer _quadVertices;

static paz::Texture _font;
static std::stringstream _msgStream;
static std::deque<std::string> _console;
static paz::ConsoleMode _consoleMode = paz::ConsoleMode::Disable;

static paz::UiDescriptor _startMenu;

static paz::Texture _defaultDiffTex;

static const paz::Object* _cameraObject; //TEMP - vector realloc breaks ptrs
static const paz::Object* _micObject; //TEMP - vector realloc breaks ptrs

static const paz::Object* _soundSrc; //TEMP - vector realloc breaks ptrs

static bool _paused;

static double _gravity;

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

void paz::App::Init(const std::string& sceneShaderPath, const std::string&
    fontPath, const UiDescriptor& startMenu)
{
    _startMenu = startMenu;

    _geometryBuffer.attach(_diffuseMap);
    // ...
    _geometryBuffer.attach(_emissMap);
    _geometryBuffer.attach(_normalMap);
    _geometryBuffer.attach(_depthMap);

    _renderBuffer.attach(_hdrRender);

#ifdef DO_FXAA
    _postBuffer.attach(_finalRender);

    _lumBuffer.attach(_finalLumMap);
#endif

    const VertexFunction geometryVert(get_builtin("geometry.vert").str());
    const VertexFunction quadVert(get_builtin("quad.vert").str());
    const VertexFunction textVert(get_builtin("text.vert").str());
    const FragmentFunction geometryFrag(get_builtin("geometry.frag").str());
    const FragmentFunction sceneFrag(get_asset(sceneShaderPath).str());
#ifdef DO_FXAA
    const FragmentFunction lumFrag(get_builtin("lum.frag").str());
    const FragmentFunction fxaaFrag(get_builtin("fxaa.frag").str());
#endif
    const FragmentFunction postFrag(get_builtin("post.frag").str());
    const FragmentFunction consoleFrag(get_builtin("console.frag").str());
    const FragmentFunction textFrag(get_builtin("text.frag").str());

    _geometryPass = RenderPass(_geometryBuffer, geometryVert, geometryFrag);
    _renderPass = RenderPass(_renderBuffer, quadVert, sceneFrag);
#ifdef DO_FXAA
    _postPass = RenderPass(_postBuffer, quadVert, postFrag);
    _lumPass = RenderPass(_lumBuffer, quadVert, lumFrag);
    _fxaaPass = RenderPass(quadVert, fxaaFrag);
#else
    _postPass = RenderPass(quadVert, postFrag);
#endif
    _consolePass = RenderPass(quadVert, consoleFrag, BlendMode::Blend);
    _textPass = RenderPass(textVert, textFrag, BlendMode::Blend);

    _quadVertices.addAttribute(2, QuadPos);

    _font = Texture(get_asset_image(fontPath)); //TEMP - note that only red channel is used

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
}

void paz::App::Run()
{
    if(true) // Has start menu
    {
        bool done = false;
        int curButton = 0;
        std::size_t maxCols = _startMenu._title.size();
        if(_startMenu._layout == UiLayout::Vertical)
        {
            for(const auto& n : _startMenu._buttons)
            {
                maxCols = std::max(maxCols, n.second.size());
            }
        }
        else
        {
            std::size_t buttonCols = _startMenu._buttons.size() > 1 ?
                _startMenu._buttons.size() - 1 : 0;
            for(const auto& n : _startMenu._buttons)
            {
                buttonCols += n.second.size();
            }
            maxCols = std::max(maxCols, buttonCols);
        }
        while(!Window::Done() && !done)
        {
            if(Window::GamepadActive())
            {
                Window::SetCursorMode(CursorMode::Disable);
            }
            else
            {
                Window::SetCursorMode(CursorMode::Normal);
            }

            if(Window::KeyPressed(Key::Space) || Window::KeyPressed(Key::Enter)
                || Window::KeyPressed(Key::KeypadEnter) || Window::
                GamepadPressed(GamepadButton::A))
            {
                switch(_startMenu._buttons[curButton].first)
                {
                    case UiAction::Start: done = true; break;
                    case UiAction::Quit: Window::Quit(); break;
                    case UiAction::ToggleFullscreen: Window::IsFullscreen() ?
                        Window::MakeWindowed() : Window::MakeFullscreen();
                        break;
                    //case UiAction::ToggleHidpi: Window::HidpiEnabled() ?
                    //    Window::DisableHidpi() : Window::EnableHidpi(); break;
                    case UiAction::None: break;
                    default: throw std::logic_error("Unrecognized UI action.");
                }
            }
            if(_startMenu._layout == UiLayout::Vertical)
            {
                if(Window::KeyPressed(Key::S) || Window::KeyPressed(Key::Down)
                    || Window::GamepadPressed(GamepadButton::Down))
                {
                    curButton = std::min(static_cast<std::size_t>(curButton) +
                        1, _startMenu._buttons.size() - 1);
                }
                if(Window::KeyPressed(Key::W) || Window::KeyPressed(Key::Up) ||
                    Window::GamepadPressed(GamepadButton::Up))
                {
                    curButton = std::max(curButton - 1, 0);
                }
            }
            else
            {
                if(Window::KeyPressed(Key::D) || Window::KeyPressed(Key::Right)
                    || Window::GamepadPressed(GamepadButton::Right))
                {
                    curButton = std::min(static_cast<std::size_t>(curButton) +
                        1, _startMenu._buttons.size() - 1);
                }
                if(Window::KeyPressed(Key::A) || Window::KeyPressed(Key::Left)
                    || Window::GamepadPressed(GamepadButton::Left))
                {
                    curButton = std::max(curButton - 1, 0);
                }
            }

            const float scale = std::round(2.f*FontScale*Window::UiScale());
            int maxVisRows = std::ceil(Window::ViewportHeight()/(scale*_font.
                height()));

            std::vector<float> highlightAttr;
            std::vector<int> characterAttr;
            std::vector<int> colAttr;
            std::vector<int> rowAttr;
            const int startRow = _startMenu._layout == UiLayout::Vertical ?
                (maxVisRows + _startMenu._buttons.size() - 1)/2 : (maxVisRows -
                1)/2;
            const int startCol = _startMenu._alignment == UiAlignment::Center ?
                std::round(0.5*(Window::ViewportWidth()/(scale*(CharWidth + 1))
                - maxCols)) : 0;
            int row = 0;
            int col = startCol;
            for(std::size_t i = 0; i < _startMenu._buttons.size() + 1; ++i)
            {
                const bool highlight = i == static_cast<std::size_t>(curButton)
                    + 1;
                for(auto n : (i ? _startMenu._buttons[i - 1].second :
                    _startMenu._title))
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
                        rowAttr.push_back(startRow - row);
                        ++col;
                    }
                }
                if(!i || _startMenu._layout == UiLayout::Vertical)
                {
                    ++row;
                    col = startCol;
                }
                else
                {
                    ++col;
                }
            }

            InstanceBuffer chars;
            chars.addAttribute(1, highlightAttr);
            chars.addAttribute(1, characterAttr);
            chars.addAttribute(1, colAttr);
            chars.addAttribute(1, rowAttr);

            _textPass.begin({paz::LoadAction::Clear});
            _textPass.read("font", _font);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", CharWidth);
            _textPass.uniform("baseWidth", _font.width());
            _textPass.uniform("baseHeight", _font.height());
            _textPass.uniform("scale", scale);
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices, chars);
            _textPass.end();

            Window::EndFrame();
        }
        if(Window::Done())
        {
            return;
        }
    }

    Window::SetCursorMode(CursorMode::Disable);

    while(!Window::Done())
    {
        if(Window::KeyPressed(Key::Escape) || Window::GamepadPressed(GamepadButton::Start))
        {
            _paused = !_paused;
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
            const double t0 = (vs + vlos)/(800.*M_PI*dist); // assuming f = 200
            const double rxPwr = std::min(maxRxPwr, txPwr*t0*t0);
            double lPwr = rxPwr*(0.6 + 0.4*(dir.dot(lEar))); //TEMP
            double rPwr = rxPwr*(0.6 + 0.4*(dir.dot(rEar))); //TEMP
            lPwr /= M_SQRT2*(0.6 + 0.4*std::cos(0.5*M_PI - 1.));
            rPwr /= M_SQRT2*(0.6 + 0.4*std::cos(0.5*M_PI - 1.));
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
                        a[j]->onCollide(*b[n]);
                        b[n]->onCollide(*a[j]);
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
            if(o->model()._i.empty())
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
        std::vector<float> lightsData;
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
                        lightsData.push_back(p(0));
                        lightsData.push_back(p(1));
                        lightsData.push_back(p(2));
                        lightsData.push_back(m[3]);
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
                    lightsData.push_back(p(0));
                    lightsData.push_back(p(1));
                    lightsData.push_back(p(2));
                    lightsData.push_back(m[3]);
                }
            }
        }
        Texture lights(TextureFormat::RGBA32Float, 1, lightsData.size()/4,
            lightsData.data());

        // Get geometry map.
        _geometryPass.begin(std::vector<LoadAction>(4, LoadAction::Clear),
            LoadAction::Clear);
        _geometryPass.cull(CullMode::Back);
        _geometryPass.depth(DepthTestMode::Less);
        _geometryPass.uniform("projection", projection);
        _geometryPass.uniform("view", convert_mat(view));
        for(const auto& n : objectsByModel)
        {
            if(n.second.back()->model()._diffTex.width())
            {
                _geometryPass.read("diffTex", n.second.back()->model().
                    _diffTex);
            }
            else
            {
                _geometryPass.read("diffTex", _defaultDiffTex);
            }
            _geometryPass.uniform("emiss", n.second.back()->model()._emiss);
            _geometryPass.draw(PrimitiveType::Triangles, n.second.back()->
                model()._v, _instances.at(n.first), n.second.back()->model().
                _i);
        }
        _geometryPass.end();

        // Render in HDR.
        _renderPass.begin();
        _renderPass.read("diffuseMap", _diffuseMap);
        // ...
        _renderPass.read("emissMap", _emissMap);
        _renderPass.read("normalMap", _normalMap);
        _renderPass.read("depthMap", _depthMap);
        _renderPass.uniform("invProjection", convert_mat(convert_mat(
            projection).inv()));
        _renderPass.read("lights", lights);
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

        // Anti-alias.
        _fxaaPass.begin();
        _fxaaPass.read("img", _finalRender);
        _fxaaPass.read("lum", _finalLumMap);
        _fxaaPass.draw(PrimitiveType::TriangleStrip, _quadVertices);
        _fxaaPass.end();
#endif

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

            _consolePass.begin();
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

            InstanceBuffer chars;
            chars.addAttribute(1, highlightAttr);
            chars.addAttribute(1, characterAttr);
            chars.addAttribute(1, colAttr);
            chars.addAttribute(1, rowAttr);

            _textPass.begin();
            _textPass.read("font", _font);
            _textPass.uniform("width", Window::ViewportWidth());
            _textPass.uniform("height", Window::ViewportHeight());
            _textPass.uniform("charWidth", CharWidth);
            _textPass.uniform("baseWidth", _font.width());
            _textPass.uniform("baseHeight", _font.height());
            _textPass.uniform("scale", scale);
            _textPass.draw(PrimitiveType::TriangleStrip, _quadVertices, chars);
            _textPass.end();
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
