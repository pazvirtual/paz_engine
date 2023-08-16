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

static paz::Framebuffer GeometryBuffer;
static paz::Framebuffer OitAccumBuffer;
static paz::Framebuffer RenderBuffer;
static paz::Framebuffer PostBuffer;
static paz::Framebuffer LumBuffer;

static paz::RenderTarget DiffuseMap(paz::TextureFormat::RGBA16Float);
// ...
static paz::RenderTarget EmissMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget NormalMap(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget DepthMap(paz::TextureFormat::Depth32Float);
static paz::RenderTarget HdrRender(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget FinalRender(paz::TextureFormat::RGBA16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget FinalLumMap(paz::TextureFormat::R16Float, paz::MinMagFilter::Linear, paz::MinMagFilter::Linear);
static paz::RenderTarget OitAccumTex(paz::TextureFormat::RGBA16Float);
static paz::RenderTarget OitRevealTex(paz::TextureFormat::RGBA8UNorm);

static std::unordered_map<void*, paz::InstanceBuffer> Instances;

static paz::RenderPass GeometryPass;
static paz::RenderPass RenderPass0;
static paz::RenderPass RenderPass1;
static paz::RenderPass OitAccumPass;
static paz::RenderPass OitCompositePass;
static paz::RenderPass PostPass0;
static paz::RenderPass LumPass;
static paz::RenderPass FxaaPass;
static paz::RenderPass PostPass1;
static paz::RenderPass ConsolePass;
static paz::RenderPass TextPass;
static paz::RenderPass CursorPass;

static paz::VertexBuffer QuadVertices;
static paz::VertexBuffer SphereVertices;
static paz::IndexBuffer SphereIndices;

static paz::Texture FontTex;
static std::stringstream MsgStream;
static std::deque<std::string> Console;
static paz::ConsoleMode CurConsoleMode = paz::ConsoleMode::Disable;

static std::string Title;

static paz::Texture DefaultDiffTex;

static paz::Texture Cursor;

static paz::ObjectPtr CameraObject;
static paz::ObjectPtr MicObject;

static paz::ObjectPtr SoundSrc;

static bool Paused;
static double AccumTime;
static constexpr double Timestep = 1./60.; //TEMP
static int LookSensitivity = 5;

static double GravAcc;
static paz::Threadpool Threads;

static bool FxaaEnabled = true;

static paz::Vec SunDir = paz::Vec::Zero(4);
static std::array<float, 4> SunIll;

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
        [](){ return FxaaEnabled ? "FXAA:        ON" : "FXAA:        OFF"; },
        [](paz::Menu&)
        {
            FxaaEnabled = !FxaaEnabled;
            paz::save_setting("fxaa", FxaaEnabled ? "1" : "0");
        }
    },
    {
        [](){ return "Sensitivity: " + std::to_string(LookSensitivity); },
        [](paz::Menu&)
        {
            ++LookSensitivity;
            if(LookSensitivity > 10)
            {
                LookSensitivity = 1;
            }
            paz::save_setting("sensitivity", std::to_string(LookSensitivity));
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
        FxaaEnabled = false;
    }
    {
        const auto str = load_setting("sensitivity");
        if(!str.empty())
        {
            LookSensitivity = std::stoi(str);
        }
    }

    Window::EnableDithering();

    Title = title;

    GeometryBuffer.attach(DiffuseMap);
    // ...
    GeometryBuffer.attach(EmissMap);
    GeometryBuffer.attach(NormalMap);
    GeometryBuffer.attach(DepthMap);

    RenderBuffer.attach(HdrRender);

    OitAccumBuffer.attach(OitAccumTex);
    OitAccumBuffer.attach(OitRevealTex);
    OitAccumBuffer.attach(DepthMap);

    PostBuffer.attach(FinalRender);

    LumBuffer.attach(FinalLumMap);

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

    GeometryPass = RenderPass(GeometryBuffer, geometryVert, geometryFrag);
    RenderPass0 = RenderPass(RenderBuffer, sceneVert0, sceneFrag0);
    RenderPass1 = RenderPass(RenderBuffer, sceneVert1, sceneFrag1, {BlendMode::
        One_One});
    OitAccumPass = RenderPass(OitAccumBuffer, oitVert, oitFrag, {BlendMode::
        One_One, BlendMode::Zero_InvSrcAlpha});
    OitCompositePass = RenderPass(RenderBuffer, quadVert, compositeFrag,
        {BlendMode::InvSrcAlpha_SrcAlpha});
    PostPass0 = RenderPass(PostBuffer, quadVert, postFrag);
    LumPass = RenderPass(LumBuffer, quadVert, lumFrag);
    FxaaPass = RenderPass(quadVert, fxaaFrag);
    PostPass1 = RenderPass(quadVert, postFrag);
    ConsolePass = RenderPass(quadVert, consoleFrag, {BlendMode::
        SrcAlpha_InvSrcAlpha});
    TextPass = RenderPass(textVert, textFrag, {BlendMode::
        SrcAlpha_InvSrcAlpha});
    CursorPass = RenderPass(cursorVert, cursorFrag, {BlendMode::
        SrcAlpha_InvSrcAlpha});

    QuadVertices.addAttribute(2, QuadPos);

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
        SphereVertices.addAttribute(4, positions[0]);
        SphereIndices = IndexBuffer(indices[0]);
    }

    FontTex = Texture(get_asset_image("font.pbm")); //TEMP - note that only red channel is used

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
    DefaultDiffTex = Texture(TextureFormat::RGBA8UNorm_sRGB, 512, 512, temp.
        data(), MinMagFilter::Linear, MinMagFilter::Linear, MipmapFilter::
        Linear, WrapMode::Repeat, WrapMode::Repeat);

    Cursor = Texture(get_asset_image("cursor.pbm")); //TEMP - note that only red channel is used
}

void paz::App::Run()
{
    Font menuFont(FontTex, 2*FontScale, CharWidth);

    {
        bool done = false;
        Menu startMenu(menuFont, Title,
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

            TextPass.begin({LoadAction::Clear});
            TextPass.read("font", startMenu.font().tex());
            TextPass.uniform("width", Window::ViewportWidth());
            TextPass.uniform("height", Window::ViewportHeight());
            TextPass.uniform("charWidth", startMenu.font().charWidth());
            TextPass.uniform("baseWidth", startMenu.font().tex().width());
            TextPass.uniform("baseHeight", startMenu.font().tex().height());
            TextPass.uniform("scale", startMenu.font().curScale());
            TextPass.draw(PrimitiveType::TriangleStrip, QuadVertices, startMenu.
                chars());
            TextPass.end();

            if(Window::MouseActive())
            {
                CursorPass.begin({LoadAction::Load});
                CursorPass.uniform("x", static_cast<float>(Window::MousePos().
                    first)/(Window::Width() - 1));
                CursorPass.uniform("y", static_cast<float>(Window::MousePos().
                    second)/(Window::Height() - 1));
                CursorPass.uniform("h", static_cast<float>(startMenu.
                    curButtonEnabled()));
                CursorPass.uniform("scale", static_cast<int>(std::round(20.*
                    Window::UiScale()))); //TEMP
                CursorPass.uniform("width", Window::ViewportWidth());
                CursorPass.uniform("height", Window::ViewportHeight());
                CursorPass.read("tex", Cursor);
                CursorPass.draw(PrimitiveType::TriangleStrip, QuadVertices);
                CursorPass.end();
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
                    Paused = false;
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
        if(!Paused && (Window::KeyPressed(Key::Escape) || Window::
            GamepadPressed(GamepadButton::Start)))
        {
            justPaused = true;
            Paused = true;
            pauseMenu.setState(0, 0);
        }

        if(!Paused && MicObject && SoundSrc)
        {
            const Vec relPos{{SoundSrc->x() - MicObject->x(), SoundSrc->y() -
                MicObject->y(), SoundSrc->z() - MicObject->z()}};
            const Vec relVel{{SoundSrc->xVel() - MicObject->xVel(), SoundSrc->
                yVel() - MicObject->yVel(), SoundSrc->zVel() - MicObject->
                zVel()}};
            const double dist = relPos.norm();
            const Vec dir = relPos/dist;
            const Vec micAtt{{MicObject->xAtt(), MicObject->yAtt(), MicObject->
                zAtt(), std::sqrt(1. - MicObject->xAtt()*MicObject->xAtt() -
                MicObject->yAtt()*MicObject->yAtt() - MicObject->zAtt()*
                MicObject->zAtt())}};
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

        if(Paused)
        {
            AccumTime = 0.;
            input.resetEvents();
        }
        else
        {
            AccumTime += Window::FrameTime();
            input.copyEvents(Timestep, 0.2*LookSensitivity);

            while(AccumTime > 0.)
            {
                if(CurConsoleMode == ConsoleMode::LatestStep)
                {
                    std::stringstream{}.swap(::MsgStream);
                }

                do_physics(GravAcc, Timestep);
                do_collisions(Threads, Timestep);
                const auto tempObjects = objects(); //TEMP - this prevents missed or multiple updates when `objects()` changes, but is not ideal
                for(const auto& n : tempObjects)
                {
                    if(objects().count(n.first))
                    {
                        reinterpret_cast<Object*>(n.first)->update(input);
                    }
                }

                input.resetEvents();
                AccumTime -= Timestep;
            }
        }

        const double cameraWAtt = std::sqrt(1. - CameraObject->xAtt()*
            CameraObject->xAtt() - CameraObject->yAtt()*CameraObject->yAtt() -
            CameraObject->zAtt()*CameraObject->zAtt());
        const Vec cameraAtt{{CameraObject->xAtt(), CameraObject->yAtt(),
            CameraObject->zAtt(), cameraWAtt}};

        const Vec cameraPos{{CameraObject->x(), CameraObject->y(),
            CameraObject->z()}};

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
            if(!Instances.count(n.first) || Instances.at(n.first).size() != n.
                second.size()) //TEMP - should only need to change _numInstances if larger, not reinit whole buffer
            {
                Instances[n.first] = InstanceBuffer(n.second.size());
                Instances.at(n.first).addAttribute(4, DataType::Float);
                Instances.at(n.first).addAttribute(2, DataType::Float);
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
            Instances.at(n.first).subAttribute(0, modelMatData[0]);
            Instances.at(n.first).subAttribute(1, modelMatData[1]);
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
        GeometryPass.begin(std::vector<LoadAction>(4, LoadAction::Clear),
            LoadAction::Clear);
        GeometryPass.cull(CullMode::Back);
        GeometryPass.depth(DepthTestMode::Less);
        GeometryPass.uniform("projection", projection);
        GeometryPass.uniform("view", convert_mat(view));
        for(const auto& n : objectsByModel)
        {
            if(n.second.back()->model()._i.empty())
            {
                continue;
            }
            if(n.second.back()->model()._diffTex.width())
            {
                GeometryPass.read("diffTex", n.second.back()->model()._diffTex);
            }
            else
            {
                GeometryPass.read("diffTex", DefaultDiffTex);
            }
            std::array<float, 4> emiss;
            std::copy(n.second.back()->model()._emiss.begin(), n.second.back()->
                model()._emiss.end(), emiss.begin());
            emiss[3] = 1.f;
            GeometryPass.uniform("emiss", emiss);
            GeometryPass.draw(PrimitiveType::Triangles, n.second.back()->
                model()._v, Instances.at(n.first), n.second.back()->model()._i);
        }
        GeometryPass.end();

        // Render in HDR.
        RenderPass0.begin({LoadAction::Clear});//, LoadAction::Load);
        RenderPass0.read("diffuseMap", DiffuseMap);
        // ...
        RenderPass0.read("emissMap", EmissMap);
        RenderPass0.read("normalMap", NormalMap);
        RenderPass0.read("depthMap", DepthMap);
        RenderPass0.uniform("invProjection", convert_mat(convert_mat(
            projection).inv()));
        RenderPass0.uniform("lightDir", convert_vec(view*SunDir));
        RenderPass0.uniform("ill", SunIll);
        RenderPass0.draw(PrimitiveType::TriangleStrip, QuadVertices);
        RenderPass0.end();

        RenderPass1.begin({LoadAction::Load});//, LoadAction::Load);
        RenderPass1.cull(CullMode::Front);
//        RenderPass1.depth(DepthTestMode::GreaterNoMask);
        RenderPass1.read("diffuseMap", DiffuseMap);
        // ...
        RenderPass1.read("emissMap", EmissMap);
        RenderPass1.read("normalMap", NormalMap);
        RenderPass1.read("depthMap", DepthMap);
        RenderPass1.uniform("projection", projection);
        RenderPass1.uniform("invProjection", convert_mat(convert_mat(
            projection).inv()));
        RenderPass1.draw(PrimitiveType::Triangles, SphereVertices, lights,
            SphereIndices);
        RenderPass1.end();

        // Render transparent objects. //TEMP - skip if possible
        unsigned int numLights = lightsData0.size()/4;
        if(numLights > 16)
        {
            throw std::runtime_error("Too many lights (" + std::to_string(numLights) + " > 16).");
        }
        OitAccumPass.begin({LoadAction::Clear, LoadAction::Clear}, LoadAction::
            Load);
        OitAccumPass.depth(DepthTestMode::LessNoMask);
        OitAccumPass.uniform("numLights", numLights);
        OitAccumPass.uniform("light0", lightsData0);
        OitAccumPass.uniform("light1", lightsData1);
        OitAccumPass.uniform("projection", projection);
        OitAccumPass.uniform("invProjection", convert_mat(convert_mat(
            projection).inv())); //TEMP - precompute - used elsewhere
        OitAccumPass.uniform("sunDir", convert_vec(view*SunDir)); //TEMP - precompute
        OitAccumPass.uniform("sunIll", SunIll);
        OitAccumPass.uniform("view", convert_mat(view));
        for(const auto& n : objectsByModel)
        {
            if(!n.second.back()->model()._transp.empty())
            {
                OitAccumPass.draw(PrimitiveType::Triangles, n.second.back()->
                    model()._transp, Instances.at(n.first));
            }
        }
        OitAccumPass.end();

        OitCompositePass.begin({LoadAction::Load});
        OitCompositePass.read("accumTex", OitAccumTex);
        OitCompositePass.read("revealTex", OitRevealTex);
        OitCompositePass.draw(PrimitiveType::TriangleStrip, QuadVertices);
        OitCompositePass.end();

        if(FxaaEnabled)
        {
            // Tonemap to linear LDR.
            PostPass0.begin();
            PostPass0.read("hdrRender", HdrRender);
            PostPass0.uniform("whitePoint", 1.f);
            PostPass0.draw(PrimitiveType::TriangleStrip, QuadVertices);
            PostPass0.end();

            // Get luminance map.
            LumPass.begin();
            LumPass.read("img", FinalRender);
            LumPass.draw(PrimitiveType::TriangleStrip, QuadVertices);
            LumPass.end();

            // Anti-alias.
            FxaaPass.begin();
            FxaaPass.read("img", FinalRender);
            FxaaPass.read("lum", FinalLumMap);
            FxaaPass.draw(PrimitiveType::TriangleStrip, QuadVertices);
            FxaaPass.end();
        }
        else
        {
            // Tonemap to linear LDR.
            PostPass1.begin();
            PostPass1.read("hdrRender", HdrRender);
            PostPass1.uniform("whitePoint", 1.f);
            PostPass1.draw(PrimitiveType::TriangleStrip, QuadVertices);
            PostPass1.end();
        }

        std::string line;
        std::size_t numNewRows = 0;
        while(std::getline(::MsgStream, line))
        {
            Console.push_back(line);
            ++numNewRows;
            if(Console.size() > MaxConsoleLines)
            {
                Console.pop_front();
            }
        }
        if(CurConsoleMode == ConsoleMode::LatestStep)
        {
            ::MsgStream.clear();
            ::MsgStream.seekg(0);
        }
        else
        {
            std::stringstream{}.swap(::MsgStream);
        }
        InstanceBuffer consoleChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(CurConsoleMode != ConsoleMode::Disable && !Console.empty())
        {
            const float scale = std::round(FontScale*Window::UiScale());
            int maxVisRows = std::ceil(Window::ViewportHeight()/(scale*FontTex.
                height()));
            maxVisRows = std::min(static_cast<std::size_t>(maxVisRows),
                Console.size());
            if(CurConsoleMode == ConsoleMode::CurrentFrame || CurConsoleMode ==
                ConsoleMode::LatestStep)
            {
                maxVisRows = std::min(static_cast<std::size_t>(maxVisRows),
                    numNewRows);
            }

            const float consoleHeight = std::min(1.f, (maxVisRows + 1.f/FontTex.
                height())*scale*FontTex.height()/Window::ViewportHeight());
            std::size_t maxCols = 0;
            for(int i = 0; i < maxVisRows; ++i)
            {
                maxCols = std::max(maxCols, Console.rbegin()[i].size());
            }
            const float consoleWidth = std::min(1.f, (1 + maxCols*(CharWidth +
                1))*scale/Window::ViewportWidth());

            ConsolePass.begin({LoadAction::Load});
            ConsolePass.uniform("width", consoleWidth);
            ConsolePass.uniform("height", consoleHeight);
            ConsolePass.draw(PrimitiveType::TriangleStrip, QuadVertices);
            ConsolePass.end();

            std::vector<float> highlightAttr;
            std::vector<int> characterAttr;
            std::vector<int> colAttr;
            std::vector<int> rowAttr;
            bool highlight = false;
            for(int row = 0; row < maxVisRows; ++row)
            {
                int col = 0;
                for(auto n : Console.rbegin()[row])
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

            TextPass.begin({LoadAction::Load});
            TextPass.read("font", FontTex);
            TextPass.uniform("width", Window::ViewportWidth());
            TextPass.uniform("height", Window::ViewportHeight());
            TextPass.uniform("charWidth", CharWidth);
            TextPass.uniform("baseWidth", FontTex.width());
            TextPass.uniform("baseHeight", FontTex.height());
            TextPass.uniform("scale", static_cast<int>(scale));
            TextPass.draw(PrimitiveType::TriangleStrip, QuadVertices,
                consoleChars);
            TextPass.end();
        }

        InstanceBuffer menuChars; //TEMP - needs wider scope because of MTLBuffer purge issue
        if(Paused)
        {
            ConsolePass.begin({LoadAction::Load});
            ConsolePass.uniform("width", 1.f);
            ConsolePass.uniform("height", 1.f);
            ConsolePass.draw(PrimitiveType::TriangleStrip, QuadVertices);
            ConsolePass.end();

            if(!justPaused && (Window::KeyPressed(Key::Escape) || Window::
                GamepadPressed(GamepadButton::Start)))
            {
                if(pauseMenu.curPage() == 0)
                {
                    Paused = false;
                    Window::SetCursorMode(CursorMode::Disable);
                }
                else if(pauseMenu.curPage() == 1)
                {
                    pauseMenu.setState(0, 1);
                }
            }

            pauseMenu.update();

            TextPass.begin({LoadAction::Load});
            TextPass.read("font", pauseMenu.font().tex());
            TextPass.uniform("width", Window::ViewportWidth());
            TextPass.uniform("height", Window::ViewportHeight());
            TextPass.uniform("charWidth", pauseMenu.font().charWidth());
            TextPass.uniform("baseWidth", pauseMenu.font().tex().width());
            TextPass.uniform("baseHeight", pauseMenu.font().tex().height());
            TextPass.uniform("scale", pauseMenu.font().curScale());
            TextPass.draw(PrimitiveType::TriangleStrip, QuadVertices, pauseMenu.
                chars());
            TextPass.end();

            if(Window::MouseActive())
            {
                CursorPass.begin({LoadAction::Load});
                CursorPass.uniform("x", static_cast<float>(Window::MousePos().
                    first)/(Window::Width() - 1));
                CursorPass.uniform("y", static_cast<float>(Window::MousePos().
                    second)/(Window::Height() - 1));
                CursorPass.uniform("h", static_cast<float>(pauseMenu.
                    curButtonEnabled()));
                CursorPass.uniform("scale", static_cast<int>(std::round(20.*
                    Window::UiScale()))); //TEMP
                CursorPass.uniform("width", Window::ViewportWidth());
                CursorPass.uniform("height", Window::ViewportHeight());
                CursorPass.read("tex", Cursor);
                CursorPass.draw(PrimitiveType::TriangleStrip, QuadVertices);
                CursorPass.end();
            }
        }

        Window::EndFrame();
    }
}

void paz::App::AttachCamera(const Object& o)
{
    CameraObject.reset(o);
}

void paz::App::AttachMic(const Object& o)
{
    MicObject.reset(o);
}

std::stringstream& paz::App::MsgStream()
{
    return ::MsgStream;
}

paz::Bytes paz::App::GetAsset(const std::string& path)
{
    return get_asset(path);
}

void paz::App::SetConsole(ConsoleMode mode)
{
    CurConsoleMode = mode;
}

void paz::App::SetGravity(double acc)
{
    GravAcc = acc;
}

void paz::App::SetSound(const Object& o, const AudioTrack& sound, bool loop)
{
    SoundSrc.reset(o);
    AudioEngine::Play(sound, loop);
}

void paz::App::SetSun(const Vec& dir, const Vec& ill)
{
    std::copy(dir.begin(), dir.end(), SunDir.begin());
    std::copy(ill.begin(), ill.end(), SunIll.begin());
}
