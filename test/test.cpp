#include "player.hpp"
#include "npc.hpp"
#include "paintball.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

static paz::Model Sphere50;
static paz::Model Sphere10;
static paz::Model FancyBox;

static constexpr double Radius = 50.;

class Droplet : public Paintball
{
    bool _stuck = false;
    double _timer = 0.;

public:
    Droplet(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec& dir) :
        Paintball(pos, vel, dir) {}
    void update() override
    {
        if(_stuck)
        {
            _timer += paz::App::PhysTime();
            if(_timer >= 30.)
            {
                addTag("done");
            }
        }
    }
    void onCollide(const Object& o) override
    {
        _stuck = true;
        Paintball::onCollide(o);
    }
};

class Fountain : public paz::Object
{
    double _timer = 0.;
    std::vector<Droplet> _droplets;
    static constexpr int NumPerLaunch = 5;

public:
    Fountain()
    {
        model() = FancyBox;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None;
        lights().push_back({0., 0., 1.3});
    }
    void update() final
    {
        _timer += paz::App::PhysTime();
        if(_timer > 0.1)
        {
            const paz::Vec pos{{x(), y(), z()}};
            const paz::Vec vel{{xVel(), yVel(), zVel()}};
            for(int i = 0; i < NumPerLaunch; ++i)
            {
                const double wAtt = std::sqrt(1. - xAtt()*xAtt() - yAtt()*yAtt()
                    - zAtt()*zAtt());
                const paz::Vec att{{xAtt(), yAtt(), zAtt(), wAtt}};
                paz::Vec dir = paz::to_mat(paz::qinv(att)).col(2);
                dir += 0.01*paz::Vec{{paz::randn(), paz::randn(), paz::
                    randn()}};
                const double offset = 1.3 + i*0.1*Droplet::LaunchSpeed/
                    NumPerLaunch;
                _droplets.emplace_back(pos + offset*dir, vel, dir);
            }
            _timer = 0.;
        }
        std::size_t i = 0;
        while(i < _droplets.size())
        {
            if(_droplets[i].isTagged("done"))
            {
                std::swap(_droplets[i], _droplets.back());
                _droplets.pop_back();
            }
            else
            {
                ++i;
            }
        }
    }
    void setDir(double xDir, double yDir, double zDir)
    {
        const paz::Vec up{{xDir, yDir, zDir}};
        paz::Vec right{{1., 0., 0.}};
        if(std::abs(up.dot(right)) > 0.99)
        {
            right = paz::Vec{{0., 1., 0.}};
        }
        const paz::Vec forward = up.cross(right).normalized();
        right = forward.cross(up).normalized();
        paz::Mat rot(3);
        rot.setCol(0, right);
        rot.setCol(1, forward);
        rot.setCol(2, up);
        paz::Vec att = paz::to_quat(rot);
        if(att(3) < 0.)
        {
            att = -att;
        }
        xAtt() = -att(0);
        yAtt() = -att(1);
        zAtt() = -att(2);
    }
};

class World : public paz::Object
{
public:
    World()
    {
        model() = Sphere50;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None; //TEMP
        stdGravParam() = 0.6*9.81*Radius*Radius;
    }
};

class World1 : public paz::Object
{
public:
    World1()
    {
        model() = Sphere10;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None; //TEMP
        stdGravParam() = 0.1*9.81*Radius*Radius;
    }
};

class World2 : public paz::Object
{
    static constexpr double AngRate = 0.1;

    double _angle;
    Fountain _f;

public:
    World2(double initialAngle) : paz::Object(), _angle(initialAngle)
    {
        model() = Sphere10;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None; //TEMP
        stdGravParam() = 9.81*10.*10.;
        _f.z() = z() + 9.5;
        _f.setDir(0., 0., 1.);
    }
    void update() override
    {
        _angle = paz::normalize_angle(_angle + AngRate*paz::App::PhysTime());
        x() = std::cos(_angle)*2.*Radius;
        y() = std::sin(_angle)*2.*Radius;
        z() = 0.;
        xVel() = AngRate*-std::sin(_angle)*2.*Radius;
        yVel() = AngRate*std::cos(_angle)*2.*Radius;
        zVel() = 0.;
        _f.x() = x();
        _f.y() = y();
        _f.z() = z() + 9.5;
        _f.xVel() = xVel();
        _f.yVel() = yVel();
        _f.zVel() = zVel();
    }
};

int main()
{
    paz::UiDescriptor startMenu;
    startMenu.setTitle("PAZ Engine Test Program");
    startMenu.alignText(paz::UiAlignment::Center);
    startMenu.setLayout(paz::UiLayout::Horizontal);
    startMenu.addButton(paz::UiAction::Start, "Start");
    startMenu.addButton(paz::UiAction::ToggleFullscreen, "Toggle Fullscreen");
    startMenu.addButton(paz::UiAction::Quit, "Quit");

    paz::App::Init("scene.frag", "font.pbm", startMenu);
    Sphere50 = paz::Model("icosphere5.obj", 0, 0., Radius, "earth-day.bmp");
    Sphere10 = paz::Model("icosphere5.obj", 0, 0., 10., "moon.bmp");
    FancyBox = paz::Model("fancybox.obj");
    Player player;
    player.z() = Radius + 10.;
    World w;
    std::array<World1, 4> w1;
    w1[0].x() = 0.9*Radius;
    w1[1].x() = -0.9*Radius;
    w1[2].y() = 0.9*Radius;
    w1[3].y() = -0.9*Radius;
    World2 w2(0.);
    Npc npc0;
    npc0.x() = 10.;
    npc0.z() = Radius;
    Npc npc1;
    npc1.x() = 2.*Radius;
    npc1.z() = 10.;
    Npc npc2 = npc0;
    npc2.x() += 2.;
    Npc npc3;
    npc3 = npc0;
    npc3.x() += 4.;
    std::array<Fountain, 4> f;
    f[0].x() = -Radius;
    f[1].x() = Radius;
    f[2].y() = -Radius;
    f[3].y() = Radius;
    for(auto& n : f)
    {
        n.z() = Radius - 0.5;
        const paz::Vec up = (0.5*paz::Vec{{0., 0., 1.}}.cross(paz::Vec{{n.x(),
            n.y(), n.z()}}).normalized() + paz::Vec{{0., 0., 0.5}}).
            normalized();
        n.setDir(up(0), up(1), up(2));
    }
    w.lights().push_back({Radius + 5. + 10., 0., 0.});
    w.lights().push_back({-Radius - 5. - 10., 0., 0.});
    w.lights().push_back({0., Radius + 5. + 10., 0.});
    w.lights().push_back({0., -Radius - 5. - 10., 0.});
    w.lights().push_back({0., 0., Radius + 10.});
    w.lights().push_back({0., 0., -Radius - 10.});
    paz::App::AttachCamera(player.head());
    paz::App::AttachMic(player.head());
    paz::App::SetConsole(paz::ConsoleMode::CurrentFrame);
    paz::App::Run();
}
