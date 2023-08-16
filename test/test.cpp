#include "player.hpp"
#include "npc.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

static paz::Model Sphere50;
static paz::Model Sphere10;
static paz::Model FancyBox;

static constexpr double Radius = 50.;

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

public:
    World2(double initialAngle) : paz::Object(), _angle(initialAngle)
    {
        model() = Sphere10;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None; //TEMP
        stdGravParam() = 9.81*10.*10.;
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
    }
};

int main()
{
    paz::App::Init("PAZ Engine Test Program");
    Sphere50 = paz::Model("icosphere5.obj", 0, 0., Radius, "earth-day.bmp");
    Sphere10 = paz::Model("icosphere5.obj", 0, 0., 10., "moon.bmp");
    FancyBox = paz::Model("fancybox.obj");
    Player player;
    player.z() = Radius + 10.;
    World w;
    w.transp().push_back({10., 5., Radius, 0., 0., Radius + 10., -10., -5.,
        Radius});
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
    w.lights().push_back({ Radius + 5. + 10., 0., 0., 1., 0., 0., 0.1});
    w.lights().push_back({-Radius - 5. - 10., 0., 0., 1., 1., 0., 0.1});
    w.lights().push_back({0.,  Radius + 5. + 10., 0., 1., 1., 1., 0.1});
    w.lights().push_back({0., -Radius - 5. - 10., 0., 0., 1., 1., 0.1});
    w.lights().push_back({0., 0.,  Radius + 10.,      0., 1., 0., 0.1});
    w.lights().push_back({0., 0., -Radius - 10.,      0., 0., 1., 0.1});
    for(auto& n : w.lights())
    {
        const float lum = 0.2126*n[3] + 0.7152*n[4] + 0.0722*n[5];
        n[3] /= lum;
        n[4] /= lum;
        n[5] /= lum;
    }
    paz::App::SetSun(paz::Vec{{1., 0., 0.}}, paz::Vec{{1., 1., 0.9}});
    paz::App::AttachCamera(player.head());
    paz::App::AttachMic(player.head());
    paz::App::SetConsole(paz::ConsoleMode::CurrentFrame);
    paz::App::Run();
}
