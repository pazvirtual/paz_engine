#include "player.hpp"
#include "npc.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

static constexpr double Radius = 50.;

static paz::Model _sphere50;
static paz::Model _sphere10;
static paz::Model _fancyBox;

class Droplet : public Paintball
{
    bool _stuck = false;
    double _timer = 0.;

public:
    Droplet(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec& dir) :
        Paintball()
    {
        lights().clear(); //TEMP
        launch(pos, vel, dir);
    }

    void update(const paz::InputData& input) override
    {
        if(_stuck)
        {
            _timer += input.timestep();
        }
        Paintball::update(input);
    }

    void onCollide(const Object& o, double xNor, double yNor, double zNor,
        double xB, double yB, double zB) override
    {
        _stuck = true;
        Paintball::onCollide(o, xNor, yNor, zNor, xB, yB, zB);
    }

    bool done() const
    {
        return _timer > 1.;
    }
};

class Fountain : public paz::Object
{
    double _timer = 0.;
    std::vector<Droplet> _droplets;
    static constexpr int NumPerLaunch = 1;

public:
    Fountain()
    {
        model() = _fancyBox;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None;
        lights().push_back({0., 0., 1.3, 0.5, 2., 2., 0.1});
    }

    void update(const paz::InputData& input) final
    {
        _timer += input.timestep();
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
            if(_droplets[i].done())
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
        model() = _sphere50;
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
        model() = _sphere10;
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

    void updateInternal()
    {
        x() = std::cos(_angle)*2.*Radius;
        y() = std::sin(_angle)*2.*Radius;
        z() = 0.;
        xVel() = AngRate*-std::sin(_angle)*2.*Radius;
        yVel() = AngRate*std::cos(_angle)*2.*Radius;
        zVel() = 0.;
        _f.x() = x();
        _f.y() = y();
        _f.z() = z() + 10.;
        _f.xVel() = xVel();
        _f.yVel() = yVel();
        _f.zVel() = zVel();
        auto att = paz::axis_angle(paz::Vec{{0., 0., 1.}}, _angle);
        if(att(3) < 0.)
        {
            att = -att;
        }
        xAtt() = att(0);
        yAtt() = att(1);
        zAtt() = att(2);
        zAngRate() = AngRate;
        _f.xAtt() = xAtt();
        _f.yAtt() = yAtt();
        _f.zAtt() = zAtt();
        _f.zAngRate() = zAngRate();
    }

public:
    World2(double initialAngle) : paz::Object(), _angle(initialAngle)
    {
        model() = _sphere10;
        collisionType() = paz::CollisionType::World;
        gravityType() = paz::GravityType::None; //TEMP
        stdGravParam() = 9.81*10.*10.;
        updateInternal();
    }

    void update(const paz::InputData& input) override
    {
        _angle = paz::normalize_angle(_angle + AngRate*input.timestep());
        updateInternal();
    }
};

int main()
{
    paz::App::Init("PAZ Engine Test Program");
    _sphere50 = paz::Model("icosphere5.pazmodel", 0, 0., Radius,
        "earth-day.bmp", {}, {{10., 5., Radius, 0., 0., Radius + 10., -10., -5.,
        Radius}});
    _sphere10 = paz::Model("icosphere5.pazmodel", 0, 0., 10., "moon.bmp");
    _fancyBox = paz::Model("fancybox.pazmodel");
    Player player;
    player.y() = 1.;
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
    std::swap(npc0, npc0);
    std::swap(npc0, npc1);
    npc0.x() = w2.x();
    npc0.y() = w2.y();
    npc0.z() = w2.z() + 10.5;
    npc0.xVel() = w2.xVel();
    npc0.yVel() = w2.yVel();
    npc0.zVel() = w2.zVel();
    Npc npc2 = npc1;
    npc2.x() += 2.;
    Npc npc3;
    npc3 = npc2;
    npc3.x() += 2.;
    Npc npc4 = std::move(npc3);
    npc3 = npc4;
    npc3.x() += 2.;
    std::vector<Npc> moreNpcs(10);
    moreNpcs[0] = npc3;
    for(std::size_t i = 0; i + 1 < moreNpcs.size(); ++i)
    {
        moreNpcs[i + 1] = moreNpcs[i];
        moreNpcs[i + 1].x() += 2.;
    }
    Npc npc5 = npc0;
    npc5.z() = -npc5.z();
    npc0.setName("NPC 0");
    npc1.setName("NPC 1");
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
    paz::App::SetConsole(paz::ConsoleMode::LatestStep);
    paz::App::Run();
}
