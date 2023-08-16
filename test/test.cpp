#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

static paz::Model Sphere50;
static paz::Model Sphere;
static paz::Model Ground;
static paz::Model Tri;
static paz::Model PlatformModel;
static paz::Model Body;
static paz::Model Head;

static constexpr double Radius = 50.;

class Player : public paz::Object
{
public:
    Player()
    {
        /*collisionHeight()*/collisionRadius() = 1.5;
    }
};

class Npc : public paz::Object
{
    double _destYaw = 0.;
    double _walkTime = 0.;

public:
    paz::Object _head;
    Npc()
    {
        _head.collisionType() = paz::CollisionType::None;
        _head.model() = Head;
        model() = Body;
    }
    void update() override
    {
        const paz::Vec up = paz::Vec{{x(), y(), z()}}.normalized();
        const paz::Vec initialAtt{{xAtt(), yAtt(), zAtt(), std::sqrt(1. - xAtt()
            *xAtt() - yAtt()*yAtt() - zAtt()*zAtt())}};
        const paz::Vec initialRight = paz::to_mat(paz::qinv(initialAtt)).col(0);
        const paz::Vec forward = up.cross(initialRight).normalized();
        const paz::Vec right = forward.cross(up);
        paz::Mat rot(3);
        rot.setCol(0, right);
        rot.setCol(1, forward);
        rot.setCol(2, up);
        const double h = 1.5 - collisionRadius();
        _head.x() = x() + h*up(0);
        _head.y() = y() + h*up(1);
        _head.z() = z() + h*up(2);
        paz::Vec att = to_quat(rot);
        if(att(3) < 0.)
        {
            att = -att;
        }
        xAtt() = -att(0);
        yAtt() = -att(1);
        zAtt() = -att(2);
        _head.xAtt() = -att(0);
        _head.yAtt() = -att(1);
        _head.zAtt() = -att(2);
    }
    void onCollide(const Object&) override
    {
        const paz::Vec att{{xAtt(), yAtt(), zAtt(), std::sqrt(1. - xAtt()*xAtt()
            - yAtt()*yAtt() - zAtt()*zAtt())}};
        const paz::Mat rot = paz::to_mat(paz::qinv(att));
        const paz::Vec right = rot.col(0);
        const paz::Vec forward = rot.col(1);
        xVel() = forward(0);
        yVel() = forward(1);
        zVel() = forward(2);
        const paz::Vec up = rot.col(2);
        const paz::Vec east = paz::Vec{{0, 0, 1}}.cross(up).normalized();
        const paz::Vec north = up.cross(east);
        const double yaw = std::atan2(forward.dot(north), forward.dot(east));
        const double deltaYaw = paz::normalize_angle(_destYaw - yaw + M_PI) -
            M_PI;
        zAngRate() = deltaYaw;
        _walkTime += paz::Window::FrameTime();
        if(_walkTime > 10.)
        {
            _walkTime = 0.;
            _destYaw = paz::uniform(0., 2.*M_PI);
        }
    }
};

class Ball : public paz::Object
{
public:
    Ball()
    {
        model() = Sphere;
        collisionHeight() = 1.;
        collisionRadius() = 1.;
    }
    void onInteract(const Object& o) override
    {
        double xDir = x() - o.x();
        double yDir = y() - o.y();
        double zDir = z() - o.z();
        const double norm = std::sqrt(xDir*xDir + yDir*yDir + zDir*zDir);
        xDir /= norm;
        yDir /= norm;
        zDir /= norm;
        xVel() = xDir;
        yVel() = yDir;
        zVel() = zDir;
    }
};

class World : public paz::Object
{
public:
    World()
    {
        model() = Sphere50;
        collisionType() = paz::CollisionType::World;
    }
};

class World1 : public paz::Object
{
public:
    World1()
    {
        model() = Ground;
        z() = Radius;
        collisionType() = paz::CollisionType::World;
    }
};

class Platform : public paz::Object
{
    double _timeAtTop = 0.1;

public:
    Platform()
    {
        model() = PlatformModel;
        z() = Radius;
        collisionType() = paz::CollisionType::World;
    }
    void update() override
    {
        if(_timeAtTop < 0.1)
        {
            if(z() < Radius + 10.)
            {
                zVel() = std::min(1., 5.2 - std::abs(z() - Radius - 5.));
            }
            else
            {
                z() = Radius + 10.;
                zVel() = 0.;
                _timeAtTop += paz::Window::FrameTime();
            }
        }
        else
        {
            if(z() > Radius)
            {
                zVel() = -std::min(1., 5.2 - std::abs(z() - Radius - 5.));
            }
            else
            {
                z() = Radius;
                zVel() = 0.;
            }
        }
    }
    void onCollide(const Object& o) override
    {
        if(o.z() > z())
        {
            _timeAtTop = 0.;
        }
    }
    void onInteract(const Object&) override
    {
        _timeAtTop = 0.;
    }
};

int main()
{
    paz::App::Init("scene.frag", "font.pbm");
//    Ground = paz::Model("plane.obj");
//    Tri = paz::Model("tri.obj");
//    PlatformModel = paz::Model("platform.obj");
    Sphere50 = paz::Model("sphere50.obj");
//    Sphere = paz::Model("unitsphere.obj");
    Body = paz::Model("persontest.obj", 0, -0.2);
    Head = paz::Model("persontest.obj", 1, -0.1);
    Player player;
    player.z() = Radius + 10.;
//    Ball b;
//    b.z() = Radius;
//    b.x() = 5.;
    std::vector<World> w(5);
    w[1].x() = Radius;
    w[2].x() = -Radius;
    w[3].y() = Radius;
    w[4].y() = -Radius;
//    std::vector<World1> w1(10);
//    for(std::size_t i = 0; i < w1.size(); ++i)
//    {
//        w1[i].y() = 20.*i;
//        w1[i].z() = Radius - 2.;
//    }
//    Platform p;
    Npc npc0;
    npc0.x() = 10.;
    npc0.z() = Radius;
    Npc npc1;
    npc1.x() = 2.*Radius;
    npc1.z() = 10.;
    paz::App::AttachCamera(player);
    paz::App::Run();
}
