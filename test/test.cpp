#include "PAZ_Engine"
#include <cmath>

static paz::Model Sphere;
static paz::Model Ground;
static paz::Model PlatformModel;

class Player : public paz::Object
{
public:
    Player()
    {
        height() = 1.5;
    }
};

class Ball : public paz::Object
{
public:
    Ball()
    {
        model() = Sphere;
        height() = 1.;
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
        model() = Ground; // wraps {VertexBuffer, IndexBuffer, vec<tri>}
        collisionType() = paz::CollisionType::World; // is collided with
    }
};

class Platform : public paz::Object
{
    double _timeAtTop = 0.1;

public:
    Platform()
    {
        model() = PlatformModel;
        collisionType() = paz::CollisionType::World;
    }
    void update() override
    {
        if(_timeAtTop < 0.1)
        {
            if(z() < 3.)
            {
                zVel() = 1.;
            }
            else
            {
                z() = 3.;
                zVel() = 0.;
                _timeAtTop += paz::Window::FrameTime();
            }
        }
        else
        {
            if(z() > 0.)
            {
                zVel() = -1.;
            }
            else
            {
                z() = 0.;
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

int main(int, char** argv)
{
    const std::string appDir = paz::split_path(argv[0])[0];
    paz::App::Init(appDir + "/scene.frag");
    Ground = paz::Model(appDir + "/plane.obj");
    PlatformModel = paz::Model(appDir + "/platform.obj");
    Sphere = paz::Model(appDir + "/unitsphere.obj");
    Player player;
    Ball b;
    World w;
    Platform p;
    p.x() = 2.;
    paz::App::AttachCamera(player);
    paz::App::Run();
}
