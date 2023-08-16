#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

#define NO_FRICTION

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
    double _pitch = 0.;
    double _prevGravPitch = 0.;
    paz::Vec _mousePos = paz::Vec::Zero(2);
    paz::Object _head; // camera

public:
    Player()
    {
        _head.collisionType() = paz::CollisionType::None;
        /*collisionHeight()*/collisionRadius() = 1.5;
    }
    void update()
    {
        _head.x() = x();
        _head.y() = y();
        _head.z() = z();

        const double wAtt = std::sqrt(1. - xAtt()*xAtt() - yAtt()*yAtt() -
            zAtt()*zAtt());
        paz::Vec cameraAtt{{xAtt(), yAtt(), zAtt(), wAtt}};

        // Get basis vectors (rows of rotation matrix).
        paz::Mat cameraRot = paz::to_mat(cameraAtt);
        paz::Vec cameraForward = cameraRot.row(1).trans();
        const paz::Vec gravDir{{xDown(), yDown(), zDown()}};
        paz::Vec cameraRight = gravDir.cross(cameraForward).normalized();

        yAngRate() = 0.;

        // Want to keep `gravPitch + _pitch` (head pitch wrt gravity)
        // constant as much as possible when grounded.
        const paz::Vec baseForward = cameraRight.cross(gravDir).normalized();
        const double gravPitch = std::acos(std::max(0., std::min(1.,
            baseForward.dot(cameraForward))))*(cameraForward.dot(gravDir) > 0. ?
            -1. : 1.);
paz::App::MsgStream() << std::fixed << std::setprecision(2) << std::setw(6) << (gravPitch + _pitch)*180./M_PI << std::endl;
        if(grounded())
        {
            _mousePos = paz::Vec::Zero(2);
            xAngRate() = 0.;
            zAngRate() = -0.1*paz::Window::MousePos().first;
            cameraRot.setRow(0, cameraRight.trans());
            cameraRot.setRow(1, baseForward.trans());
            cameraRot.setRow(2, -gravDir.trans());
            paz::Vec cameraBaseAtt = paz::to_quat(cameraRot);
            if(cameraAtt.dot(cameraBaseAtt) < 0.)
            {
                cameraBaseAtt = -cameraBaseAtt;
            }
            if((cameraBaseAtt - cameraAtt).normSq() > 1e-6)
            {
                cameraBaseAtt = (0.9*cameraAtt + 0.1*cameraBaseAtt).
                    normalized();
                if(cameraBaseAtt(3) < 0.)
                {
                    cameraBaseAtt = -cameraBaseAtt;
                }
                xAtt() = cameraBaseAtt(0);
                yAtt() = cameraBaseAtt(1);
                zAtt() = cameraBaseAtt(2);
            }

            const double deltaGravPitch = gravPitch - _prevGravPitch;
            const double deltaPitch = -deltaGravPitch + 0.1*paz::Window::
                MousePos().second*paz::Window::FrameTime();
            _pitch = std::max(-0.45*M_PI, std::min(0.45*M_PI, _pitch +
                deltaPitch));
        }
        else
        {
            if(std::abs(_pitch) > 1e-6)
            {
                cameraAtt = paz::qmult(axis_angle(paz::Vec{{1, 0, 0}}, 0.1*
                    _pitch), cameraAtt);
                xAtt() = cameraAtt(0);
                yAtt() = cameraAtt(1);
                zAtt() = cameraAtt(2);
                cameraRot = paz::to_mat(cameraAtt);
                _pitch *= 0.9;
            }

            _mousePos(0) += paz::Window::MousePos().first;
            _mousePos(1) += paz::Window::MousePos().second;
            const double norm = _mousePos.norm();
            if(norm > 100.)
            {
                _mousePos *= 100./norm;
            }
            else if(norm > 0.1)
            {
                _mousePos -= 50.0/norm*paz::Window::FrameTime()*_mousePos;
            }
            else
            {
                _mousePos = paz::Vec::Zero(2);
            }
            xAngRate() = 0.006*_mousePos(1);
            zAngRate() = -0.006*_mousePos(0);
        }
        _prevGravPitch = gravPitch;
        cameraAtt = paz::qmult(axis_angle(paz::Vec{{1, 0, 0}}, _pitch + 0.5*
            M_PI), cameraAtt);
        if(cameraAtt(3) < 0.)
        {
            cameraAtt = -cameraAtt;
        }
        _head.xAtt() = cameraAtt(0);
        _head.yAtt() = cameraAtt(1);
        _head.zAtt() = cameraAtt(2);

        cameraRight = cameraRot.row(0).trans();
        cameraForward = cameraRot.row(1).trans();
        const paz::Vec cameraUp = cameraRot.row(2).trans();

        if(grounded())
        {
            double u = 0.;
            double v = 0.;
            double w = 0.;
            if(paz::Window::KeyDown(paz::Key::A))
            {
                u -= cameraRight(0);
                v -= cameraRight(1);
                w -= cameraRight(2);
            }
            if(paz::Window::KeyDown(paz::Key::D))
            {
                u += cameraRight(0);
                v += cameraRight(1);
                w += cameraRight(2);
            }
            if(paz::Window::KeyDown(paz::Key::S))
            {
                u -= cameraForward(0);
                v -= cameraForward(1);
                w -= cameraForward(2);
            }
            if(paz::Window::KeyDown(paz::Key::W))
            {
                u += cameraForward(0);
                v += cameraForward(1);
                w += cameraForward(2);
            }
            const double norm = std::sqrt(u*u + v*v + w*w);
            if(norm)
            {
                u *= 3./norm;
                v *= 3./norm;
                w *= 3./norm;
            }
#ifdef NO_FRICTION
            xVel() = u;
            yVel() = v;
            zVel() = w;
#else
            xVel() += u;
            yVel() += v;
            zVel() += w;
#endif
        }
        else
        {
            if(paz::Window::KeyDown(paz::Key::A))
            {
                xVel() -= 12.*cameraRight(0)*paz::Window::FrameTime();
                yVel() -= 12.*cameraRight(1)*paz::Window::FrameTime();
                zVel() -= 12.*cameraRight(2)*paz::Window::FrameTime();
            }
            if(paz::Window::KeyDown(paz::Key::D))
            {
                xVel() += 12.*cameraRight(0)*paz::Window::FrameTime();
                yVel() += 12.*cameraRight(1)*paz::Window::FrameTime();
                zVel() += 12.*cameraRight(2)*paz::Window::FrameTime();
            }
            if(paz::Window::KeyDown(paz::Key::W))
            {
                xVel() += 12.*cameraForward(0)*paz::Window::FrameTime();
                yVel() += 12.*cameraForward(1)*paz::Window::FrameTime();
                zVel() += 12.*cameraForward(2)*paz::Window::FrameTime();
            }
            if(paz::Window::KeyDown(paz::Key::S))
            {
                xVel() -= 12.*cameraForward(0)*paz::Window::FrameTime();
                yVel() -= 12.*cameraForward(1)*paz::Window::FrameTime();
                zVel() -= 12.*cameraForward(2)*paz::Window::FrameTime();
            }
        }
        if(paz::Window::KeyDown(paz::Key::LeftShift))
        {
            xVel() += 12.*cameraUp(0)*paz::Window::FrameTime();
            yVel() += 12.*cameraUp(1)*paz::Window::FrameTime();
            zVel() += 12.*cameraUp(2)*paz::Window::FrameTime();
        }
        if(paz::Window::KeyDown(paz::Key::LeftControl))
        {
            xVel() -= 12.*cameraUp(0)*paz::Window::FrameTime();
            yVel() -= 12.*cameraUp(1)*paz::Window::FrameTime();
            zVel() -= 12.*cameraUp(2)*paz::Window::FrameTime();
        }
    }
    const paz::Object& head() const
    {
        return _head;
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
        paz::Vec att = paz::to_quat(rot);
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
    paz::App::AttachCamera(player.head());
    paz::App::Run();
}
