#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

#define NO_FRICTION

static paz::Model Sphere50;
static paz::Model Body;
static paz::Model Head;
static paz::Model PaintballModel;

static constexpr double Radius = 50.;

class Paintball : public paz::Object
{
public:
    Paintball(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec& dir) :
        paz::Object()
    {
        x() = pos(0);
        y() = pos(1);
        z() = pos(2);
        xVel() = vel(0) + 10.*dir(0);
        yVel() = vel(1) + 10.*dir(1);
        zVel() = vel(2) + 10.*dir(2);
        collisionRadius() = 0.01;
        model() = PaintballModel;
    }
    void update() final
    {
        if(altitude() < 0.01)
        {
            xVel() = 0.;
            yVel() = 0.;
            zVel() = 0.;
            collisionType() = paz::CollisionType::None;
        }
    }
};

class Player : public paz::Object
{
    enum class Regime
    {
        Grounded, Low, Floating
    };

    double _pitch = 0.;
    double _prevGravPitch = 0.;
    paz::Vec _mousePos = paz::Vec::Zero(2);
    paz::Object _head; // camera

    std::vector<std::shared_ptr<Paintball>> _paintballs; //TEMP - should these really be owned by `this` ?

public:
    Player()
    {
        _head.collisionType() = paz::CollisionType::None;
    }
    void update() final
    {
        const double wAtt = std::sqrt(1. - xAtt()*xAtt() - yAtt()*yAtt() -
            zAtt()*zAtt());
        paz::Vec att{{xAtt(), yAtt(), zAtt(), wAtt}};

        // Get basis vectors (rows of rotation matrix).
        paz::Mat rot = paz::to_mat(att);
        paz::Vec forward = rot.row(1).trans();
        const paz::Vec gravDir{{xDown(), yDown(), zDown()}};
        paz::Vec right = gravDir.cross(forward).normalized();

        // Check altitude regime.
        const paz::Vec nor{{localNorX(), localNorY(), localNorZ()}};
        Regime reg = Regime::Floating;
        if(altitude() < 1.)
        {
            reg = Regime::Low;
            if(altitude() < 0.01 && nor.dot(gravDir) < -0.7)
            {
                reg = Regime::Grounded;
            }
        }
const double r = std::sqrt(x()*x() + y()*y() + z()*z());
const double lat = std::asin(z()/r);
const double lon = std::atan2(y(), x());
paz::App::MsgStream() << std::fixed << std::setprecision(4) << std::setw(8) << r << " " << std::setw(9) << lat*180./M_PI << " " << std::setw(9) << lon*180./M_PI << " | " << std::setw(8) << std::sqrt(xVel()*xVel() + yVel()*yVel() + zVel()*zVel()) << std::endl;
switch(reg)
{
case Regime::Grounded: paz::App::MsgStream() << "Grounded" << std::endl; break;
case Regime::Low: paz::App::MsgStream() << "Low" << std::endl; break;
case Regime::Floating: paz::App::MsgStream() << "Floating" << std::endl; break;
}

        // Kill all roll.
        yAngRate() = 0.;

        // Want to keep `gravPitch + _pitch` (head pitch wrt gravity)
        // constant as much as possible when grounded.
        const paz::Vec baseForward = right.cross(gravDir).normalized();
        const double gravPitch = std::acos(std::max(0., std::min(1.,
            baseForward.dot(forward))))*(forward.dot(gravDir) > 0. ? -1. : 1.);
paz::App::MsgStream() << std::fixed << std::setprecision(2) << std::setw(6) << (gravPitch + _pitch)*180./M_PI << std::endl;
        if(reg != Regime::Floating)
        {
            _mousePos = paz::Vec::Zero(2);
            xAngRate() = 0.;
            zAngRate() = -0.1*(paz::Window::GamepadActive() ? 15.*paz::Window::
                GamepadRightStick().first : paz::Window::MousePos().first);
            rot.setRow(0, right.trans());
            rot.setRow(1, baseForward.trans());
            rot.setRow(2, -gravDir.trans());
            paz::Vec baseAtt = paz::to_quat(rot);
            if(att.dot(baseAtt) < 0.)
            {
                baseAtt = -baseAtt;
            }
            if((baseAtt - att).normSq() > 1e-6)
            {
                baseAtt = (0.9*att + 0.1*baseAtt).normalized();
                if(baseAtt(3) < 0.)
                {
                    baseAtt = -baseAtt;
                }
                xAtt() = baseAtt(0);
                yAtt() = baseAtt(1);
                zAtt() = baseAtt(2);
                rot = paz::to_mat(baseAtt);
            }

            const double deltaGravPitch = gravPitch - _prevGravPitch;
            const double deltaPitch = -deltaGravPitch + 0.1*(paz::Window::
                GamepadActive() ? 15.*-paz::Window::GamepadRightStick().second :
                paz::Window::MousePos().second)*paz::Window::FrameTime();
            _pitch = std::max(-0.45*M_PI, std::min(0.45*M_PI, _pitch +
                deltaPitch));
        }
        else
        {
            if(std::abs(_pitch) > 1e-6)
            {
                att = paz::qmult(paz::axis_angle(paz::Vec{{1, 0, 0}}, 0.1*
                    _pitch), att);
                xAtt() = att(0);
                yAtt() = att(1);
                zAtt() = att(2);
                rot = paz::to_mat(att);
                _pitch *= 0.9;
            }

            _mousePos(0) += (paz::Window::GamepadActive() ? 15.*paz::Window::
                GamepadRightStick().first : paz::Window::MousePos().first);
            _mousePos(1) += (paz::Window::GamepadActive() ? 15.*-paz::Window::
                GamepadRightStick().second : paz::Window::MousePos().second);
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

        right = rot.row(0).trans();
        forward = rot.row(1).trans();
        const paz::Vec up = rot.row(2).trans();

        const double h = 1.5 - collisionRadius();
        _head.x() = x() + h*up(0);
        _head.y() = y() + h*up(1);
        _head.z() = z() + h*up(2);
        paz::Vec cameraAtt = paz::qmult(paz::axis_angle(paz::Vec{{1, 0, 0}},
            _pitch), att);
        if(cameraAtt(3) < 0.)
        {
            cameraAtt = -cameraAtt;
        }
        _head.xAtt() = cameraAtt(0);
        _head.yAtt() = cameraAtt(1);
        _head.zAtt() = cameraAtt(2);

        const paz::Vec groundRight = forward.cross(nor);
        const paz::Vec groundForward = nor.cross(right);
        if(reg == Regime::Grounded)
        {
            double u = 0.;
            double v = 0.;
            double w = 0.;
            if(paz::Window::GamepadActive())
            {
                const paz::Vec net = paz::Window::GamepadLeftStick().first*
                    groundRight - paz::Window::GamepadLeftStick().second*
                    groundForward;
                const double norm = std::max(1., net.norm());
                u = 3.*net(0)/norm;
                v = 3.*net(1)/norm;
                w = 3.*net(2)/norm;
            }
            else
            {
                if(paz::Window::KeyDown(paz::Key::A))
                {
                    u -= groundRight(0);
                    v -= groundRight(1);
                    w -= groundRight(2);
                }
                if(paz::Window::KeyDown(paz::Key::D))
                {
                    u += groundRight(0);
                    v += groundRight(1);
                    w += groundRight(2);
                }
                if(paz::Window::KeyDown(paz::Key::S))
                {
                    u -= groundForward(0);
                    v -= groundForward(1);
                    w -= groundForward(2);
                }
                if(paz::Window::KeyDown(paz::Key::W))
                {
                    u += groundForward(0);
                    v += groundForward(1);
                    w += groundForward(2);
                }
                const double norm = std::sqrt(u*u + v*v + w*w);
                if(norm)
                {
                    u *= 3./norm;
                    v *= 3./norm;
                    w *= 3./norm;
                }
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
        else if(reg == Regime::Low)
        {
            if(paz::Window::GamepadActive())
            {
                const paz::Vec net = paz::Window::GamepadLeftStick().first*
                    groundRight - paz::Window::GamepadLeftStick().second*
                    groundForward;
                xVel() += 12.*net(0)*paz::Window::FrameTime();
                yVel() += 12.*net(1)*paz::Window::FrameTime();
                zVel() += 12.*net(2)*paz::Window::FrameTime();
            }
            else
            {
                if(paz::Window::KeyDown(paz::Key::A))
                {
                    xVel() -= 12.*groundRight(0)*paz::Window::FrameTime();
                    yVel() -= 12.*groundRight(1)*paz::Window::FrameTime();
                    zVel() -= 12.*groundRight(2)*paz::Window::FrameTime();
                }
                if(paz::Window::KeyDown(paz::Key::D))
                {
                    xVel() += 12.*groundRight(0)*paz::Window::FrameTime();
                    yVel() += 12.*groundRight(1)*paz::Window::FrameTime();
                    zVel() += 12.*groundRight(2)*paz::Window::FrameTime();
                }
                if(paz::Window::KeyDown(paz::Key::W))
                {
                    xVel() += 12.*groundForward(0)*paz::Window::FrameTime();
                    yVel() += 12.*groundForward(1)*paz::Window::FrameTime();
                    zVel() += 12.*groundForward(2)*paz::Window::FrameTime();
                }
                if(paz::Window::KeyDown(paz::Key::S))
                {
                    xVel() -= 12.*groundForward(0)*paz::Window::FrameTime();
                    yVel() -= 12.*groundForward(1)*paz::Window::FrameTime();
                    zVel() -= 12.*groundForward(2)*paz::Window::FrameTime();
                }
            }
        }
        else
        {
            if(paz::Window::GamepadActive())
            {
                const paz::Vec net = paz::Window::GamepadLeftStick().first*
                    right - paz::Window::GamepadLeftStick().second*forward;
                xVel() += 12.*net(0)*paz::Window::FrameTime();
                yVel() += 12.*net(1)*paz::Window::FrameTime();
                zVel() += 12.*net(2)*paz::Window::FrameTime();
            }
            else
            {
                if(paz::Window::KeyDown(paz::Key::A))
                {
                    xVel() -= 12.*right(0)*paz::Window::FrameTime();
                    yVel() -= 12.*right(1)*paz::Window::FrameTime();
                    zVel() -= 12.*right(2)*paz::Window::FrameTime();
                }
                if(paz::Window::KeyDown(paz::Key::D))
                {
                    xVel() += 12.*right(0)*paz::Window::FrameTime();
                    yVel() += 12.*right(1)*paz::Window::FrameTime();
                    zVel() += 12.*right(2)*paz::Window::FrameTime();
                }
                if(paz::Window::KeyDown(paz::Key::W))
                {
                    xVel() += 12.*forward(0)*paz::Window::FrameTime();
                    yVel() += 12.*forward(1)*paz::Window::FrameTime();
                    zVel() += 12.*forward(2)*paz::Window::FrameTime();
                }
                if(paz::Window::KeyDown(paz::Key::S))
                {
                    xVel() -= 12.*forward(0)*paz::Window::FrameTime();
                    yVel() -= 12.*forward(1)*paz::Window::FrameTime();
                    zVel() -= 12.*forward(2)*paz::Window::FrameTime();
                }
            }
        }
        if(paz::Window::GamepadActive())
        {
            const double net = paz::Window::GamepadRightTrigger() - paz::
                Window::GamepadLeftTrigger();
            xVel() += 12.*up(0)*net*paz::Window::FrameTime();
            yVel() += 12.*up(1)*net*paz::Window::FrameTime();
            zVel() += 12.*up(2)*net*paz::Window::FrameTime();
        }
        else
        {
            if(paz::Window::KeyDown(paz::Key::LeftShift))
            {
                xVel() += 12.*up(0)*paz::Window::FrameTime();
                yVel() += 12.*up(1)*paz::Window::FrameTime();
                zVel() += 12.*up(2)*paz::Window::FrameTime();
            }
            if(paz::Window::KeyDown(paz::Key::LeftControl))
            {
                xVel() -= 12.*up(0)*paz::Window::FrameTime();
                yVel() -= 12.*up(1)*paz::Window::FrameTime();
                zVel() -= 12.*up(2)*paz::Window::FrameTime();
            }
        }

        if(paz::Window::MousePressed(0) || paz::Window::GamepadPressed(paz::
            GamepadButton::A))
        {
            const paz::Vec dir = paz::to_mat(cameraAtt).row(1).trans();
            const paz::Vec pos{{_head.x(), _head.y(), _head.z()}};
            const paz::Vec vel{{xVel(), yVel(), zVel()}};
            _paintballs.push_back(std::make_shared<Paintball>(pos, vel, dir));
        }
    }
    const paz::Object& head() const
    {
        return _head;
    }
};

class Npc : public paz::Object
{
    double _destYaw = paz::uniform(0., 2.*M_PI);
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
        const paz::Vec up = -paz::Vec{{xDown(), yDown(), zDown()}}.normalized();
        const paz::Vec initialAtt{{xAtt(), yAtt(), zAtt(), std::sqrt(1. - xAtt()
            *xAtt() - yAtt()*yAtt() - zAtt()*zAtt())}};
        const paz::Vec initialRight = paz::to_mat(initialAtt).row(0).trans();
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

class World : public paz::Object
{
public:
    World()
    {
        model() = Sphere50;
        collisionType() = paz::CollisionType::World;
    }
};

int main()
{
    paz::App::Init("scene.frag", "font.pbm");
    Sphere50 = paz::Model("sphere50.obj");
    Body = paz::Model("persontest.obj", 0, -0.2);
    Head = paz::Model("persontest.obj", 1, -0.1);
    PaintballModel = paz::Model("unitsphere.obj", 0, 0., 0.05);
    Player player;
    player.z() = Radius + 10.;
    std::vector<World> w(5);
    w[1].x() = Radius;
    w[2].x() = -Radius;
    w[3].y() = Radius;
    w[4].y() = -Radius;
    Npc npc0;
    npc0.x() = 10.;
    npc0.z() = Radius;
    Npc npc1;
    npc1.x() = 2.*Radius;
    npc1.z() = 10.;
    Npc npc2;
    npc2 = npc0;
    npc2.x() += 2.;
    paz::App::AttachCamera(player.head());
    paz::App::Run();
}
