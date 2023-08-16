#include "PAZ_Engine"
#include "PAZ_Math"
#include <iomanip>

#define NO_FRICTION

static paz::Model Sphere50;
static paz::Model Sphere10;
static paz::Model Body;
static paz::Model Head;
static paz::Model PaintballModel;
static paz::Model FancyBox;

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
        xVel() = vel(0) + 15.*dir(0);
        yVel() = vel(1) + 15.*dir(1);
        zVel() = vel(2) + 15.*dir(2);
        collisionRadius() = 0.05;
        model() = PaintballModel;
    }
    void onCollide(const Object&) override
    {
        xVel() = 0.;
        yVel() = 0.;
        zVel() = 0.;
        collisionType() = paz::CollisionType::None;
        gravityType() = paz::GravityType::None;
    }
};

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
            _timer += paz::Window::FrameTime();
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
    }
    void update() final
    {
        _timer += paz::Window::FrameTime();
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
                const double offset = 1.3 + i*0.1*10./NumPerLaunch;
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

class Player : public paz::Object
{
    enum class Regime
    {
        Grounded, Low, Floating
    };

    static constexpr double LowAltitude = 10.;

    double _pitch = 0.;
    double _prevGravPitch = 0.;
    paz::Vec _mousePos = paz::Vec::Zero(2);
    paz::Object _head; // camera

    std::vector<Paintball> _paintballs; //TEMP - should these really be owned by `this` ?

public:
    Player()
    {
        _head.collisionType() = paz::CollisionType::None;
        _head.gravityType() = paz::GravityType::None;
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
        Regime reg = Regime::Floating;
        double alt;
        paz::Vec nor;
        paz::Vec surfVel;
        computeAltitude(alt, nor, surfVel);
        if(alt < LowAltitude)
        {
            reg = Regime::Low;
            if(alt < 0.1)//TEMP - should be function of maxangle
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
paz::App::MsgStream() << alt << " | " << nor.trans() << std::endl;

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
                const double fac = std::max(0., std::min(1., 1. - std::sqrt(alt/
                    LowAltitude)));
                baseAtt = ((1. - fac)*att + fac*baseAtt).normalized();
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
                const double fac = 0.1;
                att = paz::qmult(paz::axis_angle(paz::Vec{{1, 0, 0}}, fac*
                    _pitch), att);
                xAtt() = att(0);
                yAtt() = att(1);
                zAtt() = att(2);
                rot = paz::to_mat(att);
                _pitch *= 1. - fac;
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
            xVel() = surfVel(0) + u;
            yVel() = surfVel(1) + v;
            zVel() = surfVel(2) + w;
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
            _paintballs.emplace_back(pos, vel, dir);
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
        _head.gravityType() = paz::GravityType::None;
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
        _angle = paz::normalize_angle(_angle + AngRate*paz::Window::FrameTime());
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
    paz::App::Init("scene.frag", "font.pbm");
    Sphere50 = paz::Model("sphere50.obj");
    Sphere10 = paz::Model("sphere50.obj", 0, 0., 10./Radius);
    Body = paz::Model("persontest.obj", 0, -0.2);
    Head = paz::Model("persontest.obj", 1, -0.1);
    PaintballModel = paz::Model("unitsphere.obj", 0, 0., 0.1);
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
    paz::App::AttachCamera(player.head());
    paz::App::Run();
}
