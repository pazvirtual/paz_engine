#include "player.hpp"
#include <iomanip>

#define NO_FRICTION

static constexpr double CosMaxAngle = 0.7;

Player::Player()
{
    _head.collisionType() = paz::CollisionType::None;
    _head.gravityType() = paz::GravityType::None;
}

void Player::update()
{
    const double wAtt = std::sqrt(1. - xAtt()*xAtt() - yAtt()*yAtt() - zAtt()*
        zAtt());
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
        if(_collided && -nor.dot(gravDir) > CosMaxAngle)
        {
            reg = Regime::Grounded;
        }
    }
    if(_collided)
    {
        --_collided;
    }
    if(reg == Regime::Grounded && !_moving)
    {
        x() = _relX + _parent->x();
        y() = _relY + _parent->y();
        z() = _relZ + _parent->z();
        xVel() = _parent->xVel();
        yVel() = _parent->yVel();
        zVel() = _parent->zVel();
    }
const double r = std::sqrt(x()*x() + y()*y() + z()*z());
const double lat = std::asin(z()/r);
const double lon = std::atan2(y(), x());
paz::App::MsgStream() << std::fixed << std::setprecision(4) << std::setw(8) << r << " " << std::setw(9) << lat*180./paz::Pi << " " << std::setw(9) << lon*180./paz::Pi << " | " << std::setw(8) << std::sqrt(xVel()*xVel() + yVel()*yVel() + zVel()*zVel()) << std::endl;
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
    const double gravPitch = std::acos(std::max(0., std::min(1., baseForward.
        dot(forward))))*(forward.dot(gravDir) > 0. ? -1. : 1.);
paz::App::MsgStream() << std::fixed << std::setprecision(2) << std::setw(6) << (gravPitch + _pitch)*180./paz::Pi << std::endl;
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
            paz::Window::MousePos().second)*paz::App::PhysTime();
        _pitch = std::max(-0.45*paz::Pi, std::min(0.45*paz::Pi, _pitch +
            deltaPitch));
    }
    else
    {
        if(std::abs(_pitch) > 1e-6)
        {
            const double fac = 0.1;
            att = paz::qmult(paz::axis_angle(paz::Vec{{1, 0, 0}}, fac*_pitch),
                att);
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
            _mousePos -= 50.0/norm*paz::App::PhysTime()*_mousePos;
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

    _moving = false;

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
        _moving = u || v || w;
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
            xVel() += 12.*net(0)*paz::App::PhysTime();
            yVel() += 12.*net(1)*paz::App::PhysTime();
            zVel() += 12.*net(2)*paz::App::PhysTime();
        }
        else
        {
            if(paz::Window::KeyDown(paz::Key::A))
            {
                xVel() -= 12.*groundRight(0)*paz::App::PhysTime();
                yVel() -= 12.*groundRight(1)*paz::App::PhysTime();
                zVel() -= 12.*groundRight(2)*paz::App::PhysTime();
            }
            if(paz::Window::KeyDown(paz::Key::D))
            {
                xVel() += 12.*groundRight(0)*paz::App::PhysTime();
                yVel() += 12.*groundRight(1)*paz::App::PhysTime();
                zVel() += 12.*groundRight(2)*paz::App::PhysTime();
            }
            if(paz::Window::KeyDown(paz::Key::W))
            {
                xVel() += 12.*groundForward(0)*paz::App::PhysTime();
                yVel() += 12.*groundForward(1)*paz::App::PhysTime();
                zVel() += 12.*groundForward(2)*paz::App::PhysTime();
            }
            if(paz::Window::KeyDown(paz::Key::S))
            {
                xVel() -= 12.*groundForward(0)*paz::App::PhysTime();
                yVel() -= 12.*groundForward(1)*paz::App::PhysTime();
                zVel() -= 12.*groundForward(2)*paz::App::PhysTime();
            }
        }
    }
    else
    {
        if(paz::Window::GamepadActive())
        {
            const paz::Vec net = paz::Window::GamepadLeftStick().first*right -
                paz::Window::GamepadLeftStick().second*forward;
            xVel() += 12.*net(0)*paz::App::PhysTime();
            yVel() += 12.*net(1)*paz::App::PhysTime();
            zVel() += 12.*net(2)*paz::App::PhysTime();
        }
        else
        {
            if(paz::Window::KeyDown(paz::Key::A))
            {
                xVel() -= 12.*right(0)*paz::App::PhysTime();
                yVel() -= 12.*right(1)*paz::App::PhysTime();
                zVel() -= 12.*right(2)*paz::App::PhysTime();
            }
            if(paz::Window::KeyDown(paz::Key::D))
            {
                xVel() += 12.*right(0)*paz::App::PhysTime();
                yVel() += 12.*right(1)*paz::App::PhysTime();
                zVel() += 12.*right(2)*paz::App::PhysTime();
            }
            if(paz::Window::KeyDown(paz::Key::W))
            {
                xVel() += 12.*forward(0)*paz::App::PhysTime();
                yVel() += 12.*forward(1)*paz::App::PhysTime();
                zVel() += 12.*forward(2)*paz::App::PhysTime();
            }
            if(paz::Window::KeyDown(paz::Key::S))
            {
                xVel() -= 12.*forward(0)*paz::App::PhysTime();
                yVel() -= 12.*forward(1)*paz::App::PhysTime();
                zVel() -= 12.*forward(2)*paz::App::PhysTime();
            }
        }
    }
    double net = 0.;
    if(paz::Window::GamepadActive())
    {
        net = paz::Window::GamepadRightTrigger() - paz::Window::
            GamepadLeftTrigger();
    }
    else
    {
        if(paz::Window::KeyDown(paz::Key::LeftShift))
        {
            net += 1.;
        }
        if(paz::Window::KeyDown(paz::Key::LeftControl))
        {
            net -= 1.;
        }
    }
    if(net)
    {
        _moving = true;
        xVel() += 12.*up(0)*net*paz::App::PhysTime();
        yVel() += 12.*up(1)*net*paz::App::PhysTime();
        zVel() += 12.*up(2)*net*paz::App::PhysTime();
    }

    if(paz::Window::MousePressed(0) || paz::Window::GamepadPressed(paz::
        GamepadButton::A))
    {
        const paz::Vec dir = paz::to_mat(cameraAtt).row(1).trans();
        const paz::Vec pos{{_head.x(), _head.y(), _head.z()}};
        const paz::Vec vel{{xVel(), yVel(), zVel()}};
        _paintball.launch(pos, vel, dir);
    }
}

void Player::onCollide(const paz::Object& o, double xNor, double yNor, double
    zNor, double xB, double yB, double zB)
{
    const paz::Vec gravDir{{xDown(), yDown(), zDown()}};
    const paz::Vec nor{{xNor, yNor, zNor}};
    if(-nor.dot(gravDir) > CosMaxAngle)
    {
        _collided = 2;
        _parent.reset(o);
        _relX = x() - xB;
        _relY = y() - yB;
        _relZ = z() - zB;
        if(!_moving)
        {
            xVel() = _parent->xVel();
            yVel() = _parent->yVel();
            zVel() = _parent->zVel();
        }
    }
}

const paz::Object& Player::head() const
{
    return _head;
}
