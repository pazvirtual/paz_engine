#include "player.hpp"
#include <iomanip>

#define NO_FRICTION

static constexpr double CosMaxAngle = 0.7;

Player::Player()
{
    _head.collisionType() = paz::CollisionType::None;
    _head.gravityType() = paz::GravityType::None;
}

void Player::update(const paz::InputData& input)
{
    const double wAtt = std::sqrt(1. - xAtt()*xAtt() - yAtt()*yAtt() - zAtt()*
        zAtt());
    paz::Vec att{{xAtt(), yAtt(), zAtt(), wAtt}};

    // Get basis vectors (rows of rotation matrix).
    paz::Mat rot = paz::to_mat(att);
    paz::Vec forward = rot.row(0).trans();
    const paz::Vec gravDir{{xDown(), yDown(), zDown()}};
    paz::Vec left = forward.cross(gravDir).normalized();

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
        _collided = false;
    }
    if(reg == Regime::Grounded && !_moving)
    {
        const paz::Vec invParentAtt{{_parent->xAtt(), _parent->yAtt(), _parent->
            zAtt(), -std::sqrt(1. - _parent->xAtt()*_parent->xAtt() - _parent->
            yAtt()*_parent->yAtt() - _parent->zAtt()*_parent->zAtt())}};
        const auto rot = paz::to_mat(invParentAtt);
        const paz::Vec relPos = rot*_relPosPs;
        x() = relPos(0) + _parent->x();
        y() = relPos(1) + _parent->y();
        z() = relPos(2) + _parent->z();
        const paz::Vec parentAngRate{{_parent->xAngRate(), _parent->yAngRate(),
            _parent->zAngRate()}};
        const paz::Vec relVel = rot*parentAngRate.cross(_relPosPs);
        xVel() = relVel(0) + _parent->xVel();
        yVel() = relVel(1) + _parent->yVel();
        zVel() = relVel(2) + _parent->zVel();
    }
paz::App::SetReticule(reg == Regime::Grounded ? 3 : 0, reg == Regime::Grounded);
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
    const paz::Vec baseForward = gravDir.cross(left).normalized();
    const double gravPitch = std::acos(std::max(0., std::min(1., baseForward.
        dot(forward))))*(forward.dot(gravDir) > 0. ? -1. : 1.);
paz::App::MsgStream() << std::fixed << std::setprecision(2) << std::setw(6) << (gravPitch + _pitch)*180./paz::Pi << std::endl;
    if(reg != Regime::Floating)
    {
        _mousePos = paz::Vec::Zero(2);
        xAngRate() = 0.;
        zAngRate() = -0.1*(input.gamepadActive() ? 15.*input.
            gamepadRightStick().first : input.mousePos().first);
        rot.setRow(0, baseForward.trans());
        rot.setRow(1, left.trans());
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
        const double deltaPitch = -deltaGravPitch + 0.1*(input.gamepadActive() ?
            15.*-input.gamepadRightStick().second : input.mousePos().second)*
            input.timestep();
        _pitch = std::max(-0.45*paz::Pi, std::min(0.45*paz::Pi, _pitch +
            deltaPitch));
    }
    else
    {
        if(std::abs(_pitch) > 1e-6)
        {
            const double fac = 0.1;
            att = paz::qmult(paz::axis_angle(paz::Vec{{0, 1, 0}}, -fac*_pitch),
                att);
            xAtt() = att(0);
            yAtt() = att(1);
            zAtt() = att(2);
            rot = paz::to_mat(att);
            _pitch *= 1. - fac;
        }

        _mousePos(0) += (input.gamepadActive() ? 15.*input.gamepadRightStick().
            first : input.mousePos().first);
        _mousePos(1) += (input.gamepadActive() ? 15.*-input.gamepadRightStick().
            second : input.mousePos().second);
        const double norm = _mousePos.norm();
        if(norm > 100.)
        {
            _mousePos *= 100./norm;
        }
        else if(norm > 0.1)
        {
            _mousePos -= 50.0/norm*input.timestep()*_mousePos;
        }
        else
        {
            _mousePos = paz::Vec::Zero(2);
        }
        yAngRate() = -0.006*_mousePos(1);
        zAngRate() = -0.006*_mousePos(0);
    }
    _prevGravPitch = gravPitch;

    forward = rot.row(0).trans();
    left = rot.row(1).trans();
    const paz::Vec up = rot.row(2).trans();

    const double h = 1.5 - collisionRadius();
    _head.x() = x() + h*up(0);
    _head.y() = y() + h*up(1);
    _head.z() = z() + h*up(2);
    paz::Vec cameraAtt = paz::qmult(paz::axis_angle(paz::Vec{{0, 1, 0}},
        -_pitch), att);
    if(cameraAtt(3) < 0.)
    {
        cameraAtt = -cameraAtt;
    }
    _head.xAtt() = cameraAtt(0);
    _head.yAtt() = cameraAtt(1);
    _head.zAtt() = cameraAtt(2);

    _moving = false;

    const paz::Vec groundForward = left.cross(nor);
    const paz::Vec groundLeft = nor.cross(forward);
    if(reg == Regime::Grounded)
    {
        double u = 0.;
        double v = 0.;
        double w = 0.;
        if(input.gamepadActive())
        {
            const paz::Vec net = -input.gamepadLeftStick().first*groundLeft -
                input.gamepadLeftStick().second*groundForward;
            const double norm = std::max(1., net.norm());
            u = 3.*net(0)/norm;
            v = 3.*net(1)/norm;
            w = 3.*net(2)/norm;
        }
        else
        {
            if(input.keyDown(paz::Key::A))
            {
                u += groundLeft(0);
                v += groundLeft(1);
                w += groundLeft(2);
            }
            if(input.keyDown(paz::Key::D))
            {
                u -= groundLeft(0);
                v -= groundLeft(1);
                w -= groundLeft(2);
            }
            if(input.keyDown(paz::Key::S))
            {
                u -= groundForward(0);
                v -= groundForward(1);
                w -= groundForward(2);
            }
            if(input.keyDown(paz::Key::W))
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
        if(input.gamepadActive())
        {
            const paz::Vec net = -input.gamepadLeftStick().first*groundLeft -
                input.gamepadLeftStick().second*groundForward;
            xVel() += 12.*net(0)*input.timestep();
            yVel() += 12.*net(1)*input.timestep();
            zVel() += 12.*net(2)*input.timestep();
        }
        else
        {
            if(input.keyDown(paz::Key::A))
            {
                xVel() += 12.*groundLeft(0)*input.timestep();
                yVel() += 12.*groundLeft(1)*input.timestep();
                zVel() += 12.*groundLeft(2)*input.timestep();
            }
            if(input.keyDown(paz::Key::D))
            {
                xVel() -= 12.*groundLeft(0)*input.timestep();
                yVel() -= 12.*groundLeft(1)*input.timestep();
                zVel() -= 12.*groundLeft(2)*input.timestep();
            }
            if(input.keyDown(paz::Key::W))
            {
                xVel() += 12.*groundForward(0)*input.timestep();
                yVel() += 12.*groundForward(1)*input.timestep();
                zVel() += 12.*groundForward(2)*input.timestep();
            }
            if(input.keyDown(paz::Key::S))
            {
                xVel() -= 12.*groundForward(0)*input.timestep();
                yVel() -= 12.*groundForward(1)*input.timestep();
                zVel() -= 12.*groundForward(2)*input.timestep();
            }
        }
    }
    else
    {
        if(input.gamepadActive())
        {
            const paz::Vec net = -input.gamepadLeftStick().first*left - input.
                gamepadLeftStick().second*forward;
            xVel() += 12.*net(0)*input.timestep();
            yVel() += 12.*net(1)*input.timestep();
            zVel() += 12.*net(2)*input.timestep();
        }
        else
        {
            if(input.keyDown(paz::Key::A))
            {
                xVel() += 12.*left(0)*input.timestep();
                yVel() += 12.*left(1)*input.timestep();
                zVel() += 12.*left(2)*input.timestep();
            }
            if(input.keyDown(paz::Key::D))
            {
                xVel() -= 12.*left(0)*input.timestep();
                yVel() -= 12.*left(1)*input.timestep();
                zVel() -= 12.*left(2)*input.timestep();
            }
            if(input.keyDown(paz::Key::W))
            {
                xVel() += 12.*forward(0)*input.timestep();
                yVel() += 12.*forward(1)*input.timestep();
                zVel() += 12.*forward(2)*input.timestep();
            }
            if(input.keyDown(paz::Key::S))
            {
                xVel() -= 12.*forward(0)*input.timestep();
                yVel() -= 12.*forward(1)*input.timestep();
                zVel() -= 12.*forward(2)*input.timestep();
            }
        }
    }
    double net = 0.;
    if(input.gamepadActive())
    {
        net = input.gamepadRightTrigger() - input.gamepadLeftTrigger();
    }
    else
    {
        if(input.keyDown(paz::Key::LeftShift))
        {
            net += 1.;
        }
        if(input.keyDown(paz::Key::LeftControl))
        {
            net -= 1.;
        }
    }
    if(net)
    {
        _moving = true;
        xVel() += 12.*up(0)*net*input.timestep();
        yVel() += 12.*up(1)*net*input.timestep();
        zVel() += 12.*up(2)*net*input.timestep();
    }

    if(input.mousePressed(paz::MouseButton::Left) || input.gamepadPressed(paz::
        GamepadButton::A))
    {
        const paz::Vec dir = paz::to_mat(cameraAtt).row(0).trans();
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
        _collided = true;
        _parent.reset(o);
        const paz::Vec parentAtt{{_parent->xAtt(), _parent->yAtt(), _parent->
            zAtt(), std::sqrt(1. - _parent->xAtt()*_parent->xAtt() - _parent->
            yAtt()*_parent->yAtt() - _parent->zAtt()*_parent->zAtt())}};
        const paz::Vec relPos{{x() - xB, y() - yB, z() - zB}};
        const auto rot = paz::to_mat(parentAtt);
        _relPosPs = rot*relPos;
        if(!_moving)
        {
            const paz::Vec parentAngRate{{_parent->xAngRate(), _parent->
                yAngRate(), _parent->zAngRate()}};
            const paz::Vec relVel = rot.trans()*parentAngRate.cross(_relPosPs);
            xVel() = relVel(0) + _parent->xVel();
            yVel() = relVel(1) + _parent->yVel();
            zVel() = relVel(2) + _parent->zVel();
        }
    }
}

const paz::Object& Player::head() const
{
    return _head;
}
