#include "npc.hpp"
#include "PAZ_Math"

static constexpr double CosMaxAngle = 0.7;

static const paz::Model Body("persontest_body.pazmodel", 0, -0.2);
static const paz::Model Head("persontest_head.pazmodel", 0, -0.1);

Npc::Npc() : _destYaw(paz::uniform(0., paz::TwoPi)), _walkTime(0.)
{
    _head.collisionType() = paz::CollisionType::None;
    _head.gravityType() = paz::GravityType::None;
    _head.model() = Head;
    model() = Body;
}

void Npc::update(const paz::InputData& input)
{
    if(!_name.empty() && paz::uniform() < 5e-3)
    {
        paz::App::PushDialog("`" + _name + "`\nAsdfj asdf asdf asdf." + (paz::
            uniform() < 0.5 ? "\nBLAHBjLAH blah." : ""), 1.);
    }
    _walkTime += input.timestep();
    const paz::Vec up = -paz::Vec{{xDown(), yDown(), zDown()}};
    const paz::Vec initialAtt{{xAtt(), yAtt(), zAtt(), std::sqrt(1. - xAtt()
        *xAtt() - yAtt()*yAtt() - zAtt()*zAtt())}};
    const paz::Vec initialForward = paz::to_mat(initialAtt).row(0).trans();
    const paz::Vec left = up.cross(initialForward).normalized();
    const paz::Vec forward = left.cross(up);
    if(_collided)
    {
        double alt;
        paz::Vec nor;
        paz::Vec surfVel;
        computeAltitude(alt, nor, surfVel);
        if(nor.dot(up) > CosMaxAngle)
        {
            const paz::Vec invParentAtt{{_parent->xAtt(), _parent->yAtt(),
                _parent->zAtt(), -std::sqrt(1. - _parent->xAtt()*_parent->xAtt()
                - _parent->yAtt()*_parent->yAtt() - _parent->zAtt()*_parent->
                zAtt())}};
            const paz::Vec parentAngRate{{_parent->xAngRate(), _parent->
                yAngRate(), _parent->zAngRate()}};
            const paz::Vec relPos{{x() - _parent->x(), y() - _parent->y(), z() -
                _parent->z()}};
            const paz::Vec relVel = static_cast<paz::Vec>(paz::to_mat(
                invParentAtt)*parentAngRate).cross(relPos);
            xVel() = relVel(0) + _parent->xVel() + forward(0);
            yVel() = relVel(1) + _parent->yVel() + forward(1);
            zVel() = relVel(2) + _parent->zVel() + forward(2);
        }
        _collided = false;
    }
    paz::Mat rot(3);
    rot.setCol(0, forward);
    rot.setCol(1, left);
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

void Npc::onCollide(const paz::Object& o, double xNor, double yNor, double
    zNor, double xB, double yB, double zB)
{
    const paz::Vec gravDir{{xDown(), yDown(), zDown()}};
    const paz::Vec nor{{xNor, yNor, zNor}};
    if(-nor.dot(gravDir) > CosMaxAngle)
    {
        _collided = true;
        _parent.reset(o);
        const paz::Vec att{{xAtt(), yAtt(), zAtt(), std::sqrt(1. - xAtt()*xAtt()
            - yAtt()*yAtt() - zAtt()*zAtt())}};
        const paz::Mat rot = paz::to_mat(paz::qinv(att));
        const paz::Vec forward = rot.col(0);
        const paz::Vec invParentAtt{{_parent->xAtt(), _parent->yAtt(), _parent->
            zAtt(), -std::sqrt(1. - _parent->xAtt()*_parent->xAtt() - _parent->
            yAtt()*_parent->yAtt() - _parent->zAtt()*_parent->zAtt())}};
        const paz::Vec parentAngRate{{_parent->xAngRate(), _parent->yAngRate(),
            _parent->zAngRate()}};
        const paz::Vec relPos{{x() - xB, y() - yB, z() - zB}};
        const paz::Vec relVel = static_cast<paz::Vec>(paz::to_mat(invParentAtt)*
            parentAngRate).cross(relPos);
        xVel() = relVel(0) + _parent->xVel() + forward(0);
        yVel() = relVel(1) + _parent->yVel() + forward(1);
        zVel() = relVel(2) + _parent->zVel() + forward(2);
        const paz::Vec up = rot.col(2);
        const paz::Vec east = paz::Vec{{0, 0, 1}}.cross(up).normalized();
        const paz::Vec north = up.cross(east);
        const double yaw = std::atan2(forward.dot(north), forward.dot(east));
        const double deltaYaw = paz::normalize_angle(_destYaw - yaw + paz::Pi) -
            paz::Pi;
        zAngRate() = deltaYaw;
        if(_walkTime > 10.)
        {
            _walkTime = 0.;
            _destYaw = paz::uniform(0., paz::TwoPi);
        }
    }
}

void Npc::setName(const std::string& name)
{
    _name = name;
}
