#include "npc.hpp"
#include "PAZ_Math"

static const paz::Model Body("persontest.obj", 0, -0.2);
static const paz::Model Head("persontest.obj", 1, -0.1);

Npc::Npc() : _destYaw(paz::uniform(0., paz::TwoPi)), _walkTime(0.)
{
    _head.collisionType() = paz::CollisionType::None;
    _head.gravityType() = paz::GravityType::None;
    _head.model() = Head;
    model() = Body;
}

void Npc::update()
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

void Npc::onCollide(const Object&)
{
    const paz::Vec att{{xAtt(), yAtt(), zAtt(), std::sqrt(1. - xAtt()*xAtt() -
        yAtt()*yAtt() - zAtt()*zAtt())}};
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
    const double deltaYaw = paz::normalize_angle(_destYaw - yaw + paz::Pi) -
        paz::Pi;
    zAngRate() = deltaYaw;
    _walkTime += paz::App::PhysTime();
    if(_walkTime > 10.)
    {
        _walkTime = 0.;
        _destYaw = paz::uniform(0., paz::TwoPi);
    }
}
