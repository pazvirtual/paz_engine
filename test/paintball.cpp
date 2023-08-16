#include "paintball.hpp"

static const paz::Model Model("icosphere3.pazmodel", 0, 0., 0.1, "", {0.5, 2.,
    2.});
static const paz::AudioTrack SoundEffect = []()
{
    std::vector<float> samples(44'100/200);
    for(std::size_t i = 0; i < samples.size(); ++i)
    {
        samples[i] = (i*2./(samples.size() - 1) - 1.)*4e-3;
    }
    std::vector<float> moreSamples(20*samples.size());
    for(int i = 0; i < 20; ++i)
    {
        std::copy(samples.begin(), samples.end(), moreSamples.begin() + i*
            samples.size());
    }
    for(std::size_t i = 0; i < 100; ++ i)
    {
        const double fac = static_cast<double>(i)/(100 - 1);
        moreSamples[i] *= fac;
        moreSamples[moreSamples.size() - 1 - i] *= fac;
    }
    return paz::AudioTrack(moreSamples);
}();

Paintball::Paintball() : paz::Object()
{
    collisionRadius() = 0.05;
    model() = Model;
    lights().push_back({0., 0., 0., 0.5, 2., 2., 0.1});
    collisionType() = paz::CollisionType::None;
    gravityType() = paz::GravityType::None;
}

void Paintball::update(const paz::InputData&)
{
    if(_parent)
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
        xAtt() = _parent->xAtt();
        yAtt() = _parent->yAtt();
        zAtt() = _parent->zAtt();
        xAngRate() = _parent->xAngRate();
        yAngRate() = _parent->yAngRate();
        zAngRate() = _parent->zAngRate();
    }
}

void Paintball::launch(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec&
    dir)
{
    _parent.reset();
    collisionType() = paz::CollisionType::Default;
    gravityType() = paz::GravityType::Default;
    x() = pos(0);
    y() = pos(1);
    z() = pos(2);
    xVel() = vel(0) + LaunchSpeed*dir(0);
    yVel() = vel(1) + LaunchSpeed*dir(1);
    zVel() = vel(2) + LaunchSpeed*dir(2);
}

void Paintball::onCollide(const Object& o, double, double, double, double xB,
    double yB, double zB)
{
    _parent.reset(o);
    const paz::Vec parentAtt{{_parent->xAtt(), _parent->yAtt(), _parent->zAtt(),
        std::sqrt(1. - _parent->xAtt()*_parent->xAtt() - _parent->yAtt()*
        _parent->yAtt() - _parent->zAtt()*_parent->zAtt())}};
    const paz::Vec relPos{{x() - xB, y() - yB, z() - zB}};
    const auto rot = paz::to_mat(parentAtt);
    _relPosPs = rot*relPos;
    const paz::Vec parentAngRate{{_parent->xAngRate(), _parent->yAngRate(),
        _parent->zAngRate()}};
    const paz::Vec relVel = rot.trans()*parentAngRate.cross(_relPosPs);
    xVel() = relVel(0) + _parent->xVel();
    yVel() = relVel(1) + _parent->yVel();
    zVel() = relVel(2) + _parent->zVel();
    collisionType() = paz::CollisionType::None;
    gravityType() = paz::GravityType::None;
}
