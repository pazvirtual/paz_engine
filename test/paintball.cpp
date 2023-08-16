#include "paintball.hpp"

static const paz::Model Model("icosphere3.obj", 0, 0., 0.1, "", 2.f);
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

Paintball::Paintball() : paz::Object(), _parent(nullptr)
{
    collisionRadius() = 0.05;
    model() = Model;
    lights().push_back({0., 0., 0., 0.5, 2., 2., 0.1});
}

void Paintball::update()
{
    if(_parent)
    {
        x() = _relX + _parent->x();
        y() = _relY + _parent->y();
        z() = _relZ + _parent->z();
        xVel() = _parent->xVel();
        yVel() = _parent->yVel();
        zVel() = _parent->zVel();
    }
}

void Paintball::launch(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec&
    dir)
{
    _parent = nullptr;
    collisionType() = paz::CollisionType::Default;
    gravityType() = paz::GravityType::Default;
    x() = pos(0);
    y() = pos(1);
    z() = pos(2);
    xVel() = vel(0) + LaunchSpeed*dir(0);
    yVel() = vel(1) + LaunchSpeed*dir(1);
    zVel() = vel(2) + LaunchSpeed*dir(2);
}

void Paintball::onCollide(const Object& o)
{
    _parent = &o;
    _relX = x() - o.x();
    _relY = y() - o.y();
    _relZ = z() - o.z();
    xVel() = o.xVel();
    yVel() = o.yVel();
    zVel() = o.zVel();
    collisionType() = paz::CollisionType::None;
    gravityType() = paz::GravityType::None;
}
