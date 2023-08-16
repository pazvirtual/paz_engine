#include "paintball.hpp"

static const paz::Model Model("icosphere2.obj", 0, 0., 0.1);
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

Paintball::Paintball(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec&
    dir) : paz::Object()
{
    x() = pos(0);
    y() = pos(1);
    z() = pos(2);
    xVel() = vel(0) + LaunchSpeed*dir(0);
    yVel() = vel(1) + LaunchSpeed*dir(1);
    zVel() = vel(2) + LaunchSpeed*dir(2);
    collisionRadius() = 0.05;
    model() = Model;
}

void Paintball::onCollide(const Object&)
{
    paz::AudioEngine::Play(SoundEffect, false);
    xVel() = 0.;
    yVel() = 0.;
    zVel() = 0.;
    collisionType() = paz::CollisionType::None;
    gravityType() = paz::GravityType::None;
}
