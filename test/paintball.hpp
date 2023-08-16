#ifndef PAZ_ENGINE_TEST_PAINTBALL_HPP
#define PAZ_ENGINE_TEST_PAINTBALL_HPP

#include "PAZ_Engine"
#include "PAZ_Math"

class Paintball : public paz::Object
{
    const Object* _parent;
    double _relX, _relY, _relZ;

public:
    static constexpr double LaunchSpeed = 15.;

    Paintball();
    void update() override;
    void launch(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec& dir);
    void onCollide(const Object&) override;
};

#endif
