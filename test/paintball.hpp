#ifndef PAZ_ENGINE_TEST_PAINTBALL_HPP
#define PAZ_ENGINE_TEST_PAINTBALL_HPP

#include "PAZ_Engine"
#include "PAZ_Math"

class Paintball : public paz::Object
{
    paz::ObjectPtr _parent;
    paz::Vec _relPosPs;

public:
    static constexpr double LaunchSpeed = 15.;

    Paintball();
    void update(const paz::InputData&) override;
    void launch(const paz::Vec& pos, const paz::Vec& vel, const paz::Vec& dir);
    void onCollide(const Object&, double, double, double, double xB, double yB,
        double zB) override;
};

#endif
