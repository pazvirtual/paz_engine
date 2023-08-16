#ifndef PAZ_ENGINE_TEST_PLAYER_HPP
#define PAZ_ENGINE_TEST_PLAYER_HPP

#include "paintball.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"

class Player : public paz::Object
{
    enum class Regime
    {
        Grounded, Low, Floating
    };

    static constexpr double LowAltitude = 10.;

    paz::ObjectPtr _parent;
    double _relX, _relY, _relZ;
    bool _moving = false;

    double _pitch = 0.;
    double _prevGravPitch = 0.;
    paz::Vec _mousePos = paz::Vec::Zero(2);
    paz::Object _head; // camera
    int _collided = 0;

    Paintball _paintball;

public:
    Player();
    void update() final;
    void onCollide(const paz::Object& o) final;
    const paz::Object& head() const;
};

#endif
