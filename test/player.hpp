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
    paz::Vec _relPosPs;
    bool _moving = false;

    double _pitch = 0.;
    double _prevGravPitch = 0.;
    paz::Vec _mousePos = paz::Vec::Zero(2);
    paz::Object _head; // camera
    bool _collided = false;

    Paintball _paintball;

public:
    Player();
    void update(const paz::InputData& input) final;
    void onCollide(const paz::Object& o, double xNor, double yNor, double zNor,
        double xB, double yB, double zB) final;
    const paz::Object& head() const;
};

#endif
