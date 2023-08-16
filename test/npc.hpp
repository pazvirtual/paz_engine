#ifndef PAZ_ENGINE_TEST_NPC_HPP
#define PAZ_ENGINE_TEST_NPC_HPP

#include "PAZ_Engine"

class Npc : public paz::Object
{
    bool _collided;
    double _destYaw;
    double _walkTime;
    paz::Object _head;

public:
    Npc();
    void update(const paz::InputData& input) override;
    void onCollide(const Object&, double, double, double, double, double,
        double) override;
};

#endif
