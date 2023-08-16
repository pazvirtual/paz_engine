#ifndef PAZ_ENGINE_TEST_NPC_HPP
#define PAZ_ENGINE_TEST_NPC_HPP

#include "PAZ_Engine"

class Npc : public paz::Object
{
    double _destYaw;
    double _walkTime;
    paz::Object _head;
    std::string _name;

public:
    Npc();
    void update(const paz::InputData& input) override;
    void onCollide(const Object&, double, double, double, double, double,
        double) override;
    void setName(const std::string& name);
};

#endif
