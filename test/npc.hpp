#ifndef PAZ_ENGINE_TEST_NPC_HPP
#define PAZ_ENGINE_TEST_NPC_HPP

#include "PAZ_Engine"

class Npc : public paz::Object
{
    paz::ObjectPtr _parent;
    double _destYaw;
    double _walkTime;
    paz::Object _head;
    int _collided = 0;
    std::string _name;

public:
    Npc();
    void update(const paz::InputData& input) override;
    void onCollide(const paz::Object& o, double xNor, double yNor, double zNor,
        double, double, double) override;
    void setName(const std::string& name);
};

#endif
