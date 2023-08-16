#include "object.hpp"
#include "PAZ_Engine"

#define SWAP_AND_POP(x) std::swap(x[idx], x.back()); x.pop_back();

std::unordered_map<std::uintptr_t, std::size_t>& paz::objects()
{
    static std::unordered_map<std::uintptr_t, std::size_t> o;
    return o;
}

static std::vector<double> X;
static std::vector<double> Y;
static std::vector<double> Z;
static std::vector<double> XVel;
static std::vector<double> YVel;
static std::vector<double> ZVel;
static std::vector<paz::Model> Mod;
static std::vector<paz::CollisionType> CType;
static std::vector<double> LocalNorX;
static std::vector<double> LocalNorY;
static std::vector<double> LocalNorZ;
static std::vector<char> Grounded;
static std::vector<double> Height;
static std::vector<double> CRadius;

void paz::physics()
{
    const std::size_t n = X.size();
    for(std::size_t i = 0; i < n; ++i)
    {
        X[i] += paz::Window::FrameTime()*XVel[i];
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        Y[i] += paz::Window::FrameTime()*YVel[i];
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        Z[i] += paz::Window::FrameTime()*ZVel[i];
    }
}

paz::Object::Object() : _id(reinterpret_cast<std::uintptr_t>(this))
{
    objects()[_id] = X.size();
    X.emplace_back();
    Y.emplace_back();
    Z.emplace_back();
    XVel.emplace_back();
    YVel.emplace_back();
    ZVel.emplace_back();
    Mod.emplace_back();
    CType.push_back(CollisionType::Default);
    LocalNorX.push_back(0);
    LocalNorY.push_back(0);
    LocalNorZ.push_back(1);
    Grounded.push_back(false);
    Height.push_back(0.);
    CRadius.push_back(0.2);
}

paz::Object::~Object()
{
    const std::size_t idx = objects().at(_id);
    objects().erase(_id);
    SWAP_AND_POP(X);
    SWAP_AND_POP(Y);
    SWAP_AND_POP(Z);
    SWAP_AND_POP(XVel);
    SWAP_AND_POP(YVel);
    SWAP_AND_POP(ZVel);
    SWAP_AND_POP(Mod);
    SWAP_AND_POP(CType);
    SWAP_AND_POP(LocalNorX);
    SWAP_AND_POP(LocalNorY);
    SWAP_AND_POP(LocalNorZ);
    SWAP_AND_POP(Grounded);
    SWAP_AND_POP(Height);
    SWAP_AND_POP(CRadius);
}

void paz::Object::update() {}

void paz::Object::onCollide(const Object& /* o */) {}

void paz::Object::onInteract(const Object& /* o */) {}

double& paz::Object::x()
{
    return X[objects().at(_id)];
}

double paz::Object::x() const
{
    return X[objects().at(_id)];
}

double& paz::Object::y()
{
    return Y[objects().at(_id)];
}

double paz::Object::y() const
{
    return Y[objects().at(_id)];
}

double& paz::Object::z()
{
    return Z[objects().at(_id)];
}

double paz::Object::z() const
{
    return Z[objects().at(_id)];
}

double& paz::Object::xVel()
{
    return XVel[objects().at(_id)];
}

double paz::Object::xVel() const
{
    return XVel[objects().at(_id)];
}

double& paz::Object::yVel()
{
    return YVel[objects().at(_id)];
}

double paz::Object::yVel() const
{
    return YVel[objects().at(_id)];
}

double& paz::Object::zVel()
{
    return ZVel[objects().at(_id)];
}

double paz::Object::zVel() const
{
    return ZVel[objects().at(_id)];
}

paz::Model& paz::Object::model()
{
    return Mod[objects().at(_id)];
}

const paz::Model& paz::Object::model() const
{
    return Mod[objects().at(_id)];
}

paz::CollisionType& paz::Object::collisionType()
{
    return CType[objects().at(_id)];
}

const paz::CollisionType& paz::Object::collisionType() const
{
    return CType[objects().at(_id)];
}

double& paz::Object::localNorX()
{
    return LocalNorX[objects().at(_id)];
}

double paz::Object::localNorX() const
{
    return LocalNorX[objects().at(_id)];
}

double& paz::Object::localNorY()
{
    return LocalNorY[objects().at(_id)];
}

double paz::Object::localNorY() const
{
    return LocalNorY[objects().at(_id)];
}

double& paz::Object::localNorZ()
{
    return LocalNorZ[objects().at(_id)];
}

double paz::Object::localNorZ() const
{
    return LocalNorZ[objects().at(_id)];
}

char& paz::Object::grounded()
{
    return Grounded[objects().at(_id)];
}

bool paz::Object::grounded() const
{
    return Grounded[objects().at(_id)];
}

double& paz::Object::height()
{
    return Height[objects().at(_id)];
}

double paz::Object::height() const
{
    return Height[objects().at(_id)];
}

double& paz::Object::collisionRadius()
{
    return CRadius[objects().at(_id)];
}

double paz::Object::collisionRadius() const
{
    return CRadius[objects().at(_id)];
}
