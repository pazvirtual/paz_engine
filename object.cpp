#include "object.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"

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
static std::vector<double> XAtt;
static std::vector<double> YAtt;
static std::vector<double> ZAtt;
static std::vector<double> XAngRate;
static std::vector<double> YAngRate;
static std::vector<double> ZAngRate;
static std::vector<paz::Model> Mod;
static std::vector<paz::CollisionType> CType;
static std::vector<double> LocalNorX;
static std::vector<double> LocalNorY;
static std::vector<double> LocalNorZ;
static std::vector<double> XPrev;
static std::vector<double> YPrev;
static std::vector<double> ZPrev;
static std::vector<char> Grounded;
static std::vector<double> Height;
static std::vector<double> CRadius;

void paz::physics()
{
    const std::size_t n = X.size();
    XPrev.resize(n);
    YPrev.resize(n);
    ZPrev.resize(n);
    for(std::size_t i = 0; i < n; ++i)
    {
        XPrev[i] = X[i];
        X[i] += Window::FrameTime()*XVel[i];
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        YPrev[i] = Y[i];
        Y[i] += Window::FrameTime()*YVel[i];
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        ZPrev[i] = Z[i];
        Z[i] += Window::FrameTime()*ZVel[i];
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        double WAtt = std::sqrt(1. - XAtt[i]*XAtt[i] - YAtt[i]*YAtt[i] - ZAtt[i]
            *ZAtt[i]);
        const double deltaX = normalize_angle(0.5*Window::FrameTime()*XAngRate[i] + M_PI) - M_PI;
        const double deltaY = normalize_angle(0.5*Window::FrameTime()*YAngRate[i] + M_PI) - M_PI;
        const double deltaZ = normalize_angle(0.5*Window::FrameTime()*ZAngRate[i] + M_PI) - M_PI;
        XAtt[i] +=  WAtt   *deltaX - ZAtt[i]*deltaY + YAtt[i]*deltaZ;
        YAtt[i] +=  ZAtt[i]*deltaX + WAtt   *deltaY - XAtt[i]*deltaZ;
        ZAtt[i] += -YAtt[i]*deltaX + XAtt[i]*deltaY + WAtt   *deltaZ;
        WAtt    += -XAtt[i]*deltaX - YAtt[i]*deltaY - ZAtt[i]*deltaZ;
        const double signNorm = (WAtt < 0. ? -1. : 1.)*std::sqrt(XAtt[i]*XAtt[i]
            + YAtt[i]*YAtt[i] + ZAtt[i]*ZAtt[i] + WAtt*WAtt);
        XAtt[i] /= signNorm;
        YAtt[i] /= signNorm;
        ZAtt[i] /= signNorm;
    }
}

void paz::gravity()
{
    const std::size_t n = X.size();
    for(std::size_t i = 0; i < n; ++i)
    {
        if(CType[i] == CollisionType::Default)
        {
            const double radius = std::sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i]);
            XVel[i] -= 9.81*Window::FrameTime()*X[i]/radius;
            YVel[i] -= 9.81*Window::FrameTime()*Y[i]/radius;
            ZVel[i] -= 9.81*Window::FrameTime()*Z[i]/radius;
        }
    }
}

void paz::update()
{
    for(const auto& n : objects())
    {
        reinterpret_cast<Object*>(n.first)->update();
    }
}


paz::Object::Object() : _id(reinterpret_cast<std::uintptr_t>(this))
{
    objects()[_id] = X.size();
    X.push_back(0.);
    Y.push_back(0.);
    Z.push_back(0.);
    XVel.push_back(0.);
    YVel.push_back(0.);
    ZVel.push_back(0.);
    XAtt.push_back(0.);
    YAtt.push_back(0.);
    ZAtt.push_back(0.);
    XAngRate.push_back(0.);
    YAngRate.push_back(0.);
    ZAngRate.push_back(0.);
    Mod.emplace_back();
    CType.push_back(CollisionType::Default);
    LocalNorX.push_back(0.);
    LocalNorY.push_back(0.);
    LocalNorZ.push_back(1.);
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
    SWAP_AND_POP(XAtt);
    SWAP_AND_POP(YAtt);
    SWAP_AND_POP(ZAtt);
    SWAP_AND_POP(XAngRate);
    SWAP_AND_POP(YAngRate);
    SWAP_AND_POP(ZAngRate);
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

double& paz::Object::xAtt()
{
    return XAtt[objects().at(_id)];
}

double paz::Object::xAtt() const
{
    return XAtt[objects().at(_id)];
}

double& paz::Object::yAtt()
{
    return YAtt[objects().at(_id)];
}

double paz::Object::yAtt() const
{
    return YAtt[objects().at(_id)];
}

double& paz::Object::zAtt()
{
    return ZAtt[objects().at(_id)];
}

double paz::Object::zAtt() const
{
    return ZAtt[objects().at(_id)];
}

double& paz::Object::xAngRate()
{
    return XAngRate[objects().at(_id)];
}

double paz::Object::xAngRate() const
{
    return XAngRate[objects().at(_id)];
}

double& paz::Object::yAngRate()
{
    return YAngRate[objects().at(_id)];
}

double paz::Object::yAngRate() const
{
    return YAngRate[objects().at(_id)];
}

double& paz::Object::zAngRate()
{
    return ZAngRate[objects().at(_id)];
}

double paz::Object::zAngRate() const
{
    return ZAngRate[objects().at(_id)];
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

double paz::Object::xPrev() const
{
    return XPrev[objects().at(_id)];
}

double paz::Object::yPrev() const
{
    return YPrev[objects().at(_id)];
}

double paz::Object::zPrev() const
{
    return ZPrev[objects().at(_id)];
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
