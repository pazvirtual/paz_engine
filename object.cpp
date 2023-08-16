#include "object.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>
#include <unordered_map>
#include <unordered_set>

#define SWAP_AND_POP(x) std::swap(x[idx], x.back()); x.pop_back();
#define PUSH_COPY(x) x.push_back(x[otherIdx]);
#define COPY(x) x[idx] = x[otherIdx];

std::unordered_map<std::uintptr_t, std::size_t>& paz::objects()
{
    static std::unordered_map<std::uintptr_t, std::size_t> o;
    return o;
}

static std::vector<std::uintptr_t> Ids;
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
static std::vector<double> Altitude;
static std::vector<double> CRadius;
static std::vector<double> XDown;
static std::vector<double> YDown;
static std::vector<double> ZDown;
static std::unordered_map<std::string, std::unordered_set<paz::Object*>>
    ObjectsByTag;

static void grav_ode(double x, double y, double z, double u, double v, double w,
    double& dx, double& dy, double& dz, double& du, double& dv, double& dw)
{
    static constexpr double StdGravParam = 0.2*9.81*50.*50.;
    const double radius0 = std::sqrt(x*x + y*y + z*z);
    const double radius1 = std::sqrt((x - 50.)*(x - 50.) + y*y + z*z);
    const double radius2 = std::sqrt((x + 50.)*(x + 50.) + y*y + z*z);
    const double radius3 = std::sqrt(x*x + (y - 50.)*(y - 50.) + z*z);
    const double radius4 = std::sqrt(x*x + (y + 50.)*(y + 50.) + z*z);
    const double c0 = -StdGravParam/(radius0*radius0*radius0);
    const double c1 = -StdGravParam/(radius1*radius1*radius1);
    const double c2 = -StdGravParam/(radius2*radius2*radius2);
    const double c3 = -StdGravParam/(radius3*radius3*radius3);
    const double c4 = -StdGravParam/(radius4*radius4*radius4);
    const double sum = c0 + c1 + c2 + c3 + c4;
    dx = u;
    dy = v;
    dz = w;
    du = sum*x - c1*50. + c2*50.;
    dv = sum*y - c3*50. + c4*50.;
    dw = sum*z;
}

void paz::physics()
{
    XPrev = X;
    YPrev = Y;
    ZPrev = Z;
    const std::size_t n = X.size();
    for(std::size_t i = 0; i < n; ++i)
    {
        if(CType[i] == CollisionType::Default)
        {
            double dx, dy, dz, du, dv, dw;
            grav_ode(X[i], Y[i], Z[i], XVel[i], YVel[i], ZVel[i], dx, dy, dz,
                du, dv, dw);
            const double k1_1 = Window::FrameTime()*dx;
            const double k1_2 = Window::FrameTime()*dy;
            const double k1_3 = Window::FrameTime()*dz;
            const double k1_4 = Window::FrameTime()*du;
            const double k1_5 = Window::FrameTime()*dv;
            const double k1_6 = Window::FrameTime()*dw;
            grav_ode(X[i] + 0.5*k1_1, Y[i] + 0.5*k1_2, Z[i] + 0.5*k1_3, XVel[i]
                + 0.5*k1_4, YVel[i] + 0.5*k1_5, ZVel[i] + 0.5*k1_6, dx, dy, dz,
                du, dv, dw);
            const double k2_1 = Window::FrameTime()*dx;
            const double k2_2 = Window::FrameTime()*dy;
            const double k2_3 = Window::FrameTime()*dz;
            const double k2_4 = Window::FrameTime()*du;
            const double k2_5 = Window::FrameTime()*dv;
            const double k2_6 = Window::FrameTime()*dw;
            grav_ode(X[i] + 0.5*k2_1, Y[i] + 0.5*k2_2, Z[i] + 0.5*k2_3, XVel[i]
                + 0.5*k2_4, YVel[i] + 0.5*k2_5, ZVel[i] + 0.5*k2_6, dx, dy, dz,
                du, dv, dw);
            const double k3_1 = Window::FrameTime()*dx;
            const double k3_2 = Window::FrameTime()*dy;
            const double k3_3 = Window::FrameTime()*dz;
            const double k3_4 = Window::FrameTime()*du;
            const double k3_5 = Window::FrameTime()*dv;
            const double k3_6 = Window::FrameTime()*dw;
            grav_ode(X[i] + k3_1, Y[i] + k3_2, Z[i] + k3_3, XVel[i] + k3_4,
                YVel[i] + k3_5, ZVel[i] + k3_6, dx, dy, dz, du, dv, dw);
            const double k4_1 = Window::FrameTime()*dx;
            const double k4_2 = Window::FrameTime()*dy;
            const double k4_3 = Window::FrameTime()*dz;
            const double k4_4 = Window::FrameTime()*du;
            const double k4_5 = Window::FrameTime()*dv;
            const double k4_6 = Window::FrameTime()*dw;
            static constexpr double c = 1./6.;
            X[i] += c*(k1_1 + 2.*(k2_1 + k3_1) + k4_1);
            Y[i] += c*(k1_2 + 2.*(k2_2 + k3_2) + k4_2);
            Z[i] += c*(k1_3 + 2.*(k2_3 + k3_3) + k4_3);
            XDown[i] = c*(k1_4 + 2.*(k2_4 + k3_4) + k4_4);
            YDown[i] = c*(k1_5 + 2.*(k2_5 + k3_5) + k4_5);
            ZDown[i] = c*(k1_6 + 2.*(k2_6 + k3_6) + k4_6);
            XVel[i] += XDown[i];
            YVel[i] += YDown[i];
            ZVel[i] += ZDown[i];
            const double invNorm = 1./std::sqrt(XDown[i]*XDown[i] + YDown[i]*
                YDown[i] + ZDown[i]*ZDown[i]);
            if(std::isfinite(invNorm))
            {
                XDown[i] *= invNorm;
                YDown[i] *= invNorm;
                ZDown[i] *= invNorm;
            }
        }
        else
        {
            X[i] += Window::FrameTime()*XVel[i];
            Y[i] += Window::FrameTime()*YVel[i];
            Z[i] += Window::FrameTime()*ZVel[i];
        }
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        double WAtt = std::sqrt(1. - XAtt[i]*XAtt[i] - YAtt[i]*YAtt[i] - ZAtt[
            i]*ZAtt[i]);
        const double deltaX = normalize_angle(0.5*Window::FrameTime()*XAngRate[
            i] + M_PI) - M_PI;
        const double deltaY = normalize_angle(0.5*Window::FrameTime()*YAngRate[
            i] + M_PI) - M_PI;
        const double deltaZ = normalize_angle(0.5*Window::FrameTime()*ZAngRate[
            i] + M_PI) - M_PI;
        XAtt[i] +=  WAtt   *deltaX - ZAtt[i]*deltaY + YAtt[i]*deltaZ;
        YAtt[i] +=  ZAtt[i]*deltaX + WAtt   *deltaY - XAtt[i]*deltaZ;
        ZAtt[i] += -YAtt[i]*deltaX + XAtt[i]*deltaY + WAtt   *deltaZ;
        WAtt    += -XAtt[i]*deltaX - YAtt[i]*deltaY - ZAtt[i]*deltaZ;
        const double invSignNorm = (WAtt < 0. ? -1. : 1.)/std::sqrt(XAtt[i]*
            XAtt[i] + YAtt[i]*YAtt[i] + ZAtt[i]*ZAtt[i] + WAtt*WAtt);
        XAtt[i] *= invSignNorm;
        YAtt[i] *= invSignNorm;
        ZAtt[i] *= invSignNorm;
    }
    std::fill(Altitude.begin(), Altitude.end(), std::numeric_limits<double>::
        infinity());
}

paz::Object::Object() : _id(reinterpret_cast<std::uintptr_t>(this))
{
    objects()[_id] = X.size();
    Ids.push_back(_id);
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
    Altitude.push_back(std::numeric_limits<double>::infinity());
    CRadius.push_back(0.2);
    XDown.push_back(0.);
    YDown.push_back(0.);
    ZDown.push_back(0.);
}

paz::Object::Object(const Object& o) : _id(reinterpret_cast<std::uintptr_t>(
    this))
{
    const std::size_t otherIdx = objects().at(o._id);
    objects()[_id] = X.size();
    Ids.push_back(_id);
    PUSH_COPY(X)
    PUSH_COPY(Y)
    PUSH_COPY(Z)
    PUSH_COPY(XVel)
    PUSH_COPY(YVel)
    PUSH_COPY(ZVel)
    PUSH_COPY(XAtt)
    PUSH_COPY(YAtt)
    PUSH_COPY(ZAtt)
    PUSH_COPY(XAngRate)
    PUSH_COPY(YAngRate)
    PUSH_COPY(ZAngRate)
    PUSH_COPY(Mod)
    PUSH_COPY(CType)
    PUSH_COPY(LocalNorX)
    PUSH_COPY(LocalNorY)
    PUSH_COPY(LocalNorZ)
    PUSH_COPY(Altitude)
    PUSH_COPY(CRadius)
    PUSH_COPY(XDown)
    PUSH_COPY(YDown)
    PUSH_COPY(ZDown)
}

paz::Object& paz::Object::operator=(const Object& o)
{
    const std::size_t idx = objects().at(_id);
    const std::size_t otherIdx = objects().at(o._id);
    COPY(X)
    COPY(Y)
    COPY(Z)
    COPY(XVel)
    COPY(YVel)
    COPY(ZVel)
    COPY(XAtt)
    COPY(YAtt)
    COPY(ZAtt)
    COPY(XAngRate)
    COPY(YAngRate)
    COPY(ZAngRate)
    COPY(Mod)
    COPY(CType)
    COPY(LocalNorX)
    COPY(LocalNorY)
    COPY(LocalNorZ)
    COPY(Altitude)
    COPY(CRadius)
    COPY(XDown)
    COPY(YDown)
    COPY(ZDown)
    return *this;
}

paz::Object::Object(Object&& o) : _id(reinterpret_cast<std::uintptr_t>(this))
{
    const std::size_t idx = objects().at(o._id);
    objects().erase(o._id);
    objects()[_id] = idx;
    Ids[idx] = _id;
    o._moved = true;
}

paz::Object::~Object()
{
    if(_moved)
    {
        return;
    }
    const std::size_t idx = objects().at(_id);
    if(Ids.size() > 1)
    {
        objects().at(Ids.back()) = idx;
    }
    objects().erase(_id);
    SWAP_AND_POP(Ids)
    SWAP_AND_POP(X)
    SWAP_AND_POP(Y)
    SWAP_AND_POP(Z)
    SWAP_AND_POP(XVel)
    SWAP_AND_POP(YVel)
    SWAP_AND_POP(ZVel)
    SWAP_AND_POP(XAtt)
    SWAP_AND_POP(YAtt)
    SWAP_AND_POP(ZAtt)
    SWAP_AND_POP(XAngRate)
    SWAP_AND_POP(YAngRate)
    SWAP_AND_POP(ZAngRate)
    SWAP_AND_POP(Mod)
    SWAP_AND_POP(CType)
    SWAP_AND_POP(LocalNorX)
    SWAP_AND_POP(LocalNorY)
    SWAP_AND_POP(LocalNorZ)
    SWAP_AND_POP(Altitude)
    SWAP_AND_POP(CRadius)
    SWAP_AND_POP(XDown)
    SWAP_AND_POP(YDown)
    SWAP_AND_POP(ZDown)
//    for(auto& n : ObjectsByTag)
//    {
//        n.second.erase(this);
//    }
}

void paz::Object::update() {}

void paz::Object::onCollide(const Object& /* o */) {}

void paz::Object::onInteract(const Object& /* o */) {}

void paz::Object::onNotify(const Object& /* o */, const Bytes& /* data */) {}

void paz::Object::notify(Object& o, const Bytes& data) const
{
    o.onNotify(*this, data);
}

void paz::Object::notifyTagged(const std::string& tag, const Bytes& data) const
{
    if(!ObjectsByTag.count(tag))
    {
        return;
    }
    for(auto& n : ObjectsByTag.at(tag))
    {
        n->onNotify(*this, data);
    }
}

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

void paz::Object::setLocalNorX(double n)
{
    LocalNorX[objects().at(_id)] = n;
}

double paz::Object::localNorX() const
{
    return LocalNorX[objects().at(_id)];
}

void paz::Object::setLocalNorY(double n)
{
    LocalNorY[objects().at(_id)] = n;
}

double paz::Object::localNorY() const
{
    return LocalNorY[objects().at(_id)];
}

void paz::Object::setLocalNorZ(double n)
{
    LocalNorZ[objects().at(_id)] = n;
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

void paz::Object::setAltitude(double n)
{
    Altitude[objects().at(_id)] = n;
}

double paz::Object::altitude() const
{
    return Altitude[objects().at(_id)];
}

double& paz::Object::collisionRadius()
{
    return CRadius[objects().at(_id)];
}

double paz::Object::collisionRadius() const
{
    return CRadius[objects().at(_id)];
}

double paz::Object::xDown() const
{
    return XDown[objects().at(_id)];
}

double paz::Object::yDown() const
{
    return YDown[objects().at(_id)];
}

double paz::Object::zDown() const
{
    return ZDown[objects().at(_id)];
}

void paz::Object::addTag(const std::string& tag)
{
    throw std::logic_error("NOT IMPLEMENTED");
//    ObjectsByTag[tag].insert(this);
}
