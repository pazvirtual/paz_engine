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
static std::vector<paz::GravityType> GType;
static std::vector<double> XPrev;
static std::vector<double> YPrev;
static std::vector<double> ZPrev;
static std::vector<double> CRadius;
static std::vector<double> XDown;
static std::vector<double> YDown;
static std::vector<double> ZDown;
static std::vector<double> StdGravParam;
static std::unordered_map<std::string, std::unordered_set<std::uintptr_t>>
    ObjectsByTag;

static void grav_ode(double x, double y, double z, double u, double v, double w,
    double& dx, double& dy, double& dz, double& du, double& dv, double& dw,
    const std::vector<std::size_t>& massiveIds, double gravity)
{
    dx = u;
    dy = v;
    dz = w;
    du = 0;
    dv = 0;
    dw = gravity;
    for(auto n : massiveIds)
    {
        const double deltaX = x - X[n];
        const double deltaY = y - Y[n];
        const double deltaZ = z - Z[n];
        const double radius = std::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*
            deltaZ);
        const double t0 = -StdGravParam[n]/(radius*radius*radius);
        du += t0*deltaX;
        dv += t0*deltaY;
        dw += t0*deltaZ;
    }
}

void paz::physics(double gravity)
{
    std::vector<std::size_t> massiveIds;
    const auto n = X.size();
    for(std::size_t i = 0; i < n; ++i)
    {
        if(CType[i] == CollisionType::World && StdGravParam[i] > 0.) //TEMP - should go off of `GType[i]`
        {
            massiveIds.push_back(i);
        }
    }
    XPrev = X;
    YPrev = Y;
    ZPrev = Z;
    for(std::size_t i = 0; i < n; ++i)
    {
        if(GType[i] == GravityType::Default)
        {
            double dx, dy, dz, du, dv, dw;
            grav_ode(X[i], Y[i], Z[i], XVel[i], YVel[i], ZVel[i], dx, dy, dz,
                du, dv, dw, massiveIds, gravity);
            const double k1_1 = App::PhysTime()*dx;
            const double k1_2 = App::PhysTime()*dy;
            const double k1_3 = App::PhysTime()*dz;
            const double k1_4 = App::PhysTime()*du;
            const double k1_5 = App::PhysTime()*dv;
            const double k1_6 = App::PhysTime()*dw;
            grav_ode(X[i] + 0.5*k1_1, Y[i] + 0.5*k1_2, Z[i] + 0.5*k1_3, XVel[i]
                + 0.5*k1_4, YVel[i] + 0.5*k1_5, ZVel[i] + 0.5*k1_6, dx, dy, dz,
                du, dv, dw, massiveIds, gravity);
            const double k2_1 = App::PhysTime()*dx;
            const double k2_2 = App::PhysTime()*dy;
            const double k2_3 = App::PhysTime()*dz;
            const double k2_4 = App::PhysTime()*du;
            const double k2_5 = App::PhysTime()*dv;
            const double k2_6 = App::PhysTime()*dw;
            grav_ode(X[i] + 0.5*k2_1, Y[i] + 0.5*k2_2, Z[i] + 0.5*k2_3, XVel[i]
                + 0.5*k2_4, YVel[i] + 0.5*k2_5, ZVel[i] + 0.5*k2_6, dx, dy, dz,
                du, dv, dw, massiveIds, gravity);
            const double k3_1 = App::PhysTime()*dx;
            const double k3_2 = App::PhysTime()*dy;
            const double k3_3 = App::PhysTime()*dz;
            const double k3_4 = App::PhysTime()*du;
            const double k3_5 = App::PhysTime()*dv;
            const double k3_6 = App::PhysTime()*dw;
            grav_ode(X[i] + k3_1, Y[i] + k3_2, Z[i] + k3_3, XVel[i] + k3_4,
                YVel[i] + k3_5, ZVel[i] + k3_6, dx, dy, dz, du, dv, dw,
                massiveIds, gravity);
            const double k4_1 = App::PhysTime()*dx;
            const double k4_2 = App::PhysTime()*dy;
            const double k4_3 = App::PhysTime()*dz;
            const double k4_4 = App::PhysTime()*du;
            const double k4_5 = App::PhysTime()*dv;
            const double k4_6 = App::PhysTime()*dw;
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
        }
        else
        {
            X[i] += App::PhysTime()*XVel[i];
            Y[i] += App::PhysTime()*YVel[i];
            Z[i] += App::PhysTime()*ZVel[i];
            double dx, dy, dz, du, dv, dw;
            grav_ode(X[i], Y[i], Z[i], XVel[i], YVel[i], ZVel[i], dx, dy, dz,
                du, dv, dw, massiveIds, gravity);
            XDown[i] = du;
            YDown[i] = dv;
            ZDown[i] = dw;
        }
        const double normSq = XDown[i]*XDown[i] + YDown[i]*YDown[i] + ZDown[i]*
            ZDown[i];
        if(normSq > 1e-6*1e-6)
        {
            const double invNorm = 1./std::sqrt(normSq);
            XDown[i] *= invNorm;
            YDown[i] *= invNorm;
            ZDown[i] *= invNorm;
        }
        else
        {
            XDown[i] = 0.;
            YDown[i] = 0.;
            ZDown[i] = 1.;
        }
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        double WAtt = std::sqrt(1. - XAtt[i]*XAtt[i] - YAtt[i]*YAtt[i] - ZAtt[
            i]*ZAtt[i]);
        const double deltaX = normalize_angle(0.5*App::PhysTime()*XAngRate[i] +
            M_PI) - M_PI;
        const double deltaY = normalize_angle(0.5*App::PhysTime()*YAngRate[i] +
            M_PI) - M_PI;
        const double deltaZ = normalize_angle(0.5*App::PhysTime()*ZAngRate[i] +
            M_PI) - M_PI;
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
    GType.push_back(GravityType::Default);
    CRadius.push_back(0.2);
    XDown.push_back(0.);
    YDown.push_back(0.);
    ZDown.push_back(0.);
    StdGravParam.push_back(0.);
}

paz::Object::Object(const Object& o) : _id(reinterpret_cast<std::uintptr_t>(
    this))
{
    const auto otherIdx = objects().at(o._id);
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
    PUSH_COPY(GType)
    PUSH_COPY(CRadius)
    PUSH_COPY(XDown)
    PUSH_COPY(YDown)
    PUSH_COPY(ZDown)
    PUSH_COPY(StdGravParam)
    for(auto& m : ObjectsByTag) //TEMP - inefficient
    {
        if(m.second.count(o._id))
        {
            m.second.insert(_id);
        }
    }
}

paz::Object& paz::Object::operator=(const Object& o)
{
    const auto idx = objects().at(_id);
    const auto otherIdx = objects().at(o._id);
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
    COPY(GType)
    COPY(CRadius)
    COPY(XDown)
    COPY(YDown)
    COPY(ZDown)
    COPY(StdGravParam)
    for(auto& m : ObjectsByTag) //TEMP - inefficient
    {
        if(m.second.count(o._id))
        {
            m.second.insert(_id);
        }
        else if(m.second.count(_id))
        {
            m.second.erase(_id);
        }
    }
    return *this;
}

paz::Object::~Object()
{
    const auto idx = objects().at(_id);
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
    SWAP_AND_POP(GType)
    SWAP_AND_POP(CRadius)
    SWAP_AND_POP(XDown)
    SWAP_AND_POP(YDown)
    SWAP_AND_POP(ZDown)
    SWAP_AND_POP(StdGravParam)
    for(auto& n : ObjectsByTag)
    {
        n.second.erase(_id);
    }
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
    for(auto n : ObjectsByTag.at(tag))
    {
        reinterpret_cast<Object*>(n)->onNotify(*this, data);
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

paz::GravityType& paz::Object::gravityType()
{
    return GType[objects().at(_id)];
}

const paz::GravityType& paz::Object::gravityType() const
{
    return GType[objects().at(_id)];
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

double& paz::Object::stdGravParam()
{
    return StdGravParam[objects().at(_id)];
}

double paz::Object::stdGravParam() const
{
    return StdGravParam[objects().at(_id)];
}

void paz::Object::addTag(const std::string& tag)
{
    ObjectsByTag[tag].insert(_id);
}

bool paz::Object::isTagged(const std::string& tag) const
{
    return ObjectsByTag.count(tag) && ObjectsByTag.at(tag).count(_id);
}

void paz::Object::computeAltitude(double& alt, Vec& nor, Vec& vel) const
{
    nor = {{0, 0, 1}};
    vel = Vec::Zero(3);
    const auto idx = objects().at(_id);
    alt = std::numeric_limits<double>::infinity();
    for(auto n : objects())
    {
        if(n.first == _id || CType[n.second] != CollisionType::World)
        {
            continue;
        }
        const double wAtt = std::sqrt(1. - XAtt[n.second]*XAtt[n.second] - YAtt[
            n.second]*YAtt[n.second] - ZAtt[n.second]*ZAtt[n.second]);
        const Mat rot = to_mat(Vec{{XAtt[n.second], YAtt[n.second], ZAtt[n.
            second], wAtt}});
        const Vec relPos = rot*Vec{{X[idx] - X[n.second], Y[idx] - Y[n.second],
            Z[idx] - Z[n.second]}};
        const Vec dir = rot*Vec{{XDown[idx], YDown[idx], ZDown[idx]}};
        double dist;
        Vec tempNor(3);
        Mod[n.second].castRay(relPos(0), relPos(1), relPos(2), dir(0), dir(1),
            dir(2), tempNor(0), tempNor(1), tempNor(2), dist);
        if(dist < alt)
        {
            alt = dist;
            nor = rot.trans()*tempNor;
            vel(0) = XVel[n.second];
            vel(1) = YVel[n.second];
            vel(2) = ZVel[n.second];
        }
    }
    alt -= CRadius[idx];
}
