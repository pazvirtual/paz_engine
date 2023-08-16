#include "object.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>
#include <unordered_map>
#include <unordered_set>

static constexpr std::size_t NumSteps = 100;

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
static std::vector<std::vector<std::array<double, 7>>> Lights;
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

void paz::do_physics(double gravity, double timestep)
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
            const double k1_1 = timestep*dx;
            const double k1_2 = timestep*dy;
            const double k1_3 = timestep*dz;
            const double k1_4 = timestep*du;
            const double k1_5 = timestep*dv;
            const double k1_6 = timestep*dw;
            grav_ode(X[i] + 0.5*k1_1, Y[i] + 0.5*k1_2, Z[i] + 0.5*k1_3, XVel[i]
                + 0.5*k1_4, YVel[i] + 0.5*k1_5, ZVel[i] + 0.5*k1_6, dx, dy, dz,
                du, dv, dw, massiveIds, gravity);
            const double k2_1 = timestep*dx;
            const double k2_2 = timestep*dy;
            const double k2_3 = timestep*dz;
            const double k2_4 = timestep*du;
            const double k2_5 = timestep*dv;
            const double k2_6 = timestep*dw;
            grav_ode(X[i] + 0.5*k2_1, Y[i] + 0.5*k2_2, Z[i] + 0.5*k2_3, XVel[i]
                + 0.5*k2_4, YVel[i] + 0.5*k2_5, ZVel[i] + 0.5*k2_6, dx, dy, dz,
                du, dv, dw, massiveIds, gravity);
            const double k3_1 = timestep*dx;
            const double k3_2 = timestep*dy;
            const double k3_3 = timestep*dz;
            const double k3_4 = timestep*du;
            const double k3_5 = timestep*dv;
            const double k3_6 = timestep*dw;
            grav_ode(X[i] + k3_1, Y[i] + k3_2, Z[i] + k3_3, XVel[i] + k3_4,
                YVel[i] + k3_5, ZVel[i] + k3_6, dx, dy, dz, du, dv, dw,
                massiveIds, gravity);
            const double k4_1 = timestep*dx;
            const double k4_2 = timestep*dy;
            const double k4_3 = timestep*dz;
            const double k4_4 = timestep*du;
            const double k4_5 = timestep*dv;
            const double k4_6 = timestep*dw;
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
            X[i] += timestep*XVel[i];
            Y[i] += timestep*YVel[i];
            Z[i] += timestep*ZVel[i];
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
        const double deltaX = normalize_angle(0.5*timestep*XAngRate[i] + Pi) -
            Pi;
        const double deltaY = normalize_angle(0.5*timestep*YAngRate[i] + Pi) -
            Pi;
        const double deltaZ = normalize_angle(0.5*timestep*ZAngRate[i] + Pi) -
            Pi;
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

void paz::do_collisions(Threadpool& threads, double timestep)
{
    // Identify all objects that can collide and precompute as much as possible.
    std::vector<std::size_t> a;
    std::vector<std::size_t> b;
    a.reserve(objects().size());
    b.reserve(objects().size());
    for(auto& n : objects())
    {
        Object* o = reinterpret_cast<Object*>(n.first);
        if(o->collisionType() == CollisionType::Default)
        {
            a.push_back(n.second);
        }
        if(o->collisionType() == CollisionType::World)
        {
            b.push_back(n.second);
        }
    }

    // `a[i]` may collide with any `b[c[i][j]]`.
    std::vector<std::vector<std::size_t>> c(a.size());
    for(std::size_t i = 0; i < a.size(); ++i)
    {
        c[i].reserve(b.size());
        for(std::size_t j = 0; j < b.size(); ++j)
        {
            if(Mod[b[j]].sweepVol(XPrev[a[i]], YPrev[a[i]], ZPrev[a[i]], X[a[
                i]], Y[a[i]], Z[a[i]], XPrev[b[j]], YPrev[b[j]], ZPrev[b[j]], X[
                b[j]], Y[b[j]], Z[b[j]], CRadius[a[i]]))
            {
                c[i].push_back(j);
            }
        }
    }

    std::vector<double> times(NumSteps);
    for(std::size_t i = 0; i < NumSteps; ++i)
    {
        times[i] = static_cast<double>(i)/(NumSteps - 1);
    }

    std::vector<Mat> bRot(b.size());
    std::vector<std::vector<double>> bX(b.size(), std::vector<double>(
        NumSteps));
    std::vector<std::vector<double>> bY(b.size(), std::vector<double>(
        NumSteps));
    std::vector<std::vector<double>> bZ(b.size(), std::vector<double>(
        NumSteps));
    for(std::size_t i = 0; i < b.size(); ++i)
    {
        const double wAtt = std::sqrt(1. - XAtt[b[i]]*XAtt[b[i]] - YAtt[b[i]]*
            YAtt[b[i]] - ZAtt[b[i]]*ZAtt[b[i]]);
        const Vec att{{XAtt[b[i]], YAtt[b[i]], ZAtt[b[i]], wAtt}};
        bRot[i] = to_mat(att);
        for(std::size_t j = 0; j < NumSteps; ++j)
        {
            bX[i][j] = XPrev[b[i]] + times[j]*(X[b[i]] - XPrev[b[i]]);
            bY[i][j] = YPrev[b[i]] + times[j]*(Y[b[i]] - YPrev[b[i]]);
            bZ[i][j] = ZPrev[b[i]] + times[j]*(Z[b[i]] - ZPrev[b[i]]);
        }
    }

    // Find and handle collisions.
    std::vector<LockedCv> lcvs;
    lcvs.reserve(a.size());
    for(std::size_t i = 0; i < a.size(); ++i)
    {
        lcvs.push_back(threads.enqueue([=]()
        {
            for(std::size_t j = 0; j < NumSteps; ++j)
            {
                double x = XPrev[a[i]] + times[j]*(X[a[i]] - XPrev[a[i]]);
                double y = YPrev[a[i]] + times[j]*(Y[a[i]] - YPrev[a[i]]);
                double z = ZPrev[a[i]] + times[j]*(Z[a[i]] - ZPrev[a[i]]);
                double minDist = std::numeric_limits<double>::infinity();
                std::size_t idx = 0;
                double xNor = 0.;
                double yNor = 0.;
                double zNor = 1.;
                std::vector<std::pair<std::size_t, std::array<double, 3>>>
                    collisions;
                for(auto n : c[i])
                {
                    const double x1 = x - bX[n][j];
                    const double y1 = y - bY[n][j];
                    const double z1 = z - bZ[n][j];
                    const Vec relPos = bRot[n]*Vec{{x1, y1, z1}};
                    const double x2 = relPos(0);
                    const double y2 = relPos(1);
                    const double z2 = relPos(2);
                    double xNew, yNew, zNew, xNorTemp, yNorTemp, zNorTemp;
                    const double dist = Mod[b[n]].collide(x2, y2, z2, CRadius[a[
                        i]], xNew, yNew, zNew, xNorTemp, yNorTemp, zNorTemp);
                    if(dist < CRadius[a[i]] - 1e-9)
                    {
                        collisions.emplace_back(n, std::array<double, 3>{
                            xNorTemp, yNorTemp, zNorTemp});
                        const Vec newPos = bRot[n].trans()*Vec{{xNew, yNew,
                            zNew}};
                        xNew = newPos(0);
                        yNew = newPos(1);
                        zNew = newPos(2);
                        x += xNew - x1;
                        y += yNew - y1;
                        z += zNew - z1;
                        if(dist < minDist)
                        {
                            minDist = dist;
                            idx = n;
                            const Vec nor = bRot[n].trans()*Vec{{xNorTemp,
                                yNorTemp, zNorTemp}};
                            xNor = nor(0);
                            yNor = nor(1);
                            zNor = nor(2);
                        }
                    }
                }
                if(std::isfinite(minDist))
                {
                    const double xVel = XVel[a[i]] - XVel[b[idx]];
                    const double yVel = YVel[a[i]] - YVel[b[idx]];
                    const double zVel = ZVel[a[i]] - ZVel[b[idx]];
                    const double norVel = xVel*xNor + yVel*yNor + zVel*zNor;
                    if(norVel < 0.)
                    {
                        XVel[a[i]] -= norVel*xNor;
                        YVel[a[i]] -= norVel*yNor;
                        ZVel[a[i]] -= norVel*zNor;
                    }

                    // Adjust positions to fit new trajectory.
                    const double time0 = timestep*(times[j] - times[0]);
                    XPrev[a[i]] = x - time0*XVel[a[i]];
                    YPrev[a[i]] = y - time0*YVel[a[i]];
                    ZPrev[a[i]] = z - time0*ZVel[a[i]];
                    const double time1 = timestep*(times.back() - times[j]);
                    X[a[i]] = x + time1*XVel[a[i]];
                    Y[a[i]] = y + time1*YVel[a[i]];
                    Z[a[i]] = z + time1*ZVel[a[i]];

                    // Collision response.
                    std::swap(X[a[i]], x);
                    std::swap(Y[a[i]], y);
                    std::swap(Z[a[i]], z);
                    Object& aObj = *reinterpret_cast<Object*>(Ids[a[i]]);
                    for(const auto& n : collisions)
                    {
                        Object& bObj = *reinterpret_cast<Object*>(Ids[b[n.
                            first]]);
                        aObj.onCollide(bObj, n.second[0], n.second[1], n.second[
                            2], bX[n.first][j], bY[n.first][j], bZ[n.first][j]);
                    }
                    std::swap(X[a[i]], x);
                    std::swap(Y[a[i]], y);
                    std::swap(Z[a[i]], z);

                    // Adjust positions again.
                    XPrev[a[i]] = x - time0*XVel[a[i]];
                    YPrev[a[i]] = y - time0*YVel[a[i]];
                    ZPrev[a[i]] = z - time0*ZVel[a[i]];
                    X[a[i]] = x + time1*XVel[a[i]];
                    Y[a[i]] = y + time1*YVel[a[i]];
                    Z[a[i]] = z + time1*ZVel[a[i]];
                }
            }
        }
        ));
    }
    for(auto& n : lcvs)
    {
        n.wait();
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
    Lights.emplace_back();
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
    PUSH_COPY(Lights)
    for(auto& n : ObjectsByTag) //TEMP - inefficient
    {
        if(n.second.count(o._id))
        {
            n.second.insert(_id);
        }
    }
}

paz::Object& paz::Object::operator=(const Object& o)
{
    // If destination is in a valid state, copy data, otherwise, copy construct.
    const auto otherIdx = objects().at(o._id);
    if(objects().count(_id))
    {
        const auto idx = objects().at(_id);
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
        COPY(Lights)
        for(auto& n : ObjectsByTag) //TEMP - inefficient
        {
            if(n.second.count(o._id))
            {
                n.second.insert(_id);
            }
            else if(n.second.count(_id))
            {
                n.second.erase(_id);
            }
        }
    }
    else
    {
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
        PUSH_COPY(Lights)
        for(auto& n : ObjectsByTag) //TEMP - inefficient
        {
            if(n.second.count(o._id))
            {
                n.second.insert(_id);
            }
        }
    }
    return *this;
}

paz::Object::Object(Object&& o) noexcept : _id(reinterpret_cast<std::uintptr_t>(
    this))
{
    const auto idx = objects().at(o._id);
    objects()[_id] = idx;
    objects().erase(o._id);
    Ids[idx] = _id;
    for(auto& n : ObjectsByTag) //TEMP - inefficient
    {
        if(n.second.count(o._id))
        {
            n.second.erase(o._id);
            n.second.insert(_id);
        }
    }
}

paz::Object& paz::Object::operator=(Object&& o) noexcept
{
    // If destination is in a valid state, just swap, otherwise, copy construct.
    const auto otherIdx = objects().at(o._id);
    if(objects().count(_id))
    {
        const auto idx = objects().at(_id);
        objects().at(_id) = otherIdx;
        objects().at(o._id) = idx;
        std::swap(Ids[idx], Ids[otherIdx]);
    }
    else
    {
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
        PUSH_COPY(Lights)
        for(auto& n : ObjectsByTag) //TEMP - inefficient
        {
            if(n.second.count(o._id))
            {
                n.second.insert(_id);
            }
        }
    }
    return *this;
}

paz::Object::~Object()
{
    // If this object has been moved from, nothing to do.
    if(!objects().count(_id))
    {
        return;
    }
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
    SWAP_AND_POP(Lights)
    for(auto& n : ObjectsByTag) //TEMP - inefficient
    {
        n.second.erase(_id);
    }
}

void paz::Object::update(const InputData&) {}

void paz::Object::onCollide(const Object&, double, double, double, double,
    double, double) {}

void paz::Object::onInteract(const Object&) {}

void paz::Object::onNotify(const Object&, const Bytes&) {}

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

std::vector<std::array<double, 7>>& paz::Object::lights()
{
    return Lights[objects().at(_id)];
}

const std::vector<std::array<double, 7>>& paz::Object::lights() const
{
    return Lights[objects().at(_id)];
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
