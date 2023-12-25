#include "object.hpp"
#include "PAZ_Engine"
#include "PAZ_Math"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <optional>

#define SWAP_AND_POP(x) std::swap(x[idx], x.back()); x.pop_back();
#define PUSH_COPY(x) x.push_back(x[otherIdx]);
#define COPY(x) x[idx] = x[otherIdx];
#define GET_CMESH(idx) (_cMesh[idx]._t ? _cMesh[idx] : _mod[idx])

static constexpr std::size_t NumSteps = 20;

static std::vector<std::uintptr_t> _ids;
static std::vector<double> _x;
static std::vector<double> _y;
static std::vector<double> _z;
static std::vector<double> _xVel;
static std::vector<double> _yVel;
static std::vector<double> _zVel;
static std::vector<double> _xAtt;
static std::vector<double> _yAtt;
static std::vector<double> _zAtt;
static std::vector<double> _xAngRate;
static std::vector<double> _yAngRate;
static std::vector<double> _zAngRate;
static std::vector<paz::Model> _mod;
static std::vector<paz::CollisionMesh> _cMesh;
static std::vector<paz::CollisionType> _cType;
static std::vector<paz::GravityType> _gType;
static std::vector<double> _xPrev;
static std::vector<double> _yPrev;
static std::vector<double> _zPrev;
static std::vector<double> _xAttPrev;
static std::vector<double> _yAttPrev;
static std::vector<double> _zAttPrev;
static std::vector<double> _cRadius;
static std::vector<double> _xDown;
static std::vector<double> _yDown;
static std::vector<double> _zDown;
static std::vector<double> _stdGravParam;
static std::vector<std::vector<std::array<double, 7>>> _lights;
static std::unordered_map<std::string, std::unordered_set<std::uintptr_t>>
    _objectsByTag;

std::unordered_map<std::uintptr_t, std::size_t>& paz::objects()
{
    static std::unordered_map<std::uintptr_t, std::size_t> o;
    return o;
}

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
        const double deltaX = x - _x[n];
        const double deltaY = y - _y[n];
        const double deltaZ = z - _z[n];
        const double radius = std::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*
            deltaZ);
        const double t0 = -_stdGravParam[n]/(radius*radius*radius);
        du += t0*deltaX;
        dv += t0*deltaY;
        dw += t0*deltaZ;
    }
}

void paz::do_physics(double gravity, double timestep)
{
    std::vector<std::size_t> massiveIds;
    const auto n = _x.size();
    for(std::size_t i = 0; i < n; ++i)
    {
        if(_cType[i] == CollisionType::World && _stdGravParam[i] > 0.) //TEMP - should go off of `_gType[i]`
        {
            massiveIds.push_back(i);
        }
    }
    _xPrev = _x;
    _yPrev = _y;
    _zPrev = _z;
    _xAttPrev = _xAtt;
    _yAttPrev = _yAtt;
    _zAttPrev = _zAtt;
    for(std::size_t i = 0; i < n; ++i)
    {
        if(_gType[i] == GravityType::Default)
        {
            double dx, dy, dz, du, dv, dw;
            grav_ode(_x[i], _y[i], _z[i], _xVel[i], _yVel[i], _zVel[i], dx, dy,
                dz, du, dv, dw, massiveIds, gravity);
            const double k1_1 = timestep*dx;
            const double k1_2 = timestep*dy;
            const double k1_3 = timestep*dz;
            const double k1_4 = timestep*du;
            const double k1_5 = timestep*dv;
            const double k1_6 = timestep*dw;
            grav_ode(_x[i] + 0.5*k1_1, _y[i] + 0.5*k1_2, _z[i] + 0.5*k1_3,
                _xVel[i] + 0.5*k1_4, _yVel[i] + 0.5*k1_5, _zVel[i] + 0.5*k1_6,
                dx, dy, dz, du, dv, dw, massiveIds, gravity);
            const double k2_1 = timestep*dx;
            const double k2_2 = timestep*dy;
            const double k2_3 = timestep*dz;
            const double k2_4 = timestep*du;
            const double k2_5 = timestep*dv;
            const double k2_6 = timestep*dw;
            grav_ode(_x[i] + 0.5*k2_1, _y[i] + 0.5*k2_2, _z[i] + 0.5*k2_3,
                _xVel[i] + 0.5*k2_4, _yVel[i] + 0.5*k2_5, _zVel[i] + 0.5*k2_6,
                dx, dy, dz, du, dv, dw, massiveIds, gravity);
            const double k3_1 = timestep*dx;
            const double k3_2 = timestep*dy;
            const double k3_3 = timestep*dz;
            const double k3_4 = timestep*du;
            const double k3_5 = timestep*dv;
            const double k3_6 = timestep*dw;
            grav_ode(_x[i] + k3_1, _y[i] + k3_2, _z[i] + k3_3, _xVel[i] + k3_4,
                _yVel[i] + k3_5, _zVel[i] + k3_6, dx, dy, dz, du, dv, dw,
                massiveIds, gravity);
            const double k4_1 = timestep*dx;
            const double k4_2 = timestep*dy;
            const double k4_3 = timestep*dz;
            const double k4_4 = timestep*du;
            const double k4_5 = timestep*dv;
            const double k4_6 = timestep*dw;
            static constexpr double c = 1./6.;
            _x[i] += c*(k1_1 + 2.*(k2_1 + k3_1) + k4_1);
            _y[i] += c*(k1_2 + 2.*(k2_2 + k3_2) + k4_2);
            _z[i] += c*(k1_3 + 2.*(k2_3 + k3_3) + k4_3);
            _xDown[i] = c*(k1_4 + 2.*(k2_4 + k3_4) + k4_4);
            _yDown[i] = c*(k1_5 + 2.*(k2_5 + k3_5) + k4_5);
            _zDown[i] = c*(k1_6 + 2.*(k2_6 + k3_6) + k4_6);
            _xVel[i] += _xDown[i];
            _yVel[i] += _yDown[i];
            _zVel[i] += _zDown[i];
        }
        else
        {
            _x[i] += timestep*_xVel[i];
            _y[i] += timestep*_yVel[i];
            _z[i] += timestep*_zVel[i];
            double dx, dy, dz, du, dv, dw;
            grav_ode(_x[i], _y[i], _z[i], _xVel[i], _yVel[i], _zVel[i], dx, dy,
                dz, du, dv, dw, massiveIds, gravity);
            _xDown[i] = du;
            _yDown[i] = dv;
            _zDown[i] = dw;
        }
        const double normSq = _xDown[i]*_xDown[i] + _yDown[i]*_yDown[i] +
            _zDown[i]*_zDown[i];
        if(normSq > 1e-6*1e-6)
        {
            const double invNorm = 1./std::sqrt(normSq);
            _xDown[i] *= invNorm;
            _yDown[i] *= invNorm;
            _zDown[i] *= invNorm;
        }
        else
        {
            _xDown[i] = 0.;
            _yDown[i] = 0.;
            _zDown[i] = 1.;
        }
    }
    for(std::size_t i = 0; i < n; ++i)
    {
        double WAtt = std::sqrt(1. - _xAtt[i]*_xAtt[i] - _yAtt[i]*_yAtt[i] -
            _zAtt[i]*_zAtt[i]);
        const double deltaX = normalize_angle(0.5*timestep*_xAngRate[i] + Pi) -
            Pi;
        const double deltaY = normalize_angle(0.5*timestep*_yAngRate[i] + Pi) -
            Pi;
        const double deltaZ = normalize_angle(0.5*timestep*_zAngRate[i] + Pi) -
            Pi;
        _xAtt[i] +=  WAtt   *deltaX - _zAtt[i]*deltaY + _yAtt[i]*deltaZ;
        _yAtt[i] +=  _zAtt[i]*deltaX + WAtt   *deltaY - _xAtt[i]*deltaZ;
        _zAtt[i] += -_yAtt[i]*deltaX + _xAtt[i]*deltaY + WAtt   *deltaZ;
        WAtt    += -_xAtt[i]*deltaX - _yAtt[i]*deltaY - _zAtt[i]*deltaZ;
        const double invSignNorm = (WAtt < 0. ? -1. : 1.)/std::sqrt(_xAtt[i]*
            _xAtt[i] + _yAtt[i]*_yAtt[i] + _zAtt[i]*_zAtt[i] + WAtt*WAtt);
        _xAtt[i] *= invSignNorm;
        _yAtt[i] *= invSignNorm;
        _zAtt[i] *= invSignNorm;
    }
}

void paz::do_collisions(Threadpool& threads, double timestep)
{
    // Identify all objects that can collide and precompute as much as possible.
    std::vector<std::size_t> a;
    std::vector<std::size_t> b;
    a.reserve(objects().size());
    b.reserve(objects().size());
    for(const auto& n : objects())
    {
        Object* o = reinterpret_cast<Object*>(n.first);
        if(o->collisionType() == CollisionType::Default)
        {
            a.push_back(n.second);
        }
        else if(o->collisionType() == CollisionType::World)
        {
            b.push_back(n.second);
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
        const double wAtt = std::sqrt(1. - _xAtt[b[i]]*_xAtt[b[i]] - _yAtt[b[
            i]]*_yAtt[b[i]] - _zAtt[b[i]]*_zAtt[b[i]]);
        const Vec att{{_xAtt[b[i]], _yAtt[b[i]], _zAtt[b[i]], wAtt}};
        bRot[i] = to_mat(att);
        for(std::size_t j = 0; j < NumSteps; ++j)
        {
            bX[i][j] = _xPrev[b[i]] + times[j]*(_x[b[i]] - _xPrev[b[i]]);
            bY[i][j] = _yPrev[b[i]] + times[j]*(_y[b[i]] - _yPrev[b[i]]);
            bZ[i][j] = _zPrev[b[i]] + times[j]*(_z[b[i]] - _zPrev[b[i]]);
        }
    }

    // Find and handle collisions.
    std::vector<LockedCv> lcvs;
    lcvs.reserve(a.size());
    for(std::size_t i = 0; i < a.size(); ++i)
    {
        lcvs.push_back(threads.pushTask([=]()
        {
            // `a[i]` may collide with any `b[c[j].first]`'s triangles in
            // `c[j].second`.
            std::unordered_map<std::size_t, std::vector<std::size_t>> c;
            for(std::size_t j = 0; j < b.size(); ++j)
            {
                Vec relPosPrev{{_xPrev[a[i]] - _xPrev[b[j]], _yPrev[a[i]] -
                    _yPrev[b[j]], _zPrev[a[i]] - _zPrev[b[j]]}};
                relPosPrev = bRot[j]*relPosPrev;
                Vec relPos{{_x[a[i]] - _x[b[j]], _y[a[i]] - _y[b[j]], _z[a[i]] -
                    _z[b[j]]}};
                relPos = bRot[j]*relPos;
                const auto temp = GET_CMESH(b[j]).sweepVol(relPosPrev, relPos,
                    _cRadius[a[i]]);
                if(!temp.empty())
                {
                    c[j] = std::move(temp);
                }
            }

            for(std::size_t j = 0; j < NumSteps; ++j)
            {
                double x = _xPrev[a[i]] + times[j]*(_x[a[i]] - _xPrev[a[i]]);
                double y = _yPrev[a[i]] + times[j]*(_y[a[i]] - _yPrev[a[i]]);
                double z = _zPrev[a[i]] + times[j]*(_z[a[i]] - _zPrev[a[i]]);
                double minDist = std::numeric_limits<double>::infinity();
                std::size_t idx = 0;
                double xNor = 0.;
                double yNor = 0.;
                double zNor = 1.;
                std::vector<std::pair<std::size_t, std::array<double, 3>>>
                    collisions;
                for(auto n : c)
                {
                    const double x1 = x - bX[n.first][j];
                    const double y1 = y - bY[n.first][j];
                    const double z1 = z - bZ[n.first][j];
                    const Vec relPos = bRot[n.first]*Vec{{x1, y1, z1}};
                    const double x2 = relPos(0);
                    const double y2 = relPos(1);
                    const double z2 = relPos(2);
                    double xNew, yNew, zNew, xNorTemp, yNorTemp, zNorTemp;
                    const double dist = GET_CMESH(b[n.first]).collide(x2, y2,
                        z2, _cRadius[a[i]], xNew, yNew, zNew, xNorTemp,
                        yNorTemp, zNorTemp, n.second);
                    if(dist < _cRadius[a[i]])
                    {
                        collisions.emplace_back(n.first, std::array<double, 3>{
                            xNorTemp, yNorTemp, zNorTemp});
                        const Vec newPos = bRot[n.first].trans()*Vec{{xNew,
                            yNew, zNew}};
                        xNew = newPos(0);
                        yNew = newPos(1);
                        zNew = newPos(2);
                        x += xNew - x1;
                        y += yNew - y1;
                        z += zNew - z1;
                        if(dist < minDist)
                        {
                            minDist = dist;
                            idx = n.first;
                            const Vec nor = bRot[n.first].trans()*Vec{{xNorTemp,
                                yNorTemp, zNorTemp}};
                            xNor = nor(0);
                            yNor = nor(1);
                            zNor = nor(2);
                        }
                    }
                }
                if(std::isfinite(minDist))
                {
                    const double xVel = _xVel[a[i]] - _xVel[b[idx]];
                    const double yVel = _yVel[a[i]] - _yVel[b[idx]];
                    const double zVel = _zVel[a[i]] - _zVel[b[idx]];
                    const double norVel = xVel*xNor + yVel*yNor + zVel*zNor;
                    if(norVel < 0.)
                    {
                        _xVel[a[i]] -= norVel*xNor;
                        _yVel[a[i]] -= norVel*yNor;
                        _zVel[a[i]] -= norVel*zNor;
                    }

                    // Adjust positions to fit new trajectory.
                    const double time0 = timestep*(times[j] - times[0]);
                    _xPrev[a[i]] = x - time0*_xVel[a[i]];
                    _yPrev[a[i]] = y - time0*_yVel[a[i]];
                    _zPrev[a[i]] = z - time0*_zVel[a[i]];
                    const double time1 = timestep*(times.back() - times[j]);
                    _x[a[i]] = x + time1*_xVel[a[i]];
                    _y[a[i]] = y + time1*_yVel[a[i]];
                    _z[a[i]] = z + time1*_zVel[a[i]];

                    // Collision response.
                    std::swap(_x[a[i]], x);
                    std::swap(_y[a[i]], y);
                    std::swap(_z[a[i]], z);
                    Object& aObj = *reinterpret_cast<Object*>(_ids[a[i]]);
                    for(const auto& n : collisions)
                    {
                        Object& bObj = *reinterpret_cast<Object*>(_ids[b[n.
                            first]]);
                        aObj.onCollide(bObj, n.second[0], n.second[1], n.second[
                            2], bX[n.first][j], bY[n.first][j], bZ[n.first][j]);
                    }
                    std::swap(_x[a[i]], x);
                    std::swap(_y[a[i]], y);
                    std::swap(_z[a[i]], z);

                    // Adjust positions again.
                    _xPrev[a[i]] = x - time0*_xVel[a[i]];
                    _yPrev[a[i]] = y - time0*_yVel[a[i]];
                    _zPrev[a[i]] = z - time0*_zVel[a[i]];
                    _x[a[i]] = x + time1*_xVel[a[i]];
                    _y[a[i]] = y + time1*_yVel[a[i]];
                    _z[a[i]] = z + time1*_zVel[a[i]];
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
    objects()[_id] = _x.size();
    _ids.push_back(_id);
    _x.push_back(0.);
    _y.push_back(0.);
    _z.push_back(0.);
    _xVel.push_back(0.);
    _yVel.push_back(0.);
    _zVel.push_back(0.);
    _xAtt.push_back(0.);
    _yAtt.push_back(0.);
    _zAtt.push_back(0.);
    _xAngRate.push_back(0.);
    _yAngRate.push_back(0.);
    _zAngRate.push_back(0.);
    _mod.emplace_back();
    _cMesh.emplace_back();
    _cType.push_back(CollisionType::Default);
    _gType.push_back(GravityType::Default);
    _xPrev.push_back(std::nan(""));
    _yPrev.push_back(std::nan(""));
    _zPrev.push_back(std::nan(""));
    _xAttPrev.push_back(std::nan(""));
    _yAttPrev.push_back(std::nan(""));
    _zAttPrev.push_back(std::nan(""));
    _cRadius.push_back(0.2);
    _xDown.push_back(0.);
    _yDown.push_back(0.);
    _zDown.push_back(0.);
    _stdGravParam.push_back(0.);
    _lights.emplace_back();
}

paz::Object::Object(const Object& o) : _id(reinterpret_cast<std::uintptr_t>(
    this))
{
    const auto otherIdx = objects().at(o._id);
    objects()[_id] = _x.size();
    _ids.push_back(_id);
    PUSH_COPY(_x)
    PUSH_COPY(_y)
    PUSH_COPY(_z)
    PUSH_COPY(_xVel)
    PUSH_COPY(_yVel)
    PUSH_COPY(_zVel)
    PUSH_COPY(_xAtt)
    PUSH_COPY(_yAtt)
    PUSH_COPY(_zAtt)
    PUSH_COPY(_xAngRate)
    PUSH_COPY(_yAngRate)
    PUSH_COPY(_zAngRate)
    PUSH_COPY(_mod)
    PUSH_COPY(_cMesh)
    PUSH_COPY(_cType)
    PUSH_COPY(_gType)
    PUSH_COPY(_xPrev);
    PUSH_COPY(_yPrev);
    PUSH_COPY(_zPrev);
    PUSH_COPY(_xAttPrev);
    PUSH_COPY(_yAttPrev);
    PUSH_COPY(_zAttPrev);
    PUSH_COPY(_cRadius)
    PUSH_COPY(_xDown)
    PUSH_COPY(_yDown)
    PUSH_COPY(_zDown)
    PUSH_COPY(_stdGravParam)
    PUSH_COPY(_lights)
    for(auto& n : _objectsByTag) //TEMP - inefficient
    {
        if(n.second.count(o._id))
        {
            n.second.insert(_id);
        }
    }
}

paz::Object& paz::Object::operator=(const Object& o)
{
    // If source and destination are the same, do nothing.
    if(_id == o._id)
    {
        return *this;
    }

    // If destination is in a valid state, copy data, otherwise, copy construct.
    const auto otherIdx = objects().at(o._id);
    if(objects().count(_id))
    {
        const auto idx = objects().at(_id);
        COPY(_x)
        COPY(_y)
        COPY(_z)
        COPY(_xVel)
        COPY(_yVel)
        COPY(_zVel)
        COPY(_xAtt)
        COPY(_yAtt)
        COPY(_zAtt)
        COPY(_xAngRate)
        COPY(_yAngRate)
        COPY(_zAngRate)
        COPY(_mod)
        COPY(_cMesh)
        COPY(_cType)
        COPY(_gType)
        COPY(_xPrev);
        COPY(_yPrev);
        COPY(_zPrev);
        COPY(_xAttPrev);
        COPY(_yAttPrev);
        COPY(_zAttPrev);
        COPY(_cRadius)
        COPY(_xDown)
        COPY(_yDown)
        COPY(_zDown)
        COPY(_stdGravParam)
        COPY(_lights)
        for(auto& n : _objectsByTag) //TEMP - inefficient
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
        objects()[_id] = _x.size();
        _ids.push_back(_id);
        PUSH_COPY(_x)
        PUSH_COPY(_y)
        PUSH_COPY(_z)
        PUSH_COPY(_xVel)
        PUSH_COPY(_yVel)
        PUSH_COPY(_zVel)
        PUSH_COPY(_xAtt)
        PUSH_COPY(_yAtt)
        PUSH_COPY(_zAtt)
        PUSH_COPY(_xAngRate)
        PUSH_COPY(_yAngRate)
        PUSH_COPY(_zAngRate)
        PUSH_COPY(_mod)
        PUSH_COPY(_cMesh)
        PUSH_COPY(_cType)
        PUSH_COPY(_gType)
        PUSH_COPY(_xPrev);
        PUSH_COPY(_yPrev);
        PUSH_COPY(_zPrev);
        PUSH_COPY(_xAttPrev);
        PUSH_COPY(_yAttPrev);
        PUSH_COPY(_zAttPrev);
        PUSH_COPY(_cRadius)
        PUSH_COPY(_xDown)
        PUSH_COPY(_yDown)
        PUSH_COPY(_zDown)
        PUSH_COPY(_stdGravParam)
        PUSH_COPY(_lights)
        for(auto& n : _objectsByTag) //TEMP - inefficient
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
    _ids[idx] = _id;
    for(auto& n : _objectsByTag) //TEMP - inefficient
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
    // If source and destination are the same, do nothing.
    if(_id == o._id)
    {
        return *this;
    }

    // If destination is in a valid state, just swap, otherwise, copy construct.
    const auto otherIdx = objects().at(o._id);
    if(objects().count(_id))
    {
        const auto idx = objects().at(_id);
        objects().at(_id) = otherIdx;
        objects().at(o._id) = idx;
        std::swap(_ids[idx], _ids[otherIdx]);
    }
    else
    {
        objects()[_id] = _x.size();
        _ids.push_back(_id);
        PUSH_COPY(_x)
        PUSH_COPY(_y)
        PUSH_COPY(_z)
        PUSH_COPY(_xVel)
        PUSH_COPY(_yVel)
        PUSH_COPY(_zVel)
        PUSH_COPY(_xAtt)
        PUSH_COPY(_yAtt)
        PUSH_COPY(_zAtt)
        PUSH_COPY(_xAngRate)
        PUSH_COPY(_yAngRate)
        PUSH_COPY(_zAngRate)
        PUSH_COPY(_mod)
        PUSH_COPY(_cMesh)
        PUSH_COPY(_cType)
        PUSH_COPY(_gType)
        PUSH_COPY(_xPrev);
        PUSH_COPY(_yPrev);
        PUSH_COPY(_zPrev);
        PUSH_COPY(_xAttPrev);
        PUSH_COPY(_yAttPrev);
        PUSH_COPY(_zAttPrev);
        PUSH_COPY(_cRadius)
        PUSH_COPY(_xDown)
        PUSH_COPY(_yDown)
        PUSH_COPY(_zDown)
        PUSH_COPY(_stdGravParam)
        PUSH_COPY(_lights)
        for(auto& n : _objectsByTag) //TEMP - inefficient
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
    if(_ids.size() > 1)
    {
        objects().at(_ids.back()) = idx;
    }
    objects().erase(_id);
    SWAP_AND_POP(_ids)
    SWAP_AND_POP(_x)
    SWAP_AND_POP(_y)
    SWAP_AND_POP(_z)
    SWAP_AND_POP(_xVel)
    SWAP_AND_POP(_yVel)
    SWAP_AND_POP(_zVel)
    SWAP_AND_POP(_xAtt)
    SWAP_AND_POP(_yAtt)
    SWAP_AND_POP(_zAtt)
    SWAP_AND_POP(_xAngRate)
    SWAP_AND_POP(_yAngRate)
    SWAP_AND_POP(_zAngRate)
    SWAP_AND_POP(_mod)
    SWAP_AND_POP(_cMesh)
    SWAP_AND_POP(_cType)
    SWAP_AND_POP(_gType)
    SWAP_AND_POP(_xPrev);
    SWAP_AND_POP(_yPrev);
    SWAP_AND_POP(_zPrev);
    SWAP_AND_POP(_xAttPrev);
    SWAP_AND_POP(_yAttPrev);
    SWAP_AND_POP(_zAttPrev);
    SWAP_AND_POP(_cRadius)
    SWAP_AND_POP(_xDown)
    SWAP_AND_POP(_yDown)
    SWAP_AND_POP(_zDown)
    SWAP_AND_POP(_stdGravParam)
    SWAP_AND_POP(_lights)
    for(auto& n : _objectsByTag) //TEMP - inefficient
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
    if(!_objectsByTag.count(tag))
    {
        return;
    }
    for(auto n : _objectsByTag.at(tag))
    {
        reinterpret_cast<Object*>(n)->onNotify(*this, data);
    }
}

double& paz::Object::x()
{
    return _x[objects().at(_id)];
}

double paz::Object::x() const
{
    return _x[objects().at(_id)];
}

double& paz::Object::y()
{
    return _y[objects().at(_id)];
}

double paz::Object::y() const
{
    return _y[objects().at(_id)];
}

double& paz::Object::z()
{
    return _z[objects().at(_id)];
}

double paz::Object::z() const
{
    return _z[objects().at(_id)];
}

double& paz::Object::xVel()
{
    return _xVel[objects().at(_id)];
}

double paz::Object::xVel() const
{
    return _xVel[objects().at(_id)];
}

double& paz::Object::yVel()
{
    return _yVel[objects().at(_id)];
}

double paz::Object::yVel() const
{
    return _yVel[objects().at(_id)];
}

double& paz::Object::zVel()
{
    return _zVel[objects().at(_id)];
}

double paz::Object::zVel() const
{
    return _zVel[objects().at(_id)];
}

double& paz::Object::xAtt()
{
    return _xAtt[objects().at(_id)];
}

double paz::Object::xAtt() const
{
    return _xAtt[objects().at(_id)];
}

double& paz::Object::yAtt()
{
    return _yAtt[objects().at(_id)];
}

double paz::Object::yAtt() const
{
    return _yAtt[objects().at(_id)];
}

double& paz::Object::zAtt()
{
    return _zAtt[objects().at(_id)];
}

double paz::Object::zAtt() const
{
    return _zAtt[objects().at(_id)];
}

double& paz::Object::xAngRate()
{
    return _xAngRate[objects().at(_id)];
}

double paz::Object::xAngRate() const
{
    return _xAngRate[objects().at(_id)];
}

double& paz::Object::yAngRate()
{
    return _yAngRate[objects().at(_id)];
}

double paz::Object::yAngRate() const
{
    return _yAngRate[objects().at(_id)];
}

double& paz::Object::zAngRate()
{
    return _zAngRate[objects().at(_id)];
}

double paz::Object::zAngRate() const
{
    return _zAngRate[objects().at(_id)];
}

paz::Model& paz::Object::model()
{
    return _mod[objects().at(_id)];
}

const paz::Model& paz::Object::model() const
{
    return _mod[objects().at(_id)];
}

paz::CollisionMesh& paz::Object::collisionMesh()
{
    return _cMesh[objects().at(_id)];
}

const paz::CollisionMesh& paz::Object::collisionMesh() const
{
    return _cMesh[objects().at(_id)];
}

paz::CollisionType& paz::Object::collisionType()
{
    return _cType[objects().at(_id)];
}

const paz::CollisionType& paz::Object::collisionType() const
{
    return _cType[objects().at(_id)];
}

paz::GravityType& paz::Object::gravityType()
{
    return _gType[objects().at(_id)];
}

const paz::GravityType& paz::Object::gravityType() const
{
    return _gType[objects().at(_id)];
}

double paz::Object::xPrev() const
{
    return _xPrev[objects().at(_id)];
}

double paz::Object::yPrev() const
{
    return _yPrev[objects().at(_id)];
}

double paz::Object::zPrev() const
{
    return _zPrev[objects().at(_id)];
}

double paz::Object::xAttPrev() const
{
    return _xAttPrev[objects().at(_id)];
}

double paz::Object::yAttPrev() const
{
    return _yAttPrev[objects().at(_id)];
}

double paz::Object::zAttPrev() const
{
    return _zAttPrev[objects().at(_id)];
}

double& paz::Object::collisionRadius()
{
    return _cRadius[objects().at(_id)];
}

double paz::Object::collisionRadius() const
{
    return _cRadius[objects().at(_id)];
}

double paz::Object::xDown() const
{
    return _xDown[objects().at(_id)];
}

double paz::Object::yDown() const
{
    return _yDown[objects().at(_id)];
}

double paz::Object::zDown() const
{
    return _zDown[objects().at(_id)];
}

double& paz::Object::stdGravParam()
{
    return _stdGravParam[objects().at(_id)];
}

double paz::Object::stdGravParam() const
{
    return _stdGravParam[objects().at(_id)];
}

std::vector<std::array<double, 7>>& paz::Object::lights()
{
    return _lights[objects().at(_id)];
}

const std::vector<std::array<double, 7>>& paz::Object::lights() const
{
    return _lights[objects().at(_id)];
}

void paz::Object::addTag(const std::string& tag)
{
    _objectsByTag[tag].insert(_id);
}

bool paz::Object::isTagged(const std::string& tag) const
{
    return _objectsByTag.count(tag) && _objectsByTag.at(tag).count(_id);
}

void paz::Object::computeAltitude(double& alt, Vec& nor, Vec& vel) const
{
    nor = {{0, 0, 1}};
    vel = Vec::Zero(3);
    const auto idx = objects().at(_id);
    alt = std::numeric_limits<double>::infinity();
    for(auto n : objects())
    {
        if(n.first == _id || _cType[n.second] != CollisionType::World)
        {
            continue;
        }
        const double wAtt = std::sqrt(1. - _xAtt[n.second]*_xAtt[n.second] -
            _yAtt[n.second]*_yAtt[n.second] - _zAtt[n.second]*_zAtt[n.second]);
        const Mat rot = to_mat(Vec{{_xAtt[n.second], _yAtt[n.second], _zAtt[n.
            second], wAtt}});
        const Vec relPos = rot*Vec{{_x[idx] - _x[n.second], _y[idx] - _y[n.
            second], _z[idx] - _z[n.second]}};
        const Vec dir = rot*Vec{{_xDown[idx], _yDown[idx], _zDown[idx]}};
        double dist;
        Vec tempNor(3);
        GET_CMESH(n.second).castRay(relPos(0), relPos(1), relPos(2), dir(0),
            dir(1), dir(2), tempNor(0), tempNor(1), tempNor(2), dist);
        if(dist < alt)
        {
            alt = dist;
            nor = rot.trans()*tempNor;
            vel(0) = _xVel[n.second];
            vel(1) = _yVel[n.second];
            vel(2) = _zVel[n.second];
        }
    }
    alt -= _cRadius[idx];
}
