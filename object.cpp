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
static std::vector<paz::VertexBuffer> V;
static std::vector<paz::IndexBuffer> I;
static std::vector<char> Vis;

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
    V.emplace_back();
    I.emplace_back();
    Vis.push_back(false);
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
    SWAP_AND_POP(V);
    SWAP_AND_POP(I);
    SWAP_AND_POP(Vis);
}

void paz::Object::update() {}

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

paz::VertexBuffer& paz::Object::v()
{
    return V[objects().at(_id)];
}

const paz::VertexBuffer& paz::Object::v() const
{
    return V[objects().at(_id)];
}

paz::IndexBuffer& paz::Object::i()
{
    return I[objects().at(_id)];
}

const paz::IndexBuffer& paz::Object::i() const
{
    return I[objects().at(_id)];
}

char& paz::Object::vis()
{
    return Vis[objects().at(_id)];
}

bool paz::Object::vis() const
{
    return Vis[objects().at(_id)];
}
