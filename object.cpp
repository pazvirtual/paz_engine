#include "object.hpp"
#include "PAZ_Engine"

std::unordered_map<std::uintptr_t, std::size_t>& paz::objects()
{
    static std::unordered_map<std::uintptr_t, std::size_t> o;
    return o;
}

static std::vector<double> X;
static std::vector<double> Y;
static std::vector<double> Z;
static std::vector<paz::VertexBuffer> V;
static std::vector<paz::IndexBuffer> I;

paz::Object::Object() : _id(reinterpret_cast<std::uintptr_t>(this))
{
    objects()[_id] = X.size();
    X.push_back(0);
    Y.push_back(0);
    Z.push_back(0);
    V.emplace_back();
    I.emplace_back();
}

paz::Object::~Object()
{
    const std::size_t idx = objects().at(_id);
    objects().erase(_id);
    std::swap(X[idx], X.back());
    X.pop_back();
    std::swap(Y[idx], Y.back());
    Y.pop_back();
    std::swap(Z[idx], Z.back());
    Z.pop_back();
    std::swap(V[idx], V.back());
    V.pop_back();
    std::swap(I[idx], I.back());
    I.pop_back();
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

paz::VertexBuffer paz::Object::v()
{
    return V[objects().at(_id)];
}

const paz::VertexBuffer paz::Object::v() const
{
    return V[objects().at(_id)];
}

paz::IndexBuffer paz::Object::i()
{
    return I[objects().at(_id)];
}

const paz::IndexBuffer paz::Object::i() const
{
    return I[objects().at(_id)];
}
