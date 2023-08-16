#include "PAZ_Engine"
#include "object.hpp"

//TEMP - VECTOR REALLOC CHANGES IDS/PTRS (NEED ADD `OBJECT` MOVE CTOR TO FIX)

#define OBJ_EXISTS objects().count(_id) //TEMP - fails if new object has same ID

paz::ObjectPtr::ObjectPtr() : _set(false), _id(0) {}

paz::ObjectPtr::ObjectPtr(const Object& o) : _set(true), _id(reinterpret_cast<
    std::uintptr_t>(&o)) {}

paz::ObjectPtr::ObjectPtr(std::nullptr_t) : ObjectPtr() {}

paz::ObjectPtr& paz::ObjectPtr::operator=(std::nullptr_t)
{
    ObjectPtr().swap(*this);
    return *this;
}

const paz::Object& paz::ObjectPtr::operator*() const
{
    if(!_set)
    {
        throw std::runtime_error("Object pointer is null.");
    }
    if(!OBJ_EXISTS)
    {
        throw std::runtime_error("Object has been deleted.");
    }
    return *reinterpret_cast<const Object*>(_id);
}

paz::Object& paz::ObjectPtr::operator*()
{
    if(!_set)
    {
        throw std::runtime_error("Object pointer is null.");
    }
    if(!OBJ_EXISTS)
    {
        throw std::runtime_error("Object has been deleted.");
    }
    return *reinterpret_cast<Object*>(_id);
}

const paz::Object* paz::ObjectPtr::operator->() const
{
    if(!_set)
    {
        throw std::runtime_error("Object pointer is null.");
    }
    if(!OBJ_EXISTS)
    {
        throw std::runtime_error("Object has been deleted.");
    }
    return reinterpret_cast<const Object*>(_id);
}

paz::Object* paz::ObjectPtr::operator->()
{
    if(!_set)
    {
        throw std::runtime_error("Object pointer is null.");
    }
    if(!OBJ_EXISTS)
    {
        throw std::runtime_error("Object has been deleted.");
    }
    return reinterpret_cast<Object*>(_id);
}

paz::ObjectPtr::operator bool() const
{
    return _set && OBJ_EXISTS;
}

void paz::ObjectPtr::reset() noexcept
{
    _set = false;
}

void paz::ObjectPtr::reset(const Object& o) noexcept
{
    _set = true;
    _id = reinterpret_cast<std::uintptr_t>(&o);
}

void paz::ObjectPtr::swap(ObjectPtr& p) noexcept
{
    std::swap(_set, p._set);
    std::swap(_id, p._id);
}

std::ostream& paz::operator<<(std::ostream& stream, const paz::ObjectPtr& p)
{
    if(p._set)
    {
        std::ostringstream oss;
        oss << std::hex << p._id;
        stream << oss.str();
    }
    else
    {
        stream << 0;
    }
    return stream;
}
