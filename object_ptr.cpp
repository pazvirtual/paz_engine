#include "PAZ_Engine"
#include "object.hpp"

//TEMP - VECTOR REALLOC CHANGES IDS/PTRS (NEED MODIFY OBJ'S MOVE CTOR & DELETE OBJ'S COPY CTOR TO FIX)

#define OBJ_EXISTS objects().count(_id) //TEMP - wrong if new object has same ID

static const auto NullId = reinterpret_cast<std::uintptr_t>(nullptr);

paz::ObjectPtr::ObjectPtr() : _id(NullId) {}

paz::ObjectPtr::ObjectPtr(const Object& o) : _id(reinterpret_cast<std::
    uintptr_t>(&o)) {}

paz::ObjectPtr::ObjectPtr(std::nullptr_t) : ObjectPtr() {}

void paz::ObjectPtr::swap(ObjectPtr& p) noexcept
{
    std::swap(_id, p._id);
}

const paz::Object* paz::ObjectPtr::get() const
{
    if(_id == NullId)
    {
        throw std::runtime_error("Object pointer is null.");
    }
    if(!OBJ_EXISTS)
    {
        throw std::runtime_error("Object has been deleted.");
    }
    return reinterpret_cast<const Object*>(_id);
}

paz::Object* paz::ObjectPtr::get()
{
    if(_id == NullId)
    {
        throw std::runtime_error("Object pointer is null.");
    }
    if(!OBJ_EXISTS)
    {
        throw std::runtime_error("Object has been deleted.");
    }
    return reinterpret_cast<Object*>(_id);
}

paz::ObjectPtr& paz::ObjectPtr::operator=(std::nullptr_t)
{
    _id = NullId;
    return *this;
}

const paz::Object& paz::ObjectPtr::operator*() const
{
    return *get();
}

paz::Object& paz::ObjectPtr::operator*()
{
    return *get();
}

const paz::Object* paz::ObjectPtr::operator->() const
{
    return get();
}

paz::Object* paz::ObjectPtr::operator->()
{
    return get();
}

paz::ObjectPtr::operator bool() const
{
    return _id != NullId && OBJ_EXISTS;
}

bool paz::ObjectPtr::operator==(const ObjectPtr& p) const
{
    return _id == p._id;
}

bool paz::ObjectPtr::operator!=(const ObjectPtr& p) const
{
    return !(*this == p);
}

void paz::ObjectPtr::reset() noexcept
{
    _id = NullId;
}

void paz::ObjectPtr::reset(const Object& o) noexcept
{
    _id = reinterpret_cast<std::uintptr_t>(&o);
}

std::ostream& paz::operator<<(std::ostream& stream, const paz::ObjectPtr& p)
{
    stream << p.get();
    return stream;
}
