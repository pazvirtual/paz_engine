#ifndef PAZ_ENGINE_OBJECT_HPP
#define PAZ_ENGINE_OBJECT_HPP

#include "PAZ_Engine"
#include <unordered_map>

namespace paz
{
    std::unordered_map<std::uintptr_t, std::size_t>& objects();
}

#endif
