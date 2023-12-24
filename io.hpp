#ifndef PAZ_ENGINE_IO_HPP
#define PAZ_ENGINE_IO_HPP

#include <string>

namespace paz
{
    void save_setting(const std::string& name, const std::string& val);
    std::string load_setting(const std::string& name);
}

#endif
