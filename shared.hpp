#ifndef PAZ_ENGINE_SHARED_HPP
#define PAZ_ENGINE_SHARED_HPP

#include "PAZ_IO"

namespace paz
{
    paz::Bytes get_asset(const std::string& path);
    paz::Image get_asset_image(const std::string& path);
}

#endif
