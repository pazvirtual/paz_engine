#include "shared.hpp"

extern const unsigned char _binary_assets_paz_start[];
extern const unsigned char _binary_assets_paz_end[];

paz::Bytes paz::getAsset(const std::string& name)
{
    static const paz::Archive assets(paz::Bytes(_binary_assets_paz_start,
        _binary_assets_paz_end));
    return assets.get(name);
}
