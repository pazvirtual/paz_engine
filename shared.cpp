#include "shared.hpp"

extern const unsigned char binary_assets_paz_start[];
extern const unsigned char binary_assets_paz_end[];

paz::Bytes paz::getAsset(const std::string& name)
{
    static const paz::Archive assets(paz::Bytes(binary_assets_paz_start,
        binary_assets_paz_end));
    return assets.get(name);
}
