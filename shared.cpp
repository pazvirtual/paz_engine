#include "shared.hpp"

extern const unsigned char assetsStart[] asm("_binary_assets_paz_start");
extern const unsigned char assetsEnd[] asm("_binary_assets_paz_end");

paz::Bytes paz::getAsset(const std::string& name)
{
    static const paz::Archive assets(paz::Bytes(assetsStart, assetsEnd));
    return assets.get(name);
}
