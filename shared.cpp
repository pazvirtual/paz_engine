#include "shared.hpp"
#include <regex>

extern const unsigned char assetsStart[] asm("_binary_assets_paz_start");
extern const unsigned char assetsEnd[] asm("_binary_assets_paz_end");

paz::Bytes paz::get_asset(const std::string& name)
{
    static const Archive assets(Bytes(assetsStart, assetsEnd));
    return assets.get(name);
}

paz::Image paz::get_asset_image(const std::string& path)
{
    static const std::regex bmp("bmp", std::regex_constants::icase);
    static const std::regex pbm("pbm", std::regex_constants::icase);
    //static const std::regex jpg("jpe?g", std::regex_constants::icase);
    //static const std::regex png("png", std::regex_constants::icase);
    const std::string ext = split_path(path)[2];
    if(std::regex_match(ext, bmp))
    {
        return parse_bmp(get_asset(path));
    }
    if(std::regex_match(ext, pbm))
    {
        return parse_pbm(get_asset(path));
    }
    //if(std::regex_match(ext, jpg))
    //{
    //    return parse_jpg(get_asset(path));
    //}
    //if(std::regex_match(ext, png))
    //{
    //    return parse_png(get_asset(path));
    //}
    throw std::runtime_error("Unrecognized image extension \"" + ext + "\".");
}
