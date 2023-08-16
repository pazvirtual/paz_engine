#include "shared.hpp"
#include <regex>

extern const unsigned char pazAssetsStart[] asm("_paz_binary_assets_paz_start");
extern const unsigned char pazAssetsEnd[] asm("_paz_binary_assets_paz_end");
extern const unsigned char assetsStart[] asm("_binary_assets_paz_start");
extern const unsigned char assetsEnd[] asm("_binary_assets_paz_end");

static const paz::Archive& builtin_assets()
{
    static const paz::Archive assets(paz::Bytes(pazAssetsStart, pazAssetsEnd));
    return assets;
}

static const paz::Archive& user_assets()
{
    static const paz::Archive assets(paz::Bytes(assetsStart, assetsEnd));
    return assets;
}

static paz::Image get_img_internal(const paz::Archive& archive, const std::
    string& path)
{
    static const std::regex bmp("bmp", std::regex_constants::icase);
    static const std::regex pbm("pbm", std::regex_constants::icase);
    //static const std::regex jpg("jpe?g", std::regex_constants::icase);
    //static const std::regex png("png", std::regex_constants::icase);
    const std::string ext = paz::split_path(path)[2];
    if(std::regex_match(ext, bmp))
    {
        return paz::parse_bmp(archive.get(path));
    }
    if(std::regex_match(ext, pbm))
    {
        return paz::parse_pbm(archive.get(path));
    }
    //if(std::regex_match(ext, jpg))
    //{
    //    return paz::parse_jpg(archive.get(path));
    //}
    //if(std::regex_match(ext, png))
    //{
    //    return paz::parse_png(archive.get(path));
    //}
    throw std::runtime_error("Unrecognized image extension \"" + ext + "\".");
}

paz::Bytes paz::get_builtin(const std::string& name)
{
    return builtin_assets().get(name);
}

paz::Bytes paz::get_asset(const std::string& name)
{
    return user_assets().get(name);
}

paz::Image paz::get_builtin_image(const std::string& path)
{
    return get_img_internal(builtin_assets(), path);
}

paz::Image paz::get_asset_image(const std::string& path)
{
    return get_img_internal(user_assets(), path);
}
