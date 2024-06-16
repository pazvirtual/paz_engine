#include "detect_os.hpp"

#ifndef PAZ_MACOS

#include "io.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
#include <regex>

static const std::string CompanyName = "PAZ Virtual"; //TEMP
static const std::string AppName = "temp-demo"; //TEMP

static void ensure_dir(const std::filesystem::path& path)
{
    std::filesystem::file_status status;
    try
    {
        status = std::filesystem::status(path);
    }
    catch(const std::exception& e)
    {
        throw std::runtime_error("Failed to get status for \"" + path.string() +
            "\": " + e.what());
    }
    if(status.type() == std::filesystem::file_type::not_found)
    {
        try
        {
            std::filesystem::create_directory(path);
        }
        catch(const std::exception& e)
        {
            throw std::runtime_error("Failed to create directory \"" + path.
                string() + "\": " + e.what());
        }
    }
    else if(status.type() != std::filesystem::file_type::directory)
    {
        throw std::runtime_error("Path \"" + path.string() + "\" is not a direc"
            "tory.");
    }
}

static std::filesystem::path get_home()
{
#ifdef PAZ_LINUX
    static const char* homeDir = std::getenv("HOME");
#else
    static const char* homeDir = std::getenv("APPDATA");
#endif
    if(!homeDir)
    {
        throw std::runtime_error("Failed to find home directory.");
    }
    return homeDir;
}

static std::filesystem::path ensure_settings_file(const std::string& company,
    const std::string& app)
{
    auto path = get_home();
#ifdef LINUX
    path /= std::filesystem::path(".local");
    ensure_dir(path);
    path /= "share";
    ensure_dir(path);
#endif
    path /= company;
    ensure_dir(path);
    path /= app;
    ensure_dir(path);
    path /= "settings.txt";
    return path;
}

void paz::save_setting(const std::string& name, const std::string& val)
{
    const auto path = ensure_settings_file(CompanyName, AppName);

    std::map<std::string, std::string> settings = {{name, val}};

    // Get all other settings.
    {
        std::ifstream in(path);
        if(in)
        {
            std::string line;
            while(std::getline(in, line))
            {
                if(std::regex_match(line, std::regex("\\S+\\s+\\S+")))
                {
                    const auto name = regex_replace(line, std::regex(
                        "\\s+\\S+"), "");
                    if(settings.count(name))
                    {
                        continue;
                    }
                    const auto val = regex_replace(line, std::regex("\\S+\\s+"),
                        "");
                    settings[name] = val;
                }
            }
        }
    }

    // Write all settings.
    std::ofstream out(path);
    if(!out)
    {
        throw std::runtime_error("Failed to open \"" + path.string() + "\".");
    }
    for(const auto& n : settings)
    {
        out << n.first << ' ' << n.second << '\n';
    }
}

std::string paz::load_setting(const std::string& name) noexcept
{
    std::filesystem::path path;
    try
    {
        path = get_home();
#ifdef PAZ_LINUX
        path /= ".local";
        path /= "share";
#endif
        path /= CompanyName;
        path /= AppName;
        path /= "settings.txt";
    }
    catch(...)
    {
        return "";
    }
    std::ifstream in(path);
    if(!in)
    {
        return "";
    }
    std::string line;
    while(std::getline(in, line))
    {
        if(line.size() >= name.size() + 2 && line.substr(0, name.size() + 1) ==
            name + ' ')
        {
            return line.substr(name.size() + 1);
        }
    }
    return "";
}

#endif
