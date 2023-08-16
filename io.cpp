#ifdef __linux__

#include "io.hpp"
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
#include <regex>

static void ensure_dir(const std::string& path)
{
    std::filesystem::file_status status;
    try
    {
        status = std::filesystem::status(path);
    }
    catch(const std::exception& e)
    {
        throw std::runtime_error("Failed to get status for \"" + path + "\": " +
            e.what());
    }
    if(status.type() == std::filesystem::file_type::not_found)
    {
        try
        {
            std::filesystem::create_directory(path);
        }
        catch(const std::exception& e)
        {
            throw std::runtime_error("Failed to create directory \"" + path +
                "\": " + e.what());
        }
    }
    else if(status.type() != std::filesystem::file_type::directory)
    {
        throw std::runtime_error("Path \"" + path + "\" is not a directory.");
    }
}

static std::string get_home()
{
    static const char* homeDir = getenv("HOME");
    if(!homeDir)
    {
        throw std::runtime_error("Failed to find home directory.");
    }
    return homeDir;
}

static std::string ensure_settings_file(const std::string& company, const std::
    string& app)
{
    std::string path = get_home() + "/.local";
    ensure_dir(path);
    path += "/share";
    ensure_dir(path);
    path += "/" + company;
    ensure_dir(path);
    path += "/" + app;
    ensure_dir(path);
    return path + "/settings.txt";
}

void paz::save_setting(const std::string& name, const std::string& val)
{
    const auto path = ensure_settings_file("PAZ Virtual", "temp-demo");

    std::map<std::string, std::string> settings = {{name, val}};

    // Get all other settings.
    {
        std::ifstream in(path);
        if(!in)
        {
            throw std::runtime_error("Failed to open \"" + path + "\".");
        }
        std::string line;
        while(std::getline(in, line))
        {
            if(std::regex_match(line, std::regex("\\S+\\s+\\S+")))
            {
                const auto name = regex_replace(line, std::regex("\\s+\\S+"),
                    "");
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

    // Write all settings.
    std::ofstream out(path);
    if(!out)
    {
        throw std::runtime_error("Failed to open \"" + path + "\".");
    }
    for(const auto& n : settings)
    {
        out << n.first << ' ' << n.second << '\n';
    }
}

std::string paz::load_setting(const std::string& name)
{
    std::string path;
    try
    {
        path = get_home() + "/.local/share/PAZ Virtual/temp-demo/settings.txt";
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
