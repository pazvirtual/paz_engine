#include "detect_os.hpp"

#ifdef PAZ_MACOS

#include "io.hpp"
#import <Foundation/Foundation.h>
#include <sstream>
#include <fstream>
#include <regex>
#include <map>

static const bool IsSandboxed = [[[NSProcessInfo processInfo] environment]
    objectForKey:@"APP_SANDBOX_CONTAINER_ID"] != nil;

static std::string get_home()
{
    if(!IsSandboxed)
    {
        throw std::runtime_error("App is not sandboxed.");
    }

    return [NSHomeDirectory() UTF8String];
}

void paz::save_setting(const std::string& name, const std::string& val)
{
    const std::string path = get_home() + "/settings.txt";

    std::map<std::string, std::string> settings = {{name, val}};

    // Get all other settings.
    {
        std::ifstream inFile(path);
        if(inFile)
        {
            std::string line;
            while(std::getline(inFile, line))
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
    std::ofstream outFile(path);
    if(!outFile)
    {
        throw std::runtime_error("Failed to open \"" + path + "\".");
    }
    for(const auto& n : settings)
    {
        outFile << n.first << ' ' << n.second << '\n';
    }
}

std::string paz::load_setting(const std::string& name)
{
    std::string path;
    try
    {
        path = get_home() + "/settings.txt";
    }
    catch(...)
    {
        return "";
    }
    std::ifstream inFile(path);
    if(!inFile)
    {
        return "";
    }
    std::string line;
    while(std::getline(inFile, line))
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
