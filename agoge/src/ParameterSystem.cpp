#include "ParameterSystem.hpp"
#include <sstream>
#include <cctype>
#include <iostream>

namespace agoge {

ParameterSystem::ParameterSystem()
{
    setDefaults();
}

void ParameterSystem::setDefaults()
{
    defaults_["nx"] = "64";
    defaults_["ny"] = "64";
    defaults_["nz"] = "64";

    defaults_["cfl"] = "0.5";
    defaults_["use_gravity"] = "false";
    defaults_["domain"] = "[1.0, 1.0, 1.0]";
    // Remove old nsteps
    // Add new "sound_crossings"
    defaults_["sound_crossings"] = "1.0";  // 1 crossing by default

    defaults_["problem_name"] = "\"sod\"";
}

bool ParameterSystem::readYAML(const std::string &filename)
{
    bool ok = parser_.parseFile(filename);
    if(!ok) {
        std::cerr << "[ParameterSystem] Could not parse YAML: " << filename << "\n";
    }
    return ok;
}

// Everything else unchanged: getRaw, getInt, getDouble, etc.

std::string ParameterSystem::getRaw(const std::string &key) const
{
    if(parser_.hasKey(key)) {
        return parser_.getString(key);
    }
    auto it = defaults_.find(key);
    if(it != defaults_.end()) {
        return it->second;
    }
    return "";
}

int ParameterSystem::getInt(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return 0;
    }
    return std::stoi(raw);
}

double ParameterSystem::getDouble(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return 0.0;
    }
    return std::stod(raw);
}

bool ParameterSystem::getBool(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return false;
    }
    std::string lower;
    for(char c: raw) lower.push_back(std::tolower(static_cast<unsigned char>(c)));
    if(lower=="true" || lower=="on" || lower=="1" || lower=="yes") {
        return true;
    }
    return false;
}

std::string ParameterSystem::getString(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(!raw.empty() && raw.front()=='"' && raw.back()=='"' && raw.size()>1) {
        return raw.substr(1, raw.size()-2);
    }
    return raw;
}

std::vector<double> ParameterSystem::getDoubleList(const std::string &key) const
{
    std::vector<double> result;
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return result;
    }
    auto tokens = parseList(raw);
    for(const auto &t: tokens) {
        try {
            double val = std::stod(t);
            result.push_back(val);
        } catch(...) {
            // skip error
        }
    }
    return result;
}

std::vector<std::string> ParameterSystem::parseList(const std::string &raw) const
{
    std::vector<std::string> tokens;
    std::string s = trim(raw);
    if(!s.empty() && s.front()=='[') {
        s.erase(s.begin());
    }
    if(!s.empty() && s.back()==']') {
        s.pop_back();
    }
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, ',')) {
        tokens.push_back(trim(item));
    }
    return tokens;
}

std::string ParameterSystem::trim(const std::string &str) const
{
    if(str.empty()) return str;
    size_t start=0;
    while(start<str.size() && std::isspace(static_cast<unsigned char>(str[start]))) {
        start++;
    }
    size_t end=str.size()-1;
    while(end>start && std::isspace(static_cast<unsigned char>(str[end]))) {
        end--;
    }
    return str.substr(start, end-start+1);
}

} // namespace agoge
