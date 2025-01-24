#include "ParameterSystem.hpp"
#include <sstream>
#include <cctype>
#include <iostream>

namespace agoge {

ParameterSystem::ParameterSystem()
{
    // define all built-in defaults here
    setDefaults();
}

void ParameterSystem::setDefaults()
{
    // You can add as many as you need
    // For demonstration:
    defaults_["nx"] = "64";
    defaults_["ny"] = "64";
    defaults_["nz"] = "64";

    defaults_["cfl"] = "0.5";
    defaults_["use_gravity"] = "false";
    defaults_["domain"] = "[1.0, 1.0, 1.0]";
    defaults_["nsteps"] = "500";

    // etc. e.g. problem name
    defaults_["problem_name"] = "\"sod\"";

    // Or any additional keys you want as defaults
}

bool ParameterSystem::readYAML(const std::string &filename)
{
    bool ok = parser_.parseFile(filename);
    if(!ok) {
        std::cerr << "[ParameterSystem] Could not parse YAML: " << filename << "\n";
    }
    return ok;
}

std::string ParameterSystem::getRaw(const std::string &key) const
{
    // 1) check parser override
    if(parser_.hasKey(key)) {
        return parser_.getString(key);
    }
    // 2) else fallback to defaults
    auto it = defaults_.find(key);
    if(it != defaults_.end()) {
        return it->second;
    }
    // not found
    return "";
}

// Convert raw string => int
int ParameterSystem::getInt(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return 0; // or throw an exception
    }
    return std::stoi(raw);
}

// Convert raw string => double
double ParameterSystem::getDouble(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return 0.0;
    }
    return std::stod(raw);
}

// Convert raw string => bool
bool ParameterSystem::getBool(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return false;
    }
    // naive check
    std::string lower;
    for(char c: raw) lower.push_back(std::tolower(static_cast<unsigned char>(c)));
    if(lower=="true" || lower=="on" || lower=="1" || lower=="yes") {
        return true;
    }
    return false;
}

// Possibly strip quotes
std::string ParameterSystem::getString(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(!raw.empty() && raw.front()=='"' && raw.back()=='"' && raw.size()>1) {
        return raw.substr(1, raw.size()-2);
    }
    return raw;
}

// Parse bracketed list => vector<double>
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
            // skip or handle error
        }
    }
    return result;
}

// Helper: parse bracketed list [1, 2, 3]
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
    // split by comma
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, ',')) {
        tokens.push_back(trim(item));
    }
    return tokens;
}

// Helper: trim
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
