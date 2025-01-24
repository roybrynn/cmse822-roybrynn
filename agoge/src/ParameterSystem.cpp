#include "ParameterSystem.hpp"
#include <sstream>
#include <cctype>
#include <iostream>

namespace agoge {

void ParameterSystem::setDefault(const std::string &key, const std::string &rawValue)
{
    defaults_[key] = rawValue;
}

bool ParameterSystem::readYAML(const std::string &filename)
{
    bool ok = parser_.parseFile(filename);
    if(!ok) {
        std::cerr << "[ParameterSystem] Could not parse YAML file: " << filename << "\n";
    }
    return ok;
}

std::string ParameterSystem::getRaw(const std::string &key) const
{
    // 1) check parser
    if(parser_.hasKey(key)) {
        return parser_.getString(key);
    }
    // 2) else check defaults
    auto it = defaults_.find(key);
    if(it != defaults_.end()) {
        return it->second;
    }
    // not found
    return "";
}

int ParameterSystem::getInt(const std::string &key) const
{
    std::string raw = getRaw(key);
    if(raw.empty()) {
        return 0; // or throw
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
    // naive check
    // lower-case raw
    std::string lower;
    for(char c: raw) {
        lower.push_back(std::tolower(static_cast<unsigned char>(c)));
    }
    if(lower=="true" || lower=="on" || lower=="1" || lower=="yes") {
        return true;
    }
    return false;
}

std::string ParameterSystem::getString(const std::string &key) const
{
    std::string raw = getRaw(key);
    // optionally remove surrounding quotes if present
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
    // parse a bracketed list: e.g. "[1, 2.5, 3.0]"
    // we do a naive approach
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

std::vector<std::string> ParameterSystem::parseList(const std::string &raw) const
{
    // e.g. raw="[1.0,2,3.5]" => tokens: "1.0","2","3.5"
    std::vector<std::string> tokens;
    std::string s = raw;
    // remove leading '[' and trailing ']'
    // also remove whitespace
    s = trim(s);
    if(!s.empty() && s.front()=='[') {
        s.erase(s.begin());
    }
    if(!s.empty() && s.back()==']') {
        s.pop_back();
    }
    // now split by comma
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, ',')) {
        item = trim(item);
        tokens.push_back(item);
    }
    return tokens;
}

std::string ParameterSystem::trim(const std::string &str) const
{
    if(str.empty()) return str;
    size_t start = 0;
    while(start<str.size() && std::isspace(static_cast<unsigned char>(str[start]))) {
        start++;
    }
    size_t end = str.size()-1;
    while(end>start && std::isspace(static_cast<unsigned char>(str[end]))) {
        end--;
    }
    return str.substr(start, end-start+1);
}

} // namespace agoge
