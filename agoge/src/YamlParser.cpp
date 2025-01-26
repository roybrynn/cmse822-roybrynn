#include "YamlParser.hpp"
#include <fstream>
#include <sstream>
#include <cctype>

namespace agoge {

bool YamlParser::parseFile(const std::string &filename)
{
    std::ifstream fin(filename);
    if(!fin.is_open()) {
        return false; // file not found or unable to open
    }

    std::string line;
    while(std::getline(fin, line)) {
        // Trim whitespace
        line = trim(line);
        // Skip empty or comment lines
        if(line.empty() || line[0] == '#') {
            continue;
        }

        // Parse "key: value"
        auto kv = parseLine(line);
        if(!kv.first.empty()) {
            data_[kv.first] = kv.second;
        }
    }

    return true;
}

bool YamlParser::hasKey(const std::string &key) const
{
    return (data_.find(key) != data_.end());
}

std::string YamlParser::getString(const std::string &key) const
{
    auto it = data_.find(key);
    if(it != data_.end()) {
        return it->second;
    }
    return "";
}

std::pair<std::string, std::string> YamlParser::parseLine(const std::string &line) const
{
    // Find first colon
    auto pos = line.find(':');
    if(pos == std::string::npos) {
        // Not a valid "key: val" line
        return std::make_pair("", "");
    }
    // split
    std::string key = line.substr(0, pos);
    std::string val = line.substr(pos+1);

    // trim
    key = trim(key);
    val = trim(val);

    return std::make_pair(key, val);
}

std::string YamlParser::trim(const std::string &s) const
{
    if(s.empty()) return s;
    size_t start = 0;
    while(start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        start++;
    }
    size_t end = s.size() - 1;
    while(end > start && std::isspace(static_cast<unsigned char>(s[end]))) {
        end--;
    }
    return s.substr(start, end - start + 1);
}

} // namespace agoge
