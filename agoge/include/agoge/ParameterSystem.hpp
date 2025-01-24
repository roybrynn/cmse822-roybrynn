#pragma once

#include "YamlParser.hpp"
#include <string>
#include <vector>
#include <unordered_map>

/**
 * @file ParameterSystem.hpp
 * @brief A system that merges default parameters with optional YAML file overrides,
 *        providing typed getters (int, double, bool, string, list).
 */

namespace agoge {

class ParameterSystem
{
public:
    ParameterSystem();

    // Load YAML file and override defaults
    bool readYAML(const std::string &filename);

    // Typed getters
    int getInt(const std::string &key) const;
    double getDouble(const std::string &key) const;
    bool getBool(const std::string &key) const;
    std::string getString(const std::string &key) const;
    std::vector<double> getDoubleList(const std::string &key) const;

private:
    YamlParser parser_;
    std::unordered_map<std::string, std::string> defaults_;

    // Internal
    std::string getRaw(const std::string &key) const;
    std::vector<std::string> parseList(const std::string &raw) const;
    std::string trim(const std::string &s) const;
    void setDefaults();
};

} // namespace agoge
