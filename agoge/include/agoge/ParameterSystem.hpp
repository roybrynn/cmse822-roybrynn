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

/**
 * @class ParameterSystem
 * @brief Holds default parameters (internally) and merges them with YAML overrides.
 *
 * Usage:
 *   ParameterSystem params;
 *   params.readYAML("input.yaml"); // optional
 *   int Nx = params.getInt("nx");
 *   double cflVal = params.getDouble("cfl");
 */
class ParameterSystem
{
public:
    /**
     * @brief Constructor that sets default parameters internally
     */
    ParameterSystem();

    /**
     * @brief Load YAML file (optional) and override defaults
     */
    bool readYAML(const std::string &filename);

    /**
     * @brief Retrieve as int/double/bool/string/list
     */
    int getInt(const std::string &key) const;
    double getDouble(const std::string &key) const;
    bool getBool(const std::string &key) const;
    std::string getString(const std::string &key) const;
    std::vector<double> getDoubleList(const std::string &key) const;

private:
    // The minimal parser for overrides
    YamlParser parser_;

    // default map: key -> raw string
    std::unordered_map<std::string, std::string> defaults_;

    // Internal getters
    std::string getRaw(const std::string &key) const;
    std::vector<std::string> parseList(const std::string &raw) const;
    std::string trim(const std::string &s) const;

    /**
     * @brief Helper to set all built-in defaults.
     *        Called by the constructor.
     */
    void setDefaults();
};

} // namespace agoge
