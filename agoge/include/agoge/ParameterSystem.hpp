#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "YamlParser.hpp"

/**
 * @file ParameterSystem.hpp
 * @brief A system that merges default parameters with optional YAML file overrides,
 *        providing typed getters (int, double, bool, string, list).
 */

namespace agoge {

/**
 * @class ParameterSystem
 * @brief Holds default parameters and merges them with YAML overrides.
 *
 * Example usage:
 *   ParameterSystem params;
 *   // set defaults
 *   params.setDefault("nx", "64");
 *   params.setDefault("cfl", "0.5");
 *   // parse YAML
 *   params.readYAML("input.yaml");
 *   // retrieve as int/double/bool/etc.
 *   int Nx = params.getInt("nx");
 *   double cflVal = params.getDouble("cfl");
 */
class ParameterSystem
{
public:
    /**
     * @brief Assign a default raw string value for a parameter key
     */
    void setDefault(const std::string &key, const std::string &rawValue);

    /**
     * @brief Load YAML file and override defaults with the YAML content
     */
    bool readYAML(const std::string &filename);

    /**
     * @brief Retrieve parameter as raw string.
     *        If not found in either YAML or defaults, return "".
     */
    std::string getRaw(const std::string &key) const;

    /**
     * @brief Retrieve as int
     */
    int getInt(const std::string &key) const;

    /**
     * @brief Retrieve as double
     */
    double getDouble(const std::string &key) const;

    /**
     * @brief Retrieve as bool
     */
    bool getBool(const std::string &key) const;

    /**
     * @brief Retrieve as string
     */
    std::string getString(const std::string &key) const;

    /**
     * @brief Retrieve as a list of double (for YAML lines like: [1.0, 2, 3.5])
     */
    std::vector<double> getDoubleList(const std::string &key) const;

private:
    YamlParser parser_; ///< minimal parser for overrides
    std::unordered_map<std::string, std::string> defaults_;

    /**
     * @brief Helper to parse a list from raw string
     *        e.g. "[1, 2.0, 3.5]"
     */
    std::vector<std::string> parseList(const std::string &raw) const;

    /**
     * @brief trim utility
     */
    std::string trim(const std::string &s) const;
};

} // namespace agoge
