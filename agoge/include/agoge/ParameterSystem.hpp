// ParameterSystem.hpp
#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "Config.hpp"  // for BoundaryCondition

namespace agoge {

/**
 * @class ParameterSystem
 * @brief Centralized parameter handling system for the Agoge application.
 * 
 * This class provides parameter handling with defaults, YAML configuration file support,
 * and strongly-typed parameter access. It encapsulates the YAML parsing functionality
 * that was previously in a separate YamlParser class.
 */
class ParameterSystem {
public:
    /**
     * @brief Constructor initializes with built-in simulation defaults
     */
    ParameterSystem();

    /**
     * @brief Parse a YAML configuration file
     * @param filename Path to the YAML file
     * @return true if parsing was successful, false otherwise
     */
    bool readYAML(const std::string &filename);

    // Strongly-typed parameter getters
    /**
     * @brief Get parameter as integer
     * @param key Parameter name
     * @return Integer value, 0 if not found or invalid
     */
    int getInt(const std::string &key) const;

    /**
     * @brief Get parameter as floating-point value
     * @param key Parameter name
     * @return Double value, 0.0 if not found or invalid
     */
    double getDouble(const std::string &key) const;

    /**
     * @brief Get parameter as boolean
     * @param key Parameter name
     * @return Boolean value (true for "true", "yes", "on", "1"; false otherwise)
     */
    bool getBool(const std::string &key) const;

    /**
     * @brief Get parameter as string with optional quote stripping
     * @param key Parameter name
     * @return String value, empty if not found
     */
    std::string getString(const std::string &key) const;

    /**
     * @brief Get parameter as list of doubles
     * @param key Parameter name
     * @return Vector of doubles parsed from comma-separated list
     */
    std::vector<double> getDoubleList(const std::string &key) const;

    /**
     * @brief Convert a string parameter to BoundaryCondition enum
     * @param key Parameter name (e.g., "bc_xmin")
     * @return BoundaryCondition enum value (PERIODIC default if not recognized)
     */
    config::BoundaryCondition getBoundaryCondition(const std::string &key) const;

    /**
     * @brief Add or override a default parameter
     * @param key Parameter name
     * @param value Parameter value as string
     */
    void addDefault(const std::string &key, const std::string &value);

    /**
     * @brief Check if a parameter exists (either from YAML or defaults)
     * @param key Parameter name
     * @return true if parameter exists, false otherwise
     */
    bool hasParameter(const std::string &key) const;

private:
    std::unordered_map<std::string, std::string> yamlData_;
    std::unordered_map<std::string, std::string> defaults_;

    // Get raw parameter value string
    std::string getRaw(const std::string &key) const;
    
    // Helper functions for parsing
    std::vector<std::string> parseList(const std::string &raw) const;
    std::string trim(const std::string &s) const;
    std::pair<std::string, std::string> parseLine(const std::string &line) const;
    void setDefaults();
};

}  // namespace agoge
