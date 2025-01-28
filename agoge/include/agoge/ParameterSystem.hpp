// ParameterSystem.hpp
#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "Config.hpp"  // for BoundaryCondition
#include "YamlParser.hpp"

namespace agoge {

class ParameterSystem {
   public:
    ParameterSystem();

    bool readYAML(const std::string &filename);

    // typed getters
    int getInt(const std::string &key) const;
    double getDouble(const std::string &key) const;
    bool getBool(const std::string &key) const;
    std::string getString(const std::string &key) const;
    std::vector<double> getDoubleList(const std::string &key) const;

    /**
     * @brief Convert a string in ParameterSystem to a config::BoundaryCondition
     * enum. e.g. "periodic" -> BoundaryCondition::PERIODIC "outflow" ->
     * BoundaryCondition::OUTFLOW If not recognized, defaults to PERIODIC for
     * safety.
     */
    config::BoundaryCondition getBoundaryCondition(
        const std::string &key) const;

    /**
     * @brief Add a default parameter key-value pair.
     *
     * @param key The parameter key.
     * @param value The default value as a string.
     */
    void addDefault(const std::string &key, const std::string &value);

   private:
    YamlParser parser_;
    std::unordered_map<std::string, std::string> defaults_;

    std::string getRaw(const std::string &key) const;
    std::vector<std::string> parseList(const std::string &raw) const;
    std::string trim(const std::string &s) const;
    void setDefaults();
};

}  // namespace agoge
