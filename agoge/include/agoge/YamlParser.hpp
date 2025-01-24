#pragma once

#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file YamlParser.hpp
 * @brief A minimal YAML parser that reads key-value pairs (and lists),
 *        storing them in a string-based map. No external YAML libraries.
 *
 * Example YAML subset supported:
 * ----------------------------------
 * # This is a comment
 * nx: 64
 * cfl: 0.5
 * use_gravity: true
 * problem_name: "sod"
 * domain: [1.0, 1.0, 1.0]
 * ----------------------------------
 *
 * The parser will:
 *  - Ignore lines starting with '#' or empty lines
 *  - Split on the first ':' to get key and value
 *  - Trim whitespace
 *  - Store the raw value as a string
 */

namespace agoge {

/**
 * @class YamlParser
 * @brief Minimal parser for a subset of YAML, storing key->string_value.
 */
class YamlParser
{
public:
    /**
     * @brief Parse a YAML file (subset) from disk. 
     * @param filename path to the YAML file
     * @return true if the file was parsed successfully, false otherwise
     */
    bool parseFile(const std::string &filename);

    /**
     * @brief Check if a key exists in the map
     */
    bool hasKey(const std::string &key) const;

    /**
     * @brief Get the raw string value for a key
     *        (or empty string if not found)
     */
    std::string getString(const std::string &key) const;

private:
    /**
     * @brief Internal storage: key -> raw string
     */
    std::unordered_map<std::string, std::string> data_;

    /**
     * @brief Helper: parse a line "key: value"
     *        Return (key, value) or empty if invalid
     */
    std::pair<std::string, std::string> parseLine(const std::string &line) const;

    /**
     * @brief Helper: trim leading/trailing whitespace
     */
    std::string trim(const std::string &s) const;
};

} // namespace agoge
