// ParameterSystem.cpp
#include "ParameterSystem.hpp"

#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>

namespace agoge {

// Constructor initializes the internal defaults.
ParameterSystem::ParameterSystem() {
    // Calls helper to populate defaults map with typical simulation parameters.
    setDefaults();
}

// Sets baseline simulation parameters for domain sizes, CFL, boundary conditions, etc.
void ParameterSystem::setDefaults() {
    defaults_["nx"] = "64";
    defaults_["ny"] = "64";
    defaults_["nz"] = "64";

    defaults_["nghost"] = "1";

    defaults_["cfl"] = "0.5";
    defaults_["use_gravity"] = "false";
    defaults_["xmin"] = "0.0";
    defaults_["xmax"] = "1.0";
    defaults_["ymin"] = "0.0";
    defaults_["ymax"] = "1.0";
    defaults_["zmin"] = "0.0";
    defaults_["zmax"] = "1.0";
    defaults_["sound_crossings"] = "1.0";
    defaults_["dt_init"] = "0.01";
    defaults_["t_max"] = "1.0";

    defaults_["problem_name"] = "\"sod\"";

    // possible values: "periodic", "outflow"
    defaults_["bc_xmin"] = "periodic";
    defaults_["bc_xmax"] = "periodic";
    defaults_["bc_ymin"] = "periodic";
    defaults_["bc_ymax"] = "periodic";
    defaults_["bc_zmin"] = "periodic";
    defaults_["bc_zmax"] = "periodic";

    // new default
    defaults_["do_euler_update"] =
        "true";  // or "true" if you want skip by default
    
    // Turn off IO         
    defaults_["do_io"] = "true";

    // new default for wall clock time: 14400 seconds = 4 hours
    defaults_["max_wallclock_time"] = "14400";

    // number of steps between screen output
    defaults_["screen_out_interval"] = "2";

    // output directory
    defaults_["output_dir"] = "\"output\"";
}

// Adds or overrides a default parameter with a given key-value pair.
void ParameterSystem::addDefault(const std::string &key,
                                 const std::string &value) {
    defaults_[key] = value;
}

// Check if parameter exists
bool ParameterSystem::hasParameter(const std::string &key) const {
    return !getRaw(key).empty();
}

// Parse a YAML file and populate yamlData_ map
bool ParameterSystem::readYAML(const std::string &filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "[ParameterSystem] Could not open YAML file: " << filename << "\n";
        return false;
    }

    std::string line;
    while (std::getline(fin, line)) {
        // Trim whitespace
        line = trim(line);
        // Skip empty or comment lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse "key: value"
        auto kv = parseLine(line);
        if (!kv.first.empty()) {
            yamlData_[kv.first] = kv.second;
        }
    }

    return true;
}

// Parse a line of the format "key: value"
std::pair<std::string, std::string> ParameterSystem::parseLine(const std::string &line) const {
    // Find first colon
    auto pos = line.find(':');
    if (pos == std::string::npos) {
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

// Returns the parameter value as a raw string; checks both parsed YAML and defaults.
std::string ParameterSystem::getRaw(const std::string &key) const {
    auto it = yamlData_.find(key);
    if (it != yamlData_.end()) {
        return it->second;
    }
    
    it = defaults_.find(key);
    if (it != defaults_.end()) {
        return it->second;
    }
    return "";
}

// Converts the parameter string to an integer, with fallback behavior if missing.
int ParameterSystem::getInt(const std::string &key) const {
    std::string raw = getRaw(key);
    if (raw.empty()) {
        return 0;
    }
    return std::stoi(raw);
}

// Converts the parameter string to a double, returning 0.0 if missing.
double ParameterSystem::getDouble(const std::string &key) const {
    std::string raw = getRaw(key);
    if (raw.empty()) {
        return 0.0;
    }
    return std::stod(raw);
}

// Interprets the parameter string as a boolean, accepting various true/false forms.
bool ParameterSystem::getBool(const std::string &key) const {
    std::string raw = getRaw(key);
    if (raw.empty()) {
        return false;
    }
    std::string lower;
    for (char c : raw) {
        lower.push_back(std::tolower(static_cast<unsigned char>(c)));
    }
    if (lower == "true" || lower == "on" || lower == "1" || lower == "yes") {
        return true;
    }
    return false;
}

// Strips surrounding quotes if present.
std::string ParameterSystem::getString(const std::string &key) const {
    std::string raw = getRaw(key);
    if (!raw.empty() && raw.front() == '"' && raw.back() == '"' &&
        raw.size() > 1) {
        return raw.substr(1, raw.size() - 2);
    }
    return raw;
}

// Parses a comma-separated list into a vector of doubles, ignoring parsing errors.
std::vector<double> ParameterSystem::getDoubleList(
    const std::string &key) const {
    std::vector<double> result;
    std::string raw = getRaw(key);
    if (raw.empty()) {
        return result;
    }
    auto tokens = parseList(raw);
    for (const auto &t : tokens) {
        try {
            double val = std::stod(t);
            result.push_back(val);
        } catch (...) {
            // skip error
        }
    }
    return result;
}

// Translates a boundary condition string (e.g., "periodic", "outflow") into an enum type.
config::BoundaryCondition ParameterSystem::getBoundaryCondition(
    const std::string &key) const {
    // e.g. param might be "bc_xmin" -> "outflow"
    std::string bcStr = getString(key);
    // convert to lowercase
    std::string lower;
    for (char c : bcStr) {
        lower.push_back(std::tolower(static_cast<unsigned char>(c)));
    }
    if (lower == "outflow") {
        return config::BoundaryCondition::OUTFLOW;
    }
    // default
    return config::BoundaryCondition::PERIODIC;
}

// Helper to remove brackets from a list-string and split the values by commas.
std::vector<std::string> ParameterSystem::parseList(
    const std::string &raw) const {
    std::vector<std::string> tokens;
    std::string s = trim(raw);
    if (!s.empty() && s.front() == '[') {
        s.erase(s.begin());
    }
    if (!s.empty() && s.back() == ']') {
        s.pop_back();
    }
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        tokens.push_back(trim(item));
    }
    return tokens;
}

// Trims leading and trailing whitespace from a string.
std::string ParameterSystem::trim(const std::string &str) const {
    if (str.empty()) return str;
    size_t start = 0;
    while (start < str.size() &&
           std::isspace(static_cast<unsigned char>(str[start]))) {
        start++;
    }
    size_t end = str.size() - 1;
    while (end > start && std::isspace(static_cast<unsigned char>(str[end]))) {
        end--;
    }
    return str.substr(start, end - start + 1);
}

}  // namespace agoge
