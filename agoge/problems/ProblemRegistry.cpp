// problems/ProblemFactory.cpp
#include "Problem.hpp"
#include <memory>
#include <unordered_map>
#include <string>

// Forward declare or include problem classes:
#include "SodShockTube.hpp"
#include "GravityCollapse.hpp"

namespace agoge {
namespace problems {

std::unique_ptr<Problem> createProblem(const std::string &name)
{
    if(name == "sod") {
        return std::make_unique<SodShockTube>();
    } else if(name == "collapse") {
        return std::make_unique<GravityCollapse>();
    }
    // etc.
    // If unknown:
    return nullptr;
}

} // namespace problems
} // namespace agoge
