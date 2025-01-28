// problems/ProblemRegistry.cpp
#include "ProblemRegistry.hpp"

#include "GaussianPulse.hpp"
#include "GravityCollapse.hpp"
#include "Problem.hpp"  // Include definitions of Problem and its derivatives
#include "SodShockTube.hpp"

namespace agoge {
namespace problems {

std::unique_ptr<Problem> createProblem(const std::string &name) {
    if (name == "sod") {
        return std::make_unique<SodShockTube>();
    } else if (name == "collapse") {
        return std::make_unique<GravityCollapse>();
    } else if (name == "GaussianPulse") {  // Updated identifier
        return std::make_unique<GaussianPulse>();
    }
    // Add more problem types as needed

    // If unknown:
    return nullptr;
}

}  // namespace problems
}  // namespace agoge