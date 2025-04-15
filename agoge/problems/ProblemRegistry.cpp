/**
 * @file ProblemRegistry.cpp
 * @brief Implementation of the problem registry for the Agoge application.
 *
 * This file implements the createProblem factory function that instantiates
 * different types of problems based on a string identifier.
 */

#include "ProblemRegistry.hpp"
#include "GaussianPulse.hpp"
#include "GravityCollapse.hpp"
#include "Problem.hpp"  // Include definitions of Problem and its derivatives
#include "SodShockTube.hpp"
#include "Sedov.hpp"

namespace agoge {
namespace problems {

std::unique_ptr<Problem> createProblem(const std::string &name) {
    if (name == "sod") {
        return std::make_unique<SodShockTube>();
    } else if (name == "collapse") {
        return std::make_unique<GravityCollapse>();
    } else if (name == "GaussianPulse") { 
        return std::make_unique<GaussianPulse>();
    } else if (name == "Sedov") {  
        return std::make_unique<Sedov>();
    }
    // Add more problem types as needed

    // If unknown:
    return nullptr;
}

}  // namespace problems
}  // namespace agoge