/**
 * @file ProblemRegistry.hpp
 * @brief Declarations for problem registry functions and forward declarations.
 *
 * This file declares the createProblem function and forward declares various
 * problem classes used in the Agoge application.
 */

#pragma once // Alternatively, use include guards

#include <memory>
#include <string>

// Forward declarations of problem classes
namespace agoge
{
    namespace problems
    {
        class Problem;
        class SodShockTube;
        class GravityCollapse;
        class GaussianPulse; 
        class Sedov;

        /**
         * @brief Factory function to create a Problem instance by name.
         *
         * @param name Identifier of the problem.
         * @return std::unique_ptr<Problem> A pointer to the created problem, or nullptr if unknown.
         */
        std::unique_ptr<Problem> createProblem(const std::string &name);
    }
}
