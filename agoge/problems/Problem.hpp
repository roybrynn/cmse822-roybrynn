/**
 * @file Problem.hpp
 * @brief Base class for specifying problem-specific setup and physics choices.
 *
 * This file declares the abstract Problem class which requires derived classes
 * to implement initialization and parameter registration routines.
 */

#pragma once

#include <string>

#include "Field3d.hpp"
#include "ParameterSystem.hpp"

namespace agoge {
namespace problems {

/**
 * @brief Base class for problem-specific settings.
 */
class Problem {
   public:
    virtual ~Problem() = default;
    /**
     * @brief Initialize the Field3D with problem-specific initial conditions.
     *
     * @param Q The field to initialize.
     * @param params The parameters for the problem.
     */
    virtual void initialize(Field3D &Q, const ParameterSystem &params) = 0;

    /**
     * @brief Register problem-specific default parameters.
     *
     * @param params The ParameterSystem instance.
     */
    virtual void registerParameters(ParameterSystem &params) const = 0;

    /**
     * @brief Indicate whether gravity is used in this problem.
     * @return true if gravity is enabled; false otherwise.
     */
    virtual bool useGravity() const { return false; }

    /**
     * @brief Get the name of the problem.
     * @return The problem name.
     */
    virtual std::string name() const = 0;
};

}  // namespace problems
}  // namespace agoge