// problems/Problem.hpp
#pragma once

#include <string>

#include "Field3d.hpp"
#include "ParameterSystem.hpp"

namespace agoge {
namespace problems {

/**
 * @brief Base class for specifying problem-specific setup and physics choices.
 */
class Problem {
   public:
    virtual ~Problem() = default;

    /**
     * @brief Initialize the Field3D with problem-specific initial conditions.
     *
     * @param Q The field to initialize.
     * @param params The ParameterSystem instance containing parameters.
     */
    virtual void initialize(Field3D &Q, const ParameterSystem &params) = 0;

    /**
     * @brief Register problem-specific default parameters with the
     * ParameterSystem.
     *
     * @param params The ParameterSystem instance.
     */
    virtual void registerParameters(ParameterSystem &params) const = 0;

    /**
     * @brief Whether or not gravity should be enabled for this problem.
     *
     * @return True if gravity is needed, false otherwise.
     */
    virtual bool useGravity() const { return false; }

    /**
     * @brief Get the name identifier of the problem.
     *
     * @return Problem name as a string.
     */
    virtual std::string name() const = 0;
};

}  // namespace problems
}  // namespace agoge