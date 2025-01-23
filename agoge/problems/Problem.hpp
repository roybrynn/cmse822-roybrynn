// problems/Problem.hpp
#pragma once

#include <string>
#include "Field3d.hpp"

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
     */
    virtual void initialize(Field3D &Q) = 0;

    /**
     * @brief Whether or not gravity should be enabled for this problem.
     *
     * @return True if gravity is needed, false otherwise.
     */
    virtual bool useGravity() const { return false; }

    // Optionally: a name identifier or other problem-specific parameters
    virtual std::string name() const = 0;
};

} // namespace problems
} // namespace agoge
