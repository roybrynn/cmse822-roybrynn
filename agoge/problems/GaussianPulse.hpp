// problems/GaussianPulse.hpp
#pragma once

#include "../include/agoge/Field3d.hpp"
#include "Problem.hpp"

/**
 * @file GaussianPulse.hpp
 * @brief A problem setup for a Gaussian density pulse advecting in a uniform
 * velocity field.
 *
 * This can be run in 1D (Ny=Nz=1), 2D (Nz=1), or 3D. Domain is periodic,
 * ignoring gravity.
 */

namespace agoge {
namespace problems {

/**
 * @class GaussianPulse
 * @brief Implements a Gaussian pulse in density, uniform velocity in x (or any
 * direction), and constant pressure, leading to a simple advection test with
 * periodic BC.
 */
class GaussianPulse : public Problem {
   public:
    GaussianPulse() = default;
    virtual ~GaussianPulse() = default;

    /// Initialize the field with a Gaussian pulse in density
    void initialize(Field3D &Q, const ParameterSystem &params) override;

    /// Register problem-specific parameters
    void registerParameters(ParameterSystem &params) const override;

    /// No gravity for this test
    bool useGravity() const override { return false; }

    /// Problem name
    std::string name() const override { return "GaussianPulse"; }
};

}  // namespace problems
}  // namespace agoge