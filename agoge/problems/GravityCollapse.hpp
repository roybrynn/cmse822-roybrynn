#pragma once

#include <string>

#include "../include/agoge/Field3d.hpp"
#include "../include/agoge/ParameterSystem.hpp"
#include "Problem.hpp"

/**
 * @class GravityCollapse
 * @brief Implements a cold dust sphere with self-gravity.
 *
 * Initializes a uniform-density sphere (with total mass and Jeans radius
 * parameters), solves Poisson's equation for the gravitational potential using
 * the custom FFT-based solver, and computes L₁ and L₂ error norms by comparing
 * with the analytic uniform-sphere potential.
 */
namespace agoge {
namespace problems {

class GravityCollapse : public Problem {
   public:
    GravityCollapse() = default;
    virtual ~GravityCollapse() = default;

    /**
     * @brief Initialize the Field3D with problem-specific initial conditions.
     *
     * Sets up a spherical density distribution: inside the sphere a uniform
     * density, outside a low floor value. Velocities and energy are set to near
     * zero (cold dust). Then solves Poisson's equation and computes error
     * norms.
     *
     * @param Q      The Field3D object to initialize.
     * @param params The ParameterSystem containing problem parameters.
     */
    void initialize(Field3D &Q, const ParameterSystem &params) override;

    /**
     * @brief Register default parameters for the GravityCollapse problem.
     *
     * Defaults include:
     *   - mass_solar (total mass)
     *   - r_jeans (Jeans radius)
     *   - grav_method ("cooley_tukey" by default)
     *
     * @param params The ParameterSystem to which defaults are added.
     */
    void registerParameters(ParameterSystem &params) const override;

    /// Indicate that this problem uses gravity.
    bool useGravity() const override { return true; }

    /// Return the name of the problem.
    std::string name() const override { return "GravityCollapse"; }
};

}  // namespace problems
}  // namespace agoge
