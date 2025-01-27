#pragma once

#include "../include/agoge/Field3d.hpp"
#include "../include/agoge/ParameterSystem.hpp"
#include "Problem.hpp"

/**
 * @class GravityCollapse
 * @brief Implements a 3D cold dust sphere with periodic BC, checking
 * gravitational potential.
 */
namespace agoge {
namespace problems {

/**
 * @brief GravityCollapse initializes a uniform-density sphere of mass=1 (by
 * default), radius = R_jeans (by default), domain from -1.2*R to +1.2*R in each
 * dimension. We solve Poisson eqn and compare with the exact uniform-sphere
 * potential for error norms.
 */
class GravityCollapse : public Problem {
   public:
    GravityCollapse() = default;
    virtual ~GravityCollapse() = default;

    /**
     * @brief Initialize the cold dust collapse problem in Q:
     *        - Nx, Ny, Nz from user (ParameterSystem).
     *        - Mass = 1.0 (by default) => radius => R_jeans => domain =>
     * [-1.2R, 1.2R].
     *        - Inside sphere => constant rho, outside => 0.
     *        - velocities=0, E=0.
     *        Then solve Poisson, compare phi with analytic, print L1,L2.
     *
     * @param Q The field to initialize.
     * @param params The ParameterSystem instance containing parameters.
     */
    void initialize(Field3D &Q, const ParameterSystem &params) override;

    /**
     * @brief Register problem-specific default parameters with the
     * ParameterSystem.
     *
     * For instance:
     *   GravityCollapse.mass_solar => "1.0"
     *   GravityCollapse.jeans_radius => "1.0"
     *   etc.
     *
     * @param params The ParameterSystem instance.
     */
    void registerParameters(ParameterSystem &params) const override;

    /// We do use gravity
    bool useGravity() const override { return true; }

    /// Problem name
    std::string name() const override { return "GravityCollapse"; }
};

}  // namespace problems
}  // namespace agoge