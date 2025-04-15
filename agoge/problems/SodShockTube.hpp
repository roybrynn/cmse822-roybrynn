/**
 * @file SodShockTube.hpp
 * @brief Declaration of the SodShockTube problem for the Agoge application.
 *
 * This file declares the SodShockTube class, a concrete Problem that sets up
 * a Sod shock tube initial condition.
 */

#pragma once

#include "../include/agoge/Field3d.hpp"
#include "Problem.hpp"

namespace agoge {
namespace problems {

/**
 * @brief Problem implementation for the Sod shock tube.
 */
class SodShockTube : public Problem {
   public:
    SodShockTube() = default;
    virtual ~SodShockTube() = default;

    /**
     * @brief Initialize the Field3D with Sod shock tube conditions.
     * @param Q The field to initialize.
     * @param params The parameter system with configured values.
     */
    void initialize(Field3D &Q, const ParameterSystem &params) override;

    /**
     * @brief Register default parameters for the Sod shock tube.
     * @param params The ParameterSystem instance.
     */
    void registerParameters(ParameterSystem &params) const override;

    /**
     * @brief Check if gravity is used.
     * @return false as gravity is not used.
     */
    bool useGravity() const override { return false; }

    /**
     * @brief Get the problem name.
     * @return "SodShockTube" as the identifier.
     */
    std::string name() const override { return "SodShockTube"; }
};

}  // namespace problems
}  // namespace agoge