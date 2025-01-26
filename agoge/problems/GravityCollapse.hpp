// problems/GravityCollapse.hpp
#pragma once

#include "../include/agoge/Field3d.hpp"
#include "Problem.hpp"

/**
 * @class GravityCollapse
 * @brief Implements a gravity-driven collapse problem.
 */
namespace agoge {
namespace problems {

class GravityCollapse : public Problem {
   public:
    GravityCollapse() = default;
    virtual ~GravityCollapse() = default;

    /**
     * @brief Initialize the field with gravity-driven initial conditions.
     *
     * @param Q The field to initialize.
     * @param params The ParameterSystem instance containing parameters.
     */
    void initialize(Field3D &Q, const ParameterSystem &params) override;

    /**
     * @brief Register problem-specific default parameters with the
     * ParameterSystem.
     *
     * @param params The ParameterSystem instance.
     */
    void registerParameters(ParameterSystem &params) const override;

    /// Use gravity for this problem
    bool useGravity() const override { return true; }

    /// Problem name
    std::string name() const override { return "GravityCollapse"; }
};

}  // namespace problems
}  // namespace agoge