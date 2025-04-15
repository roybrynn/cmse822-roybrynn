// Example: problems/Sedov.hpp
#pragma once

#include "../include/agoge/Field3d.hpp"
#include "Problem.hpp"

namespace agoge {
namespace problems {

class Sedov : public Problem {
   public:
    Sedov() = default;
    virtual ~Sedov() = default;

    void initialize(Field3D &Q, const ParameterSystem &params) override;
    void registerParameters(ParameterSystem &params) const override;

    bool useGravity() const override { return false; }

    std::string name() const override { return "Sedov"; }
};

}  // namespace problems
}  // namespace agoge