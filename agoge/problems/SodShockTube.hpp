// Example: problems/SodShockTube.hpp
#pragma once

#include "../include/agoge/Field3d.hpp"
#include "Problem.hpp"

namespace agoge {
namespace problems {

class SodShockTube : public Problem {
   public:
    SodShockTube() = default;
    virtual ~SodShockTube() = default;

    void initialize(Field3D &Q, const ParameterSystem &params) override;
    void registerParameters(ParameterSystem &params) const override;

    bool useGravity() const override { return false; }

    std::string name() const override { return "SodShockTube"; }
};

}  // namespace problems
}  // namespace agoge