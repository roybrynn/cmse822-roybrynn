// problems/SodShockTube.hpp
#pragma once
#include "Problem.hpp"

namespace agoge {
namespace problems {

class SodShockTube : public Problem {
public:
    void initialize(Field3D &Q) override;
    bool useGravity() const override { return false; }
    std::string name() const override { return "SodShockTube"; }
};

} // namespace problems
} // namespace agoge
