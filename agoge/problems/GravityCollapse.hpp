// problems/GravityCollapse.hpp
#pragma once
#include "Problem.hpp"

namespace agoge {
namespace problems {

class GravityCollapse : public Problem {
public:
    void initialize(Field3D &Q) override;
    bool useGravity() const override { return true; }
    std::string name() const override { return "GravityCollapse"; }
};

} // namespace problems
} // namespace agoge
