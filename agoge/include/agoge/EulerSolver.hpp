#pragma once
#include "Field.hpp"

namespace agoge {

// Public function to do one RK2 step, for example
void doTimeStep(Field3D &Q, double dt);

// Could also have a function to compute the fluxes, or to set initial conditions
// ...
} // namespace agoge
