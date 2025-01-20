#pragma once
#include "Field.hpp"

namespace agoge {

// Solve Poisson's eqn for phi: nabla^2(phi) = 4*pi*G*rho (periodic BC)
void solvePoissonFFT(Field3D &Q, double G = 1.0);

} // namespace agoge
