#pragma once

#include "Field3D.hpp"

/**
 * @file GravitySolver.hpp
 * @brief Declaration of the FFT-based Poisson solver for self-gravity in the Agoge project.
 */

namespace agoge {
namespace gravity {

/**
 * @brief Solve the Poisson equation \f$\nabla^2 \phi = 4 \pi G \rho\f$ 
 * for the gravitational potential \f$\phi\f$ under periodic boundary conditions,
 * using the 3D FFTW library (real-to-complex transform).
 *
 * This function modifies the \p Q.phi array in-place, storing the result of the
 * inverse transform that represents the solution of \f$\phi\f$.
 *
 * Steps:
 *  1. Copy \p Q.rho into a temporary real buffer of size Nx*Ny*Nz.
 *  2. Perform a forward real-to-complex FFT with FFTW.
 *  3. Multiply by the factor \f$-4 \pi G / k^2\f$ in wavenumber space, setting the \f$k=0\f$ mode to zero.
 *  4. Perform an inverse FFT to obtain \f$\phi\f$ in real space.
 *  5. Normalize by \f$1 / (Nx * Ny * Nz)\f$ because FFTW does not do this automatically.
 *
 * @param [in,out] Q Field3D containing \p rho. On output, Q.phi holds the gravitational potential.
 */
void solvePoissonFFT(Field3D &Q);

} // namespace gravity
} // namespace agoge
