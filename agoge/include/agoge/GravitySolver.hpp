#pragma once

#include <complex>
#include <vector>

#include "Config.hpp"
#include "Field3d.hpp"

/**
 * @file GravitySolver.hpp
 * @brief Defines the gravitational Poisson solver for agoge.
 *
 * Provides functions to solve Poisson's equation:
 *   ∇²φ = 4π G ρ
 * using either a naive DFT (O(N⁶), for educational purposes) or a Cooley–Tukey
 * FFT (O(N³ log N)) approach on a Field3D. Periodic boundary conditions are
 * assumed.
 */

namespace agoge {
namespace gravity {

/**
 * @enum GravityMethod
 * @brief Options for the gravitational Poisson solver approach.
 */
enum class GravityMethod {
    NAIVE_DFT,    ///< O(N⁶) naive direct DFT (educational only)
    COOLEY_TUKEY  ///< O(N³ log N) Cooley–Tukey FFT (recommended)
};

/**
 * @brief Solve Poisson's equation ∇²φ = 4π G ρ on a Field3D.
 *
 * The solver reads the interior density field (rho) from Q, transforms it into
 * Fourier space, applies the Green's function (with periodic BC and
 * power-of-two grid sizes assumed), and then transforms back to obtain the
 * gravitational potential, which is stored in Q.phi.
 *
 * @param Q      The Field3D object containing the density (rho) and where φ is
 * stored.
 * @param method The method to use: NAIVE_DFT or COOLEY_TUKEY.
 */
void solvePoisson(Field3D &Q, GravityMethod method);

}  // namespace gravity
}  // namespace agoge
