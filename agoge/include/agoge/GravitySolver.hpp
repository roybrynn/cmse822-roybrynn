#pragma once

#include "Field3d.hpp"
#include <string>

namespace agoge {
namespace gravity {

/**
 * @enum GravityMethod
 * @brief Options for the gravitational Poisson solver approach.
 */
enum class GravityMethod {
    NAIVE_DFT, ///< \f$O(N^6)\f$ naive direct DFT
    COOLEY_TUKEY ///< \f$O(N^3 \log N)\f$ Cooley–Tukey FFT
};

/**
 * @brief Solve Poisson's eqn \(\nabla^2 \phi = 4\pi G\rho\) using periodic BC
 * under either a naive DFT or a Cooley–Tukey FFT approach.
 *
 * @param[in,out] Q     Field3D with \p rho. On output, \p Q.phi is set.
 * @param[in]     method GravityMethod::NAIVE_DFT or GravityMethod::COOLEY_TUKEY.
 *
 * If method = NAIVE_DFT, the solve is \f$O(N^6)\f$ (educational only).
 * If method = COOLEY_TUKEY, the solve is \f$O(N^3 \log N)\f$, which is far faster
 * but still a manually-coded 1D-FFT-based approach (power-of-two sizes assumed).
 */
void solvePoisson(Field3D &Q, GravityMethod method);

} // namespace gravity
} // namespace agoge
