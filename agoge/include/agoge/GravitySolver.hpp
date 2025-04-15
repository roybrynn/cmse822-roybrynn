#pragma once

#include <vector>
#include <complex>

#include "Field3d.hpp"

namespace agoge
{
    namespace gravity
    {

        /**
         * @enum GravityMethod
         * @brief Options for the gravitational Poisson solver approach.
         */
        enum class GravityMethod
        {
            NAIVE_DFT,   ///< O(N^6) demonstration
            COOLEY_TUKEY ///< O(N^3 log N) demonstration
        };

        /**
         * @brief Solve Poisson's eqn \(\nabla^2 \phi = 4\pi G (\rho - rhoMean)\)
         * using periodic BC with either naive DFT or Cooley–Tukey FFT approach.
         *
         * This follows the "Jeans swindle" approach used by Athena:
         *  - subtract average domain density => net mass = 0
         *  - transform
         *  - multiply by -1/k^2 in frequency domain (except k=0 => 0)
         *  - inverse transform => potential fluctuations
         *
         * The domain is purely periodic. Nx,Ny,Nz must be power-of-two for COOLEY_TUKEY,
         * but can be anything for NAIVE_DFT.
         */
        void solvePoisson(Field3D &Q, GravityMethod method);

    } // namespace gravity
} // namespace agoge