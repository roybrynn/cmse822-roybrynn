#pragma once

#include "Field3d.hpp"

/**
 * @file EulerSolver.hpp
 * @brief Declarations of functions for computing the Euler flux divergence and
 *        performing time integration using a shock-capturing finite difference
 * method.
 */

namespace agoge {
namespace euler {

/**
 * @brief Computes the 3D flux divergence term \(-\nabla \cdot F(Q)\) for the
 * compressible Euler equations using a second-order finite-difference scheme,
 * augmented by artificial viscosity for shock capture.
 *
 * Optionally adds gravitational source terms if \p gravField is provided
 * and contains a valid gravitational potential \(\phi\).
 *
 * @param[in]  Q         The current state field (density, momentum, energy).
 * @param[out] LQ        The computed right-hand side for each variable.
 * @param[in]  gravField Optional field containing gravitational potential
 *                       (phi). If non-null, source terms for gravity are added
 *                       to the momentum and energy equations.
 */
void computeL(const Field3D &Q, Field3D &LQ,
              const Field3D *gravField = nullptr);

/**
 * @brief Performs a two-stage (midpoint) Runge-Kutta time integration on the
 * compressible Euler system.
 *
 * Calls computeL(Q, ...) internally for both stages.
 *
 * @param[in,out] Q  The field to be updated in-place.
 * @param[in]     dt The time step.
 */
void runRK2(Field3D &Q, double dt);

/**
 * @brief Computes a stable time step size \(\Delta t\) based on the CFL
 * condition.
 *
 * Scans all cells in \p Q, computing local wave speed (|u|+a) to find the
 * maximum wave speed. Returns \(\Delta t = \text{cfl} \times \min(dx,dy,dz) /
 * \max(\text{wave speed})\).
 *
 * @param[in] Q   The current state field (density, momentum, energy).
 * @param[in] cfl The desired CFL number (e.g., 0.5).
 * @return The computed time step size.
 */
double computeTimeStep(const Field3D &Q, double cfl);

}  // namespace euler
}  // namespace agoge
