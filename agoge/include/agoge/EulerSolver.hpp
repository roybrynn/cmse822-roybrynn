#pragma once

#include "Field3D.hpp"

/**
 * @file EulerSolver.hpp
 * @brief Declarations of functions for computing the Euler flux divergence and
 *        performing time integration using Runge-Kutta methods.
 */

namespace agoge {
namespace euler {

/**
 * @brief Computes the 3D flux divergence term \(-\nabla \cdot F(Q)\) for the
 * compressible Euler equations using simple finite differences.
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
void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField = nullptr);

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

} // namespace euler
} // namespace agoge
