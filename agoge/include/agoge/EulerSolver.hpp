#pragma once

#include "Field3d.hpp"
#include <limits>

/**
 * @file EulerSolver.hpp
 * @brief Declarations of functions for computing the Euler flux divergence and
 *        performing time integration using a shock-capturing finite difference
 * method.
 */
#ifndef AGOGE_NTILE
#define AGOGE_NTILE 16
#endif

namespace agoge {
namespace euler {

constexpr int Ntile = AGOGE_NTILE;                       // 1D tile size
const int Nghost = config::Ng;              // Number of ghost cells
const int NtileGhost = Ntile + 2 * Nghost;  // 1D tile size with ghost cells
const int tileSize = NtileGhost * NtileGhost * NtileGhost;  // 3D tile size
constexpr int Nflux = 5;  // Number of fluxes in tile excluding phi

enum var {
    rho = 0,
    rhou = 1,
    rhov = 2,
    rhow = 3,
    E = 4
};

constexpr double huge_val = std::numeric_limits<double>::max();
const std::array<double, Nflux> floor = {1.0e-14, -huge_val, -huge_val,
                                         -huge_val, 1.0e-10};

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
void computeL(const Field3D &Q, Field3D &LQ);

/// @brief  Computes the 3D flux divergence term \(-\nabla \cdot F(Q)\) for the
/// @param tileQ 
/// @param tileL  
/// @param xExtent 
/// @param yExtent 
/// @param zExtent 
/// @param dx_inv 
/// @param dy_inv 
/// @param dz_inv 
/// @param tilePhi 
void computeLtile(std::vector<double> &tileQ,
                  std::vector<double> &tileL, int xExtent,
                  int yExtent, int zExtent, double dx_inv, double dy_inv,
                  double dz_inv, std::vector<double> &tilePhi);

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
