#pragma once

/**
 * @file Config.hpp
 * @brief Global configuration constants for the Agoge solver.
 *
 * This header provides fundamental constants and toggles used throughout
 * the Agoge codebase, such as the gas constant (gamma), gravitational constant,
 * and whether self-gravity is enabled.
 */

#ifndef AGOGE_NGHOST
#define AGOGE_NGHOST 2
#endif

namespace agoge {
namespace config {

/// @brief Number of ghost cells for boundary conditions.
constexpr int Ng = AGOGE_NGHOST;

/**
 * @brief Ratio of specific heats for an ideal gas.
 *
 * Defaults to 1.4 for diatomic gas (like air at standard conditions).
 * Modify as needed for different equations of state.
 */
constexpr double gamma_gas = 1.4;

/**
 * @brief Gravitational constant (in code units).
 *
 * Set to 1.0 by default for convenience in dimensionless or code unit systems.
 * Adjust to a physical value if using CGS or SI units.
 */
constexpr double G = 6.67e-8;

/**
 * @brief Toggle for enabling or disabling self-gravity.
 *
 * When true, the code will compute gravitational potential via Poisson's equation
 * and add corresponding source terms to the fluid equations.
 */
constexpr bool use_gravity = true;

enum class BoundaryCondition {
    PERIODIC,
    OUTFLOW
    // Add more boundary conditions if needed
};

} // namespace config
} // namespace agoge
