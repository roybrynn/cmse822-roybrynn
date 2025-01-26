// problems/GravityCollapse.cpp
#include "GravityCollapse.hpp"

#include <cmath>
#include <iostream>

namespace agoge {
namespace problems {

/**
 * @brief Initialize the Field3D with gravity-driven initial conditions.
 *
 * @param Q The field to initialize.
 * @param params The ParameterSystem instance containing parameters.
 */
void GravityCollapse::initialize(Field3D &Q, const ParameterSystem &params) {
    // Retrieve problem-specific parameters using prefix
    std::string prefix = name() + ".";

    double initialDensity = params.getDouble(prefix + "initial_density");
    double collapseStrength = params.getDouble(prefix + "collapse_strength");
    double gravity =
        params.getDouble(prefix + "gravity");  // Example additional parameter

    // Example initialization logic
    for (int k = 0; k < Q.Nz; ++k) {
        for (int j = 0; j < Q.Ny; ++j) {
            for (int i = 0; i < Q.Nx; ++i) {
                int idx = Q.index(i, j, k);

                // Set density
                Q.rho[idx] = initialDensity;

                // Set velocity (e.g., zero or based on collapse strength)
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;

                // Set energy based on pressure (example)
                double pressure = 1.0;  // Could be parameterized
                double gamma = params.getDouble(
                    prefix + "gamma");  // Assuming gamma is a parameter
                Q.E[idx] = pressure / (gamma - 1.0);

                // Set gravitational potential if applicable
                Q.phi[idx] =
                    gravity * k * Q.dz;  // Simple gravity in z-direction
            }
        }
    }

    std::cout << "[GravityCollapse] Initialized gravity-driven collapse with "
              << "initial_density=" << initialDensity
              << ", collapse_strength=" << collapseStrength
              << ", gravity=" << gravity << "\n";
}

/**
 * @brief Register problem-specific default parameters with the ParameterSystem.
 *
 * @param params The ParameterSystem instance.
 */
void GravityCollapse::registerParameters(ParameterSystem &params) const {
    // Prefix parameters with the problem name to avoid key collisions
    std::string prefix = name() + ".";

    // Define default values for problem-specific parameters
    params.addDefault(prefix + "initial_density", "1.0");
    params.addDefault(prefix + "collapse_strength", "0.5");
    params.addDefault(prefix + "gravity",
                      "9.81");  // Example gravitational constant
    params.addDefault(prefix + "gamma", "1.4");  // Assuming ideal gas
    // Add more parameters as needed
}

}  // namespace problems
}  // namespace agoge