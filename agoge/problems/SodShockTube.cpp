/**
 * @file SodShockTube.cpp
 * @brief Implementation of the SodShockTube problem for the Agoge application.
 *
 * This file implements the SodShockTube class methods to initialize a 
 * shock tube problem and register its parameters.
 */

#include "SodShockTube.hpp"

#include <cmath>
#include <iostream>

namespace agoge {
namespace problems {

/**
 * @brief Registers the parameters for the Sod Shock Tube problem.
 *
 * @param params The parameter system to register the parameters with.
 */
void SodShockTube::registerParameters(ParameterSystem &params) const {
    std::string prefix = name() + ".";

    params.addDefault(prefix + "left_density", "1.0");
    params.addDefault(prefix + "right_density", "0.125");
    params.addDefault(prefix + "left_pressure", "1.0");
    params.addDefault(prefix + "right_pressure", "0.1");
    params.addDefault(prefix + "u0", "0.0");
    params.addDefault(prefix + "gamma", "1.4");
}

/**
 * @brief Initializes the field for the Sod Shock Tube problem.
 *
 * @param Q The field to initialize.
 * @param params The parameter system to use for initialization.
 */
void SodShockTube::initialize(Field3D &Q, const ParameterSystem &params) {
    std::string prefix = name() + ".";

    double leftDensity = params.getDouble(prefix + "left_density");
    double rightDensity = params.getDouble(prefix + "right_density");
    double leftPressure = params.getDouble(prefix + "left_pressure");
    double rightPressure = params.getDouble(prefix + "right_pressure");
    double u0 = params.getDouble(prefix + "u0");
    double gamma = params.getDouble(prefix + "gamma");
    double midpoint = Q.bbox.xmin + 0.5 * (Q.bbox.xmax - Q.bbox.xmin);

    for (int k = 0; k < Q.Nz; ++k) {
        for (int j = 0; j < Q.Ny; ++j) {
            for (int i = 0; i < Q.Nx; ++i) {
                int idx = Q.interiorIndex(i, j, k);

                // Determine left or right side based on position                
                if (Q.xCenter(i) < midpoint) {
                    Q.rho[idx] = leftDensity;
                    Q.rhou[idx] = leftDensity * u0;
                    Q.rhov[idx] = 0.0;
                    Q.rhow[idx] = 0.0;
                    Q.E[idx] = leftPressure / (gamma - 1.0) +
                               0.5 * leftDensity * u0 * u0;
                } else {
                    Q.rho[idx] = rightDensity;
                    Q.rhou[idx] = rightDensity * u0;
                    Q.rhov[idx] = 0.0;
                    Q.rhow[idx] = 0.0;
                    Q.E[idx] = rightPressure / (gamma - 1.0) +
                               0.5 * rightDensity * u0 * u0;
                }
                // No gravity
                Q.phi[idx] = 0.0;
            }
        }
    }

    std::cout << "[SodShockTube] Initialized Sod Shock Tube with "
              << "left_density=" << leftDensity
              << ", right_density=" << rightDensity
              << ", left_pressure=" << leftPressure
              << ", right_pressure=" << rightPressure << ", u0=" << u0
              << ", gamma=" << gamma << "\n";
}

}  // namespace problems
}  // namespace agoge