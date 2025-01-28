// Example: problems/SodShockTube.cpp
#include "SodShockTube.hpp"

#include <cmath>
#include <iostream>

namespace agoge {
namespace problems {

void SodShockTube::registerParameters(ParameterSystem &params) const {
    std::string prefix = name() + ".";

    params.addDefault(prefix + "left_density", "1.0");
    params.addDefault(prefix + "right_density", "0.125");
    params.addDefault(prefix + "left_pressure", "1.0");
    params.addDefault(prefix + "right_pressure", "0.1");
    params.addDefault(prefix + "u0", "0.0");
    params.addDefault(prefix + "gamma", "1.4");
}

void SodShockTube::initialize(Field3D &Q, const ParameterSystem &params) {
    std::string prefix = name() + ".";

    double leftDensity = params.getDouble(prefix + "left_density");
    double rightDensity = params.getDouble(prefix + "right_density");
    double leftPressure = params.getDouble(prefix + "left_pressure");
    double rightPressure = params.getDouble(prefix + "right_pressure");
    double u0 = params.getDouble(prefix + "u0");
    double gamma = params.getDouble(prefix + "gamma");

    for (int k = 0; k < Q.Nz; ++k) {
        for (int j = 0; j < Q.Ny; ++j) {
            for (int i = 0; i < Q.Nx; ++i) {
                int idx = Q.interiorIndex(i, j, k);

                // Determine left or right side based on position
                double x = (i + 0.5) * Q.dx;
                double midpoint = 0.5 * Q.Nx * Q.dx;
                if (x < midpoint) {
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