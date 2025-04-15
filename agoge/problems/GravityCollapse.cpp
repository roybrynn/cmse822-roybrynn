#include "GravityCollapse.hpp"

#include <chrono>
#include <cmath>
#include <iostream>

#include "../include/agoge/Config.hpp"  // for agoge::config::G
#include "../include/agoge/Field3d.hpp"
#include "../include/agoge/GravitySolver.hpp"  // for solvePoisson
#include "../include/agoge/ParameterSystem.hpp"

namespace agoge {
namespace problems {

void GravityCollapse::registerParameters(ParameterSystem &params) const {
    std::string prefix = name() + ".";
    // Total mass in code units (default 1.0)
    params.addDefault(prefix + "mass_solar", "1.0");
    // Jeans radius in code units (default 1.0)
    params.addDefault(prefix + "r_jeans", "1.0");
    // Gravity solver method: "naive_dft" or "cooley_tukey" (default
    // "cooley_tukey")
    params.addDefault(prefix + "grav_method", "cooley_tukey");
}

void GravityCollapse::initialize(Field3D &Q, const ParameterSystem &params) {
    std::string prefix = name() + ".";
    // Read user-specified parameters.
    double M = params.getDouble(prefix + "mass_solar");
    double Rjean = params.getDouble(prefix + "r_jeans");
    double Gval = agoge::config::G;  // typically 1.0 in code units
    std::string gravMethod = params.getString(prefix + "grav_method");
    if (gravMethod.empty()) {
        gravMethod = "cooley_tukey";  // default if not specified
    }

    // Compute the density for a uniform sphere.
    double volumeSphere = (4.0 / 3.0) * M_PI * (Rjean * Rjean * Rjean);
    double rhoInside = M / volumeSphere;

    // Compute domain center from the bounding box.
    double xMid = (Q.bbox.xmax - Q.bbox.xmin) / 2.0;
    double yMid = (Q.bbox.ymax - Q.bbox.ymin) / 2.0;
    double zMid = (Q.bbox.zmax - Q.bbox.zmin) / 2.0;

    // Initialize density, velocities, energy, and potential.
    for (int k = 0; k < Q.Nz; k++) {
        for (int j = 0; j < Q.Ny; j++) {
            for (int i = 0; i < Q.Nx; i++) {
                int idx = Q.interiorIndex(i, j, k);
                // Compute cell center coordinates (relative to domain center).
                double xC = Q.xCenter(i) - xMid;
                double yC = Q.yCenter(j) - yMid;
                double zC = Q.zCenter(k) - zMid;
                double r2 = xC * xC + yC * yC + zC * zC;
                double R2 = Rjean * Rjean;

                // Inside sphere: uniform density; outside: a very low floor
                // value.
                if (r2 <= R2) {
                    Q.rho[idx] = rhoInside;
                } else {
                    Q.rho[idx] = 1.0e-10;
                }
                // Set velocities to zero (cold dust).
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;
                // Set energy to a near-zero floor.
                Q.E[idx] = 1.0e-10;
                // Initialize gravitational potential to zero.
                Q.phi[idx] = 0.0;
            }
        }
    }

    // Determine the gravity solver method.
    agoge::gravity::GravityMethod method =
        agoge::gravity::GravityMethod::NAIVE_DFT;
    if (gravMethod == "cooley_tukey") {
        method = agoge::gravity::GravityMethod::COOLEY_TUKEY;
    }
    std::cout << "[GravityCollapse] Using gravity solver method = "
              << gravMethod << "\n";

    // Time the Poisson solve.
    bool gravityEnabled = params.getBool("use_gravity");
    if (!gravityEnabled) {
        auto start = std::chrono::high_resolution_clock::now();
        agoge::gravity::solvePoisson(Q, method);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "[GravityCollapse] Gravity solver completed in "
                  << elapsed.count() << " seconds.\n";
    }

    // Compute L₁ and L₂ error norms comparing Q.phi with the analytic
    // uniform-sphere potential.
    double sumAbs = 0.0;
    double sumSqr = 0.0;
    long count = Q.Nx * Q.Ny * Q.Nz;
    for (int k = 0; k < Q.Nz; k++) {
        for (int j = 0; j < Q.Ny; j++) {
            for (int i = 0; i < Q.Nx; i++) {
                int idx = Q.interiorIndex(i, j, k);
                // Compute the radial distance from the domain center.
                double rr = std::sqrt(Q.xCenter(i) * Q.xCenter(i) +
                                      Q.yCenter(j) * Q.yCenter(j) +
                                      Q.zCenter(k) * Q.zCenter(k));
                double phiExact = 0.0;
                if (rr <= Rjean) {
                    // Analytic potential inside a uniform sphere.
                    double r2 = rr * rr;
                    phiExact = -(Gval * M) / (2.0 * Rjean * Rjean * Rjean) *
                               (3.0 * Rjean * Rjean - r2);
                } else {
                    // Analytic potential outside the sphere.
                    phiExact = -(Gval * M) / rr;
                }
                double phiNum = Q.phi[idx];
                double diff = phiNum - phiExact;
                sumAbs += std::fabs(diff);
                sumSqr += diff * diff;
            }
        }
    }
    double l1 = sumAbs / static_cast<double>(count);
    double l2 = std::sqrt(sumSqr / static_cast<double>(count));
    std::cout << "[GravityCollapse] Uniform-sphere potential check:\n"
              << "   L₁ = " << l1 << ", L₂ = " << l2 << " (over " << count
              << " cells)\n";
}

}  // namespace problems
}  // namespace agoge
