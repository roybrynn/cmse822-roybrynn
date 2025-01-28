// problems/GaussianPulse.cpp
#include "GaussianPulse.hpp"

#include <cmath>
#include <iostream>

namespace agoge {
namespace problems {

void GaussianPulse::registerParameters(ParameterSystem &params) const {
    // Prefix parameters with the problem name to avoid key collisions
    std::string prefix = name() + ".";

    // Define default values for problem-specific parameters
    params.addDefault(prefix + "amplitude", "0.1");
    params.addDefault(prefix + "sigma", "0.05");
    params.addDefault(prefix + "u0", "1.0");  // Advection velocity
    params.addDefault(prefix + "base_rho", "1.0");
    params.addDefault(prefix + "base_p", "1.0");
    params.addDefault(prefix + "gamma", "1.4");  // Assuming ideal gas
}

void GaussianPulse::initialize(Field3D &Q, const ParameterSystem &params) {
    // Retrieve problem-specific parameters using prefix
    std::string prefix = name() + ".";

    double amp = params.getDouble(prefix + "amplitude");
    double sigma = params.getDouble(prefix + "sigma");
    double u0 = params.getDouble(prefix + "u0");
    double baseRho = params.getDouble(prefix + "base_rho");
    double baseP = params.getDouble(prefix + "base_p");
    double gamma = params.getDouble(prefix + "gamma");

    // center
    double cx = 0.5 * Q.Nx * Q.dx;
    double cy = 0.5 * Q.Ny * Q.dy;
    double cz = 0.5 * Q.Nz * Q.dz;

    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    // The total energy E = p/(gamma-1) + 0.5 * rho * (u^2 + v^2 + w^2).
    // We'll do v=0, w=0, uniform u0 => easy.

    for (int k = 0; k < Nz; k++) {
        double zC = (k + 0.5) * dz;
        for (int j = 0; j < Ny; j++) {
            double yC = (j + 0.5) * dy;
            for (int i = 0; i < Nx; i++) {
                double xC = (i + 0.5) * dx;
                int idx = Q.index(i, j, k);

                // distance^2 from center
                double dx_ = xC - cx;
                double dy_ = yC - cy;
                double dz_ = zC - cz;
                double r2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;

                // Gaussian bump
                double localRho =
                    baseRho + amp * std::exp(-r2 / (2.0 * sigma * sigma));
                Q.rho[idx] = localRho;

                // uniform velocity
                double u = u0;
                double v = 0.0;
                double w = 0.0;

                Q.rhou[idx] = localRho * u;
                Q.rhov[idx] = localRho * v;
                Q.rhow[idx] = localRho * w;

                // uniform pressure
                double p = baseP;
                double eInt = p / (gamma - 1.0);
                double eKin = 0.5 * localRho * (u * u + v * v + w * w);
                Q.E[idx] = eInt + eKin;

                // no gravity => phi=0
                Q.phi[idx] = 0.0;
            }
        }
    }

    std::cout << "[GaussianPulse] Initialized a Gaussian density bump, "
              << "u0=" << u0 << ", p=" << baseP << ", gamma=" << gamma
              << ", amplitude=" << amp << ", sigma=" << sigma << "\n";
}

}  // namespace problems
}  // namespace agoge