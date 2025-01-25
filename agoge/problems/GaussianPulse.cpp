#include "GaussianPulse.hpp"

#include <cmath>
#include <iostream>

#include "../include/agoge/Config.hpp"

namespace agoge {
namespace problems {

void GaussianPulse::initialize(Field3D &Q) {
    // We'll set:
    //   - uniform velocity in x-direction (u0 = 1, for example)
    //   - uniform pressure = 1
    //   - base density = 1
    //   - an added Gaussian "bump" in density, e.g. amplitude=0.1, sigma=0.05
    //
    // The domain is [0..(Nx*dx)] x [0..(Ny*dy)] x [0..(Nz*dz)] with periodic
    // BC. If user sets Nx>1, or Ny=1 => 1D, etc. It's all consistent.

    double u0 = 1.0;  // advection velocity in x
    double baseRho = 1.0;
    double baseP = 1.0;
    double gamma = agoge::config::gamma_gas;

    double amp = 0.1;     // amplitude of the Gaussian
    double sigma = 0.05;  // width

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