#include <iostream>
#include <cstdlib> // for std::atoi
#include "Field3D.hpp"
#include "EulerSolver.hpp"
#include "GravitySolver.hpp"
#include "HDF5_IO.hpp"
#include "Config.hpp"

/**
 * @file main.cpp
 * @brief Main driver for the Agoge solver with adaptive time stepping.
 */

int main(int argc, char* argv[])
{
    // Parse domain parameters (optionally from cmd line)
    int Nx = 64, Ny = 64, Nz = 64;
    double Lx = 1.0, Ly = 1.0, Lz = 1.0;
    if (argc > 3) {
        Nx = std::atoi(argv[1]);
        Ny = std::atoi(argv[2]);
        Nz = std::atoi(argv[3]);
    }
    if (argc > 6) {
        Lx = std::atof(argv[4]);
        Ly = std::atof(argv[5]);
        Lz = std::atof(argv[6]);
    }

    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double dz = Lz / Nz;

    std::cout << "Initializing Field3D with domain size: "
              << Nx << "x" << Ny << "x" << Nz << "\n"
              << "Physical lengths: " << Lx << " x " << Ly << " x " << Lz << "\n"
              << "Cell sizes: " << dx << ", " << dy << ", " << dz << "\n";

    // Create Field3D
    agoge::Field3D Q(Nx, Ny, Nz, dx, dy, dz);

    // Simple initialization
    double rho0 = 1.0;
    double p0   = 1.0;
    double E0   = p0 / (agoge::config::gamma_gas - 1.0);

    for(int k = 0; k < Nz; ++k) {
        for(int j = 0; j < Ny; ++j) {
            for(int i = 0; i < Nx; ++i) {
                int idx = Q.index(i,j,k);

                double x = (i + 0.5) * dx;
                double y = (j + 0.5) * dy;
                double z = (k + 0.5) * dz;

                double cx = 0.5 * Lx;
                double cy = 0.5 * Ly;
                double cz = 0.5 * Lz;
                double r2 = (x - cx)*(x - cx)
                          + (y - cy)*(y - cy)
                          + (z - cz)*(z - cz);
                double radius2 = 0.01; // e.g., sphere radius=0.1 => 0.1^2=0.01

                double local_rho = (r2 < radius2) ? (rho0 + 0.2*rho0) : rho0;

                Q.rho [idx] = local_rho;
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;
                Q.E   [idx] = E0 * (local_rho / rho0);
                Q.phi [idx] = 0.0;
            }
        }
    }

    // Time-stepping
    int nSteps = 500;
    double cflVal = 0.5; // for example
    std::cout << "Beginning time-stepping with nSteps=" << nSteps
              << ", cfl=" << cflVal << "\n"
              << "Self-gravity is "
              << (agoge::config::use_gravity ? "ENABLED" : "DISABLED") << "\n";

    for(int step = 0; step < nSteps; ++step) {
        // If gravity is enabled, solve Poisson
        if (agoge::config::use_gravity) {
            agoge::gravity::solvePoissonFFT(Q);
        }

        // Compute adaptive dt from the CFL condition
        double dt = agoge::euler::computeTimeStep(Q, cflVal);

        // Perform a two-stage RK2 step
        agoge::euler::runRK2(Q, dt);

        if(step % 50 == 0) {
            std::cout << "Step " << step << " / " << nSteps
                      << ", dt=" << dt << "\n";
        }
    }

    // Write final state
    std::string outFile = "agoge_final.h5";
    agoge::io::writeFieldHDF5(Q, outFile);

    std::cout << "Simulation finished. Final data written to: "
              << outFile << std::endl;

    return 0;
}
