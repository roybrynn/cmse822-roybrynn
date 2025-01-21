#include <iostream>
#include <cstdlib>     // for std::atoi
#include "Field3D.hpp"
#include "EulerSolver.hpp"
#include "GravitySolver.hpp"
#include "HDF5_IO.hpp"
#include "Config.hpp"

/**
 * @file main.cpp
 * @brief Main driver for the Agoge solver.
 *
 * Demonstrates:
 *  - Initializing a 3D domain (Nx, Ny, Nz) with cell sizes (dx, dy, dz)
 *  - Simple uniform initial condition with a small density perturbation
 *  - Optional self-gravity via solvePoissonFFT
 *  - Two-stage RK2 time integration
 *  - Writing final data to agoge_final.h5
 */

int main(int argc, char* argv[])
{
    // -------------------------------------------------------------------------
    // 1. Parse or define domain parameters
    //    This example uses command-line arguments or defaults for Nx, Ny, Nz,
    //    as well as domain length Lx, Ly, Lz.
    // -------------------------------------------------------------------------
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

    // Compute cell sizes
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double dz = Lz / Nz;

    std::cout << "Initializing Field3D with domain size: "
              << Nx << "x" << Ny << "x" << Nz << "\n"
              << "Physical lengths: " << Lx << " x " << Ly << " x " << Lz << "\n"
              << "Cell sizes: " << dx << ", " << dy << ", " << dz << "\n";

    // -------------------------------------------------------------------------
    // 2. Construct the Field3D object
    // -------------------------------------------------------------------------
    agoge::Field3D Q(Nx, Ny, Nz, dx, dy, dz);

    // -------------------------------------------------------------------------
    // 3. Simple initialization:
    //    - Uniform density = 1.0
    //    - Zero velocities (rhou, rhov, rhow)
    //    - Energy from p/(gamma-1) or trivial
    //    - Small perturbation in density in a spherical region for demonstration
    // -------------------------------------------------------------------------
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

                // Small spherical perturbation around center (0.5, 0.5, 0.5)
                double cx = 0.5 * Lx;
                double cy = 0.5 * Ly;
                double cz = 0.5 * Lz;
                double r2 = (x - cx)*(x - cx)
                          + (y - cy)*(y - cy)
                          + (z - cz)*(z - cz);

                double radius2 = 0.01; // (0.1)^2 if radius=0.1 * domain length
                double local_rho = (r2 < radius2) ? (rho0 + 0.2 * rho0) : rho0;

                Q.rho [idx] = local_rho;
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;
                // Scale the energy in the region of higher density, if desired
                Q.E   [idx] = E0 * (local_rho / rho0);
                Q.phi [idx] = 0.0; // Initially zero potential
            }
        }
    }

    // -------------------------------------------------------------------------
    // 4. Time-stepping
    // -------------------------------------------------------------------------
    int nSteps = 500;
    double dt  = 0.001; // A simple fixed time step

    std::cout << "Beginning time-stepping with nSteps=" << nSteps
              << ", dt=" << dt << "\n"
              << "Self-gravity is " << (agoge::config::use_gravity ? "ENABLED" : "DISABLED") << "\n";

    for(int step = 0; step < nSteps; ++step) {
        // 4a. (Optional) Solve Poisson's eqn for self-gravity if enabled
        if (agoge::config::use_gravity) {
            agoge::gravity::solvePoissonFFT(Q);
        }

        // 4b. RK2 integrator
        agoge::euler::runRK2(Q, dt);

        // 4c. Print progress
        if(step % 50 == 0) {
            std::cout << "Step " << step << " / " << nSteps
                      << " done, dt=" << dt << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // 5. Write final state to an HDF5 file
    // -------------------------------------------------------------------------
    std::string outFile = "agoge_final.h5";
    agoge::io::writeFieldHDF5(Q, outFile);

    std::cout << "Simulation finished. Final data written to: "
              << outFile << std::endl;

    return 0;
}
