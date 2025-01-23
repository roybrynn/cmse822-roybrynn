#include <iostream>
#include <string>
#include <cstdlib>

// Agoge core headers
#include "Field3d.hpp"
#include "EulerSolver.hpp"
#include "GravitySolver.hpp"
#include "HDF5_IO.hpp"
#include "Config.hpp"

// Problem-based approach: interface & registry
#include "../problems/Problem.hpp"
#include "../problems/ProblemRegistry.cpp"

// Performance monitoring
#include "PerformanceMonitor.hpp"

// Choose the gravity solver method
agoge::gravity::GravityMethod method = agoge::gravity::GravityMethod::COOLEY_TUKEY;
// Alternatively:
// agoge::gravity::GravityMethod method = agoge::gravity::GravityMethod::NAIVE_DFT;

int main(int argc, char** argv)
{
    // Start timing the entire main program
    agoge::PerformanceMonitor::instance().startTimer("main");

    // 1. Check usage: we expect at least 1 argument = problem name
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <problemName> [Nx Ny Nz Lx Ly Lz]\n"
                  << "Example: ./agoge_run sod 64 1 1 1.0 1.0 1.0\n"
                  << "Available problems might be: sod, collapse, etc.\n";
        return 1;
    }

    // 2. Parse problem name
    std::string problemName = argv[1];

    // 3. Parse optional grid and domain arguments
    int Nx = 64, Ny = 64, Nz = 64;
    double Lx = 1.0, Ly = 1.0, Lz = 1.0;

    // If we have at least 4 additional args, interpret them as Nx, Ny, Nz
    if (argc > 4) {
        Nx = std::atoi(argv[2]);
        Ny = std::atoi(argv[3]);
        Nz = std::atoi(argv[4]);
    }
    // If we have at least 7 additional args, interpret them as Lx, Ly, Lz
    if (argc > 7) {
        Lx = std::atof(argv[5]);
        Ly = std::atof(argv[6]);
        Lz = std::atof(argv[7]);
    }

    // 4. Create a problem instance via the registry
    auto problem = agoge::problems::createProblem(problemName);
    if (!problem) {
        std::cerr << "Error: Unrecognized problem name '" << problemName << "'\n";
        return 1;
    }

    std::cout << "Selected problem: " << problem->name() << "\n";

    // 5. Allocate the Field3D
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double dz = Lz / Nz;

    agoge::Field3D Q(Nx, Ny, Nz, dx, dy, dz);

    // 6. Initialize the field with the chosen problem setup
    problem->initialize(Q);

    // 7. Decide if we use gravity
    bool gravityEnabled = problem->useGravity();
    std::cout << "Gravity is " << (gravityEnabled ? "ENABLED" : "DISABLED") << "\n";

    // 8. Basic time-stepping setup
    int nSteps = 500;
    double cflVal = 0.5;
    std::cout << "Beginning time-stepping with nSteps=" << nSteps
              << ", CFL=" << cflVal << "\n";

    // (Optional) Start a timer for the main time loop
    agoge::PerformanceMonitor::instance().startTimer("timeLoop");

    // 9. Main time loop
    for (int step = 0; step < nSteps; ++step) {
        // If gravity is on, solve Poisson's equation
        if (gravityEnabled) {
            agoge::PerformanceMonitor::instance().startTimer("solvePoisson");
            agoge::gravity::solvePoisson(Q, method);
            agoge::PerformanceMonitor::instance().stopTimer("solvePoisson");
        }

        // Compute adaptive dt from Euler solver & CFL
        double dt = agoge::euler::computeTimeStep(Q, cflVal);

        // Perform a 2-stage RK2 update
        agoge::euler::runRK2(Q, dt);

        // Optional progress output
        if (step % 50 == 0) {
            std::cout << "Step " << step << "/" << nSteps
                      << ", dt = " << dt << "\n";
        }
    }

    // (Optional) Stop the time loop timer
    agoge::PerformanceMonitor::instance().stopTimer("timeLoop");

    // 10. Write final state to HDF5
    agoge::io::writeFieldHDF5(Q, "agoge_final.h5");
    std::cout << "Simulation finished. Final data written to agoge_final.h5\n";

    // Stop overall main timer
    agoge::PerformanceMonitor::instance().stopTimer("main");

    // Print a performance report from PerformanceMonitor
    agoge::PerformanceMonitor::instance().printReport();

    return 0;
}
