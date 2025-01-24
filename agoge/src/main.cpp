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
#include "../problems/ProblemRegistry.hpp"

// Performance monitoring
#include "PerformanceMonitor.hpp"

// Our parameter system
#include "ParameterSystem.hpp"

// Choose the gravity solver method globally or pass it in:
agoge::gravity::GravityMethod method = agoge::gravity::GravityMethod::COOLEY_TUKEY;

int main(int argc, char** argv)
{
    // Start timing the entire main program
    agoge::PerformanceMonitor::instance().startTimer("main");

    // Create a ParameterSystem with built-in defaults
    agoge::ParameterSystem params;

    // We might expect 2 arguments:
    // 1) Problem name
    // 2) Optional YAML file for overrides
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <problemName> [yaml_file]\n"
                  << "Example: ./agoge_run sod input.yaml\n"
                  << "Or if no YAML file, just do: ./agoge_run sod\n";
        return 1;
    }

    // Arg1 => problem name
    std::string problemName = argv[1];

    // Arg2 => optional YAML
    if(argc >= 3) {
        std::string yamlFile = argv[2];
        bool ok = params.readYAML(yamlFile);
        if(!ok) {
            std::cerr << "[main] WARNING: Could not parse " << yamlFile << "\n";
        }
    }

    // 1) Create the problem
    auto problem = agoge::problems::createProblem(problemName);
    if(!problem) {
        std::cerr << "Error: Unrecognized problem name '" << problemName << "'\n";
        return 1;
    }
    std::cout << "Selected problem: " << problem->name() << "\n";

    // 2) Get Nx, Ny, Nz, etc. from ParameterSystem
    int Nx = params.getInt("nx");
    int Ny = params.getInt("ny");
    int Nz = params.getInt("nz");

    // For domain, we can use e.g. "domain: [1.0, 2.0, 3.0]"
    auto domainVec = params.getDoubleList("domain"); // default was [1.0,1.0,1.0]
    if(domainVec.size() < 3) {
        std::cerr << "[main] WARNING: domain list < 3 elements, using defaults.\n";
        domainVec = {1.0, 1.0, 1.0};
    }
    double Lx = domainVec[0];
    double Ly = domainVec[1];
    double Lz = domainVec[2];

    // 3) Allocate the Field3D
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double dz = Lz / Nz;
    agoge::Field3D Q(Nx, Ny, Nz, dx, dy, dz);

    // 4) Initialize with the chosen problem
    problem->initialize(Q);

    // 5) Gravity usage
    bool gravityEnabled = params.getBool("use_gravity");
    // or check problem->useGravity(), or combine them
    std::cout << "Gravity is " << (gravityEnabled ? "ENABLED" : "DISABLED") << "\n";

    // 6) Time stepping
    int nSteps = params.getInt("nsteps"); // default was 500
    double cflVal = params.getDouble("cfl"); // default 0.5

    std::cout << "Beginning time-stepping with nSteps=" << nSteps
              << ", CFL=" << cflVal << "\n";

    // Time loop
    agoge::PerformanceMonitor::instance().startTimer("timeLoop");
    for(int step = 0; step < nSteps; ++step) {
        if(gravityEnabled) {
            agoge::PerformanceMonitor::instance().startTimer("solvePoisson");
            agoge::gravity::solvePoisson(Q, method);
            agoge::PerformanceMonitor::instance().stopTimer("solvePoisson");
        }
        double dt = agoge::euler::computeTimeStep(Q, cflVal);
        agoge::euler::runRK2(Q, dt);

        if(step % 50 == 0) {
            std::cout << "Step " << step << "/" << nSteps
                      << ", dt=" << dt << "\n";
        }
    }
    agoge::PerformanceMonitor::instance().stopTimer("timeLoop");

    // Output
    agoge::io::writeFieldHDF5(Q, "agoge_final.h5");
    std::cout << "Simulation finished. Final data written to agoge_final.h5\n";

    // Stop overall main timer & print report
    agoge::PerformanceMonitor::instance().stopTimer("main");
    agoge::PerformanceMonitor::instance().printReport();

    return 0;
}
