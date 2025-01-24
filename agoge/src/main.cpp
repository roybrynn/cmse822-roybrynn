#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

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

/** 
 * @brief Helper to compute just the maximum wave speed (|u|+a) in Q 
 *        (similar to computeTimeStep, but ignoring cfl).
 */
static double findMaxWaveSpeed(const agoge::Field3D &Q)
{
    double maxSpeed = 0.0;
    double gamma_gas = agoge::config::gamma_gas;
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = Q.index(i, j, k);
                double r = Q.rho[idx];
                if (r <= 0.0) continue;

                double ru = Q.rhou[idx];
                double rv = Q.rhov[idx];
                double rw = Q.rhow[idx];
                double e  = Q.E[idx];

                double u = ru / r;
                double v = rv / r;
                double w = rw / r;
                double speed = std::sqrt(u*u + v*v + w*w);

                // pressure
                double p = (gamma_gas - 1.0)*( e - 0.5*r*(u*u + v*v + w*w) );
                if(p < 0.0) continue;
                double a = std::sqrt(gamma_gas * p / r);
                double localWave = speed + a;
                if(localWave > maxSpeed) {
                    maxSpeed = localWave;
                }
            }
        }
    }
    return maxSpeed;
}

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

    // 2) Get Nx, Ny, Nz, domain, etc.
    int Nx = params.getInt("nx");
    int Ny = params.getInt("ny");
    int Nz = params.getInt("nz");

    auto domainVec = params.getDoubleList("domain");
    if(domainVec.size() < 3) {
        std::cerr << "[main] WARNING: domain list < 3 elements, fallback [1,1,1].\n";
        domainVec = {1.0, 1.0, 1.0};
    }
    double Lx = domainVec[0];
    double Ly = domainVec[1];
    double Lz = domainVec[2];

    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double dz = Lz / Nz;
    agoge::Field3D Q(Nx, Ny, Nz, dx, dy, dz);

    // 3) Initialize with the chosen problem
    problem->initialize(Q);

    bool gravityEnabled = params.getBool("use_gravity");
    std::cout << "Gravity is " << (gravityEnabled ? "ENABLED" : "DISABLED") << "\n";

    // 4) time stepping with "sound_crossings"
    double cflVal = params.getDouble("cfl");
    double crossingCount = params.getDouble("sound_crossings");
    std::cout << "CFL=" << cflVal 
              << ", sound_crossings=" << crossingCount << "\n";

    // Compute initial max wave speed
    double initMaxSpeed = findMaxWaveSpeed(Q);
    if(initMaxSpeed < 1e-14) {
        initMaxSpeed = 1e-14; // avoid divide by zero
    }
    double Lmax = std::max({Lx, Ly, Lz});
    double crossingTime = Lmax / initMaxSpeed;
    double totalTime = crossingTime * crossingCount;

    std::cout << "Initial max wave speed= " << initMaxSpeed
              << ", crossingTime= " << crossingTime
              << ", totalTime= " << totalTime << "\n";

    // Start main time loop, but in terms of totalTime
    agoge::PerformanceMonitor::instance().startTimer("timeLoop");

    double currentTime = 0.0;
    int step = 0;
    while(currentTime < totalTime) {
        // If gravity is on, solve Poisson
        if(gravityEnabled) {
            agoge::PerformanceMonitor::instance().startTimer("solvePoisson");
            agoge::gravity::solvePoisson(Q, method);
            agoge::PerformanceMonitor::instance().stopTimer("solvePoisson");
        }

        // Compute dt from Euler solver & cfl
        double dt = agoge::euler::computeTimeStep(Q, cflVal);

        // If dt is so tiny or zero => break
        if(dt < 1e-15) {
            std::cerr << "[main] dt is extremely small, stopping.\n";
            break;
        }

        // If next step would exceed totalTime, clamp dt
        if( (currentTime + dt) > totalTime ) {
            dt = totalTime - currentTime;
        }

        // run one RK2 step
        agoge::euler::runRK2(Q, dt);
        currentTime += dt;
        step++;

        if(step % 2 == 0) {
            std::cout << "Step=" << step 
                      << ", time=" << currentTime 
                      << "/" << totalTime 
                      << ", dt=" << dt << "\n";
        }
    }

    agoge::PerformanceMonitor::instance().stopTimer("timeLoop");

    // Output
    agoge::io::writeFieldHDF5(Q, "agoge_final.h5");
    std::cout << "Simulation finished. Final time=" << currentTime
              << ", step count=" << step << "\n";
    std::cout << "Final data written to agoge_final.h5\n";

    agoge::PerformanceMonitor::instance().stopTimer("main");
    agoge::PerformanceMonitor::instance().printReport();

    return 0;
}
