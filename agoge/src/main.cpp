#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

// Agoge core headers
#include "Config.hpp"
#include "EulerSolver.hpp"
#include "Field3d.hpp"
#include "GravitySolver.hpp"
#include "HDF5_IO.hpp"

// Problem-based approach: interface & registry
#include "../problems/Problem.hpp"
#include "../problems/ProblemRegistry.hpp"

// Performance monitoring
#include "PerformanceMonitor.hpp"

// Our parameter system
#include "BoundaryManager.hpp"
#include "ParameterSystem.hpp"

// Choose the gravity solver method globally or pass it in:
// agoge::gravity::GravityMethod method =
//     agoge::gravity::GravityMethod::COOLEY_TUKEY;

/**
 * @brief Helper to compute just the maximum wave speed (|u|+a) in Q
 *        (similar to computeTimeStep, but ignoring cfl).
 */
static double findMaxWaveSpeed(const agoge::Field3D& Q) {
    double maxSpeed = 0.0;
    double gamma_gas = agoge::config::gamma_gas;
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = Q.interiorIndex(i, j, k);
                double r = Q.rho[idx];
                if (r <= 0.0) continue;

                double ru = Q.rhou[idx];
                double rv = Q.rhov[idx];
                double rw = Q.rhow[idx];
                double e = Q.E[idx];

                double u = ru / r;
                double v = rv / r;
                double w = rw / r;
                double speed = std::sqrt(u * u + v * v + w * w);

                // pressure
                double p =
                    (gamma_gas - 1.0) * (e - 0.5 * r * (u * u + v * v + w * w));
                if (p < 0.0) continue;
                double a = std::sqrt(gamma_gas * p / r);
                double localWave = speed + a;
                if (localWave > maxSpeed) {
                    maxSpeed = localWave;
                }
            }
        }
    }
    return maxSpeed;
}

int main(int argc, char** argv) {
    // Start timing the entire main program
    agoge::PerformanceMonitor::instance().startTimer("main");

    // We might expect 1 argument:
    // 1) YAML file for parameters including problem name
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [yaml_file]\n"
                  << "Example: ./agoge_run Sod.yaml\n";
        return 1;
    }

    // Create a ParameterSystem with built-in defaults
    agoge::ParameterSystem params;

    // Register global parameters are already set in ParameterSystem's
    // constructor Register problem-specific parameters will be handled after
    // problem creation

    // Read global parameters first
    if (!params.readYAML(argv[1])) {
        std::cerr << "Failed to read configuration file.\n";
        return 1;
    }

    // Get the problem name
    std::string problem_name = params.getString("problem_name");
    if (problem_name.empty()) {
        std::cerr << "Problem name not specified in configuration.\n";
        return 1;
    }

    // Create the Problem instance
    std::unique_ptr<agoge::problems::Problem> problem =
        agoge::problems::createProblem(problem_name);
    if (!problem) {
        std::cerr << "Unknown problem name: " << problem_name << "\n";
        return 1;
    }

    // Register problem-specific parameters
    problem->registerParameters(params);

    // Re-read YAML to allow overriding problem-specific defaults
    if (!params.readYAML(argv[1])) {
        std::cerr << "Failed to re-read configuration file.\n";
        return 1;
    }

    // set the boundary conditions
    // agoge::BoundaryManager::initBCsFromParameters(params);

    std::string gravMethod = params.getString("GravityCollapse.grav_method");
    agoge::gravity::GravityMethod method =
        agoge::gravity::GravityMethod::NAIVE_DFT;
    if (gravMethod == "cooley_tukey") {
        method = agoge::gravity::GravityMethod::COOLEY_TUKEY;
    }

    bool doEulerUpdate = params.getBool("do_euler_update");

    // 2) Get Nx, Ny, Nz, domain, etc.
    int Nx = params.getInt("nx");
    int Ny = params.getInt("ny");
    int Nz = params.getInt("nz");

    double xmin = params.getDouble("xmin");
    double xmax = params.getDouble("xmax");
    double ymin = params.getDouble("ymin");
    double ymax = params.getDouble("ymax");
    double zmin = params.getDouble("zmin");
    double zmax = params.getDouble("zmax");

    agoge::BoundingBox bbox = {xmin, xmax, ymin, ymax, zmin, zmax};

    // 3) Initialize with the chosen problem
    agoge::Field3D Q(Nx, Ny, Nz, bbox, 1);
    problem->initialize(Q, params);

    bool gravityEnabled = params.getBool("use_gravity");
    std::cout << "Gravity is " << (gravityEnabled ? "ENABLED" : "DISABLED")
              << ", G: " << agoge::config::G << "\n";

    // 4) Set up boundary conditions (once), reading from param
    Q.bc_xmin = params.getBoundaryCondition("bc_xmin");
    Q.bc_xmax = params.getBoundaryCondition("bc_xmax");
    Q.bc_ymin = params.getBoundaryCondition("bc_ymin");
    Q.bc_ymax = params.getBoundaryCondition("bc_ymax");
    Q.bc_zmin = params.getBoundaryCondition("bc_zmin");
    Q.bc_zmax = params.getBoundaryCondition("bc_zmax");

    // 5) time stepping with "sound_crossings"
    double cflVal = params.getDouble("cfl");
    double crossingCount = params.getDouble("sound_crossings");
    std::cout << "CFL=" << cflVal << ", sound_crossings=" << crossingCount
              << "\n";

    // Compute initial max wave speed
    double initMaxSpeed = findMaxWaveSpeed(Q);
    if (initMaxSpeed < 1e-14) {
        initMaxSpeed = 1e-14;  // avoid divide by zero
    }
    double Lmax = std::max(
        {bbox.xmax - bbox.xmin, bbox.ymax - bbox.ymin, bbox.zmax - bbox.zmin});
    double crossingTime = Lmax / initMaxSpeed;
    double totalTime = crossingTime * crossingCount;

    std::cout << "Initial max wave speed= " << initMaxSpeed
              << ", crossingTime= " << crossingTime
              << ", totalTime= " << totalTime << "\n";

    agoge::io::writeFieldHDF5(Q, "agoge_init.h5");

    // Start main time loop, but in terms of totalTime
    agoge::PerformanceMonitor::instance().startTimer("timeLoop");

    double currentTime = 0.0;
    int step = 0;
    while (currentTime < totalTime) {
        // If gravity is on, solve Poisson
        if (gravityEnabled) {
            agoge::PerformanceMonitor::instance().startTimer("solvePoisson");
            agoge::gravity::solvePoisson(Q, method);
            agoge::PerformanceMonitor::instance().stopTimer("solvePoisson");
        }

        // Compute dt from Euler solver & cfl
        double dt = agoge::euler::computeTimeStep(Q, cflVal);

        // If dt is so tiny or zero => break
        if (dt < 1e-15) {
            std::cerr << "[main] dt is extremely small, stopping.\n";
            break;
        }

        // If next step would exceed totalTime, clamp dt
        if ((currentTime + dt) > totalTime) {
            dt = totalTime - currentTime;
        }

        // run one RK2 step
        if (doEulerUpdate) {
            agoge::PerformanceMonitor::instance().startTimer("EulerSolve");
            agoge::euler::runRK2(Q, dt);
            agoge::PerformanceMonitor::instance().stopTimer("EulerSolve");
        }
        currentTime += dt;
        step++;

        if (step % 2 == 0) {
            std::cout << "Step=" << step << ", time=" << currentTime << "/"
                      << totalTime << ", dt=" << dt << "\n";
        }
    }

    agoge::PerformanceMonitor::instance().stopTimer("timeLoop");

    // Example values; replace with actual simulation data
    long totalZones = Nx * Ny * Nz;  // Total number of zones in the domain

    // Set the steps and zones in the PerformanceMonitor
    agoge::PerformanceMonitor::instance().setSteps(step);
    agoge::PerformanceMonitor::instance().setZones(totalZones);

    // Output
    agoge::io::writeFieldHDF5(Q, "agoge_final.h5");
    std::cout << "Simulation finished. Final time=" << currentTime
              << ", step count=" << step << "\n";
    std::cout << "Final data written to agoge_final.h5\n";

    agoge::PerformanceMonitor::instance().stopTimer("main");
    agoge::PerformanceMonitor::instance().printReport();

    return 0;
}