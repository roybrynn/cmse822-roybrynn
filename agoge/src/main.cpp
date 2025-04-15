#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>  // already included above if not, ensure it is

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
#include "ParameterSystem.hpp"

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

    MPI_Init(&argc, &argv);

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    std::string gravMethod = params.getString("GravityCollapse.grav_method");
    agoge::gravity::GravityMethod method =
        agoge::gravity::GravityMethod::NAIVE_DFT;
    if (gravMethod == "cooley_tukey") {
        method = agoge::gravity::GravityMethod::COOLEY_TUKEY;
    }

    bool doEulerUpdate = params.getBool("do_euler_update");
    bool doIO = params.getBool("do_io");

    const int scrn_out_freq = params.getInt("screen_out_interval");
    
    // Get the output directory parameter
    std::string outputDir = params.getString("output_dir");

    // 2) Get GLOBAL Nx, Ny, Nz, domain, etc.
    const int global_Nx = params.getInt("nx");
    const int global_Ny = params.getInt("ny");
    const int global_Nz = params.getInt("nz");
    const int nghost = params.getInt("nghost");

    double global_xmin = params.getDouble("xmin");
    double global_xmax = params.getDouble("xmax");
    double global_ymin = params.getDouble("ymin");
    double global_ymax = params.getDouble("ymax");
    double global_zmin = params.getDouble("zmin");
    double global_zmax = params.getDouble("zmax");

    // --- Begin explicit MPI domain decomposition (no MPI_Cart APIs) ---
    // Compute factors (Px, Py, Pz) for a 3D subdomain decomposition.
    // We use a simple heuristic: find greatest factor for Px s.t. (size % Px ==
    // 0), then for Py and Pz.
    int Px = 1, Py = 1, Pz = 1;
    {
        // Determine Px by scanning downward from cube root
        int cube = round(std::pow(size, 1.0 / 3.0));
        for (int i = cube; i >= 1; i--) {
            if (size % i == 0) {
                Px = i;
                break;
            }
        }
        int rem = size / Px;
        int sq = round(std::sqrt(rem));
        for (int i = sq; i >= 1; i--) {
            if (rem % i == 0) {
                Py = i;
                break;
            }
        }
        Pz = rem / Py;
    }
    // Compute explicit coordinates (myI, myJ, myK)
    int myI = rank % Px;
    int myJ = (rank / Px) % Py;
    int myK = rank / (Px * Py);

    // Compute local grid sizes with remainder distribution
    int local_Nx = global_Nx / Px;
    int remX = global_Nx % Px;
    if (myI < remX)
        local_Nx++;  // distribute extra cells in lower-index subdomains
    int local_Ny = global_Ny / Py;
    int remY = global_Ny % Py;
    if (myJ < remY) local_Ny++;
    int local_Nz = global_Nz / Pz;
    int remZ = global_Nz % Pz;
    if (myK < remZ) local_Nz++;

    // Compute physical cell sizes (global uniform grid)
    double dxGlobal = (global_xmax - global_xmin) / global_Nx;
    double dyGlobal = (global_ymax - global_ymin) / global_Ny;
    double dzGlobal = (global_zmax - global_zmin) / global_Nz;

    // Compute offsets: number of cells preceding current subdomain in each
    // direction
    int offsetX = myI * (global_Nx / Px) + std::min(myI, remX);
    int offsetY = myJ * (global_Ny / Py) + std::min(myJ, remY);
    int offsetZ = myK * (global_Nz / Pz) + std::min(myK, remZ);

    // Compute local bounding box
    double local_xmin = global_xmin + offsetX * dxGlobal;
    double local_xmax = local_xmin + local_Nx * dxGlobal;
    double local_ymin = global_ymin + offsetY * dyGlobal;
    double local_ymax = local_ymin + local_Ny * dyGlobal;
    double local_zmin = global_zmin + offsetZ * dzGlobal;
    double local_zmax = local_zmin + local_Nz * dzGlobal;
    agoge::BoundingBox localBox = {local_xmin, local_xmax, local_ymin,
                                   local_ymax, local_zmin, local_zmax};
    // --- End explicit MPI domain decomposition ---
    

    // Update Field3D initialization with local sizes and local bounding box:
    agoge::Field3D Q(local_Nx, local_Ny, local_Nz, localBox);

    // --- Set metadata for global domain reconstruction ---
    Q.global_bbox = {global_xmin, global_xmax, global_ymin,
                     global_ymax, global_zmin, global_zmax};
    Q.global_Nx = global_Nx;
    Q.global_Ny = global_Ny;
    Q.global_Nz = global_Nz;
    Q.Px = Px;
    Q.Py = Py;
    Q.Pz = Pz;
    Q.subdomain_x = myI;
    Q.subdomain_y = myJ;
    Q.subdomain_z = myK;
    Q.myRank = rank;
    Q.mpiSize = size;

    // Step 1: Standard neighbor assignment for all ranks (works for both outflow and periodic)
    // X-direction neighbors
    Q.rankMinusX = (myI > 0) ? (rank - 1) : MPI_PROC_NULL;
    Q.rankPlusX = (myI < Px - 1) ? (rank + 1) : MPI_PROC_NULL;

    // Y-direction neighbors
    Q.rankMinusY = (myJ > 0) ? (rank - Px) : MPI_PROC_NULL;
    Q.rankPlusY = (myJ < Py - 1) ? (rank + Px) : MPI_PROC_NULL;

    // Z-direction neighbors
    Q.rankMinusZ = (myK > 0) ? (rank - (Px * Py)) : MPI_PROC_NULL;
    Q.rankPlusZ = (myK < Pz - 1) ? (rank + (Px * Py)) : MPI_PROC_NULL;

    // Step 2: Override only boundary connections for periodic BCs
    // Only update the domain boundaries that are periodic, leaving outflow as MPI_PROC_NULL
    if (params.getBoundaryCondition("bc_xmin") == agoge::config::BoundaryCondition::PERIODIC && myI == 0) {
        Q.rankMinusX = myK * (Px * Py) + myJ * Px + (Px - 1);
    }
    if (params.getBoundaryCondition("bc_xmax") == agoge::config::BoundaryCondition::PERIODIC && myI == Px - 1) {
        Q.rankPlusX = myK * (Px * Py) + myJ * Px + 0;
    }
    if (params.getBoundaryCondition("bc_ymin") == agoge::config::BoundaryCondition::PERIODIC && myJ == 0) {
        Q.rankMinusY = myK * (Px * Py) + (Py - 1) * Px + myI;
    }
    if (params.getBoundaryCondition("bc_ymax") == agoge::config::BoundaryCondition::PERIODIC && myJ == Py - 1) {
        Q.rankPlusY = myK * (Px * Py) + 0 * Px + myI;
    }
    if (params.getBoundaryCondition("bc_zmin") == agoge::config::BoundaryCondition::PERIODIC && myK == 0) {
        Q.rankMinusZ = (Pz - 1) * (Px * Py) + myJ * Px + myI;
    }
    if (params.getBoundaryCondition("bc_zmax") == agoge::config::BoundaryCondition::PERIODIC && myK == Pz - 1) {
        Q.rankPlusZ = 0 * (Px * Py) + myJ * Px + myI;
    }

    Q.allocateMPIBuffers();

    // Print out neighbor ranks for debugging
    if (rank == 0) {
        std::cout << "Global domain: Nx=" << global_Nx << ", Ny=" << global_Ny
                  << ", Nz=" << global_Nz << "\n";
        std::cout << "Domain decomposition: Px=" << Px << ", Py=" << Py
                  << ", Pz=" << Pz << "\n";
        std::cout << "Rank 0 neighbors: "
                  << "rankMinusX=" << Q.rankMinusX << ", rankPlusX=" << Q.rankPlusX
                  << ", rankMinusY=" << Q.rankMinusY << ", rankPlusY=" << Q.rankPlusY
                  << ", rankMinusZ=" << Q.rankMinusZ << ", rankPlusZ=" << Q.rankPlusZ
                  << "\n";
    }

    // 3) Initialize with the chosen problem (uses local Q)
    problem->initialize(Q, params);

    bool gravityEnabled = params.getBool("use_gravity");

    // 4) Set up boundary conditions (once), reading from param
    // Replace boundary condition setting with initBCsFromParameters method
    Q.initBCsFromParameters(params);

    // 5) time stepping with "sound_crossings"
    double cflVal = params.getDouble("cfl");
    double crossingCount = params.getDouble("sound_crossings");

    // Compute initial max wave speed
    double initMaxSpeed = findMaxWaveSpeed(Q);
    if (initMaxSpeed < 1e-14) {
        initMaxSpeed = 1e-14;  // avoid divide by zero
    }
    double Lmax =
        std::max({localBox.xmax - localBox.xmin, localBox.ymax - localBox.ymin,
                  localBox.zmax - localBox.zmin});
    double crossingTime = Lmax / initMaxSpeed;
    double totalTime = params.getDouble("t_max");

    // Initial output (epoch 0) - replaced with performFieldIO
    if (doIO) {
        agoge::io::performFieldIO(Q, problem_name, rank, 0.0, outputDir);
    }

    // Start main time loop, but in terms of totalTime
    agoge::PerformanceMonitor::instance().setRank(rank);
    agoge::PerformanceMonitor::instance().setCommSize(size);
    agoge::PerformanceMonitor::instance().startTimer("timeLoop");

    double currentTime = 0.0;
    int step = 0;
    double dt = params.getDouble("dt_init");
    dt = std::min(dt, agoge::euler::computeTimeStep(Q, cflVal));
    {
        // Synchronize dt across all ranks
        double global_dt;
        MPI_Allreduce(&dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dt = global_dt;
    }

    // Record wall clock start and get maximum allowed wall clock time
    auto wallStart = std::chrono::high_resolution_clock::now();
    double maxWallTime = params.getDouble("max_wallclock_time");

    while (currentTime < totalTime) {
        // Check wall clock time in an MPI-safe manner
        auto wallNow = std::chrono::high_resolution_clock::now();
        double wallElapsed = std::chrono::duration<double>(wallNow - wallStart).count();
        int localStop = (wallElapsed >= maxWallTime) ? 1 : 0;
        int globalStop = 0;
        MPI_Allreduce(&localStop, &globalStop, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (globalStop != 0) {
            if (rank == 0) {
                std::cout << "[main] Max wall clock time reached (" << wallElapsed
                          << " sec), stopping simulation.\n";
            }
            break;
        }

        // If gravity is on, solve Poisson
        if (gravityEnabled) {
            agoge::PerformanceMonitor::instance().startTimer("solvePoisson");
            agoge::gravity::solvePoisson(Q, method);
            agoge::PerformanceMonitor::instance().stopTimer("solvePoisson");
        }

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

        // Compute new dt and synchronize across all ranks:
        dt = std::min(1.2 * dt, agoge::euler::computeTimeStep(Q, cflVal));
        {
            double global_dt;
            MPI_Allreduce(&dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            dt = global_dt;
        }

        if (step % scrn_out_freq == 0 && rank == 0) {
            std::cout << "Step=" << step << ", time=" << currentTime << "/"
                      << totalTime << ", dt=" << dt << "\n";
        }
        // Compute dt from Euler solver & cfl for the next step
        dt = std::min(1.2 * dt, agoge::euler::computeTimeStep(Q, cflVal));
    }

    agoge::PerformanceMonitor::instance().stopTimer("timeLoop");

    // Example values; replace with actual simulation data
    long totalZones =
        global_Nx * global_Ny * global_Nz;  // total zones in the domain

    // Set the steps and zones in the PerformanceMonitor
    agoge::PerformanceMonitor::instance().setSteps(step);
    agoge::PerformanceMonitor::instance().setZones(totalZones);

    // Final output (e.g., epoch 1) - replaced with performFieldIO
    if (doIO) {
        agoge::io::performFieldIO(Q, problem_name, rank, currentTime, outputDir);
    }

    if (rank == 0) {
        std::cout << "Simulation finished. Final time=" << currentTime
                  << ", step count=" << step << "\n";
    }

    agoge::PerformanceMonitor::instance().stopTimer("main");
    agoge::PerformanceMonitor::instance().compileReport();
    MPI_Finalize();
    return 0;
}