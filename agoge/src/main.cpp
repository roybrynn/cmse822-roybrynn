#include <iostream>
#include <string>
#include <agoge/Config.hpp>
#include <agoge/Field3D.hpp>
#include <agoge/EulerSolver.hpp>
#include <agoge/GravitySolver.hpp>
#include <agoge/HDF5_IO.hpp>

int main(int argc, char** argv)
{
    // Parse args or config
    // E.g., read Nx, Ny, Nz, or use defaults from Config.hpp

    // Create Field3D
    Field3D Q(Nx, Ny, Nz, dx, dy, dz);

    // Initialize Q
    initialize(Q);

    // Time-stepping
    for(int step=0; step<nSteps; step++){
        // If gravity, solve Poisson
#ifdef AGOGE_USE_GRAVITY
        solvePoissonFFT(Q, phiArray);
#endif

        // Compute L(Q) then do RK2 or other integration
        // e.g. runRK2(Q, dt);
    }

    // Write final data
    writeFieldHDF5(Q, "agoge_final.h5");
    return 0;
}

//------------------------------------------------------------
// Main function
//------------------------------------------------------------
int main()
{
    Field3D Q(Nx, Ny, Nz, dx, dy, dz);
    Field3D Qtemp(Nx, Ny, Nz, dx, dy, dz);

    // L(Q) arrays
    std::vector<double> L_rho (Nx*Ny*Nz, 0.0);
    std::vector<double> L_rhou(Nx*Ny*Nz, 0.0);
    std::vector<double> L_rhov(Nx*Ny*Nz, 0.0);
    std::vector<double> L_rhow(Nx*Ny*Nz, 0.0);
    std::vector<double> L_E   (Nx*Ny*Nz, 0.0);

    std::vector<double> L_rho_temp (Nx*Ny*Nz, 0.0);
    std::vector<double> L_rhou_temp(Nx*Ny*Nz, 0.0);
    std::vector<double> L_rhov_temp(Nx*Ny*Nz, 0.0);
    std::vector<double> L_rhow_temp(Nx*Ny*Nz, 0.0);
    std::vector<double> L_E_temp   (Nx*Ny*Nz, 0.0);

    // Initialize with a simple perturbation
    double rho0 = 1.0;
    double p0   = 1.0;
    double E0   = p0 / (gamma_gas - 1.0); // no velocity
    for(int k = 0; k < Nz; k++) {
        for(int j = 0; j < Ny; j++) {
            for(int i = 0; i < Nx; i++) {
                int idx = Q.index(i,j,k);

                double x = (i + 0.5)*dx;
                double y = (j + 0.5)*dy;
                double z = (k + 0.5)*dz;

                double r   = std::sqrt((x-0.5)*(x-0.5) + 
                                       (y-0.5)*(y-0.5) +
                                       (z-0.5)*(z-0.5));
                double rho = (r < 0.1) ? (rho0 + 0.2*rho0) : rho0;

                Q.rho [idx] = rho;
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;
                // Scale E with density in the perturbation region
                Q.E   [idx] = E0 * (rho / rho0);
            }
        }
    }

    // Time-stepping
    const int nSteps = 500;
    for(int step = 0; step < nSteps; step++)
    {
        // 1. Determine stable dt from max wave speed
        double maxSpeed = 0.0;
        // #pragma omp parallel for reduction(max:maxSpeed)
        for(int k = 0; k < Nz; k++) {
            for(int j = 0; j < Ny; j++) {
                for(int i = 0; i < Nx; i++) {
                    int idx = Q.index(i,j,k);
                    double s = maxWaveSpeed(
                        Q.rho [idx],
                        Q.rhou[idx],
                        Q.rhov[idx],
                        Q.rhow[idx],
                        Q.E   [idx]);
                    maxSpeed = std::max(maxSpeed, s);
                }
            }
        }
        double dt = CFL * std::min(std::min(dx, dy), dz) / maxSpeed;

        // 2. RK2 Stage 1: Qtemp = Q + dt * L(Q)
        computeL(Q, L_rho, L_rhou, L_rhov, L_rhow, L_E);

        // #pragma omp parallel for
        for(int k = 0; k < Nz; k++) {
            for(int j = 0; j < Ny; j++) {
                for(int i = 0; i < Nx; i++) {
                    int idx = Q.index(i,j,k);
                    Qtemp.rho [idx] = Q.rho [idx]  + dt * L_rho [idx];
                    Qtemp.rhou[idx] = Q.rhou[idx] + dt * L_rhou[idx];
                    Qtemp.rhov[idx] = Q.rhov[idx] + dt * L_rhov[idx];
                    Qtemp.rhow[idx] = Q.rhow[idx] + dt * L_rhow[idx];
                    Qtemp.E   [idx] = Q.E   [idx] + dt * L_E   [idx];
                }
            }
        }

        // 3. RK2 Stage 2: Q^{n+1} = 0.5 * [ Q + Qtemp + dt * L(Qtemp) ]
        computeL(Qtemp, L_rho_temp, L_rhou_temp, L_rhov_temp, L_rhow_temp, L_E_temp);

        // #pragma omp parallel for
        for(int k = 0; k < Nz; k++) {
            for(int j = 0; j < Ny; j++) {
                for(int i = 0; i < Nx; i++) {
                    int idx = Q.index(i,j,k);
                    Q.rho [idx]  = 0.5*( Q.rho [idx]  + Qtemp.rho [idx]
                                     + dt * L_rho_temp [idx] );
                    Q.rhou[idx]  = 0.5*( Q.rhou[idx]  + Qtemp.rhou[idx]
                                     + dt * L_rhou_temp[idx] );
                    Q.rhov[idx]  = 0.5*( Q.rhov[idx]  + Qtemp.rhov[idx]
                                     + dt * L_rhov_temp[idx] );
                    Q.rhow[idx]  = 0.5*( Q.rhow[idx]  + Qtemp.rhow[idx]
                                     + dt * L_rhow_temp[idx] );
                    Q.E   [idx]  = 0.5*( Q.E   [idx]  + Qtemp.E   [idx]
                                     + dt * L_E_temp   [idx] );
                }
            }
        }

        // Print info every 50 steps
        if(step % 50 == 0) {
            std::cout << "Step = " << step 
                      << ", dt = " << dt 
                      << ", maxSpeed = " << maxSpeed << std::endl;
        }
    }

    // Write final fields to HDF5
    writeFieldHDF5(Q, "agoge_final.h5");

    std::cout << "Simulation finished. Final data written to agoge_final.h5\n";
    return 0;
}
