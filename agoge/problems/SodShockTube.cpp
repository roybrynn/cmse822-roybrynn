// problems/SodShockTube.cpp
#include "SodShockTube.hpp"
#include "Config.hpp"
#include <cmath>

namespace agoge {
namespace problems {

void SodShockTube::initialize(Field3D &Q)
{
    // For example: half the domain is high-density,
    // half is low-density, etc.
    // Gravity disabled in useGravity().
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double gamma_gas = config::gamma_gas;
    double xMid = 0.5 * Nx * Q.dx; // domain midpoint in x

    for(int k=0; k<Nz; ++k) {
        for(int j=0; j<Ny; ++j) {
            for(int i=0; i<Nx; ++i) {
                int idx = Q.index(i,j,k);

                double x = (i + 0.5) * Q.dx; 
                // For a 1D Sod problem, we might ignore j,k or keep them uniform

                if(x < xMid) {
                    // left state
                    Q.rho [idx] = 1.0;
                    Q.rhou[idx] = 0.0;
                    Q.rhov[idx] = 0.0;
                    Q.rhow[idx] = 0.0;
                    double p0 = 1.0;
                    double E0 = p0 / (gamma_gas - 1.0);
                    Q.E   [idx] = E0;
                } else {
                    // right state
                    Q.rho [idx] = 0.125;
                    Q.rhou[idx] = 0.0;
                    Q.rhov[idx] = 0.0;
                    Q.rhow[idx] = 0.0;
                    double p0 = 0.1;
                    double E0 = p0 / (gamma_gas - 1.0);
                    Q.E   [idx] = E0;
                }
                Q.phi [idx] = 0.0; // no gravity
            }
        }
    }
}

} // namespace problems
} // namespace agoge
