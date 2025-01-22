// problems/GravityCollapse.cpp
#include "GravityCollapse.hpp"
#include "Config.hpp"
#include <cmath>

namespace agoge {
namespace problems {

void GravityCollapse::initialize(Field3D &Q)
{
    // Example: uniform sphere with a slight overdensity in center
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;
    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    double rho0 = 1.0;
    double p0   = 1.0;
    double E0   = p0 / (config::gamma_gas - 1.0);

    double Lx = Nx * dx;
    double Ly = Ny * dy;
    double Lz = Nz * dz;

    for(int k=0; k<Nz; ++k) {
        for(int j=0; j<Ny; ++j) {
            for(int i=0; i<Nx; ++i) {
                int idx = Q.index(i,j,k);

                double x = (i+0.5)*dx;
                double y = (j+0.5)*dy;
                double z = (k+0.5)*dz;

                double cx = 0.5 * Lx;
                double cy = 0.5 * Ly;
                double cz = 0.5 * Lz;

                double r2 = (x - cx)*(x - cx)
                          + (y - cy)*(y - cy)
                          + (z - cz)*(z - cz);

                double radius2 = 0.01 * Lx*Lx;  // e.g. 0.1 * domain radius
                double local_rho = (r2 < radius2) ? (rho0 * 2.0) : rho0;

                Q.rho [idx] = local_rho;
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;
                Q.E   [idx] = E0 * (local_rho / rho0);
                Q.phi [idx] = 0.0;  // solved via Poisson
            }
        }
    }
}

} // namespace problems
} // namespace agoge
