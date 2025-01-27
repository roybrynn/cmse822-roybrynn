#include "GravityCollapse.hpp"

#include <cmath>
#include <iostream>

#include "../include/agoge/Config.hpp"  // for config::G
#include "../include/agoge/Field3d.hpp"
#include "../include/agoge/GravitySolver.hpp"  // for solvePoisson
#include "../include/agoge/ParameterSystem.hpp"

namespace agoge {
namespace problems {

/**
 * @brief Register problem-specific parameters
 */
void GravityCollapse::registerParameters(ParameterSystem &params) const {
    std::string prefix = name() + ".";  // "GravityCollapse."

    // The total mass in code units (like 1.0 for 1 Msun)
    params.addDefault(prefix + "mass_solar", "1.0");
    // The Jeans radius in code units (like 1.0). We will override domain
    params.addDefault(prefix + "r_jeans", "1.0");
    // The user can supply Nx, Ny, Nz, but we'll override domain=2.4*r_jeans
    // The method used for gravity solver => "naive_dft" or "cooley_tukey"
    // We'll rely on the existing "gravity_method" or similar. If not, we
    // define:
    params.addDefault(prefix + "grav_method", "cooley_tukey");
}

/**
 * @brief Initialize Q for cold dust sphere, solve Poisson, compute error norms.
 */
void GravityCollapse::initialize(Field3D &Q, const ParameterSystem &params) {
    std::string prefix = name() + ".";  // "GravityCollapse."

    // 1) Read user param or defaults
    double M = params.getDouble(prefix + "mass_solar");   // total mass
    double Rjean = params.getDouble(prefix + "r_jeans");  // desired radius
    // We assume config::G = 1 in code units, or user can also specify some G if
    // desired
    double Gval = agoge::config::G;  // typically 1.0
    std::string gravMethod = params.getString(prefix + "grav_method");
    if (gravMethod.empty()) gravMethod = "cooley_tukey";  // default

    // 2) We fix domain => [-1.2R, +1.2R] in each dimension => length = 2.4R
    // The user presumably sets Nx,Ny,Nz in the top-level param (like
    // param["nx"], etc.)
    int Nx = Q.Nx;  // we read the Nx from the existing Q. The user must set Nx
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double sideLength = 2.4 * Rjean;  // for each dimension
    // We'll assume the user *passed* Nx,Ny,Nz but we override Q.dx, Q.dy, Q.dz
    // so that domain is exactly 2.4 R in each dimension.
    Q.dx = sideLength / Nx;
    Q.dy = sideLength / Ny;
    Q.dz = sideLength / Nz;

    // We also might want to shift the domain center to [-1.2R, +1.2R].
    // Let's define an origin shift so that cell centers range from -1.2R to
    // +1.2R
    double x0 = -1.2 * Rjean;  // left edge in x
    double y0 = -1.2 * Rjean;
    double z0 = -1.2 * Rjean;

    std::cout << "[GravityCollapse] Setting domain to "
              << "[" << x0 << "," << (x0 + sideLength) << "] x"
              << "[" << y0 << "," << (y0 + sideLength) << "] x"
              << "[" << z0 << "," << (z0 + sideLength) << "]\n"
              << "with Nx=" << Nx << ",Ny=" << Ny << ",Nz=" << Nz
              << " => dx=" << Q.dx << ",dy=" << Q.dy << ",dz=" << Q.dz << "\n";

    // 3) Fill Q:
    // inside radius => uniform density
    // outside => zero. velocities=0, E=0 => cold dust
    double volumeSphere = (4.0 / 3.0) * M_PI * (Rjean * Rjean * Rjean);
    double rhoInside = M / volumeSphere;  // total mass / volume of sphere

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = Q.index(i, j, k);

                // cell center
                double xC = x0 + (i + 0.5) * Q.dx;
                double yC = y0 + (j + 0.5) * Q.dy;
                double zC = z0 + (k + 0.5) * Q.dz;
                double r2 = xC * xC + yC * yC + zC * zC;
                double R2 = Rjean * Rjean;

                if (r2 <= R2) {
                    // inside sphere
                    Q.rho[idx] = rhoInside;
                } else {
                    // outside sphere
                    Q.rho[idx] = 1.0e-10;
                }
                // velocities=0
                Q.rhou[idx] = 0.0;
                Q.rhov[idx] = 0.0;
                Q.rhow[idx] = 0.0;
                // E=0 => truly cold dust
                Q.E[idx] = 1.0e-10;
                // gravitational potential init => 0
                Q.phi[idx] = 0.0;
            }
        }
    }

    // 4) Solve Poisson eqn => Q.phi.
    // We interpret "gravMethod" => "naive_dft" or "cooley_tukey"
    agoge::gravity::GravityMethod method =
        agoge::gravity::GravityMethod::NAIVE_DFT;
    if (gravMethod == "cooley_tukey") {
        method = agoge::gravity::GravityMethod::COOLEY_TUKEY;
    }
    std::cout << "[GravityCollapse] Using gravity solver method= " << gravMethod
              << "\n";
    agoge::gravity::solvePoisson(Q, method);

    // 5) Compare Q.phi with the analytic uniform-sphere potential
    //    Inside r <= R: phi = - G*M/(2 R^3)*(3R^2 - r^2)
    //    Outside r> R : phi = - G*M / r
    // We'll do an L1, L2 norm over the entire grid.
    double sumAbs = 0.0;
    double sumSqr = 0.0;
    long count = Nx * Ny * Nz;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = Q.index(i, j, k);

                double xC = x0 + (i + 0.5) * Q.dx;
                double yC = y0 + (j + 0.5) * Q.dy;
                double zC = z0 + (k + 0.5) * Q.dz;
                double rr = std::sqrt(xC * xC + yC * yC + zC * zC);

                double phiExact = 0.0;
                if (rr <= Rjean) {
                    // inside
                    double r2 = rr * rr;
                    phiExact = -(Gval * M) / (2.0 * Rjean * Rjean * Rjean) *
                               (3.0 * Rjean * Rjean - r2);
                } else {
                    // outside
                    phiExact = -(Gval * M) / rr;
                }

                double phiNum = Q.phi[idx];
                double diff = phiNum - phiExact;
                sumAbs += std::fabs(diff);
                sumSqr += diff * diff;
            }
        }
    }

    double l1 = sumAbs / double(count);
    double l2 = std::sqrt(sumSqr / double(count));

    std::cout << "[GravityCollapse] Uniform-sphere potential check:\n"
              << "   L1= " << l1 << ", L2= " << l2 << "  (over " << count
              << " cells)\n";
}

}  // namespace problems
}  // namespace agoge