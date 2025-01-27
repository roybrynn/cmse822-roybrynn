#include "EulerSolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

#include "Config.hpp"
#include "Field3d.hpp"
#include "PerformanceMonitor.hpp"
// We include BoundaryManager to get neighbor indices
#include "BoundaryManager.hpp"

namespace agoge {
namespace euler {

//----------------------------------------
// Utility: pressure
static inline double pressure(double rho, double rhou, double rv, double rw,
                              double E) {
    double gamma_gas = config::gamma_gas;
    double u = rhou / rho;
    double v = rv / rho;
    double w = rw / rho;
    double kinetic = 0.5 * rho * (u * u + v * v + w * w);
    double p = (gamma_gas - 1.0) * (E - kinetic);
    return p;
}

//----------------------------------------
// Gravity Acceleration Calculation
static inline void computeGravityAccel(const Field3D &gravField, int i, int j,
                                       int k, double &gx, double &gy,
                                       double &gz) {
    // Compute gravitational acceleration as the negative gradient of phi
    // g = -∇phi
    // Using central differences

    // Handle boundary conditions using BoundaryManager
    int iL = BoundaryManager::getNeighborIndexX(i, gravField.Nx, true);   // i-1
    int iR = BoundaryManager::getNeighborIndexX(i, gravField.Nx, false);  // i+1
    int jD = BoundaryManager::getNeighborIndexY(j, gravField.Ny, true);   // j-1
    int jU = BoundaryManager::getNeighborIndexY(j, gravField.Ny, false);  // j+1
    int kB = BoundaryManager::getNeighborIndexZ(k, gravField.Nz, true);   // k-1
    int kF = BoundaryManager::getNeighborIndexZ(k, gravField.Nz, false);  // k+1

    // Central differences
    double phiL = gravField.phi[gravField.index(iL, j, k)];
    double phiR = gravField.phi[gravField.index(iR, j, k)];
    double phiD = gravField.phi[gravField.index(i, jD, k)];
    double phiU = gravField.phi[gravField.index(i, jU, k)];
    double phiB = gravField.phi[gravField.index(i, j, kB)];
    double phiF = gravField.phi[gravField.index(i, j, kF)];

    // Avoid division by zero by ensuring grid spacing is positive
    double dx = gravField.dx;
    double dy = gravField.dy;
    double dz = gravField.dz;

    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        std::cerr
            << "[computeGravityAccel] Error: Non-positive grid spacing.\n";
        gx = gy = gz = 0.0;
        return;
    }

    gx = -(phiR - phiL) / (2.0 * dx);
    gy = -(phiU - phiD) / (2.0 * dy);
    gz = -(phiF - phiB) / (2.0 * dz);
}

//----------------------------------------
// Minimal artificial viscosity (unchanged)
static void applyArtificialViscosity(const Field3D &Q, Field3D &LQ,
                                     double alpha /* = 0.1 */) {

    agoge::PerformanceMonitor::instance().startTimer("applyArtificialViscosity");

    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    // 1) Compute divergence of velocity in each cell
    std::vector<double> divU(Nx * Ny * Nz, 0.0);
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int c = Q.index(i, j, k);
                double rhoC = Q.rho[c];
                if (rhoC < 1.e-14) {
                    divU[c] = 0.0;
                    continue;
                }

                // local velocity
                double uC = Q.rhou[c] / rhoC;
                double vC = Q.rhov[c] / rhoC;
                double wC = Q.rhow[c] / rhoC;

                // X neighbors
                int iL = BoundaryManager::getNeighborIndexX(i, Nx, true);
                int iR = BoundaryManager::getNeighborIndexX(i, Nx, false);
                double rhoL = Q.rho[Q.index(iL, j, k)];
                double rhoR = Q.rho[Q.index(iR, j, k)];
                double uL =
                    (Q.rhou[Q.index(iL, j, k)] + 1.e-30) / (rhoL + 1.e-30);
                double uR =
                    (Q.rhou[Q.index(iR, j, k)] + 1.e-30) / (rhoR + 1.e-30);
                double dudx = (uR - uL) / (2.0 * dx);

                // Y neighbors
                int jD = BoundaryManager::getNeighborIndexY(j, Ny, true);
                int jU = BoundaryManager::getNeighborIndexY(j, Ny, false);
                double rhoD = Q.rho[Q.index(i, jD, k)];
                double rhoU = Q.rho[Q.index(i, jU, k)];
                double vD =
                    (Q.rhov[Q.index(i, jD, k)] + 1.e-30) / (rhoD + 1.e-30);
                double vU =
                    (Q.rhov[Q.index(i, jU, k)] + 1.e-30) / (rhoU + 1.e-30);
                double dvdy = (vU - vD) / (2.0 * dy);

                // Z neighbors
                int kB = BoundaryManager::getNeighborIndexZ(k, Nz, true);
                int kF = BoundaryManager::getNeighborIndexZ(k, Nz, false);
                double rhoB = Q.rho[Q.index(i, j, kB)];
                double rhoF = Q.rho[Q.index(i, j, kF)];
                double wB =
                    (Q.rhow[Q.index(i, j, kB)] + 1.e-30) / (rhoB + 1.e-30);
                double wF =
                    (Q.rhow[Q.index(i, j, kF)] + 1.e-30) / (rhoF + 1.e-30);
                double dwdz = (wF - wB) / (2.0 * dz);

                divU[c] = dudx + dvdy + dwdz;
            }
        }
    }

    // 2) Compute the Laplacian of each field (rho, rhou, rhov, rhow, E)
    std::vector<double> lapRho(Nx * Ny * Nz, 0.0);
    std::vector<double> lapRhou(Nx * Ny * Nz, 0.0);
    std::vector<double> lapRhov(Nx * Ny * Nz, 0.0);
    std::vector<double> lapRhow(Nx * Ny * Nz, 0.0);
    std::vector<double> lapE(Nx * Ny * Nz, 0.0);

    auto laplacian = [&](const std::vector<double> &field,
                         std::vector<double> &lap) {
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int c = Q.index(i, j, k);
                    double valC = field[c];

                    // X neighbors
                    int iL = BoundaryManager::getNeighborIndexX(i, Nx, true);
                    int iR = BoundaryManager::getNeighborIndexX(i, Nx, false);
                    double valL = field[Q.index(iL, j, k)];
                    double valR = field[Q.index(iR, j, k)];
                    double d2x = (valR - 2.0 * valC + valL) / (dx * dx);

                    // Y neighbors
                    int jD = BoundaryManager::getNeighborIndexY(j, Ny, true);
                    int jU = BoundaryManager::getNeighborIndexY(j, Ny, false);
                    double valD = field[Q.index(i, jD, k)];
                    double valU = field[Q.index(i, jU, k)];
                    double d2y = (valU - 2.0 * valC + valD) / (dy * dy);

                    // Z neighbors
                    int kB = BoundaryManager::getNeighborIndexZ(k, Nz, true);
                    int kF = BoundaryManager::getNeighborIndexZ(k, Nz, false);
                    double valB = field[Q.index(i, j, kB)];
                    double valF = field[Q.index(i, j, kF)];
                    double d2z = (valF - 2.0 * valC + valB) / (dz * dz);

                    lap[c] = d2x + d2y + d2z;
                }
            }
        }
    };

    laplacian(Q.rho, lapRho);
    laplacian(Q.rhou, lapRhou);
    laplacian(Q.rhov, lapRhov);
    laplacian(Q.rhow, lapRhow);
    laplacian(Q.E, lapE);

    // 3) Apply artificial viscosity where divergence is negative
    double h = std::min({dx, dy, dz});

    for (int c = 0; c < Nx * Ny * Nz; c++) {
        double q = (divU[c] < 0.0) ? -divU[c] : 0.0;  // only if negative
        double nu = alpha * h * h * q;

        // Add to LQ
        LQ.rho[c] += nu * lapRho[c];
        LQ.rhou[c] += nu * lapRhou[c];
        LQ.rhov[c] += nu * lapRhov[c];
        LQ.rhow[c] += nu * lapRhow[c];
        LQ.E[c] += nu * lapE[c];
    }
    agoge::PerformanceMonitor::instance().stopTimer("applyArtificialViscosity");
}

//----------------------------------------
// Modified computeL to include gravitational acceleration
// Note: Removed 'static' keyword to match header declaration
void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField) {
    agoge::PerformanceMonitor::instance().startTimer("computeL");
    // Clear LQ
    std::fill(LQ.rho.begin(), LQ.rho.end(), 0.0);
    std::fill(LQ.rhou.begin(), LQ.rhou.end(), 0.0);
    std::fill(LQ.rhov.begin(), LQ.rhov.end(), 0.0);
    std::fill(LQ.rhow.begin(), LQ.rhow.end(), 0.0);
    std::fill(LQ.E.begin(), LQ.E.end(), 0.0);
    std::fill(LQ.phi.begin(), LQ.phi.end(), 0.0);

    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    auto fluxX = [&](double r, double ru, double rv, double rw, double e) {
        double p = pressure(r, ru, rv, rw, e);
        double u = ru / r;
        return std::array<double, 5>{ru, ru * u + p, ru * (rv / r),
                                     ru * (rw / r), (e + p) * u};
    };
    auto fluxY = [&](double r, double ru, double rv, double rw, double e) {
        double p = pressure(r, ru, rv, rw, e);
        double v = rv / r;
        return std::array<double, 5>{rv, rv * (ru / r), rv * v + p,
                                     rv * (rw / r), (e + p) * v};
    };
    auto fluxZ = [&](double r, double ru, double rv, double rw, double e) {
        double p = pressure(r, ru, rv, rw, e);
        double w = rw / r;
        return std::array<double, 5>{rw, rw * (ru / r), rw * (rv / r),
                                     rw * w + p, (e + p) * w};
    };

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int c = Q.index(i, j, k);

                // X neighbors
                int iL = BoundaryManager::getNeighborIndexX(i, Nx, true);
                int iR = BoundaryManager::getNeighborIndexX(i, Nx, false);
                int xm = Q.index(iL, j, k);
                int xp = Q.index(iR, j, k);

                double inv2dx = 1.0 / (2.0 * dx);
                auto FXm = fluxX(Q.rho[xm], Q.rhou[xm], Q.rhov[xm], Q.rhow[xm],
                                 Q.E[xm]);
                auto FXp = fluxX(Q.rho[xp], Q.rhou[xp], Q.rhov[xp], Q.rhow[xp],
                                 Q.E[xp]);
                double dFxdx_rho = (FXp[0] - FXm[0]) * inv2dx;
                double dFxdx_rhou = (FXp[1] - FXm[1]) * inv2dx;
                double dFxdx_rhov = (FXp[2] - FXm[2]) * inv2dx;
                double dFxdx_rhow = (FXp[3] - FXm[3]) * inv2dx;
                double dFxdx_E = (FXp[4] - FXm[4]) * inv2dx;

                // Y neighbors
                int jD = BoundaryManager::getNeighborIndexY(j, Ny, true);
                int jU = BoundaryManager::getNeighborIndexY(j, Ny, false);
                int ym = Q.index(i, jD, k);
                int yp = Q.index(i, jU, k);
                double inv2dy = 1.0 / (2.0 * dy);
                auto FYm = fluxY(Q.rho[ym], Q.rhou[ym], Q.rhov[ym], Q.rhow[ym],
                                 Q.E[ym]);
                auto FYp = fluxY(Q.rho[yp], Q.rhou[yp], Q.rhov[yp], Q.rhow[yp],
                                 Q.E[yp]);
                double dFydy_rho = (FYp[0] - FYm[0]) * inv2dy;
                double dFydy_rhou = (FYp[1] - FYm[1]) * inv2dy;
                double dFydy_rhov = (FYp[2] - FYm[2]) * inv2dy;
                double dFydy_rhow = (FYp[3] - FYm[3]) * inv2dy;
                double dFydy_E = (FYp[4] - FYm[4]) * inv2dy;

                // Z neighbors
                int kB = BoundaryManager::getNeighborIndexZ(k, Nz, true);
                int kF = BoundaryManager::getNeighborIndexZ(k, Nz, false);
                int zm = Q.index(i, j, kB);
                int zp = Q.index(i, j, kF);
                double inv2dz = 1.0 / (2.0 * dz);
                auto FZm = fluxZ(Q.rho[zm], Q.rhou[zm], Q.rhov[zm], Q.rhow[zm],
                                 Q.E[zm]);
                auto FZp = fluxZ(Q.rho[zp], Q.rhou[zp], Q.rhov[zp], Q.rhow[zp],
                                 Q.E[zp]);
                double dFzdz_rho = (FZp[0] - FZm[0]) * inv2dz;
                double dFzdz_rhou = (FZp[1] - FZm[1]) * inv2dz;
                double dFzdz_rhov = (FZp[2] - FZm[2]) * inv2dz;
                double dFzdz_rhow = (FZp[3] - FZm[3]) * inv2dz;
                double dFzdz_E = (FZp[4] - FZm[4]) * inv2dz;

                // Compute RHS without gravity
                double rhs_rho = -(dFxdx_rho + dFydy_rho + dFzdz_rho);
                double rhs_rhou = -(dFxdx_rhou + dFydy_rhou + dFzdz_rhou);
                double rhs_rhov = -(dFxdx_rhov + dFydy_rhov + dFzdz_rhov);
                double rhs_rhow = -(dFxdx_rhow + dFydy_rhow + dFzdz_rhow);
                double rhs_E = -(dFxdx_E + dFydy_E + dFzdz_E);

                // Gravity Acceleration
                double gx = 0.0, gy = 0.0, gz = 0.0;
                if (gravField != nullptr) {
                    computeGravityAccel(*gravField, i, j, k, gx, gy, gz);
                }

                // Compute velocity components
                double u = 0.0, v = 0.0, w = 0.0;
                if (Q.rho[c] > 1.e-14) {  // Avoid division by zero
                    u = Q.rhou[c] / Q.rho[c];
                    v = Q.rhov[c] / Q.rho[c];
                    w = Q.rhow[c] / Q.rho[c];
                }

                // Add gravitational body forces to momentum equations
                rhs_rhou += Q.rho[c] * gx;
                rhs_rhov += Q.rho[c] * gy;
                rhs_rhow += Q.rho[c] * gz;

                // Add gravitational work to energy equation
                rhs_E += Q.rho[c] * (u * gx + v * gy + w * gz);

                // Store in LQ
                LQ.rho[c] = rhs_rho;
                LQ.rhou[c] = rhs_rhou;
                LQ.rhov[c] = rhs_rhov;
                LQ.rhow[c] = rhs_rhow;
                LQ.E[c] = rhs_E;
            }
        }

    }
    // Artificial viscosity
    applyArtificialViscosity(Q, LQ, 1.0);

    agoge::PerformanceMonitor::instance().stopTimer("computeL");
}

//----------------------------------------
// Modified runRK2 to pass the gravitational field
void runRK2(Field3D &Q, double dt) {

    Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQ(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

    // Stage 1
    computeL(Q, LQ, &Q);  // Pass Q as gravField since phi is part of Q
    for (size_t n = 0; n < Q.rho.size(); n++) {
        Qtemp.rho[n] = std::max(Q.rho[n] + dt * LQ.rho[n], 1.0e-10);
        Qtemp.rhou[n] = Q.rhou[n] + dt * LQ.rhou[n];
        Qtemp.rhov[n] = Q.rhov[n] + dt * LQ.rhov[n];
        Qtemp.rhow[n] = Q.rhow[n] + dt * LQ.rhow[n];
        Qtemp.E[n] = std::max(Q.E[n] + dt * LQ.E[n], 1.0e-10);
    }

    // Stage 2
    computeL(Qtemp, LQtemp, &Qtemp);  // Pass Qtemp as gravField
    for (size_t n = 0; n < Q.rho.size(); n++) {
        Q.rho[n] = std::max(0.5 * (Q.rho[n] + Qtemp.rho[n] + dt * LQtemp.rho[n]), 1.0e-10);
        Q.rhou[n] = 0.5 * (Q.rhou[n] + Qtemp.rhou[n] + dt * LQtemp.rhou[n]);
        Q.rhov[n] = 0.5 * (Q.rhov[n] + Qtemp.rhov[n] + dt * LQtemp.rhov[n]);
        Q.rhow[n] = 0.5 * (Q.rhow[n] + Qtemp.rhow[n] + dt * LQtemp.rhow[n]);
        Q.E[n] = std::max(0.5 * (Q.E[n] + Qtemp.E[n] + dt * LQtemp.E[n]), 1.0e-10);
    }

}

//----------------------------------------
// Compute Time Step (unchanged)
double computeTimeStep(const Field3D &Q, double cfl) {
    double maxSpeed = 0.0;
    double gamma_gas = config::gamma_gas;

    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int c = Q.index(i, j, k);
                double r = Q.rho[c];
                if (r <= 1.e-14) continue;

                double ru = Q.rhou[c];
                double rv = Q.rhov[c];
                double rw = Q.rhow[c];
                double e = Q.E[c];

                double u = ru / r;
                double v = rv / r;
                double w = rw / r;
                double speed = std::sqrt(u * u + v * v + w * w);

                double p = pressure(r, ru, rv, rw, e);
                if (p < 0.0) continue;
                double a = std::sqrt(gamma_gas * p / r);

                double localWave = speed + a;
                if (localWave > maxSpeed) {
                    maxSpeed = localWave;
                }
            }
        }
    }

    double minDx = std::min({Q.dx, Q.dy, Q.dz});
    double dt = 1.e20;
    if (maxSpeed > 1.e-14) {
        dt = cfl * (minDx / maxSpeed);
    }

    return dt;
}

}  // namespace euler
}  // namespace agoge
