#include "EulerSolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

#include "Config.hpp"
#include "Field3d.hpp"
#include "PerformanceMonitor.hpp"

namespace agoge {
namespace euler {

static inline double pressure(double rho, double rhou, double rhov, double rhow,
                              double E) {
    double gamma_gas = agoge::config::gamma_gas;
    double u = rhou / rho;
    double v = rhov / rho;
    double w = rhow / rho;
    double kinetic = 0.5 * rho * (u * u + v * v + w * w);
    double p = (gamma_gas - 1.0) * (E - kinetic);
    return p;
}

/**
 * @brief Compute gravitational acceleration via central differences on phi
 */
static inline void computeGravityAccel(const Field3D &gravField, int i, int j,
                                       int k, double &gx, double &gy,
                                       double &gz) {
    int Nx = gravField.Nx;
    int Ny = gravField.Ny;
    int Nz = gravField.Nz;

    double dx = gravField.dx;
    double dy = gravField.dy;
    double dz = gravField.dz;

    // Periodic neighbors
    auto idx = [&](int ii, int jj, int kk) {
        return gravField.index(ii, jj, kk);
    };
    int iL = (i == 0) ? (Nx - 1) : (i - 1);
    int iR = (i == Nx - 1) ? 0 : (i + 1);
    int jL = (j == 0) ? (Ny - 1) : (j - 1);
    int jR = (j == Ny - 1) ? 0 : (j + 1);
    int kL = (k == 0) ? (Nz - 1) : (k - 1);
    int kR = (k == Nz - 1) ? 0 : (k + 1);

    double phiXp = gravField.phi[idx(iR, j, k)];
    double phiXm = gravField.phi[idx(iL, j, k)];
    double phiYp = gravField.phi[idx(i, jR, k)];
    double phiYm = gravField.phi[idx(i, jL, k)];
    double phiZp = gravField.phi[idx(i, j, kR)];
    double phiZm = gravField.phi[idx(i, j, kL)];

    double dphidx = (phiXp - phiXm) / (2.0 * dx);
    double dphidy = (phiYp - phiYm) / (2.0 * dy);
    double dphidz = (phiZp - phiZm) / (2.0 * dz);

    gx = -dphidx;
    gy = -dphidy;
    gz = -dphidz;
}

/**
 * @brief applyArtificialViscosity:
 *        We compute a small Laplacian-like term for each field, scaled by a
 * factor nu(i,j,k) that triggers near shocks (large negative divergence).
 */
static void applyArtificialViscosity(const Field3D &Q, Field3D &LQ,
                                     double alpha = 0.1) {
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    // local divergence -> compression detection
    // We'll do a naive estimate of \nabla \cdot u in each cell
    std::vector<double> divU(Nx * Ny * Nz, 0.0);

    auto idx = [&](int i, int j, int k) { return Q.index(i, j, k); };

    // compute divergence of velocity for each cell
    for (int k = 0; k < Nz; k++) {
        // periodic or out-of-bound checks...
        int kL = (k == 0) ? (Nz - 1) : (k - 1);
        int kR = (k == Nz - 1) ? 0 : (k + 1);
        for (int j = 0; j < Ny; j++) {
            int jL = (j == 0) ? (Ny - 1) : (j - 1);
            int jR = (j == Ny - 1) ? 0 : (j + 1);
            for (int i = 0; i < Nx; i++) {
                int iL = (i == 0) ? (Nx - 1) : (i - 1);
                int iR = (i == Nx - 1) ? 0 : (i + 1);

                double rhoC = Q.rho[idx(i, j, k)];
                double rhouC = Q.rhou[idx(i, j, k)];
                double rhovC = Q.rhov[idx(i, j, k)];
                double rhowC = Q.rhow[idx(i, j, k)];

                if (rhoC <= 1.e-14) {
                    divU[idx(i, j, k)] = 0.0;
                    continue;
                }
                double uC = rhouC / rhoC;
                double vC = rhovC / rhoC;
                double wC = rhowC / rhoC;

                // neighbors
                double rhoL = Q.rho[idx(iL, j, k)];
                double rhoR = Q.rho[idx(iR, j, k)];
                double rhouL = Q.rhou[idx(iL, j, k)];
                double rhouR = Q.rhou[idx(iR, j, k)];
                double rhovU = Q.rhov[idx(i, jR, k)];
                double rhovD = Q.rhov[idx(i, jL, k)];
                double rhowF = Q.rhow[idx(i, j, kR)];
                double rhowB = Q.rhow[idx(i, j, kL)];

                double uL = (rhouL / (rhoL + 1.e-30));
                double uR = (rhouR / (rhoR + 1.e-30));
                double vU = (rhovU / (Q.rho[idx(i, jR, k)] + 1.e-30));
                double vD = (rhovD / (Q.rho[idx(i, jL, k)] + 1.e-30));
                double wF = (rhowF / (Q.rho[idx(i, j, kR)] + 1.e-30));
                double wB = (rhowB / (Q.rho[idx(i, j, kL)] + 1.e-30));

                double dudx = (uR - uL) / (2.0 * dx);
                double dvdy = (vU - vD) / (2.0 * dy);
                double dwdz = (wF - wB) / (2.0 * dz);

                divU[idx(i, j, k)] = dudx + dvdy + dwdz;
            }
        }
    }

    // Now for each cell, we add \nu * Laplacian(Q) where \nu = alpha * h^2 *
    // max(0, -divU) h ~ local cell size, let's pick h = min(dx, dy, dz)
    double h = std::min({dx, dy, dz});

    // We'll do a naive Laplacian: (Q_{i+1} - 2Q_i + Q_{i-1})/dx^2 + ...
    // apply to each field: rho, rhou, rhov, rhow, E
    // store partial results in temporary arrays
    std::vector<double> lapRho(Nx * Ny * Nz, 0.0);
    std::vector<double> lapRhou(Nx * Ny * Nz, 0.0);
    std::vector<double> lapRhov(Nx * Ny * Nz, 0.0);
    std::vector<double> lapRhow(Nx * Ny * Nz, 0.0);
    std::vector<double> lapE(Nx * Ny * Nz, 0.0);

    auto laplacian = [&](const std::vector<double> &arr,
                         std::vector<double> &lap) {
        // compute in 3D with periodic bc for demonstration
        for (int k = 0; k < Nz; k++) {
            int kL = (k == 0) ? (Nz - 1) : (k - 1);
            int kR = (k == Nz - 1) ? 0 : (k + 1);
            for (int j = 0; j < Ny; j++) {
                int jL = (j == 0) ? (Ny - 1) : (j - 1);
                int jR = (j == Ny - 1) ? 0 : (j + 1);
                for (int i = 0; i < Nx; i++) {
                    int iL = (i == 0) ? (Nx - 1) : (i - 1);
                    int iR = (i == Nx - 1) ? 0 : (i + 1);

                    int idxC = i + Nx * (j + Ny * k);
                    int idxXm = iL + Nx * (j + Ny * k);
                    int idxXp = iR + Nx * (j + Ny * k);
                    int idxYm = i + Nx * (jL + Ny * k);
                    int idxYp = i + Nx * (jR + Ny * k);
                    int idxZm = i + Nx * (j + Ny * kL);
                    int idxZp = i + Nx * (j + Ny * kR);

                    double d2x =
                        (arr[idxXp] - 2.0 * arr[idxC] + arr[idxXm]) / (dx * dx);
                    double d2y =
                        (arr[idxYp] - 2.0 * arr[idxC] + arr[idxYm]) / (dy * dy);
                    double d2z =
                        (arr[idxZp] - 2.0 * arr[idxC] + arr[idxZm]) / (dz * dz);

                    lap[idxC] = d2x + d2y + d2z;
                }
            }
        }
    };

    laplacian(Q.rho, lapRho);
    laplacian(Q.rhou, lapRhou);
    laplacian(Q.rhov, lapRhov);
    laplacian(Q.rhow, lapRhow);
    laplacian(Q.E, lapE);

    // add it to LQ, scaled by \nu(i,j,k)
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idxC = i + Nx * (j + Ny * k);
                // if divergence < 0 => compression => apply artificial
                // viscosity
                double divVal = divU[idxC];
                double q =
                    (divVal < 0.0) ? (-divVal) : 0.0;  // only if negative
                double nu = alpha * h * h * q;         // e.g. alpha ~ 0.1

                LQ.rho[idxC] += nu * lapRho[idxC];
                LQ.rhou[idxC] += nu * lapRhou[idxC];
                LQ.rhov[idxC] += nu * lapRhov[idxC];
                LQ.rhow[idxC] += nu * lapRhow[idxC];
                LQ.E[idxC] += nu * lapE[idxC];
            }
        }
    }
}

/**
 * @brief computeL:
 *        1) second-order central difference for flux derivatives
 *        2) if gravField != nullptr, add gravity
 *        3) apply artificial viscosity for shock handling
 */
void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField) {
    agoge::PerformanceMonitor::instance().startTimer("computeL");

    // 1) Zero out LQ
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

    // local lambdas for fluxes
    auto fluxX = [&](double r, double ru, double rv, double rw, double e) {
        double pcell = pressure(r, ru, rv, rw, e);
        double u = ru / r;
        return std::array<double, 5>{ru, ru * u + pcell, ru * (rv / r),
                                     ru * (rw / r), (e + pcell) * u};
    };
    auto fluxY = [&](double r, double ru, double rv, double rw, double e) {
        double pcell = pressure(r, ru, rv, rw, e);
        double v = rv / r;
        return std::array<double, 5>{rv, rv * (ru / r), rv * v + pcell,
                                     rv * (rw / r), (e + pcell) * v};
    };
    auto fluxZ = [&](double r, double ru, double rv, double rw, double e) {
        double pcell = pressure(r, ru, rv, rw, e);
        double w = rw / r;
        return std::array<double, 5>{rw, rw * (ru / r), rw * (rv / r),
                                     rw * w + pcell, (e + pcell) * w};
    };

    // compute second-order derivative of flux
    auto idx = [&](int i, int j, int k) { return Q.index(i, j, k); };

    for (int k = 0; k < Nz; k++) {
        int kL = (k == 0) ? (Nz - 1) : (k - 1);
        int kR = (k == Nz - 1) ? 0 : (k + 1);
        for (int j = 0; j < Ny; j++) {
            int jL = (j == 0) ? (Ny - 1) : (j - 1);
            int jR = (j == Ny - 1) ? 0 : (j + 1);
            for (int i = 0; i < Nx; i++) {
                int iL = (i == 0) ? (Nx - 1) : (i - 1);
                int iR = (i == Nx - 1) ? 0 : (i + 1);

                int c = idx(i, j, k);
                int xm = idx(iL, j, k);
                int xp = idx(iR, j, k);
                int ym = idx(i, jL, k);
                int yp = idx(i, jR, k);
                int zm = idx(i, j, kL);
                int zp = idx(i, j, kR);

                double rhoC = Q.rho[c];
                double rhouC = Q.rhou[c];
                double rhovC = Q.rhov[c];
                double rhowC = Q.rhow[c];
                double eC = Q.E[c];

                // flux in x-direction
                auto FXm = fluxX(Q.rho[xm], Q.rhou[xm], Q.rhov[xm], Q.rhow[xm],
                                 Q.E[xm]);
                auto FXp = fluxX(Q.rho[xp], Q.rhou[xp], Q.rhov[xp], Q.rhow[xp],
                                 Q.E[xp]);
                double inv2dx = 1.0 / (2.0 * dx);

                double dFxdx_rho = (FXp[0] - FXm[0]) * inv2dx;
                double dFxdx_rhou = (FXp[1] - FXm[1]) * inv2dx;
                double dFxdx_rhov = (FXp[2] - FXm[2]) * inv2dx;
                double dFxdx_rhow = (FXp[3] - FXm[3]) * inv2dx;
                double dFxdx_E = (FXp[4] - FXm[4]) * inv2dx;

                // flux in y-direction
                auto FYm = fluxY(Q.rho[ym], Q.rhou[ym], Q.rhov[ym], Q.rhow[ym],
                                 Q.E[ym]);
                auto FYp = fluxY(Q.rho[yp], Q.rhou[yp], Q.rhov[yp], Q.rhow[yp],
                                 Q.E[yp]);
                double inv2dy = 1.0 / (2.0 * dy);

                double dFydy_rho = (FYp[0] - FYm[0]) * inv2dy;
                double dFydy_rhou = (FYp[1] - FYm[1]) * inv2dy;
                double dFydy_rhov = (FYp[2] - FYm[2]) * inv2dy;
                double dFydy_rhow = (FYp[3] - FYm[3]) * inv2dy;
                double dFydy_E = (FYp[4] - FYm[4]) * inv2dy;

                // flux in z-direction
                auto FZm = fluxZ(Q.rho[zm], Q.rhou[zm], Q.rhov[zm], Q.rhow[zm],
                                 Q.E[zm]);
                auto FZp = fluxZ(Q.rho[zp], Q.rhou[zp], Q.rhov[zp], Q.rhow[zp],
                                 Q.E[zp]);
                double inv2dz = 1.0 / (2.0 * dz);

                double dFzdz_rho = (FZp[0] - FZm[0]) * inv2dz;
                double dFzdz_rhou = (FZp[1] - FZm[1]) * inv2dz;
                double dFzdz_rhov = (FZp[2] - FZm[2]) * inv2dz;
                double dFzdz_rhow = (FZp[3] - FZm[3]) * inv2dz;
                double dFzdz_E = (FZp[4] - FZm[4]) * inv2dz;

                double rhs_rho = -(dFxdx_rho + dFydy_rho + dFzdz_rho);
                double rhs_rhou = -(dFxdx_rhou + dFydy_rhou + dFzdz_rhou);
                double rhs_rhov = -(dFxdx_rhov + dFydy_rhov + dFzdz_rhov);
                double rhs_rhow = -(dFxdx_rhow + dFydy_rhow + dFzdz_rhow);
                double rhs_E = -(dFxdx_E + dFydy_E + dFzdz_E);

                // Gravity
                if (gravField) {
                    double gx = 0.0, gy = 0.0, gz = 0.0;
                    computeGravityAccel(*gravField, i, j, k, gx, gy, gz);

                    rhs_rhou += rhoC * gx;
                    rhs_rhov += rhoC * gy;
                    rhs_rhow += rhoC * gz;

                    // add E source
                    double u = (rhouC + 1.e-30) / (rhoC + 1.e-30);
                    double v = (rhovC + 1.e-30) / (rhoC + 1.e-30);
                    double w = (rhowC + 1.e-30) / (rhoC + 1.e-30);
                    rhs_E += rhoC * (u * gx + v * gy + w * gz);
                }

                // store in LQ
                LQ.rho[c] = rhs_rho;
                LQ.rhou[c] = rhs_rhou;
                LQ.rhov[c] = rhs_rhov;
                LQ.rhow[c] = rhs_rhow;
                LQ.E[c] = rhs_E;
            }
        }
    }

    // 2) artificial viscosity for shock capture
    applyArtificialViscosity(Q, LQ, 0.);
    // alpha=0.2 or so for stronger shock damping

    agoge::PerformanceMonitor::instance().stopTimer("computeL");
}

void runRK2(Field3D &Q, double dt) {
    agoge::PerformanceMonitor::instance().startTimer("runRK2");

    Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQ(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

    // Stage 1
    computeL(Q, LQ, nullptr);
    for (size_t n = 0; n < Q.rho.size(); n++) {
        Qtemp.rho[n] = Q.rho[n] + dt * LQ.rho[n];
        Qtemp.rhou[n] = Q.rhou[n] + dt * LQ.rhou[n];
        Qtemp.rhov[n] = Q.rhov[n] + dt * LQ.rhov[n];
        Qtemp.rhow[n] = Q.rhow[n] + dt * LQ.rhow[n];
        Qtemp.E[n] = Q.E[n] + dt * LQ.E[n];
    }

    // Stage 2
    computeL(Qtemp, LQtemp, nullptr);
    for (size_t n = 0; n < Q.rho.size(); n++) {
        Q.rho[n] = 0.5 * (Q.rho[n] + Qtemp.rho[n] + dt * LQtemp.rho[n]);
        Q.rhou[n] = 0.5 * (Q.rhou[n] + Qtemp.rhou[n] + dt * LQtemp.rhou[n]);
        Q.rhov[n] = 0.5 * (Q.rhov[n] + Qtemp.rhov[n] + dt * LQtemp.rhov[n]);
        Q.rhow[n] = 0.5 * (Q.rhow[n] + Qtemp.rhow[n] + dt * LQtemp.rhow[n]);
        Q.E[n] = 0.5 * (Q.E[n] + Qtemp.E[n] + dt * LQtemp.E[n]);
    }

    agoge::PerformanceMonitor::instance().stopTimer("runRK2");
}

double computeTimeStep(const Field3D &Q, double cfl) {
    agoge::PerformanceMonitor::instance().startTimer("computeTimeStep");
    double maxSpeed = 0.0;
    double gamma_gas = agoge::config::gamma_gas;

    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = Q.index(i, j, k);
                double r = Q.rho[idx];
                if (r <= 1.e-14) continue;

                double ru = Q.rhou[idx];
                double rv = Q.rhov[idx];
                double rw = Q.rhow[idx];
                double e = Q.E[idx];

                double u = ru / r;
                double v = rv / r;
                double w = rw / r;
                double speed = std::sqrt(u * u + v * v + w * w);

                // pressure
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

    agoge::PerformanceMonitor::instance().stopTimer("computeTimeStep");
    return dt;
}

}  // namespace euler
}  // namespace agoge
