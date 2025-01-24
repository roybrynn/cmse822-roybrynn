#include "EulerSolver.hpp"

#include "Config.hpp"
#include "Field3d.hpp"

// NEW:
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

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

// Add a helper for neighbor indices
inline int neighborIndex(int i, int N, bool left,
                         config::BoundaryCondition bc) {
    if (bc == config::BoundaryCondition::PERIODIC) {
        if (left)
            return (i == 0) ? (N - 1) : (i - 1);
        else
            return (i == N - 1) ? 0 : (i + 1);
    } else if (bc == config::BoundaryCondition::OUTFLOW) {
        if (left)
            return (i == 0) ? 0 : (i - 1);
        else
            return (i == N - 1) ? (N - 1) : (i + 1);
    }
    // ...extend if more BCs...
    return i;
}

/**
 * @brief Compute the gravitational acceleration via central differences on phi
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
    int iL = (i == 0) ? (Nx - 1) : (i - 1);
    int iR = (i == Nx - 1) ? 0 : (i + 1);
    int jL = (j == 0) ? (Ny - 1) : (j - 1);
    int jR = (j == Ny - 1) ? 0 : (j + 1);
    int kL = (k == 0) ? (Nz - 1) : (k - 1);
    int kR = (k == Nz - 1) ? 0 : (k + 1);

    int idxC = gravField.index(i, j, k);
    int idxXm = gravField.index(iL, j, k);
    int idxXp = gravField.index(iR, j, k);
    int idxYm = gravField.index(i, jL, k);
    int idxYp = gravField.index(i, jR, k);
    int idxZm = gravField.index(i, j, kL);
    int idxZp = gravField.index(i, j, kR);

    double phiXm = gravField.phi[idxXm];
    double phiXp = gravField.phi[idxXp];
    double phiYm = gravField.phi[idxYm];
    double phiYp = gravField.phi[idxYp];
    double phiZm = gravField.phi[idxZm];
    double phiZp = gravField.phi[idxZp];

    double dphidx = (phiXp - phiXm) / (2.0 * dx);
    double dphidy = (phiYp - phiYm) / (2.0 * dy);
    double dphidz = (phiZp - phiZm) / (2.0 * dz);

    gx = -dphidx;
    gy = -dphidy;
    gz = -dphidz;
}

/**
 * @brief computeL: second-order central difference for flux derivatives,
 *        plus optional gravity source if gravField is not null.
 */
void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField) {
    // START TIMER for computeL
    agoge::PerformanceMonitor::instance().startTimer("computeL");

    // 1) Clear LQ
    std::fill(LQ.rho.begin(), LQ.rho.end(), 0.0);
    std::fill(LQ.rhou.begin(), LQ.rhou.end(), 0.0);
    std::fill(LQ.rhov.begin(), LQ.rhov.end(), 0.0);
    std::fill(LQ.rhow.begin(), LQ.rhow.end(), 0.0);
    std::fill(LQ.E.begin(), LQ.E.end(), 0.0);
    // LQ.phi is not strictly used for fluxes, so zero for consistency:
    std::fill(LQ.phi.begin(), LQ.phi.end(), 0.0);

    // 2) Setup
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    // 3) Define local lambdas for flux in each direction
    auto fluxX = [&](double r, double ru, double rv, double rw,
                     double e) -> std::array<double, 5> {
        double pcell = pressure(r, ru, rv, rw, e);
        double u = ru / r;
        return {ru,              // rho*u
                ru * u + pcell,  // rho*u^2 + p
                ru * (rv / r),   // rho*u*v
                ru * (rw / r),   // rho*u*w
                (e + pcell) * u};
    };
    auto fluxY = [&](double r, double ru, double rv, double rw,
                     double e) -> std::array<double, 5> {
        double pcell = pressure(r, ru, rv, rw, e);
        double v = rv / r;
        return {rv,              // rho*v
                rv * (ru / r),   // rho*v*u
                rv * v + pcell,  // rho*v^2 + p
                rv * (rw / r),   // rho*v*w
                (e + pcell) * v};
    };
    auto fluxZ = [&](double r, double ru, double rv, double rw,
                     double e) -> std::array<double, 5> {
        double pcell = pressure(r, ru, rv, rw, e);
        double w = rw / r;
        return {rw,              // rho*w
                rw * (ru / r),   // rho*w*u
                rw * (rv / r),   // rho*w*v
                rw * w + pcell,  // rho*w^2 + p
                (e + pcell) * w};
    };

    // 4) Loop over all cells to compute flux derivatives
    //    We'll do second-order central difference in x,y,z with periodic or outlfow BC.
    config::BoundaryCondition bc =
        config::BoundaryCondition::OUTFLOW;  // Or read from somewhere
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                // Set the indexing for neighbors based on the boundary conditions
                int iL = neighborIndex(i, Nx, true, bc);
                int iR = neighborIndex(i, Nx, false, bc);
                int jL = neighborIndex(j, Ny, true, bc);
                int jR = neighborIndex(j, Ny, false, bc);
                int kL = neighborIndex(k, Nz, true, bc);
                int kR = neighborIndex(k, Nz, false, bc);

                int idxC = Q.index(i, j, k);

                // X-direction neighbors
                int idxXm = Q.index(iL, j, k);
                int idxXp = Q.index(iR, j, k);

                // Y-direction neighbors
                int idxYm = Q.index(i, jL, k);
                int idxYp = Q.index(i, jR, k);

                // Z-direction neighbors
                int idxZm = Q.index(i, j, kL);
                int idxZp = Q.index(i, j, kR);

                // Current cell variables
                double rhoC = Q.rho[idxC];
                double rhouC = Q.rhou[idxC];
                double rhovC = Q.rhov[idxC];
                double rhowC = Q.rhow[idxC];
                double eC = Q.E[idxC];

                // X- neighbors
                double rhoXm = Q.rho[idxXm];
                double rhoXp = Q.rho[idxXp];
                double rhouXm = Q.rhou[idxXm];
                double rhouXp = Q.rhou[idxXp];
                double rhovXm = Q.rhov[idxXm];
                double rhovXp = Q.rhov[idxXp];
                double rhowXm = Q.rhow[idxXm];
                double rhowXp = Q.rhow[idxXp];
                double eXm = Q.E[idxXm];
                double eXp = Q.E[idxXp];

                // Y- neighbors
                double rhoYm = Q.rho[idxYm];
                double rhoYp = Q.rho[idxYp];
                double rhouYm = Q.rhou[idxYm];
                double rhouYp = Q.rhou[idxYp];
                double rhovYm = Q.rhov[idxYm];
                double rhovYp = Q.rhov[idxYp];
                double rhowYm = Q.rhow[idxYm];
                double rhowYp = Q.rhow[idxYp];
                double eYm = Q.E[idxYm];
                double eYp = Q.E[idxYp];

                // Z- neighbors
                double rhoZm = Q.rho[idxZm];
                double rhoZp = Q.rho[idxZp];
                double rhouZm = Q.rhou[idxZm];
                double rhouZp = Q.rhou[idxZp];
                double rhovZm = Q.rhov[idxZm];
                double rhovZp = Q.rhov[idxZp];
                double rhowZm = Q.rhow[idxZm];
                double rhowZp = Q.rhow[idxZp];
                double eZm = Q.E[idxZm];
                double eZp = Q.E[idxZp];

                // 4a) Compute fluxes at i-1 and i+1 in X direction
                auto FXm = fluxX(rhoXm, rhouXm, rhovXm, rhowXm, eXm);
                auto FXp = fluxX(rhoXp, rhouXp, rhovXp, rhowXp, eXp);

                double inv2dx = 1.0 / (2.0 * dx);

                // dFxdx for each conserved variable
                double dFxdx_rho = (FXp[0] - FXm[0]) * inv2dx;
                double dFxdx_rhou = (FXp[1] - FXm[1]) * inv2dx;
                double dFxdx_rhov = (FXp[2] - FXm[2]) * inv2dx;
                double dFxdx_rhow = (FXp[3] - FXm[3]) * inv2dx;
                double dFxdx_E = (FXp[4] - FXm[4]) * inv2dx;

                // 4b) Y direction
                auto FYm = fluxY(rhoYm, rhouYm, rhovYm, rhowYm, eYm);
                auto FYp = fluxY(rhoYp, rhouYp, rhovYp, rhowYp, eYp);

                double inv2dy = 1.0 / (2.0 * dy);

                double dFydy_rho = (FYp[0] - FYm[0]) * inv2dy;
                double dFydy_rhou = (FYp[1] - FYm[1]) * inv2dy;
                double dFydy_rhov = (FYp[2] - FYm[2]) * inv2dy;
                double dFydy_rhow = (FYp[3] - FYm[3]) * inv2dy;
                double dFydy_E = (FYp[4] - FYm[4]) * inv2dy;

                // 4c) Z direction
                auto FZm = fluxZ(rhoZm, rhouZm, rhovZm, rhowZm, eZm);
                auto FZp = fluxZ(rhoZp, rhouZp, rhovZp, rhowZp, eZp);

                double inv2dz = 1.0 / (2.0 * dz);

                double dFzdz_rho = (FZp[0] - FZm[0]) * inv2dz;
                double dFzdz_rhou = (FZp[1] - FZm[1]) * inv2dz;
                double dFzdz_rhov = (FZp[2] - FZm[2]) * inv2dz;
                double dFzdz_rhow = (FZp[3] - FZm[3]) * inv2dz;
                double dFzdz_E = (FZp[4] - FZm[4]) * inv2dz;

                // 5) Combine for L(Q) = - ( dF/dx + dF/dy + dF/dz )
                double rhs_rho = -(dFxdx_rho + dFydy_rho + dFzdz_rho);
                double rhs_rhou = -(dFxdx_rhou + dFydy_rhou + dFzdz_rhou);
                double rhs_rhov = -(dFxdx_rhov + dFydy_rhov + dFzdz_rhov);
                double rhs_rhow = -(dFxdx_rhow + dFydy_rhow + dFzdz_rhow);
                double rhs_E = -(dFxdx_E + dFydy_E + dFzdz_E);

                // 6) Gravity source terms if gravField != nullptr
                if (gravField) {
                    double gx = 0.0, gy = 0.0, gz = 0.0;
                    computeGravityAccel(*gravField, i, j, k, gx, gy, gz);

                    // + (rho*g) on momentum
                    rhs_rhou += rhoC * gx;
                    rhs_rhov += rhoC * gy;
                    rhs_rhow += rhoC * gz;

                    // + (rho * (u dot g)) on energy
                    double u = rhouC / rhoC;
                    double v = rhovC / rhoC;
                    double w = rhowC / rhoC;
                    rhs_E += rhoC * (u * gx + v * gy + w * gz);
                }

                // 7) Store results in LQ
                LQ.rho[idxC] = rhs_rho;
                LQ.rhou[idxC] = rhs_rhou;
                LQ.rhov[idxC] = rhs_rhov;
                LQ.rhow[idxC] = rhs_rhow;
                LQ.E[idxC] = rhs_E;
            }
        }
    }

    // STOP TIMER for computeL
    agoge::PerformanceMonitor::instance().stopTimer("computeL");
}

/**
 * @brief Runge-Kutta 2 steps remain unchanged from the partial code
 */
void runRK2(Field3D &Q, double dt) {
    agoge::PerformanceMonitor::instance().startTimer("runRK2");

    Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQ(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

    // Stage 1
    computeL(Q, LQ, nullptr);  // If gravity stored in Q, pass &Q
    for (size_t n = 0; n < Q.rho.size(); ++n) {
        Qtemp.rho[n] = Q.rho[n] + dt * LQ.rho[n];
        Qtemp.rhou[n] = Q.rhou[n] + dt * LQ.rhou[n];
        Qtemp.rhov[n] = Q.rhov[n] + dt * LQ.rhov[n];
        Qtemp.rhow[n] = Q.rhow[n] + dt * LQ.rhow[n];
        Qtemp.E[n] = Q.E[n] + dt * LQ.E[n];
    }

    // Stage 2
    computeL(Qtemp, LQtemp, nullptr);  // If gravity stored, pass &Qtemp
    for (size_t n = 0; n < Q.rho.size(); ++n) {
        Q.rho[n] = 0.5 * (Q.rho[n] + Qtemp.rho[n] + dt * LQtemp.rho[n]);
        Q.rhou[n] = 0.5 * (Q.rhou[n] + Qtemp.rhou[n] + dt * LQtemp.rhou[n]);
        Q.rhov[n] = 0.5 * (Q.rhov[n] + Qtemp.rhov[n] + dt * LQtemp.rhov[n]);
        Q.rhow[n] = 0.5 * (Q.rhow[n] + Qtemp.rhow[n] + dt * LQtemp.rhow[n]);
        Q.E[n] = 0.5 * (Q.E[n] + Qtemp.E[n] + dt * LQtemp.E[n]);
    }

    agoge::PerformanceMonitor::instance().stopTimer("runRK2");
}

/**
 * @brief computeTimeStep is unchanged from the partial code snippet
 */
double computeTimeStep(const Field3D &Q, double cfl) {
    agoge::PerformanceMonitor::instance().startTimer("computeTimeStep");
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
                double e = Q.E[idx];

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

    agoge::PerformanceMonitor::instance().stopTimer("computeTimeStep");
    return dt;
}

}  // namespace euler
}  // namespace agoge
