#include "EulerSolver.hpp"
#include "Field3D.hpp"
#include "Config.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace agoge {
namespace euler {

/**
 * @brief Compute pressure from conserved variables.
 *
 * p = (gamma - 1) * [E - 0.5 * rho * (u^2 + v^2 + w^2)]
 */
static inline double pressure(double rho, double rhou, double rhov, double rhow, double E)
{
    double gamma_gas = agoge::config::gamma_gas;
    double u = rhou / rho;
    double v = rhov / rho;
    double w = rhow / rho;
    double kinetic = 0.5 * rho * (u*u + v*v + w*w);
    double p = (gamma_gas - 1.0) * (E - kinetic);
    return p;
}

/**
 * @brief Helper to compute gravitational acceleration from phi, if gravField provided.
 *
 * g = - grad(phi). Using simple central differences with periodic indexing for demonstration.
 */
static inline void computeGravityAccel(const Field3D &gravField,
                                       int i, int j, int k,
                                       double &gx, double &gy, double &gz)
{
    int Nx = gravField.Nx;
    int Ny = gravField.Ny;
    int Nz = gravField.Nz;

    double dx = gravField.dx;
    double dy = gravField.dy;
    double dz = gravField.dz;

    // Periodic neighbors
    int iL = (i == 0)        ? (Nx - 1) : (i - 1);
    int iR = (i == Nx - 1)   ? 0        : (i + 1);
    int jL = (j == 0)        ? (Ny - 1) : (j - 1);
    int jR = (j == Ny - 1)   ? 0        : (j + 1);
    int kL = (k == 0)        ? (Nz - 1) : (k - 1);
    int kR = (k == Nz - 1)   ? 0        : (k + 1);

    int idxC = gravField.index(i,j,k);
    int idxXm = gravField.index(iL,j,k);
    int idxXp = gravField.index(iR,j,k);
    int idxYm = gravField.index(i,jL,k);
    int idxYp = gravField.index(i,jR,k);
    int idxZm = gravField.index(i,j,kL);
    int idxZp = gravField.index(i,j,kR);

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

void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField)
{
    // Zero out LQ fields
    std::fill(LQ.rho .begin(), LQ.rho .end(),  0.0);
    std::fill(LQ.rhou.begin(), LQ.rhou.end(), 0.0);
    std::fill(LQ.rhov.begin(), LQ.rhov.end(), 0.0);
    std::fill(LQ.rhow.begin(), LQ.rhow.end(), 0.0);
    std::fill(LQ.E   .begin(), LQ.E   .end(), 0.0);
    std::fill(LQ.phi .begin(), LQ.phi .end(), 0.0); // Not used in flux, but zero for safety

    // For flux calculations
    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;

    // Loop over each cell
    for(int k = 0; k < Nz; k++) {
        for(int j = 0; j < Ny; j++) {
            for(int i = 0; i < Nx; i++) {

                // Periodic neighbors
                int iL = (i == 0)        ? (Nx - 1) : (i - 1);
                int iR = (i == Nx - 1)   ? 0        : (i + 1);
                int jL = (j == 0)        ? (Ny - 1) : (j - 1);
                int jR = (j == Ny - 1)   ? 0        : (j + 1);
                int kL = (k == 0)        ? (Nz - 1) : (k - 1);
                int kR = (k == Nz - 1)   ? 0        : (k + 1);

                int idxC  = Q.index(i, j, k);
                int idxXm = Q.index(iL, j, k);
                int idxXp = Q.index(iR, j, k);
                int idxYm = Q.index(i, jL, k);
                int idxYp = Q.index(i, jR, k);
                int idxZm = Q.index(i, j, kL);
                int idxZp = Q.index(i, j, kR);

                // Load current cell
                double rhoC  = Q.rho [idxC];
                double rhouC = Q.rhou[idxC];
                double rhovC = Q.rhov[idxC];
                double rhowC = Q.rhow[idxC];
                double EC    = Q.E   [idxC];

                // X-direction neighbors
                double rhoXm  = Q.rho [idxXm];
                double rhoXp  = Q.rho [idxXp];
                double rhouXm = Q.rhou[idxXm];
                double rhouXp = Q.rhou[idxXp];
                double rhovXm = Q.rhov[idxXm];
                double rhovXp = Q.rhov[idxXp];
                double rhowXm = Q.rhow[idxXm];
                double rhowXp = Q.rhow[idxXp];
                double EXm    = Q.E   [idxXm];
                double EXp    = Q.E   [idxXp];

                // Y-direction neighbors
                double rhoYm  = Q.rho [idxYm];
                double rhoYp  = Q.rho [idxYp];
                double rhouYm = Q.rhou[idxYm];
                double rhouYp = Q.rhou[idxYp];
                double rhovYm = Q.rhov[idxYm];
                double rhovYp = Q.rhov[idxYp];
                double rhowYm = Q.rhow[idxYm];
                double rhowYp = Q.rhow[idxYp];
                double EYm    = Q.E   [idxYm];
                double EYp    = Q.E   [idxYp];

                // Z-direction neighbors
                double rhoZm  = Q.rho [idxZm];
                double rhoZp  = Q.rho [idxZp];
                double rhouZm = Q.rhou[idxZm];
                double rhouZp = Q.rhou[idxZp];
                double rhovZm = Q.rhov[idxZm];
                double rhovZp = Q.rhov[idxZp];
                double rhowZm = Q.rhow[idxZm];
                double rhowZp = Q.rhow[idxZp];
                double EZm    = Q.E   [idxZm];
                double EZp    = Q.E   [idxZp];

                // Flux in x-direction
                auto fluxX = [&](double r, double ru, double rv, double rw, double e){
                    double pcell = pressure(r, ru, rv, rw, e);
                    double u = ru / r;
                    return std::array<double,5>{
                        ru,
                        ru*u + pcell,
                        ru*(rv/r),
                        ru*(rw/r),
                        (e + pcell)*u
                    };
                };

                auto FXm = fluxX(rhoXm, rhouXm, rhovXm, rhowXm, EXm);
                auto FXp = fluxX(rhoXp, rhouXp, rhovXp, rhowXp, EXp);

                double inv2dx = 1.0/(2.0*dx);
                double dFxdx_rho  = (FXp[0] - FXm[0]) * inv2dx;
                double dFxdx_rhou = (FXp[1] - FXm[1]) * inv2dx;
                double dFxdx_rhov = (FXp[2] - FXm[2]) * inv2dx;
                double dFxdx_rhow = (FXp[3] - FXm[3]) * inv2dx;
                double dFxdx_E    = (FXp[4] - FXm[4]) * inv2dx;

                // Flux in y-direction
                auto fluxY = [&](double r, double ru, double rv, double rw, double e){
                    double pcell = pressure(r, ru, rv, rw, e);
                    double v = rv / r;
                    return std::array<double,5>{
                        rv,
                        rv*(ru/r),
                        rv*v + pcell,
                        rv*(rw/r),
                        (e + pcell)*v
                    };
                };

                auto FYm = fluxY(rhoYm, rhouYm, rhovYm, rhowYm, EYm);
                auto FYp = fluxY(rhoYp, rhouYp, rhovYp, rhowYp, EYp);

                double inv2dy = 1.0/(2.0*dy);
                double dFydy_rho  = (FYp[0] - FYm[0]) * inv2dy;
                double dFydy_rhou = (FYp[1] - FYm[1]) * inv2dy;
                double dFydy_rhov = (FYp[2] - FYm[2]) * inv2dy;
                double dFydy_rhow = (FYp[3] - FYm[3]) * inv2dy;
                double dFydy_E    = (FYp[4] - FYm[4]) * inv2dy;

                // Flux in z-direction
                auto fluxZ = [&](double r, double ru, double rv, double rw, double e){
                    double pcell = pressure(r, ru, rv, rw, e);
                    double w = rw / r;
                    return std::array<double,5>{
                        rw,
                        rw*(ru/r),
                        rw*(rv/r),
                        rw*w + pcell,
                        (e + pcell)*w
                    };
                };

                auto FZm = fluxZ(rhoZm, rhouZm, rhovZm, rhowZm, EZm);
                auto FZp = fluxZ(rhoZp, rhouZp, rhovZp, rhowZp, EZp);

                double inv2dz = 1.0/(2.0*dz);
                double dFzdz_rho  = (FZp[0] - FZm[0]) * inv2dz;
                double dFzdz_rhou = (FZp[1] - FZm[1]) * inv2dz;
                double dFzdz_rhov = (FZp[2] - FZm[2]) * inv2dz;
                double dFzdz_rhow = (FZp[3] - FZm[3]) * inv2dz;
                double dFzdz_E    = (FZp[4] - FZm[4]) * inv2dz;

                double rhs_rho  = - (dFxdx_rho  + dFydy_rho  + dFzdz_rho );
                double rhs_rhou = - (dFxdx_rhou + dFydy_rhou + dFzdz_rhou);
                double rhs_rhov = - (dFxdx_rhov + dFydy_rhov + dFzdz_rhov);
                double rhs_rhow = - (dFxdx_rhow + dFydy_rhow + dFzdz_rhow);
                double rhs_E    = - (dFxdx_E    + dFydy_E    + dFzdz_E   );

                // Gravity source
                if(gravField) {
                    double gx=0.0, gy=0.0, gz=0.0;
                    computeGravityAccel(*gravField, i, j, k, gx, gy, gz);

                    rhs_rhou += rhoC * gx;
                    rhs_rhov += rhoC * gy;
                    rhs_rhow += rhoC * gz;

                    double u = rhouC / rhoC;
                    double v = rhovC / rhoC;
                    double w = rhowC / rhoC;
                    rhs_E += rhoC * (u*gx + v*gy + w*gz);
                }

                // Store in LQ
                LQ.rho [idxC] = rhs_rho;
                LQ.rhou[idxC] = rhs_rhou;
                LQ.rhov[idxC] = rhs_rhov;
                LQ.rhow[idxC] = rhs_rhow;
                LQ.E   [idxC] = rhs_E;
            }
        }
    }
}

void runRK2(Field3D &Q, double dt)
{
    Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQ   (Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

    // Stage 1
    computeL(Q, LQ, nullptr); // or pass gravity field if you store it in Q
    for(size_t n = 0; n < Q.rho.size(); ++n) {
        Qtemp.rho [n]  = Q.rho [n]  + dt * LQ.rho [n];
        Qtemp.rhou[n]  = Q.rhou[n] + dt * LQ.rhou[n];
        Qtemp.rhov[n]  = Q.rhov[n] + dt * LQ.rhov[n];
        Qtemp.rhow[n]  = Q.rhow[n] + dt * LQ.rhow[n];
        Qtemp.E   [n]  = Q.E   [n]  + dt * LQ.E   [n];
    }

    // Stage 2
    computeL(Qtemp, LQtemp, nullptr); // or pass gravity field if needed
    for(size_t n = 0; n < Q.rho.size(); ++n) {
        Q.rho [n]  = 0.5 * ( Q.rho [n]  + Qtemp.rho [n]  + dt * LQtemp.rho [n] );
        Q.rhou[n]  = 0.5 * ( Q.rhou[n] + Qtemp.rhou[n] + dt * LQtemp.rhou[n] );
        Q.rhov[n]  = 0.5 * ( Q.rhov[n] + Qtemp.rhov[n] + dt * LQtemp.rhov[n] );
        Q.rhow[n]  = 0.5 * ( Q.rhow[n] + Qtemp.rhow[n] + dt * LQtemp.rhow[n] );
        Q.E   [n]  = 0.5 * ( Q.E   [n]  + Qtemp.E   [n]  + dt * LQtemp.E   [n] );
    }
}

double computeTimeStep(const Field3D &Q, double cfl)
{
    // Scan all cells for maximum wave speed
    double maxSpeed = 0.0;
    double gamma_gas = agoge::config::gamma_gas;

    int Nx = Q.Nx;
    int Ny = Q.Ny;
    int Nz = Q.Nz;

    for(int k = 0; k < Nz; ++k) {
        for(int j = 0; j < Ny; ++j) {
            for(int i = 0; i < Nx; ++i) {
                int idx = Q.index(i,j,k);
                double r  = Q.rho [idx];
                double ru = Q.rhou[idx];
                double rv = Q.rhov[idx];
                double rw = Q.rhow[idx];
                double e  = Q.E   [idx];

                if(r <= 0.0) {
                    // Avoid negative density or divide-by-zero.
                    // You could handle error or skip.
                    continue;
                }

                double u = ru / r;
                double v = rv / r;
                double w = rw / r;
                double speed = std::sqrt(u*u + v*v + w*w);

                double p = pressure(r, ru, rv, rw, e);
                if(p < 0.0) {
                    // Negative pressure is unphysical. Handle or skip.
                    continue;
                }

                double a = std::sqrt(gamma_gas * p / r); // Sound speed
                double localWave = speed + a;
                if(localWave > maxSpeed) {
                    maxSpeed = localWave;
                }
            }
        }
    }

    if(maxSpeed < 1.e-14) {
        // If the flow is effectively stationary, choose a small dt anyway
        return 1.e20; // or something large if there's no motion
    }

    // min cell size:
    double minDx = std::min({Q.dx, Q.dy, Q.dz});

    return cfl * (minDx / maxSpeed);
}

} // namespace euler
} // namespace agoge
