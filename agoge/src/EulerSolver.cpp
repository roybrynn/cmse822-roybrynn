#include "EulerSolver.hpp"
#include "Field3D.hpp"

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
    double gamma_gas = 1.4; // Or retrieve from config if desired.
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

    double phiC  = gravField.phi[idxC];
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
    std::fill(LQ.phi .begin(), LQ.phi .end(), 0.0); // Not necessarily used, but zero for safety

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

                // Compute flux differences via simple central difference
                // (F_{i+1} - F_{i-1}) / 2 dx, etc.

                // 1. X-direction flux difference
                //    F_x(rho)   = rho * u
                //    F_x(rhou)  = rho u^2 + p
                //    F_x(rhov)  = rho u v
                //    F_x(rhow)  = rho u w
                //    F_x(E)     = (E + p) u
                auto fluxX = [&](double r, double ru, double rv, double rw, double Ecell){
                    double pcell = pressure(r, ru, rv, rw, Ecell);
                    double u     = ru / r;
                    return std::array<double,5>{
                        ru,                 // rho u
                        ru*u + pcell,       // rho u^2 + p
                        ru*(rv/r),          // rho u v
                        ru*(rw/r),          // rho u w
                        (Ecell + pcell)*u
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

                // 2. Y-direction flux difference
                auto fluxY = [&](double r, double ru, double rv, double rw, double Ecell){
                    double pcell = pressure(r, ru, rv, rw, Ecell);
                    double v     = rv / r;
                    return std::array<double,5>{
                        rv,                 // rho v
                        rv*(ru/r),          // rho v u
                        rv*v + pcell,       // rho v^2 + p
                        rv*(rw/r),          // rho v w
                        (Ecell + pcell)*v
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

                // 3. Z-direction flux difference
                auto fluxZ = [&](double r, double ru, double rv, double rw, double Ecell){
                    double pcell = pressure(r, ru, rv, rw, Ecell);
                    double w     = rw / r;
                    return std::array<double,5>{
                        rw,                // rho w
                        rw*(ru/r),         // rho w u
                        rw*(rv/r),         // rho w v
                        rw*w + pcell,      // rho w^2 + p
                        (Ecell + pcell)*w
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

                // Accumulate in LQ
                // L(Q) = - ( dFxdx + dFydy + dFzdz )
                double rhs_rho  = - ( dFxdx_rho  + dFydy_rho  + dFzdz_rho  );
                double rhs_rhou = - ( dFxdx_rhou + dFydy_rhou + dFzdz_rhou );
                double rhs_rhov = - ( dFxdx_rhov + dFydy_rhov + dFzdz_rhov );
                double rhs_rhow = - ( dFxdx_rhow + dFydy_rhow + dFzdz_rhow );
                double rhs_E    = - ( dFxdx_E    + dFydy_E    + dFzdz_E    );

                // Gravity source terms, if gravField != nullptr
                if(gravField) {
                    double gx=0.0, gy=0.0, gz=0.0;
                    computeGravityAccel(*gravField, i, j, k, gx, gy, gz);

                    // Add + rho*g to momentum, + rho*(u dot g) to energy
                    rhs_rhou += rhoC * gx;
                    rhs_rhov += rhoC * gy;
                    rhs_rhow += rhoC * gz;

                    // Add energy source: E' = +rho * (u dot g)
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
    // Temporary field for stage 1
    Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

    // Temporary field for L(Q) computations
    Field3D LQ (Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
    Field3D LQtemp (Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

    // Stage 1: Qtemp = Q + dt * L(Q)
    computeL(Q, LQ, nullptr); // If gravity is needed, pass &Q or separate field
    for(size_t n = 0; n < Q.rho.size(); ++n) {
        Qtemp.rho [n] = Q.rho [n]  + dt * LQ.rho [n];
        Qtemp.rhou[n] = Q.rhou[n] + dt * LQ.rhou[n];
        Qtemp.rhov[n] = Q.rhov[n] + dt * LQ.rhov[n];
        Qtemp.rhow[n] = Q.rhow[n] + dt * LQ.rhow[n];
        Qtemp.E   [n] = Q.E   [n] + dt * LQ.E   [n];
    }

    // Stage 2: Q^(n+1) = 0.5*[Q + Qtemp + dt * L(Qtemp)]
    computeL(Qtemp, LQtemp, nullptr); // If gravity, pass it in
    for(size_t n = 0; n < Q.rho.size(); ++n) {
        Q.rho [n] = 0.5 * ( Q.rho [n]  + Qtemp.rho [n]
                         + dt * LQtemp.rho [n] );
        Q.rhou[n] = 0.5 * ( Q.rhou[n] + Qtemp.rhou[n]
                         + dt * LQtemp.rhou[n] );
        Q.rhov[n] = 0.5 * ( Q.rhov[n] + Qtemp.rhov[n]
                         + dt * LQtemp.rhov[n] );
        Q.rhow[n] = 0.5 * ( Q.rhow[n] + Qtemp.rhow[n]
                         + dt * LQtemp.rhow[n] );
        Q.E   [n] = 0.5 * ( Q.E   [n] + Qtemp.E   [n]
                         + dt * LQtemp.E   [n] );
    }
}

} // namespace euler
} // namespace agoge
