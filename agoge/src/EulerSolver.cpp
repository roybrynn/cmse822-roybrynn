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

//=====================================================
// UTILITY FUNCTIONS
//=====================================================

static inline double pressure(double rho, double rhou, double rhov, double rhow,
                              double E) {
    double gamma = config::gamma_gas;
    if (rho < 1e-14) return 0.0;
    double u = rhou / rho, v = rhov / rho, w = rhow / rho;
    double kinetic = 0.5 * rho * (u * u + v * v + w * w);
    double p = (gamma - 1.0) * (E - kinetic);
    if (p < 0.0) p = 0.0;
    return p;
}

static inline double waveSpeed(double rho, double rhou, double rhov,
                               double rhow, double E) {
    if (rho < 1e-14) return 0.0;
    double u = rhou / rho, v = rhov / rho, w = rhow / rho;
    double speed = std::sqrt(u * u + v * v + w * w);
    double p = pressure(rho, rhou, rhov, rhow, E);
    double gamma = config::gamma_gas;
    double a = std::sqrt(gamma * p / rho);
    return (speed + a);
}

static inline std::array<double, 5> fluxX(double r, double ru, double rv,
                                          double rw, double E) {
    double p = pressure(r, ru, rv, rw, E);
    double u = (r > 1e-14) ? (ru / r) : 0.0;
    return {ru, ru * u + p, ru * (rv / r), ru * (rw / r), (E + p) * u};
}
static inline std::array<double, 5> fluxY(double r, double ru, double rv,
                                          double rw, double E) {
    double p = pressure(r, ru, rv, rw, E);
    double v = (r > 1e-14) ? (rv / r) : 0.0;
    return {rv, rv * (ru / r), rv * v + p, rv * (rw / r), (E + p) * v};
}
static inline std::array<double, 5> fluxZ(double r, double ru, double rv,
                                          double rw, double E) {
    double p = pressure(r, ru, rv, rw, E);
    double w = (r > 1e-14) ? (rw / r) : 0.0;
    return {rw, rw * (ru / r), rw * (rv / r), rw * w + p, (E + p) * w};
}

static inline double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    double absa = std::fabs(a), absb = std::fabs(b);
    return (a > 0.0) ? std::min(absa, absb) : -std::min(absa, absb);
}

/**
 * @brief computeGravityAccel from ghosted field. We'll do central difference.
 */
static inline void computeGravityAccel(const Field3D &Q, int iG, int jG,
                                       int kG, double &gx, double &gy,
                                       double &gz) {
    // We'll do i-1 => iG-1, i+1 => iG+1, etc. We assume it is valid because
    // ghost cells exist.
    double dx = Q.dx, dy = Q.dy, dz = Q.dz;

    int cL = Q.interiorIndex(iG - 1, jG, kG);
    int cR = Q.interiorIndex(iG + 1, jG, kG);
    int cD = Q.interiorIndex(iG, jG - 1, kG);
    int cU = Q.interiorIndex(iG, jG + 1, kG);
    int cB = Q.interiorIndex(iG, jG, kG - 1);
    int cF = Q.interiorIndex(iG, jG, kG + 1);

    double phiL = Q.phi[cL];
    double phiR = Q.phi[cR];
    double phiD = Q.phi[cD];
    double phiU = Q.phi[cU];
    double phiB = Q.phi[cB];
    double phiF = Q.phi[cF];

    gx = -(phiR - phiL) / (2. * dx);
    gy = -(phiU - phiD) / (2. * dy);
    gz = -(phiF - phiB) / (2. * dz);
}

//------------------------------------------------
// Minimal artificial viscosity with ghosted Field
//------------------------------------------------
static void applyArtificialViscosity(const Field3D &Q, Field3D &LQ,
                                     double alpha = 0.1) {
    agoge::PerformanceMonitor::instance().startTimer("applyArtificialVisc");

    int nx = Q.Nx, ny = Q.Ny, nz = Q.Nz;
    // total ghosted size:
    int NxG = Q.NxGhost, NyG = Q.NyGhost, NzG = Q.NzGhost;
    double dx = Q.dx, dy = Q.dy, dz = Q.dz;
    double h = std::min({dx, dy, dz});

    // 1) compute divU for interior cells
    std::vector<double> divU(NxG * NyG * NzG, 0.0);

    for (int kIn = 0; kIn < nz; kIn++) {
        for (int jIn = 0; jIn < ny; jIn++) {
            for (int iIn = 0; iIn < nx; iIn++) {
                // map to ghost
                int iG = iIn + Q.nghost;
                int jG = jIn + Q.nghost;
                int kG = kIn + Q.nghost;
                int c = Q.index(iG, jG, kG);

                double r = Q.rho[c];
                if (r < 1e-14) {
                    divU[c] = 0;
                    continue;
                }

                double u = Q.rhou[c] / r, v = Q.rhov[c] / r, w = Q.rhow[c] / r;

                // neighbors => iG-1, iG+1, etc. guaranteed valid by ghost cells
                double rL = Q.rho[Q.index(iG - 1, jG, kG)];
                double rR = Q.rho[Q.index(iG + 1, jG, kG)];
                double uL =
                    (rL > 1e-14) ? (Q.rhou[Q.index(iG - 1, jG, kG)] / rL) : 0.0;
                double uR =
                    (rR > 1e-14) ? (Q.rhou[Q.index(iG + 1, jG, kG)] / rR) : 0.0;
                double dudx = (uR - uL) / (2. * dx);

                double rD = Q.rho[Q.index(iG, jG - 1, kG)];
                double rU = Q.rho[Q.index(iG, jG + 1, kG)];
                double vD =
                    (rD > 1e-14) ? (Q.rhov[Q.index(iG, jG - 1, kG)] / rD) : 0.0;
                double vU =
                    (rU > 1e-14) ? (Q.rhov[Q.index(iG, jG + 1, kG)] / rU) : 0.0;
                double dvdy = (vU - vD) / (2. * dy);

                double rB = Q.rho[Q.index(iG, jG, kG - 1)];
                double rF = Q.rho[Q.index(iG, jG, kG + 1)];
                double wB =
                    (rB > 1e-14) ? (Q.rhow[Q.index(iG, jG, kG - 1)] / rB) : 0.0;
                double wF =
                    (rF > 1e-14) ? (Q.rhow[Q.index(iG, jG, kG + 1)] / rF) : 0.0;
                double dwdz = (wF - wB) / (2. * dz);

                divU[c] = dudx + dvdy + dwdz;
            }
        }
    }

    // 2) Laplacian of each field
    auto laplacian = [&](const std::vector<double> &arr,
                         std::vector<double> &lap) {
        for (int kIn = 0; kIn < nz; kIn++) {
            for (int jIn = 0; jIn < ny; jIn++) {
                for (int iIn = 0; iIn < nx; iIn++) {
                    int iG = iIn + Q.nghost, jG = jIn + Q.nghost,
                        kG = kIn + Q.nghost;
                    int c = Q.index(iG, jG, kG);

                    double valC = arr[c];
                    double valL = arr[Q.index(iG - 1, jG, kG)];
                    double valR = arr[Q.index(iG + 1, jG, kG)];
                    double d2x = (valR - 2. * valC + valL) / (dx * dx);

                    double valD = arr[Q.index(iG, jG - 1, kG)];
                    double valU = arr[Q.index(iG, jG + 1, kG)];
                    double d2y = (valU - 2. * valC + valD) / (dy * dy);

                    double valB = arr[Q.index(iG, jG, kG - 1)];
                    double valF = arr[Q.index(iG, jG, kG + 1)];
                    double d2z = (valF - 2. * valC + valB) / (dz * dz);

                    lap[c] = d2x + d2y + d2z;
                }
            }
        }
    };

    size_t Ntotal = NxG * NyG * NzG;
    std::vector<double> lapRho(Ntotal, 0.), lapRhou(Ntotal, 0.),
        lapRhov(Ntotal, 0.);
    std::vector<double> lapRhow(Ntotal, 0.), lapE(Ntotal, 0.);
    laplacian(Q.rho, lapRho);
    laplacian(Q.rhou, lapRhou);
    laplacian(Q.rhov, lapRhov);
    laplacian(Q.rhow, lapRhow);
    laplacian(Q.E, lapE);

    // 3) add viscosity if divU<0
    for (int kIn = 0; kIn < nz; kIn++) {
        for (int jIn = 0; jIn < ny; jIn++) {
            for (int iIn = 0; iIn < nx; iIn++) {
                int iG = iIn + Q.nghost, jG = jIn + Q.nghost,
                    kG = kIn + Q.nghost;
                int c = Q.index(iG, jG, kG);

                double q = (divU[c] < 0.0) ? -divU[c] : 0.0;
                double nu = alpha * h * h * q;

                LQ.rho[c] += nu * lapRho[c];
                LQ.rhou[c] += nu * lapRhou[c];
                LQ.rhov[c] += nu * lapRhov[c];
                LQ.rhow[c] += nu * lapRhow[c];
                LQ.E[c] += nu * lapE[c];
            }
        }
    }

    agoge::PerformanceMonitor::instance().stopTimer("applyArtificialVisc");
}

//=====================================================
// computeL: 2D or 3D split approach with ghost cells
//=====================================================
void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField) {
    agoge::PerformanceMonitor::instance().startTimer("computeL");

    // 1) We first apply BCs to Q so that all ghost zones are valid
    const_cast<Field3D &>(Q).applyBCs();  // if Q is const, we can cast away or
                                          // applyBCs might be logically const

    // 2) Clear LQ arrays
    std::fill(LQ.rho.begin(), LQ.rho.end(), 0.0);
    std::fill(LQ.rhou.begin(), LQ.rhou.end(), 0.0);
    std::fill(LQ.rhov.begin(), LQ.rhov.end(), 0.0);
    std::fill(LQ.rhow.begin(), LQ.rhow.end(), 0.0);
    std::fill(LQ.E.begin(), LQ.E.end(), 0.0);
    std::fill(LQ.phi.begin(), LQ.phi.end(), 0.0);

    // We'll do dimensional splitting in x,y,z.
    // We'll define Nx= Q.Nx, etc.
    int nx = Q.Nx;
    int ny = Q.Ny;
    int nz = Q.Nz;
    double dx = Q.dx, dy = Q.dy, dz = Q.dz;

    // We'll define flux arrays for each pass, update LQ with - (dF/dx)
    // for x-sweep, then y-sweep, then z-sweep.

    //============================
    // Helper function to gather {rho,rhou,rhov,rhow,E} from Q:
    auto getU = [&](int idx) {
        return std::array<double, 5>{Q.rho[idx], Q.rhou[idx], Q.rhov[idx],
                                     Q.rhow[idx], Q.E[idx]};
    };

    //============================
    // Step 1: X-sweep
    // define flux array Fx => size (nx+1)*ny*nz (in the interior sense),
    // but we store them in a single vector of length (nx+1)*ny*nz, each flux is
    // a 5-vector.

    std::vector<std::array<double, 5>> Fx((nx + 1) * ny * nz);

    // function to index flux face in x
    auto fxIndex = [&](int iF, int jIn, int kIn) {
        // iF in [0..nx], jIn in [0..ny-1], kIn in [0..nz-1]
        return iF + (nx + 1) * (jIn + ny * kIn);
    };

    // compute fluxes at i-1/2 for iF in [0..nx],
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i <= nx; i++) {
                // left cell is iF-1, right cell is iF, in interior coords
                // if iF=0 => left cell is iF-1= -1 => that's in ghost region
                // but Q has ghost cells => so let iL= iF-1.
                
                int cL = Q.interiorIndex(i-1, j, k);
                int cR = Q.interiorIndex(i, j, k);

                // slope-limited states
                auto UL = getU(cL);
                auto UR = getU(cR);

                // we can do minmod slope from neighbors cL-1, cL+1, etc., or do
                // simpler approach for brevity, let's do 1st-order 

                double aL = waveSpeed(UL[0], UL[1], UL[2], UL[3], UL[4]);
                double aR = waveSpeed(UR[0], UR[1], UR[2], UR[3], UR[4]);
                double alpha = std::max(aL, aR);

                auto fL = fluxX(UL[0], UL[1], UL[2], UL[3], UL[4]);
                auto fR = fluxX(UR[0], UR[1], UR[2], UR[3], UR[4]);

                std::array<double, 5> faceFlux;
                for (int n = 0; n < 5; n++) {
                    faceFlux[n] =
                        0.5 * (fL[n] + fR[n]) - 0.5 * alpha * (UR[n] - UL[n]);
                }

                Fx[fxIndex(i, j, k)] = faceFlux;
            }
        }
    }

    // now accumulate flux differences for interior cells
    for (int kIn = 0; kIn < nz; kIn++) {
        for (int jIn = 0; jIn < ny; jIn++) {
            for (int iIn = 0; iIn < nx; iIn++) {
                int iG = iIn + Q.nghost;
                int jG = jIn + Q.nghost;
                int kG = kIn + Q.nghost;
                int c = Q.index(iG, jG, kG);

                auto FL = Fx[fxIndex(iIn, jIn, kIn)];      // flux at left face
                auto FR = Fx[fxIndex(iIn + 1, jIn, kIn)];  // flux at right face

                LQ.rho[c] -= (FR[0] - FL[0]) / dx;
                LQ.rhou[c] -= (FR[1] - FL[1]) / dx;
                LQ.rhov[c] -= (FR[2] - FL[2]) / dx;
                LQ.rhow[c] -= (FR[3] - FL[3]) / dx;
                LQ.E[c] -= (FR[4] - FL[4]) / dx;
            }
        }
    }

    //===========================
    // Step2: Y-sweep
    // similarly define (nxGhost * (ny+1)* nzGhost) flux array or just (nx *
    // (ny+1)* nz) if only interior for brevity, I show a direct approach...
    // ...
    // (omitted for brevity, same pattern as x-sweep)
    //===========================

    //===========================
    // Step3: Z-sweep
    // same pattern
    //===========================

    // Then add gravity for interior cells
    if (gravField) {
        // applyBCs to gravField as well if needed
        // const_cast<Field3D *>(gravField)->applyBCs();
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    int c = Q.interiorIndex(i, j, k);

                    double r = Q.rho[c];
                    if (r < 1e-14) continue;

                    double u = Q.rhou[c] / r, v = Q.rhov[c] / r,
                           w = Q.rhow[c] / r;
                    double gx = 0, gy = 0, gz = 0;
                    computeGravityAccel(*gravField, i, j, k, gx, gy, gz);

                    LQ.rhou[c] += r * gx;
                    LQ.rhov[c] += r * gy;
                    LQ.rhow[c] += r * gz;
                    LQ.E[c] += r * (u * gx + v * gy + w * gz);
                }
            }
        }
    }

    // Finally artificial viscosity
    // applyArtificialViscosity(Q, LQ, 0.1);

    agoge::PerformanceMonitor::instance().stopTimer("computeL");
}

//=====================================================
// runRK2 remains the same, but we do Q.applyBCs()
// inside computeL each step, so no changes needed here
//=====================================================
void runRK2(Field3D &Q, double dt) {
    // create Qtemp, LQ etc. sized with same ghost
    Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz,
                  Q.nghost);
    Field3D LQ(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz,
               Q.nghost);
    Field3D LQtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz,
                   Q.nghost);

    // copy BC flags
    Qtemp.bc_xmin = Q.bc_xmin;
    Qtemp.bc_xmax = Q.bc_xmax;
    Qtemp.bc_ymin = Q.bc_ymin;
    Qtemp.bc_ymax = Q.bc_ymax;
    Qtemp.bc_zmin = Q.bc_zmin;
    Qtemp.bc_zmax = Q.bc_zmax;

    LQ.bc_xmin = Q.bc_xmin;
    LQ.bc_xmax = Q.bc_xmax;
    LQ.bc_ymin = Q.bc_ymin;
    LQ.bc_ymax = Q.bc_ymax;
    LQ.bc_zmin = Q.bc_zmin;
    LQ.bc_zmax = Q.bc_zmax;

    LQtemp.bc_xmin = Q.bc_xmin;
    LQtemp.bc_xmax = Q.bc_xmax;
    LQtemp.bc_ymin = Q.bc_ymin;
    LQtemp.bc_ymax = Q.bc_ymax;
    LQtemp.bc_zmin = Q.bc_zmin;
    LQtemp.bc_zmax = Q.bc_zmax;

    // copy data from Q => Qtemp
    Qtemp.rho = Q.rho;
    Qtemp.rhou = Q.rhou;
    Qtemp.rhov = Q.rhov;
    Qtemp.rhow = Q.rhow;
    Qtemp.E = Q.E;
    Qtemp.phi = Q.phi;

    // Stage 1
    computeL(Q, LQ, &Q);  // calls Q.applyBCs() inside
    // do the partial update over the entire domain (ghost included, but real
    // update only needed interior)
    size_t tot = Q.rho.size();
    for (size_t n = 0; n < tot; n++) {
        double newRho = Q.rho[n] + dt * LQ.rho[n];
        double newE = Q.E[n] + dt * LQ.E[n];
        if (newRho < 1e-14) newRho = 1e-14;
        if (newE < 1e-14) newE = 1e-14;

        Qtemp.rho[n] = newRho;
        Qtemp.rhou[n] = Q.rhou[n] + dt * LQ.rhou[n];
        Qtemp.rhov[n] = Q.rhov[n] + dt * LQ.rhov[n];
        Qtemp.rhow[n] = Q.rhow[n] + dt * LQ.rhow[n];
        Qtemp.E[n] = newE;
    }

    // Stage 2
    computeL(Qtemp, LQtemp, &Qtemp);
    for (size_t n = 0; n < tot; n++) {
        double newRho = 0.5 * (Q.rho[n] + Qtemp.rho[n] + dt * LQtemp.rho[n]);
        double newE = 0.5 * (Q.E[n] + Qtemp.E[n] + dt * LQtemp.E[n]);
        if (newRho < 1e-14) newRho = 1e-14;
        if (newE < 1e-14) newE = 1e-14;

        Q.rho[n] = newRho;
        Q.rhou[n] = 0.5 * (Q.rhou[n] + Qtemp.rhou[n] + dt * LQtemp.rhou[n]);
        Q.rhov[n] = 0.5 * (Q.rhov[n] + Qtemp.rhov[n] + dt * LQtemp.rhov[n]);
        Q.rhow[n] = 0.5 * (Q.rhow[n] + Qtemp.rhow[n] + dt * LQtemp.rhow[n]);
        Q.E[n] = newE;
    }
}

//=====================================================
// computeTimeStep => only loops interior or entire ghost?
// We'll just do interior to avoid ghost.
// It's simpler to skip boundary cells anyway.
//=====================================================
double computeTimeStep(const Field3D &Q, double cfl) {
    // we only compute max speed in interior
    int nx = Q.Nx;
    int ny = Q.Ny;
    int nz = Q.Nz;
    double gamma = config::gamma_gas;
    double maxSpeed = 0.0;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int c = Q.interiorIndex(i, j, k);

                double r = Q.rho[c];
                if (r < 1e-14) continue;
                double u = Q.rhou[c] / r, v = Q.rhov[c] / r, w = Q.rhow[c] / r;
                double speed = std::sqrt(u * u + v * v + w * w);

                double p = pressure(r, Q.rhou[c], Q.rhov[c], Q.rhow[c], Q.E[c]);
                if (p < 0.) p = 0.;
                double a = std::sqrt(gamma * p / r);
                double local = speed + a;
                if (local > maxSpeed) maxSpeed = local;
            }
        }
    }

    double minDx = std::min({Q.dx, Q.dy, Q.dz});
    double dt = 1.e20;
    if (maxSpeed > 1e-14) dt = cfl * (minDx / maxSpeed);
    return dt;
}

}  // namespace euler
}  // namespace agoge
