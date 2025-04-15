/**
 * @file EulerSolver.cpp
 * @brief Implementation of Euler solver routines for the Agoge application.
 *
 * This file contains functions to compute numerical fluxes, update solution
 * using RK2 time stepping, and compute the time step.
 */

#include "EulerSolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

#include "Config.hpp"
#include "Field3d.hpp"
#include "ParameterSystem.hpp"
#include "PerformanceMonitor.hpp"

namespace agoge {
namespace euler {

//=====================================================
// UTILITY FUNCTIONS
//=====================================================

static inline double pressure(double rho, double rhou, double rhov, double rhow,
                              double E) {
    double gamma = config::gamma_gas;
    double u = rhou / rho, v = rhov / rho, w = rhow / rho;
    double kinetic = 0.5 * rho * (u * u + v * v + w * w);
    double p = (gamma - 1.0) * (E - kinetic);
    return p;
}

static inline double waveSpeed(double rho, double rhou, double rhov,
                               double rhow, double E) {
    double u = rhou / rho, v = rhov / rho, w = rhow / rho;
    double speed = std::sqrt(u * u + v * v + w * w);
    double p = pressure(rho, rhou, rhov, rhow, E);
    double gamma = config::gamma_gas;
    double a = std::sqrt(gamma * p / rho);
    return (speed + a);
}

/// @brief Computes the minmod limited slope
/// @param a First slope
/// @param b Second slope
/// @return Limited slope using minmod limiter
static inline double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

/// @brief Computes the minmod limited slope of three values
/// @param a First slope
/// @param b Second slope
/// @param c Third slope
/// @return Limited slope using minmod limiter
static inline double minmod(double a, double b, double c) {
    return minmod(a, minmod(b, c));
}

static inline std::array<double, 5> fluxX(double r, double ru, double rv,
                                          double rw, double E) {
    double p = pressure(r, ru, rv, rw, E);
    double u = (ru / r);
    return {ru, ru * u + p, ru * (rv / r), ru * (rw / r), (E + p) * u};
}
static inline std::array<double, 5> fluxY(double r, double ru, double rv,
                                          double rw, double E) {
    double p = pressure(r, ru, rv, rw, E);
    double v = (rv / r);
    return {rv, rv * (ru / r), rv * v + p, rv * (rw / r), (E + p) * v};
}
static inline std::array<double, 5> fluxZ(double r, double ru, double rv,
                                          double rw, double E) {
    double p = pressure(r, ru, rv, rw, E);
    double w = (rw / r);
    return {rw, rw * (ru / r), rw * (rv / r), rw * w + p, (E + p) * w};
}

// Move computeFaceFlux here at namespace scope:
template <typename FluxFunction>
inline std::array<double, 5> computeFaceFlux(const std::array<double, 5> &UL,
                                             const std::array<double, 5> &UR,
                                             FluxFunction &&fluxFunc) {
    double aL = waveSpeed(UL[0], UL[1], UL[2], UL[3], UL[4]);
    double aR = waveSpeed(UR[0], UR[1], UR[2], UR[3], UR[4]);
    double alpha = std::max(aL, aR);

    auto fL = fluxFunc(UL[0], UL[1], UL[2], UL[3], UL[4]);
    auto fR = fluxFunc(UR[0], UR[1], UR[2], UR[3], UR[4]);

    std::array<double, 5> faceFlux;
    for (auto v : {var::rho, var::rhou, var::rhov, var::rhow, var::E}) {
        faceFlux[v] = 0.5 * (fL[v] + fR[v]) - 0.5 * alpha * (UR[v] - UL[v]);
    }
    return faceFlux;
}

// Thread-local storage for slopes to avoid repeated allocations
std::vector<double> slopesScratch(tileSize * Nflux, 0.0);

template <typename FluxFunction>
inline void performSweep(std::vector<double> &tileQ, std::vector<double> &tileL,
                         int start1, int stop1, int start2, int stop2,
                         int start3, int stop3, int N1, int N2,
                         double delta_inv, std::vector<double> &fluxN,
                         FluxFunction &&fluxFunc) {
    auto indTile = [](int i1, int i2, int i3, int n1, int n2) {
        return i1 + i2 * n1 + i3 * n1 * n2;
    };
    
    // Compute slopes using minmod limiter
    for (int i3 = start3; i3 < stop3; ++i3) {
        for (int i2 = start2; i2 < stop2; ++i2) {
#pragma omp simd
            for (int i1 = start1 - 1; i1 <= stop1; ++i1) {
                int cM = indTile(i1 - 1, i2, i3, N1, N2);
                int c  = indTile(i1, i2, i3, N1, N2);
                int cP = indTile(i1 + 1, i2, i3, N1, N2);
                
                for (auto v : {var::rho, var::rhou, var::rhov, var::rhow, var::E}) {
                    double vm = tileQ[cM + v * tileSize];
                    double v0 = tileQ[c + v * tileSize];
                    double vp = tileQ[cP + v * tileSize];
                    
                    // Compute differences and limited slope
                    double dL = v0 - vm;
                    double dR = vp - v0;
                    slopesScratch[c + v * tileSize] = minmod(dL, dR);
                }
            }
        }
    }
    
    //============================
    // Compute fluxes at i-1/2 using higher-order reconstruction
    for (int i3 = start3; i3 < stop3; ++i3) {
        for (int i2 = start2; i2 < stop2; ++i2) {
            for (int i1 = start1; i1 <= stop1; ++i1) {
                int cL = indTile(i1 - 1, i2, i3, N1, N2);
                int cR = indTile(i1, i2, i3, N1, N2);

                // Reconstruct left and right states at interface
                std::array<double, 5> UL, UR;
                for (auto v : {var::rho, var::rhou, var::rhov, var::rhow, var::E}) {
                    // Right state of left cell (i-1)
                    UL[v] = tileQ[cL + v * tileSize] + 0.5 * slopesScratch[cL + v * tileSize];
                    
                    // Left state of right cell (i)
                    UR[v] = tileQ[cR + v * tileSize] - 0.5 * slopesScratch[cR + v * tileSize];
                }

                // Compute flux with reconstructed states
                auto flux = computeFaceFlux(UL, UR, fluxFunc);
                for (auto v : {var::rho, var::rhou, var::rhov, var::rhow, var::E}) {
                    fluxN[cR + v * tileSize] = flux[v];
                }
            }
        }
    }
    
    // now accumulate flux differences for interior cells
    for (int i3 = start3; i3 < stop3; ++i3) {
        for (int i2 = start2; i2 < stop2; ++i2) {
#pragma omp simd
            for (int i1 = start1; i1 < stop1; ++i1) {
                int cP = indTile(i1 + 1, i2, i3, N1, N2);
                int cM = indTile(i1, i2, i3, N1, N2);

                for (auto v : {var::rho, var::rhou, var::rhov, var::rhow, var::E}) {
                    auto FP = fluxN[cP + v * tileSize];
                    auto FM = fluxN[cM + v * tileSize];
                    tileL[cM + v * tileSize] -= (FP - FM) * delta_inv;
                }
            }
        }
    }
}

/// @brief Performs a swap of first and second indices (X-Y) of 3D array
/// @param src Array of 5 components (fluid state variables)
/// @param n1 X dimension size
/// @param n2 Y dimension size
/// @param n3 Z dimension size
inline void transpose12(std::vector<double> &src, int n1, int n2, int n3) {
    // Create a complete copy of the source data to work with
    const std::vector<double> copy = src;

    for (int n = 0; n < Nflux; ++n) {
        // Process all planes
        for (int i3 = 0; i3 < n3; ++i3) {
            for (int i2 = 0; i2 < n2; ++i2) {
#pragma omp simd
                for (int i1 = 0; i1 < n1; ++i1) {
                    // Original index (i1, i2, i3)
                    int srcIdx = i1 + i2 * n1 + i3 * (n1 * n2);

                    // Transposed index (i2, i1, i3) - swap dimensions 1
                    // and 2
                    int dstIdx = i2 + i1 * n2 + i3 * (n1 * n2);

                    src[dstIdx + n * tileSize] = copy[srcIdx + n * tileSize];
                }
            }
        }
    }
}

/// @brief Performs a swap of first and third indices (X-Z) of 3D array
/// @param src Array of 5 components (fluid state variables)
/// @param n1 X dimension size
/// @param n2 Y dimension size
/// @param n3 Z dimension size
inline void transpose13(std::vector<double> &src, int n1, int n2, int n3) {
    // Create a complete copy of the source data
    std::vector<double> copy = src;

    // Process each component separately
    for (int n = 0; n < Nflux; ++n) {
        for (int i2 = 0; i2 < n2; ++i2) {
            for (int i1 = 0; i1 < n1; ++i1) {
#pragma clang loop vectorize(enable) interleave(enable)
                for (int i3 = 0; i3 < n3; ++i3) {
                    // Original index (i1, i2, i3)
                    int srcIdx = i1 + i2 * n1 + i3 * (n1 * n2);

                    // Transposed index (i3, i2, i1)
                    int dstIdx = i3 + i2 * n3 + i1 * (n3 * n2);

                    src[dstIdx + n * tileSize] = copy[srcIdx + n * tileSize];
                }
            }
        }
    }
}

//=====================================================
// computeL: 2D or 3D split approach with ghost cells
//=====================================================

/**
 * @brief Compute the residual fluxes over the domain.
 *
 * @param Q The current field.
 * @param LQ The computed residual flux.
 * @param gravField Optional pointer to a gravity field.
 */
void computeL(const Field3D &Q, Field3D &LQ) {
    agoge::PerformanceMonitor::instance().startTimer("computeL");

    // Initialize LQ arrays to zero
    std::fill(LQ.rho.begin(), LQ.rho.end(), 0.0);
    std::fill(LQ.rhou.begin(), LQ.rhou.end(), 0.0);
    std::fill(LQ.rhov.begin(), LQ.rhov.end(), 0.0);
    std::fill(LQ.rhow.begin(), LQ.rhow.end(), 0.0);
    std::fill(LQ.E.begin(), LQ.E.end(), 0.0);

    const int nx = Q.Nx;
    const int ny = Q.Ny;
    const int nz = Q.Nz;
    const double dx = Q.dx, dy = Q.dy, dz = Q.dz;
    const double dx_inv = 1.0 / dx, dy_inv = 1.0 / dy, dz_inv = 1.0 / dz;

    // Setup the tiles for update.
    // Tiles will always be cubic.
    const int nTileX = (Q.Nx + Ntile - 1) / Ntile;
    const int nTileY = (Q.Ny + Ntile - 1) / Ntile;
    const int nTileZ = (Q.Nz + Ntile - 1) / Ntile;

    // std::cout << "Tile sizes: " << nTileX << " " << nTileY << " " << nTileZ
    //           << std::endl;
    // And storage for the tile data and deltas.
    std::vector<double> tileQ(tileSize * Nflux);
    std::vector<double> tileL(tileSize * Nflux);
    std::vector<double> tilePhi(tileSize);

    // Begin processing tiles. We will perform X- Y- Z-sweeps together
    // Indexing here will be global, so we need to adjust for the ghost cells.
    // But we must account for the fact that there is overlap of ghost cells.
    for (int tz = 0; tz < nTileZ; ++tz) {
        const int zStart = tz * Ntile;
        const int zExtent = std::min(NtileGhost, Q.NzGhost - zStart);
        for (int ty = 0; ty < nTileY; ++ty) {
            const int yStart = ty * Ntile;
            const int yExtent = std::min(NtileGhost, Q.NyGhost - yStart);
            for (int tx = 0; tx < nTileX; ++tx) {
                const int xStart = tx * Ntile;
                const int xExtent = std::min(NtileGhost, Q.NxGhost - xStart);

                // std::cout << "Tile: " << tx << " " << ty << " " << tz
                //           << std::endl;
                // std::cout << "Tile extent: " << xExtent << " " << yExtent <<
                // " "
                //           << zExtent << std::endl;
                // std::cout << "Tile start: " << xStart << " " << yStart << " "
                //           << zStart << std::endl;
                // Copy data from Q to tileQ
                // following i,j,k reference the _tile_ frame
                for (int k = 0; k < zExtent; ++k) {
                    for (int j = 0; j < yExtent; ++j) {
                        for (int i = 0; i < xExtent; ++i) {
                            // shift indices to _Field_ global frame
                            int idx =
                                Q.index(i + xStart, j + yStart, k + zStart);
                            int idxTile =
                                i + j * xExtent + k * xExtent * yExtent;
                            tileQ[idxTile + var::rho * tileSize] = Q.rho[idx];
                            tileQ[idxTile + var::rhou * tileSize] = Q.rhou[idx];
                            tileQ[idxTile + var::rhov * tileSize] = Q.rhov[idx];
                            tileQ[idxTile + var::rhow * tileSize] = Q.rhow[idx];
                            tileQ[idxTile + var::E * tileSize] = Q.E[idx];
                            tilePhi[idxTile] =
                                Q.phi[idx];  // Store phi separately
                        }
                    }
                }
                computeLtile(tileQ, tileL, xExtent, yExtent, zExtent, dx_inv,
                             dy_inv, dz_inv, tilePhi);
                // Copy data from tileL to LQ
                for (int k = Nghost; k < zExtent - Nghost; ++k) {
                    for (int j = Nghost; j < yExtent - Nghost; ++j) {
                        for (int i = Nghost; i < xExtent - Nghost; ++i) {
                            int idx =
                                Q.index(i + xStart, j + yStart, k + zStart);
                            int idxTile =
                                i + j * xExtent + k * xExtent * yExtent;
                            LQ.rho[idx] += tileL[idxTile + var::rho * tileSize];
                            LQ.rhou[idx] +=
                                tileL[idxTile + var::rhou * tileSize];
                            LQ.rhov[idx] +=
                                tileL[idxTile + var::rhov * tileSize];
                            LQ.rhow[idx] +=
                                tileL[idxTile + var::rhow * tileSize];
                            LQ.E[idx] += tileL[idxTile + var::E * tileSize];
                        }
                    }
                }
            }
        }
    }
    agoge::PerformanceMonitor::instance().stopTimer("computeL");
}

/// @brief Computes the 3D flux divergence term for the compressible Euler
/// equations in a single tile
/// @param tileQ Vector containing tile state variables (density, momentum,
/// energy)
/// @param tileL Vector containing computed residual fluxes for the tile
/// @param xExtent Extent of tile in x-direction
/// @param yExtent Extent of tile in y-direction
/// @param zExtent Extent of tile in z-direction
/// @param dx_inv Inverse of grid spacing in x-direction
/// @param dy_inv Inverse of grid spacing in y-direction
/// @param dz_inv Inverse of grid spacing in z-direction
/// @param tilePhi Vector containing gravitational potential in the tile
void computeLtile(std::vector<double> &tileQ, std::vector<double> &tileL,
                  int xExtent, int yExtent, int zExtent, double dx_inv,
                  double dy_inv, double dz_inv, std::vector<double> &tilePhi) {
    // clear the tileL
    tileL.assign(tileSize * Nflux, 0.0);

    auto indTile = [](int i1, int i2, int i3, int n1, int n2) {
        return i1 + i2 * n1 + i3 * n1 * n2;
    };

    // Helper for accessing component values in the linearized vectors
    auto getQ = [&tileQ](int n, int idx) -> double & {
        return tileQ[idx + n * tileSize];
    };

    auto getL = [&tileL](int n, int idx) -> double & {
        return tileL[idx + n * tileSize];
    };

    // Start with adding gravity for interior cells before any transposes
    for (int k = Nghost; k < zExtent - Nghost; k++) {
        for (int j = Nghost; j < yExtent - Nghost; j++) {
            for (int i = Nghost; i < xExtent - Nghost; i++) {
                int c = indTile(i, j, k, xExtent, yExtent);
                int cL = indTile(i - 1, j, k, xExtent, yExtent);
                int cR = indTile(i + 1, j, k, xExtent, yExtent);
                int cD = indTile(i, j - 1, k, xExtent, yExtent);
                int cU = indTile(i, j + 1, k, xExtent, yExtent);
                int cB = indTile(i, j, k - 1, xExtent, yExtent);
                int cF = indTile(i, j, k + 1, xExtent, yExtent);

                double gx = -(tilePhi[cR] - tilePhi[cL]) * (0.5 * dx_inv);
                double gy = -(tilePhi[cU] - tilePhi[cD]) * (0.5 * dy_inv);
                double gz = -(tilePhi[cF] - tilePhi[cB]) * (0.5 * dz_inv);

                double rho = getQ(var::rho, c);
                double rhou = getQ(var::rhou, c);
                double rhov = getQ(var::rhov, c);
                double rhow = getQ(var::rhow, c);

                getL(var::rhou, c) = rho * gx;
                getL(var::rhov, c) = rho * gy;
                getL(var::rhow, c) = rho * gz;
                getL(var::E, c) = rhou * gx + rhov * gy + rhow * gz;
            }
        }
    }

    // storage for normal-direction fluxes
    std::vector<double> fluxN(tileSize * Nflux);

    // starting and stopping indices
    int startX = Nghost;
    int stopX = xExtent - Nghost;
    int startY = Nghost;
    int stopY = yExtent - Nghost;
    int startZ = Nghost;
    int stopZ = zExtent - Nghost;

    agoge::PerformanceMonitor::instance().startTimer("computeLtileX");
    //============================
    // Step 1: X-sweep
    // compute fluxes at i-1/2 for iF in [0..nx],
    //============================
    performSweep(tileQ, tileL, startX, stopX, startY, stopY, startZ, stopZ,
                 xExtent, yExtent, dx_inv, fluxN, fluxX);
    agoge::PerformanceMonitor::instance().stopTimer("computeLtileX");

    //===========================
    // Step2: Y-sweep
    // compute fluxes at j-1/2 for jF in [0..ny],
    //===========================
    agoge::PerformanceMonitor::instance().startTimer("transposeY");
    // First, we need to transpose the data; order will be YXZ
    transpose12(tileQ, xExtent, yExtent, zExtent);
    // Now transpose the deltas
    transpose12(tileL, xExtent, yExtent, zExtent);
    agoge::PerformanceMonitor::instance().stopTimer("transposeY");

    agoge::PerformanceMonitor::instance().startTimer("computeLtileY");
    performSweep(tileQ, tileL, startY, stopY, startX, stopX, startZ, stopZ,
                 yExtent, xExtent, dy_inv, fluxN, fluxY);
    agoge::PerformanceMonitor::instance().stopTimer("computeLtileY");

    //===========================
    // Step3: Z-sweep
    // compute fluxes at k-1/2 for kF in [0..nz],
    //===========================

    agoge::PerformanceMonitor::instance().startTimer("transposeZ");
    // First, we need to transpose the data; order will be ZXY
    transpose13(tileQ, yExtent, xExtent, zExtent);
    // transpose the deltas
    transpose13(tileL, yExtent, xExtent, zExtent);
    agoge::PerformanceMonitor::instance().stopTimer("transposeZ");

    agoge::PerformanceMonitor::instance().startTimer("computeLtileZ");
    performSweep(tileQ, tileL, startZ, stopZ, startX, stopX, startY, stopY,
                 zExtent, xExtent, dz_inv, fluxN, fluxZ);
    agoge::PerformanceMonitor::instance().stopTimer("computeLtileZ");

    //===========================
    // Step4: Transpose back to original order
    // transpose the deltas ZXY -> YXZ -> XYZ
    //===========================
    agoge::PerformanceMonitor::instance().startTimer("transposeX");
    transpose13(tileL, zExtent, xExtent, yExtent);
    transpose12(tileL, yExtent, xExtent, zExtent);
    agoge::PerformanceMonitor::instance().stopTimer("transposeX");
}

//=====================================================
// runRK2 remains the same, but we do Q.applyBCs()
// inside computeL each step, so no changes needed here
//=====================================================

/**
 * @brief Advances the solution using a second-order Runge-Kutta scheme.
 *
 * @param Q The field to update.
 * @param dt The time-step.
 */
void runRK2(Field3D &Q, double dt) {
    // create Qtemp, LQ etc. sized with same ghost
    Field3D Qtemp = Q;
    Field3D LQ(Q.Nx, Q.Ny, Q.Nz, Q.bbox);

    // array of references to the fields in Q
    std::array<std::reference_wrapper<std::vector<double>>, 5> fieldsQ = {
        std::ref(Q.rho), std::ref(Q.rhou), std::ref(Q.rhov), std::ref(Q.rhow),
        std::ref(Q.E)};

    // array of references to the fields in Qtemp
    std::array<std::reference_wrapper<std::vector<double>>, 5> fieldsQtemp = {
        std::ref(Qtemp.rho), std::ref(Qtemp.rhou), std::ref(Qtemp.rhov),
        std::ref(Qtemp.rhow), std::ref(Qtemp.E)};

    // array of references to the fields in LQ
    std::array<std::reference_wrapper<std::vector<double>>, 5> fieldsLQ = {
        std::ref(LQ.rho), std::ref(LQ.rhou), std::ref(LQ.rhov),
        std::ref(LQ.rhow), std::ref(LQ.E)};

    // Stage 1
    Q.applyBCs();
    computeL(Q, LQ);  // calls Q.applyBCs() inside
    // do the partial update over the entire domain (ghost included,
    // but real update only needed interior)
    agoge::PerformanceMonitor::instance().startTimer("update1");
    for (int var = 0; var < Nflux; ++var) {
        auto &thisQ = fieldsQ[var].get();
        auto &thisQtemp = fieldsQtemp[var].get();
        auto &thisLQ = fieldsLQ[var].get();

        for (int k = 0; k < Q.Nz; k++) {
            for (int j = 0; j < Q.Ny; ++j) {
                for (int i = 0; i < Q.Nx; ++i) {
                    int n = Q.interiorIndex(i, j, k);
                    thisQtemp[n] =
                        std::max(thisQ[n] + dt * thisLQ[n], floor[var]);
                }
            }
        }
    }
    agoge::PerformanceMonitor::instance().stopTimer("update1");

    // Stage 2
    Qtemp.applyBCs();
    computeL(Qtemp, LQ);
    agoge::PerformanceMonitor::instance().startTimer("update2");
    for (int var = 0; var < Nflux; ++var) {
        auto &thisQ = fieldsQ[var].get();
        auto &thisQtemp = fieldsQtemp[var].get();
        auto &thisLQ = fieldsLQ[var].get();

        for (int k = 0; k < Q.Nz; k++) {
            for (int j = 0; j < Q.Ny; ++j) {
                for (int i = 0; i < Q.Nx; ++i) {
                    int n = Q.interiorIndex(i, j, k);
                    thisQ[n] = std::max(
                        0.5 * (thisQ[n] + thisQtemp[n] + dt * thisLQ[n]),
                        floor[var]);
                }
            }
        }
    }
    agoge::PerformanceMonitor::instance().stopTimer("update2");
}

//=====================================================
// computeTimeStep => only loops interior or entire ghost?
// We'll just do interior to avoid ghost.
// It's simpler to skip boundary cells anyway.
//=====================================================

/**
 * @brief Computes the time step based on the CFL condition.
 *
 * @param Q The current field.
 * @param cfl The CFL number.
 * @return double The computed time step.
 */
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
