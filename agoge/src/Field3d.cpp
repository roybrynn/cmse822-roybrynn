/// @file Field3d.cpp
/// @brief Implementation of the Field3D class for the Agoge application.
///
/// This file provides method definitions for spatial coordinate calculations,
/// boundary conditions application, and MPI halo exchange for 3D fields.

// Add preprocessor definition to switch between blocking and non-blocking MPI
#ifndef AGOGE_USE_NONBLOCKING_MPI
#define AGOGE_USE_NONBLOCKING_MPI 1  // Default to blocking MPI (0=blocking, 1=non-blocking)
#endif

#include "Field3d.hpp"
#include "ParameterSystem.hpp"

#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <array>

namespace agoge {

Field3D::Field3D(int nx, int ny, int nz, const BoundingBox &bbox_in, int ghost)
    : Nx(nx),
      Ny(ny),
      Nz(nz),
      nghost(ghost),
      dx((bbox_in.xmax - bbox_in.xmin) / nx),
      dy((bbox_in.ymax - bbox_in.ymin) / ny),
      dz((bbox_in.zmax - bbox_in.zmin) / nz),
      bbox(bbox_in),  // Initialize BoundingBox member
      bc_xmin(config::BoundaryCondition::PERIODIC),
      bc_xmax(config::BoundaryCondition::PERIODIC),
      bc_ymin(config::BoundaryCondition::PERIODIC),
      bc_ymax(config::BoundaryCondition::PERIODIC),
      bc_zmin(config::BoundaryCondition::PERIODIC),
      bc_zmax(config::BoundaryCondition::PERIODIC),
      reqCountX(0),
      reqCountY(0),
      reqCountZ(0) {
    NxGhost = Nx + 2 * nghost;
    NyGhost = Ny + 2 * nghost;
    NzGhost = Nz + 2 * nghost;

    size_t totalSize = NxGhost * NyGhost * NzGhost;
    rho.resize(totalSize, 0.0);
    rhou.resize(totalSize, 0.0);
    rhov.resize(totalSize, 0.0);
    rhow.resize(totalSize, 0.0);
    E.resize(totalSize, 0.0);
    phi.resize(totalSize, 0.0);
    
    // Initialize MPI request arrays to MPI_REQUEST_NULL
    std::fill(requestsX.begin(), requestsX.end(), MPI_REQUEST_NULL);
    std::fill(requestsY.begin(), requestsY.end(), MPI_REQUEST_NULL);
    std::fill(requestsZ.begin(), requestsZ.end(), MPI_REQUEST_NULL);
}

//======================================================
// Spatial coordinate methods
//======================================================
// The left edge of the interior domain is at the start of cell i=0 =>
// if you want a global offset X0, define it as a member.
// For now, assume X0=0.
// Cell iIn's center => (iIn+0.5)*dx.
// We do NOT add nghost here because that's an array index offset only.
// If the physical domain truly starts at X0, then:
//   x= X0 + (iIn+0.5)*dx
// That ensures iIn in [0..Nx-1] => x in [0.5*dx.. Nx-0.5]*dx.

double Field3D::xCenter(int iIn) const { return bbox.xmin + (iIn + 0.5) * dx; }
double Field3D::xLeftEdge(int iIn) const { return bbox.xmin + iIn * dx; }
double Field3D::xRightEdge(int iIn) const { return bbox.xmin + (iIn + 1) * dx; }

double Field3D::yCenter(int jIn) const { return bbox.ymin + (jIn + 0.5) * dy; }
double Field3D::yLeftEdge(int jIn) const { return bbox.ymin + jIn * dy; }
double Field3D::yRightEdge(int jIn) const { return bbox.ymin + (jIn + 1) * dy; }

double Field3D::zCenter(int kIn) const { return bbox.zmin + (kIn + 0.5) * dz; }
double Field3D::zLeftEdge(int kIn) const { return bbox.zmin + kIn * dz; }
double Field3D::zRightEdge(int kIn) const { return bbox.zmin + (kIn + 1) * dz; }

/// @brief Determine MPI tags for sending and receiving data between neighboring
/// ranks.
///
/// This function helps to avoid tag collisions by assigning unique tags based
/// on the ranks of the current process (myRank) and its neighbor
/// (neighborRank). The base value is used to differentiate between different
/// types of communications (e.g., x, y, z directions).
///
/// @param myRank The rank of the current process.
/// @param neighborRank The rank of the neighboring process.
/// @param base The base value for the tag to differentiate communication types.
/// @return A pair of integers representing the send and receive tags.
// inline std::pair<int, int> getTags(int myRank, int neighborRank, int base) {
//     // ...obsolete code, remove...
// }

// Simple symmetric tag computation:
//   base = direction_offset + phase_offset + (min(myRank,neighborRank)*10 +
//   max(myRank,neighborRank))
// where for example, for 'Z': direction_offset = 20000, phase0 -> phase_offset
// = 0, phase1 -> phase_offset = 100. Then if (myRank < neighborRank) return
// {base+1, base+2}; else return {base+2, base+1};
inline std::pair<int, int> computeSymmetricTagPair(int myRank, int neighborRank,
                                                   char direction, int phase) {
    int directionOffset = 0;
    if (direction == 'X') {
        directionOffset = 0;
    } else if (direction == 'Y') {
        directionOffset = 10000;
    } else if (direction == 'Z') {
        directionOffset = 20000;
    }
    int phaseOffset = (phase == 0) ? 0 : 100;
    int base =
        directionOffset + phaseOffset +
        (std::min(myRank, neighborRank) * 10 + std::max(myRank, neighborRank));
    if (myRank < neighborRank) {
        return {base + 1, base + 2};  // lower rank sends with (base+1), higher
                                      // receives with (base+2)
    } else if (myRank > neighborRank) {
        return {base + 2, base + 1};  // higher rank sends with (base+2), lower
                                      // receives with (base+1)
    } else {
        return {base + 1, base + 1};  // self-communication (unlikely)
    }
}

void Field3D::applyBCs() {
    // Helper lambdas for periodic boundary conditions - simplified for clarity
    auto periodicCoordMinus = [this](int ghost, int interiorSize) -> int {
        // For ghost cell on the "minus" side, map to the opposite end of the domain
        // Ghost 0 -> Last interior cell, Ghost 1 -> Second to last interior cell, etc.
        return nghost + interiorSize - 1 - ghost;
    };

    auto periodicCoordPlus = [this](int physicalGhost, int interiorSize) -> int {
        // For ghost cell on the "plus" side, map to the beginning of the interior domain
        // We take the modulo to handle the wrapping directly
        return nghost + (physicalGhost % interiorSize);
    };

    auto fillBC_x = [&](std::vector<double> &arr) {
        // x-minus physical BC: iGhost in [0, nghost)
        for (int k = 0; k < NzGhost; k++) {
            for (int j = 0; j < NyGhost; j++) {
                for (int iGhost = 0; iGhost < nghost; iGhost++) {
                    int ghostIdx = index(iGhost, j, k);
                    if (bc_xmin == config::BoundaryCondition::PERIODIC) {
                        int iMirror = periodicCoordMinus(iGhost, Nx);
                        int src = index(iMirror, j, k);
                        arr[ghostIdx] = arr[src];
                    } else if (bc_xmin == config::BoundaryCondition::OUTFLOW) {
                        int src = index(nghost, j, k);
                        arr[ghostIdx] = arr[src];
                    }
                }
            }
        }

        // x-plus physical BC: iGhost in [nghost+Nx, NxGhost)
        for (int k = 0; k < NzGhost; k++) {
            for (int j = 0; j < NyGhost; j++) {
                for (int iGhost = nghost + Nx; iGhost < NxGhost; iGhost++) {
                    int ghostIdx = index(iGhost, j, k);
                    if (bc_xmax == config::BoundaryCondition::PERIODIC) {
                        // Compute physical ghost index and then find mirror coordinate
                        int physicalGhost = iGhost - nghost;
                        int iMirror = periodicCoordPlus(physicalGhost, Nx);
                        int src = index(iMirror, j, k);
                        arr[ghostIdx] = arr[src];
                    } else if (bc_xmax == config::BoundaryCondition::OUTFLOW) {
                        int src = index(NxGhost - nghost - 1, j, k);
                        arr[ghostIdx] = arr[src];
                    }
                }
            }
        }
    };

    auto fillBC_y = [&](std::vector<double> &arr) {
        // y-min physical BC: jGhost in [0, nghost).
        for (int k = 0; k < NzGhost; k++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int jGhost = 0; jGhost < nghost; jGhost++) {
                    int ghostIdx = index(i, jGhost, k);
                    if (bc_ymin == config::BoundaryCondition::PERIODIC) {
                        int jMirror = periodicCoordMinus(jGhost, Ny);
                        int src = index(i, jMirror, k);
                        arr[ghostIdx] = arr[src];
                    } else if (bc_ymin == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, nghost, k);
                        arr[ghostIdx] = arr[src];
                    }
                }
            }
        }
        
        // y-max physical BC: jGhost in [nghost+Ny, NyGhost).
        for (int k = 0; k < NzGhost; k++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int jGhost = nghost + Ny; jGhost < NyGhost; jGhost++) {
                    int ghostIdx = index(i, jGhost, k);
                    if (bc_ymax == config::BoundaryCondition::PERIODIC) {
                        int physicalGhost = jGhost - nghost;
                        int jMirror = periodicCoordPlus(physicalGhost, Ny);
                        int src = index(i, jMirror, k);
                        arr[ghostIdx] = arr[src];
                    } else if (bc_ymax == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, Ny + nghost - 1, k);
                        arr[ghostIdx] = arr[src];
                    }
                }
            }
        }
    };

    auto fillBC_z = [&](std::vector<double> &arr) {
        // z-min physical BC: kGhost in [0, nghost).
        for (int j = 0; j < NyGhost; j++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int kGhost = 0; kGhost < nghost; kGhost++) {
                    int ghostIdx = index(i, j, kGhost);
                    if (bc_zmin == config::BoundaryCondition::PERIODIC) {
                        int kMirror = periodicCoordMinus(kGhost, Nz);
                        int src = index(i, j, kMirror);
                        arr[ghostIdx] = arr[src];
                    } else if (bc_zmin == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, j, nghost);
                        arr[ghostIdx] = arr[src];
                    }
                }
            }
        }
        
        // z-max physical BC: kGhost in [nghost+Nz, NzGhost).
        for (int j = 0; j < NyGhost; j++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int kGhost = nghost + Nz; kGhost < NzGhost; kGhost++) {
                    int ghostIdx = index(i, j, kGhost);
                    if (bc_zmax == config::BoundaryCondition::PERIODIC) {
                        int physicalGhost = kGhost - nghost;
                        int kMirror = periodicCoordPlus(physicalGhost, Nz);
                        int src = index(i, j, kMirror);
                        arr[ghostIdx] = arr[src];
                    } else if (bc_zmax == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, j, Nz + nghost - 1);
                        arr[ghostIdx] = arr[src];
                    }
                }
            }
        }
    };

    fillBC_x(rho);
    fillBC_x(rhou);
    fillBC_x(rhov);
    fillBC_x(rhow);
    fillBC_x(E);
    fillBC_x(phi);
    fillBC_y(rho);
    fillBC_y(rhou);
    fillBC_y(rhov);
    fillBC_y(rhow);
    fillBC_y(E);
    fillBC_y(phi);
    fillBC_z(rho);
    fillBC_z(rhou);
    fillBC_z(rhov);
    fillBC_z(rhow);
    fillBC_z(E);
    fillBC_z(phi);

    // --- Begin MPI halo exchange ---
    // The following applies if there is a valid MPI neighbor,
    // while for global boundaries (neighbor rank == -1) physical BCs remain.

    // Unified dimension specifier
    enum class Dimension { X, Y, Z };

#if AGOGE_USE_NONBLOCKING_MPI
    // --- Non-blocking MPI communication implementation ---
    const int numArrays = 6;  // rho, rhou, rhov, rhow, E, phi
    
    // X-direction exchange using non-blocking MPI
    if (rankMinusX != MPI_PROC_NULL || rankPlusX != MPI_PROC_NULL) {
        int slabSize = NyGhost * NzGhost * nghost;
        int totalCount = slabSize * numArrays;
        
        std::fill(requestsX.begin(), requestsX.end(), MPI_REQUEST_NULL);
        reqCountX = 0;

        // Post receives first
        if (rankMinusX != MPI_PROC_NULL) {
            MPI_Irecv(recvBufXMinus.data(), totalCount, MPI_DOUBLE,
                      rankMinusX, 0, MPI_COMM_WORLD,
                      &requestsX[reqCountX++]);
        }
        if (rankPlusX != MPI_PROC_NULL) {
            MPI_Irecv(recvBufXPlus.data(), totalCount, MPI_DOUBLE,
                      rankPlusX, 0, MPI_COMM_WORLD,
                      &requestsX[reqCountX++]);
        }

        // Pack & send
        if (rankPlusX != MPI_PROC_NULL) {
            packGhosts(sendBufXPlus, 0, true);
            MPI_Isend(sendBufXPlus.data(), totalCount, MPI_DOUBLE,
                      rankPlusX, 0, MPI_COMM_WORLD,
                      &requestsX[reqCountX++]);
        }
        if (rankMinusX != MPI_PROC_NULL) {
            packGhosts(sendBufXMinus, 0, false);
            MPI_Isend(sendBufXMinus.data(), totalCount, MPI_DOUBLE,
                      rankMinusX, 0, MPI_COMM_WORLD,
                      &requestsX[reqCountX++]);
        }

        // Wait and unpack
        MPI_Waitall(reqCountX, requestsX.data(), MPI_STATUSES_IGNORE);
        if (rankMinusX != MPI_PROC_NULL) {
            unpackGhosts(recvBufXMinus, 0, false);
        }
        if (rankPlusX != MPI_PROC_NULL) {
            unpackGhosts(recvBufXPlus, 0, true);
        }
    }
    
    // Y-direction exchange using non-blocking MPI
    if (rankMinusY != MPI_PROC_NULL || rankPlusY != MPI_PROC_NULL) {
        int slabSize = NxGhost * NzGhost * nghost;
        int totalCount = slabSize * numArrays;
        
        std::fill(requestsY.begin(), requestsY.end(), MPI_REQUEST_NULL);
        reqCountY = 0;

        // Post receives
        if (rankMinusY != MPI_PROC_NULL) {
            MPI_Irecv(recvBufYMinus.data(), totalCount, MPI_DOUBLE,
                      rankMinusY, 1, MPI_COMM_WORLD,
                      &requestsY[reqCountY++]);
        }
        if (rankPlusY != MPI_PROC_NULL) {
            MPI_Irecv(recvBufYPlus.data(), totalCount, MPI_DOUBLE,
                      rankPlusY, 1, MPI_COMM_WORLD,
                      &requestsY[reqCountY++]);
        }

        // Pack & send
        if (rankPlusY != MPI_PROC_NULL) {
            packGhosts(sendBufYPlus, 1, true);
            MPI_Isend(sendBufYPlus.data(), totalCount, MPI_DOUBLE,
                      rankPlusY, 1, MPI_COMM_WORLD,
                      &requestsY[reqCountY++]);
        }
        if (rankMinusY != MPI_PROC_NULL) {
            packGhosts(sendBufYMinus, 1, false);
            MPI_Isend(sendBufYMinus.data(), totalCount, MPI_DOUBLE,
                      rankMinusY, 1, MPI_COMM_WORLD,
                      &requestsY[reqCountY++]);
        }

        MPI_Waitall(reqCountY, requestsY.data(), MPI_STATUSES_IGNORE);
        if (rankMinusY != MPI_PROC_NULL) {
            unpackGhosts(recvBufYMinus, 1, false);
        }
        if (rankPlusY != MPI_PROC_NULL) {
            unpackGhosts(recvBufYPlus, 1, true);
        }
    }
    
    // Z-direction exchange using non-blocking MPI
    if (rankMinusZ != MPI_PROC_NULL || rankPlusZ != MPI_PROC_NULL) {
        int slabSize = NxGhost * NyGhost * nghost;
        int totalCount = slabSize * numArrays;
        
        std::fill(requestsZ.begin(), requestsZ.end(), MPI_REQUEST_NULL);
        reqCountZ = 0;

        // Post receives
        if (rankMinusZ != MPI_PROC_NULL) {
            MPI_Irecv(recvBufZMinus.data(), totalCount, MPI_DOUBLE,
                      rankMinusZ, 2, MPI_COMM_WORLD,
                      &requestsZ[reqCountZ++]);
        }
        if (rankPlusZ != MPI_PROC_NULL) {
            MPI_Irecv(recvBufZPlus.data(), totalCount, MPI_DOUBLE,
                      rankPlusZ, 2, MPI_COMM_WORLD,
                      &requestsZ[reqCountZ++]);
        }

        // Pack & send
        if (rankPlusZ != MPI_PROC_NULL) {
            packGhosts(sendBufZPlus, 2, true);
            MPI_Isend(sendBufZPlus.data(), totalCount, MPI_DOUBLE,
                      rankPlusZ, 2, MPI_COMM_WORLD,
                      &requestsZ[reqCountZ++]);
        }
        if (rankMinusZ != MPI_PROC_NULL) {
            packGhosts(sendBufZMinus, 2, false);
            MPI_Isend(sendBufZMinus.data(), totalCount, MPI_DOUBLE,
                      rankMinusZ, 2, MPI_COMM_WORLD,
                      &requestsZ[reqCountZ++]);
        }

        MPI_Waitall(reqCountZ, requestsZ.data(), MPI_STATUSES_IGNORE);
        if (rankMinusZ != MPI_PROC_NULL) {
            unpackGhosts(recvBufZMinus, 2, false);
        }
        if (rankPlusZ != MPI_PROC_NULL) {
            unpackGhosts(recvBufZPlus, 2, true);
        }
    }

#else
    // --- Blocking MPI communication implementation ---

    const int numArrays = 6;
    
    // Calculate the maximum buffer size needed for any direction
    int slabSizeX = NyGhost * NzGhost * nghost;
    int slabSizeY = NxGhost * NzGhost * nghost;
    int slabSizeZ = NxGhost * NyGhost * nghost;
    int maxslabSize = std::max({slabSizeX, slabSizeY, slabSizeZ});
    int totalCount = maxslabSize * numArrays;
    
    // Reuse the same buffer for all directions
    std::vector<double> sendBuf(totalCount), recvBuf(totalCount);
    int tag = 0;

    // X-direction using MPI_Sendrecv with two buffers
    if (rankMinusX != MPI_PROC_NULL || rankPlusX != MPI_PROC_NULL) {
        // Phase 1: send plus-side data, receive from minus-side
        if (rankPlusX != MPI_PROC_NULL) packGhosts(sendBuf, 0, true);
        
        MPI_Sendrecv(sendBuf.data(), slabSizeX * numArrays, MPI_DOUBLE, rankPlusX, tag,
                     recvBuf.data(), slabSizeX * numArrays, MPI_DOUBLE, rankMinusX, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rankMinusX != MPI_PROC_NULL) unpackGhosts(recvBuf, 0, false);

        // Phase 2: send minus-side data, receive from plus-side
        if (rankMinusX != MPI_PROC_NULL) packGhosts(sendBuf, 0, false);
        
        MPI_Sendrecv(sendBuf.data(), slabSizeX * numArrays, MPI_DOUBLE, rankMinusX, tag,
                     recvBuf.data(), slabSizeX * numArrays, MPI_DOUBLE, rankPlusX, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rankPlusX != MPI_PROC_NULL) unpackGhosts(recvBuf, 0, true);
    }

    // --- Y-DIRECTION HALO EXCHANGE USING MPI_Sendrecv ---
    if (rankMinusY != MPI_PROC_NULL || rankPlusY != MPI_PROC_NULL) {
        // Phase 1: send plus-side data, receive minus-side
        if (rankPlusY != MPI_PROC_NULL) packGhosts(sendBuf, 1, true);
        
        MPI_Sendrecv(sendBuf.data(), slabSizeY * numArrays, MPI_DOUBLE, rankPlusY, tag,
                     recvBuf.data(), slabSizeY * numArrays, MPI_DOUBLE, rankMinusY, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rankMinusY != MPI_PROC_NULL) unpackGhosts(recvBuf, 1, false);

        // Phase 2: send minus-side data, receive plus-side
        if (rankMinusY != MPI_PROC_NULL) packGhosts(sendBuf, 1, false);
        
        MPI_Sendrecv(sendBuf.data(), slabSizeY * numArrays, MPI_DOUBLE, rankMinusY, tag,
                     recvBuf.data(), slabSizeY * numArrays, MPI_DOUBLE, rankPlusY, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rankPlusY != MPI_PROC_NULL) unpackGhosts(recvBuf, 1, true);
    }

    // --- Z-DIRECTION HALO EXCHANGE USING MPI_Sendrecv ---
    if (rankMinusZ != MPI_PROC_NULL || rankPlusZ != MPI_PROC_NULL) {
        // Phase 1: send to plus-side, receive from minus-side
        if (rankPlusZ != MPI_PROC_NULL) packGhosts(sendBuf, 2, true);
        
        MPI_Sendrecv(sendBuf.data(), slabSizeZ * numArrays, MPI_DOUBLE, rankPlusZ, tag,
                     recvBuf.data(), slabSizeZ * numArrays, MPI_DOUBLE, rankMinusZ, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rankMinusZ != MPI_PROC_NULL) unpackGhosts(recvBuf, 2, false);

        // Phase 2: send minus-side data, receive plus-side data
        if (rankMinusZ != MPI_PROC_NULL) packGhosts(sendBuf, 2, false);
        
        MPI_Sendrecv(sendBuf.data(), slabSizeZ * numArrays, MPI_DOUBLE, rankMinusZ, tag,
                     recvBuf.data(), slabSizeZ * numArrays, MPI_DOUBLE, rankPlusZ, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rankPlusZ != MPI_PROC_NULL) unpackGhosts(recvBuf, 2, true);
    }
#endif
}

void Field3D::allocateMPIBuffers() {
    const int numArrays = 6; // rho, rhou, rhov, rhow, E, phi

    // X-direction buffers
    int slabSizeX = NyGhost * NzGhost * nghost;
    int totalCountX = slabSizeX * numArrays;
    if (rankMinusX != MPI_PROC_NULL) {
        sendBufXMinus.resize(totalCountX);
        recvBufXMinus.resize(totalCountX);
    } else {
        std::vector<double>().swap(sendBufXMinus);
        std::vector<double>().swap(recvBufXMinus);
    }
    if (rankPlusX != MPI_PROC_NULL) {
        sendBufXPlus.resize(totalCountX);
        recvBufXPlus.resize(totalCountX);
    } else {
        std::vector<double>().swap(sendBufXPlus);
        std::vector<double>().swap(recvBufXPlus);
    }

    // Y-direction buffers
    int slabSizeY = NxGhost * NzGhost * nghost;
    int totalCountY = slabSizeY * numArrays;
    if (rankMinusY != MPI_PROC_NULL) {
        sendBufYMinus.resize(totalCountY);
        recvBufYMinus.resize(totalCountY);
    } else {
        std::vector<double>().swap(sendBufYMinus);
        std::vector<double>().swap(recvBufYMinus);
    }
    if (rankPlusY != MPI_PROC_NULL) {
        sendBufYPlus.resize(totalCountY);
        recvBufYPlus.resize(totalCountY);
    } else {
        std::vector<double>().swap(sendBufYPlus);
        std::vector<double>().swap(recvBufYPlus);
    }

    // Z-direction buffers
    int slabSizeZ = NxGhost * NyGhost * nghost;
    int totalCountZ = slabSizeZ * numArrays;
    if (rankMinusZ != MPI_PROC_NULL) {
        sendBufZMinus.resize(totalCountZ);
        recvBufZMinus.resize(totalCountZ);
    } else {
        std::vector<double>().swap(sendBufZMinus);
        std::vector<double>().swap(recvBufZMinus);
    }
    if (rankPlusZ != MPI_PROC_NULL) {
        sendBufZPlus.resize(totalCountZ);
        recvBufZPlus.resize(totalCountZ);
    } else {
        std::vector<double>().swap(sendBufZPlus);
        std::vector<double>().swap(recvBufZPlus);
    }
}

// Promote local pack to a public function:
void Field3D::packGhosts(std::vector<double> &buffer, int dim, bool sendPlus) const {
    // Set loop boundaries for all dimensions
    int iStart = 0, iEnd = NxGhost;
    int jStart = 0, jEnd = NyGhost;
    int kStart = 0, kEnd = NzGhost;
    
    // Adjust the loop boundaries for the target dimension
    switch (dim) {
        case 0: // X
            iStart = sendPlus ? (NxGhost - 2 * nghost) : nghost;
            iEnd = sendPlus ? (NxGhost - nghost) : (2 * nghost);
            break;
        case 1: // Y
            jStart = sendPlus ? (NyGhost - 2 * nghost) : nghost;
            jEnd = sendPlus ? (NyGhost - nghost) : (2 * nghost);
            break;
        case 2: // Z
            kStart = sendPlus ? (NzGhost - 2 * nghost) : nghost;
            kEnd = sendPlus ? (NzGhost - nghost) : (2 * nghost);
            break;
    }
    
    int cnt = 0;
    // Single triple loop with dynamic boundaries
    for (int k = kStart; k < kEnd; k++) {
        for (int j = jStart; j < jEnd; j++) {
            for (int i = iStart; i < iEnd; i++) {
                buffer[cnt++] = rho[index(i, j, k)];
                buffer[cnt++] = rhou[index(i, j, k)];
                buffer[cnt++] = rhov[index(i, j, k)];
                buffer[cnt++] = rhow[index(i, j, k)];
                buffer[cnt++] = E[index(i, j, k)];
                buffer[cnt++] = phi[index(i, j, k)];
            }
        }
    }
}

// Promote local unpack to a public function:
void Field3D::unpackGhosts(const std::vector<double> &buffer, int dim, bool plusSide) {
    // Set loop boundaries for all dimensions
    int iStart = 0, iEnd = NxGhost;
    int jStart = 0, jEnd = NyGhost;
    int kStart = 0, kEnd = NzGhost;
    
    // Adjust the loop boundaries for the target dimension
    switch (dim) {
        case 0: // X
            iStart = plusSide ? (NxGhost - nghost) : 0;
            iEnd = plusSide ? NxGhost : nghost;
            break;
        case 1: // Y
            jStart = plusSide ? (NyGhost - nghost) : 0;
            jEnd = plusSide ? NyGhost : nghost;
            break;
        case 2: // Z
            kStart = plusSide ? (NzGhost - nghost) : 0;
            kEnd = plusSide ? NzGhost : nghost;
            break;
    }
    
    int cnt = 0;
    // Single triple loop with dynamic boundaries
    for (int k = kStart; k < kEnd; k++) {
        for (int j = jStart; j < jEnd; j++) {
            for (int i = iStart; i < iEnd; i++) {
                rho[index(i, j, k)] = buffer[cnt++];
                rhou[index(i, j, k)] = buffer[cnt++];
                rhov[index(i, j, k)] = buffer[cnt++];
                rhow[index(i, j, k)] = buffer[cnt++];
                E[index(i, j, k)] = buffer[cnt++];
                phi[index(i, j, k)] = buffer[cnt++];
            }
        }
    }
}

// Methods migrated from BoundaryManager
void Field3D::initBCsFromParameters(const ParameterSystem &params) {
    bc_xmin = params.getBoundaryCondition("bc_xmin");
    bc_xmax = params.getBoundaryCondition("bc_xmax");
    bc_ymin = params.getBoundaryCondition("bc_ymin");
    bc_ymax = params.getBoundaryCondition("bc_ymax");
    bc_zmin = params.getBoundaryCondition("bc_zmin");
    bc_zmax = params.getBoundaryCondition("bc_zmax");

    std::cout << "[Field3D] Initialized BCs:\n"
              << "  X: (" << int(bc_xmin) << "," << int(bc_xmax) << ")\n"
              << "  Y: (" << int(bc_ymin) << "," << int(bc_ymax) << ")\n"
              << "  Z: (" << int(bc_zmin) << "," << int(bc_zmax) << ")\n";
}

int Field3D::getNeighborIndexX(int i, bool leftSide) const {
    if (leftSide) {
        // x-min
        if (i == 0) {
            if (bc_xmin == config::BoundaryCondition::PERIODIC) {
                return Nx - 1;
            } else if (bc_xmin == config::BoundaryCondition::OUTFLOW) {
                return 0;  // clamp
            }
        } else {
            return i - 1;
        }
    } else {
        // x-max
        if (i == Nx - 1) {
            if (bc_xmax == config::BoundaryCondition::PERIODIC) {
                return 0;
            } else if (bc_xmax == config::BoundaryCondition::OUTFLOW) {
                return Nx - 1;
            }
        } else {
            return i + 1;
        }
    }
    return i;  // fallback
}

int Field3D::getNeighborIndexY(int j, bool lowerSide) const {
    if (lowerSide) {
        // y-min
        if (j == 0) {
            if (bc_ymin == config::BoundaryCondition::PERIODIC) {
                return Ny - 1;
            } else if (bc_ymin == config::BoundaryCondition::OUTFLOW) {
                return 0;  // clamp
            }
        } else {
            return j - 1;
        }
    } else {
        // y-max
        if (j == Ny - 1) {
            if (bc_ymax == config::BoundaryCondition::PERIODIC) {
                return 0;
            } else if (bc_ymax == config::BoundaryCondition::OUTFLOW) {
                return Ny - 1;
            }
        } else {
            return j + 1;
        }
    }
    return j;
}

int Field3D::getNeighborIndexZ(int k, bool lowerSide) const {
    if (lowerSide) {
        // z-min
        if (k == 0) {
            if (bc_zmin == config::BoundaryCondition::PERIODIC) {
                return Nz - 1;
            } else if (bc_zmin == config::BoundaryCondition::OUTFLOW) {
                return 0;  // clamp
            }
        } else {
            return k - 1;
        }
    } else {
        // z-max
        if (k == Nz - 1) {
            if (bc_zmax == config::BoundaryCondition::PERIODIC) {
                return 0;
            } else if (bc_zmax == config::BoundaryCondition::OUTFLOW) {
                return Nz - 1;
            }
        } else {
            return k + 1;
        }
    }
    return k;
}

}  // namespace agoge
