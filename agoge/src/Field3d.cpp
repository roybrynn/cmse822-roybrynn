#include "Field3d.hpp"

#include <algorithm>
#include <iostream>

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
      bc_zmax(config::BoundaryCondition::PERIODIC) {
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

double Field3D::xCenter(int iIn) const { return (iIn + 0.5) * dx; }
double Field3D::xLeftEdge(int iIn) const { return (double)(iIn)*dx; }
double Field3D::xRightEdge(int iIn) const { return (double)(iIn + 1) * dx; }

double Field3D::yCenter(int jIn) const { return (jIn + 0.5) * dy; }
double Field3D::yLeftEdge(int jIn) const { return jIn * dy; }
double Field3D::yRightEdge(int jIn) const { return (jIn + 1) * dy; }

double Field3D::zCenter(int kIn) const { return (kIn + 0.5) * dz; }
double Field3D::zLeftEdge(int kIn) const { return kIn * dz; }
double Field3D::zRightEdge(int kIn) const { return (kIn + 1) * dz; }

/**
 * @brief copy boundary or wrap for periodic
 */
static inline void copyCell(std::vector<double> &arr, int dst, int src) {
    arr[dst] = arr[src];
}

void Field3D::applyBCs() {
    // We'll fill ghost zones in x, y, z
    // Example for x-min side: i in [0..nghost-1], we do one of:
    //   if PERIODIC => copy from i + Nx
    //   if OUTFLOW => copy from i=nghost
    //
    // We'll do it for each array: rho, rhou, rhov, rhow, E, phi
    // We'll define a small lambda to handle 1 variable at a time:

    auto fillBC_x = [&](std::vector<double> &arr) {
        // x-min
        for (int k = 0; k < NzGhost; k++) {
            for (int j = 0; j < NyGhost; j++) {
                for (int iGhost = 0; iGhost < nghost; iGhost++) {
                    int iInterior = nghost;  // the interior boundary
                    int ghostIdx = index(iGhost, j, k);

                    if (bc_xmin == config::BoundaryCondition::PERIODIC) {
                        // map iGhost => iInteriorLeft + Nx - nghost
                        // effectively i = Nx + iGhost
                        int iMirror = Nx + iGhost;  // wrap
                        if (iMirror >= Nx)
                            iMirror -= Nx;  // if needed
                        int src = index(iMirror + nghost, j, k);
                        copyCell(arr, ghostIdx, src);
                    } else if (bc_xmin == config::BoundaryCondition::OUTFLOW) {
                        // clamp to iInteriorLeft
                        int src = index(iInterior, j, k);
                        copyCell(arr, ghostIdx, src);
                    }
                }
            }
        }

        // x-max
        for (int k = 0; k < NzGhost; k++) {
            for (int j = 0; j < NyGhost; j++) {
                for (int iGhost = NxGhost - nghost; iGhost < NxGhost; iGhost++) {
                    int iInterior = NxGhost - nghost - 1;  // interior boundary
                    int ghostIdx = index(iGhost, j, k);

                    if (bc_xmax == config::BoundaryCondition::PERIODIC) {
                        // OLD code:
                        // int iMirror = iGhost - Nx;  <-- ends up referencing
                        // ghost cells Fix: offset by nghost so iMirror points
                        // into the real interior
                        int iMirror = (iGhost - Nx) + nghost;
                        int src = index(iMirror, j, k);
                        copyCell(arr, ghostIdx, src);
                    } else if (bc_xmax == config::BoundaryCondition::OUTFLOW) {
                        int src = index(iInterior, j, k);
                        copyCell(arr, ghostIdx, src);
                    }
                }
            }
        }
    };

    // Do similarly for fillBC_y, fillBC_z or unify logic:
    auto fillBC_y = [&](std::vector<double> &arr) {
        // y-min
        for (int k = 0; k < NzGhost; k++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int jGhost = 0; jGhost < nghost; jGhost++) {
                    int jInterior = nghost;
                    int ghostIdx = index(i, jGhost, k);

                    if (bc_ymin == config::BoundaryCondition::PERIODIC) {
                        int jMirror = Ny + jGhost;
                        if (jMirror >= Ny) jMirror -= Ny;
                        int src = index(i, jMirror + nghost, k);
                        copyCell(arr, ghostIdx, src);
                    } else if (bc_ymin == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, jInterior, k);
                        copyCell(arr, ghostIdx, src);
                    }
                }
            }
        }
        // y-max
        for (int k = 0; k < NzGhost; k++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int jGhost = nghost + Ny; jGhost < NyGhost;
                     jGhost++) {
                    int jInterior = nghost + Ny - 1;
                    int ghostIdx = index(i, jGhost, k);
                    if (bc_ymax == config::BoundaryCondition::PERIODIC) {
                        int jMirror = jGhost - Ny;
                        int src = index(i, jMirror, k);
                        copyCell(arr, ghostIdx, src);
                    } else if (bc_ymax == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, jInterior, k);
                        copyCell(arr, ghostIdx, src);
                    }
                }
            }
        }
    };

    auto fillBC_z = [&](std::vector<double> &arr) {
        // z-min
        for (int j = 0; j < NyGhost; j++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int kGhost = 0; kGhost < nghost; kGhost++) {
                    int kInterior = nghost;
                    int ghostIdx = index(i, j, kGhost);
                    if (bc_zmin == config::BoundaryCondition::PERIODIC) {
                        int kMirror = Nz + kGhost;
                        if (kMirror >= Nz) kMirror -= Nz;
                        int src = index(i, j, kMirror + nghost);
                        copyCell(arr, ghostIdx, src);
                    } else if (bc_zmin == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, j, kInterior);
                        copyCell(arr, ghostIdx, src);
                    }
                }
            }
        }
        // z-max
        for (int j = 0; j < NyGhost; j++) {
            for (int i = 0; i < NxGhost; i++) {
                for (int kGhost = nghost + Nz; kGhost < NzGhost;
                     kGhost++) {
                    int kInterior = nghost + Nz - 1;
                    int ghostIdx = index(i, j, kGhost);
                    if (bc_zmax == config::BoundaryCondition::PERIODIC) {
                        int kMirror = kGhost - Nz;
                        int src = index(i, j, kMirror);
                        copyCell(arr, ghostIdx, src);
                    } else if (bc_zmax == config::BoundaryCondition::OUTFLOW) {
                        int src = index(i, j, kInterior);
                        copyCell(arr, ghostIdx, src);
                    }
                }
            }
        }
    };

    // Now call them for each array:
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
}

}  // namespace agoge
