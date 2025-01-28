#pragma once

#include <vector>

#include "Config.hpp"  // for BoundaryCondition

namespace agoge {

/**
 * @class Field3D
 * @brief Data structure storing 3D fields with ghost zones for boundary
 * conditions.
 */
class Field3D {
   public:
    /**
     * @brief Constructor
     * @param nx,ny,nz interior cells
     * @param dx,dy,dz cell sizes
     * @param nghost number of ghost cells on each side
     */
    Field3D(int nx, int ny, int nz, double dx, double dy, double dz,
            int nghost = 0);

    /**
     * @brief The total domain's interior size
     */
    int Nx, Ny, Nz;

    /**
     * @brief The total including ghost
     */
    int NxGhost, NyGhost, NzGhost;

    /**
     * @brief how many ghost cells per side
     */
    int nghost;

    /**
     * @brief Physical cell sizes
     */
    double dx, dy, dz;

    /**
     * @brief The arrays: size NxGhost*NyGhost*NzGhost
     */
    std::vector<double> rho;
    std::vector<double> rhou;
    std::vector<double> rhov;
    std::vector<double> rhow;
    std::vector<double> E;
    std::vector<double> phi;

    /**
     * @brief Indexing for the full array (including ghosts).
     * i in [0..NxGhost-1], j in [0..NyGhost-1], k in [0..NzGhost-1].
     */
    int index(int i, int j, int k) const {
        return i + NxGhost * (j + NyGhost * k);
    }

    /**
     * @brief Convert interior cell (iIn, jIn, kIn) in [0..Nx-1] =>
     * actual ghosted index
     */
    int interiorIndex(int iIn, int jIn, int kIn) const {
        int i = iIn + nghost;
        int j = jIn + nghost;
        int k = kIn + nghost;
        return index(i, j, k);
    }

    /**
     * @brief Apply boundary conditions to fill ghost zones.
     * We store BC types in bcXmin, bcXmax, etc. for demonstration.
     */
    void applyBCs();

    // boundary condition flags
    config::BoundaryCondition bc_xmin, bc_xmax;
    config::BoundaryCondition bc_ymin, bc_ymax;
    config::BoundaryCondition bc_zmin, bc_zmax;
};

}  // namespace agoge
