#pragma once

#include <vector>

#include "Config.hpp"  // for BoundaryCondition

namespace agoge {

// Define a BoundingBox structure
struct BoundingBox {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
};

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
    // Updated constructor with BoundingBox and default ghost value
    Field3D(int nx, int ny, int nz, const BoundingBox &bbox, int nghost = 0);

    // Accessor for BoundingBox
    const BoundingBox& getBoundingBox() const { return bbox; }

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
     * @brief The bounding box for the field
     */
    BoundingBox bbox;

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

    /**
     * @brief Return x-center coordinate for interior cell iIn in [0..Nx-1].
     *        We assume origin=0 for interior cell i=0 is at left edge
     *        plus ghost offset => so the cell-center is (iIn + 0.5)*dx,
     *        plus possibly an offset for the ghost region: iTotal= iIn + nghost
     */
    double xCenter(int iIn) const;

    /**
     * @brief xLeftEdge for the interior cell iIn
     */
    double xLeftEdge(int iIn) const;

    /**
     * @brief xRightEdge for the interior cell iIn
     */
    double xRightEdge(int iIn) const;

    /**
     * @brief Similarly for y
     */
    double yCenter(int jIn) const;
    double yLeftEdge(int jIn) const;
    double yRightEdge(int jIn) const;

    /**
     * @brief Similarly for z
     */
    double zCenter(int kIn) const;
    double zLeftEdge(int kIn) const;
    double zRightEdge(int kIn) const;
};

}  // namespace agoge
