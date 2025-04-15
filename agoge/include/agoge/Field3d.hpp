#pragma once

#include <mpi.h>

#include <vector>
#include <array>

#include "Config.hpp"  // for BoundaryCondition

namespace agoge {

// Define a BoundingBox structure
struct BoundingBox {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
};

// A 3D field data structure with ghost zones for boundary conditions
class Field3D {
   public:
    // Constructor: nx,ny,nz are interior cells, dx,dy,dz are cell sizes
    // nghost specifies number of ghost cells on each side
    Field3D(int nx, int ny, int nz, const BoundingBox& bbox, int nghost = 0);

    // Explicit copy and assignment constructors
    Field3D(const Field3D&) = default;
    Field3D& operator=(const Field3D&) = default;

    // Accessor for BoundingBox
    const BoundingBox& getBoundingBox() const { return bbox; }

    // The total domain's interior size
    int Nx, Ny, Nz;

    // The total including ghost cells
    int NxGhost, NyGhost, NzGhost;

    // Number of ghost cells per side
    int nghost;

    // Physical cell sizes
    double dx, dy, dz;

    // The bounding box for the field
    BoundingBox bbox;

    // New metadata for global domain reconstruction:
    // Global bounding box for the entire domain
    BoundingBox global_bbox;
    // Global interior grid dimensions
    int global_Nx, global_Ny, global_Nz;
    // Total subdomains in each direction
    int Px, Py, Pz;
    // rank id 
    int myRank = 0; // Default to rank 0
    // Comm world size
    int mpiSize;
    // This rank's subdomain indices
    int subdomain_x, subdomain_y, subdomain_z;

    /**
     * @brief The arrays: size NxGhost*NyGhost*NzGhost
     */
    std::vector<double> rho;
    std::vector<double> rhou;
    std::vector<double> rhov;
    std::vector<double> rhow;
    std::vector<double> E;
    std::vector<double> phi;

    // Indexing for the full array (including ghosts)
    // i in [0..NxGhost-1], j in [0..NyGhost-1], k in [0..NzGhost-1]
    int index(int i, int j, int k) const {
        return i + NxGhost * (j + NyGhost * k);
    }

    // Convert interior cell (iIn, jIn, kIn) in [0..Nx-1] to actual ghosted
    // index
    int interiorIndex(int iIn, int jIn, int kIn) const {
        int i = iIn + nghost;
        int j = jIn + nghost;
        int k = kIn + nghost;
        return index(i, j, k);
    }

    // Apply boundary conditions to fill ghost zones
    void applyBCs();

    // Boundary condition flags
    config::BoundaryCondition bc_xmin, bc_xmax;
    config::BoundaryCondition bc_ymin, bc_ymax;
    config::BoundaryCondition bc_zmin, bc_zmax;

    // Return x-center coordinate for interior cell iIn in [0..Nx-1]
    // Assumes origin=0 for interior cell i=0 at left edge plus ghost offset
    double xCenter(int iIn) const;

    // Edge getters for x dimension
    double xLeftEdge(int iIn) const;
    double xRightEdge(int iIn) const;

    // Edge and center getters for y dimension
    double yCenter(int jIn) const;
    double yLeftEdge(int jIn) const;
    double yRightEdge(int jIn) const;

    // Edge and center getters for z dimension
    double zCenter(int kIn) const;
    double zLeftEdge(int kIn) const;
    double zRightEdge(int kIn) const;

    // MPI neighbor ranks
    int rankMinusX, rankPlusX;
    int rankMinusY, rankPlusY;
    int rankMinusZ, rankPlusZ;

    // MPI communication buffers
    std::vector<double> sendBufXMinus, sendBufXPlus;
    std::vector<double> recvBufXMinus, recvBufXPlus;
    std::vector<double> sendBufYMinus, sendBufYPlus;
    std::vector<double> recvBufYMinus, recvBufYPlus;
    std::vector<double> sendBufZMinus, sendBufZPlus;
    std::vector<double> recvBufZMinus, recvBufZPlus;

    // MPI request arrays for non-blocking communication
    std::array<MPI_Request, 4> requestsX;  // e.g. recvMinus, recvPlus, sendPlus, sendMinus
    std::array<MPI_Request, 4> requestsY;
    std::array<MPI_Request, 4> requestsZ;

    int reqCountX = 0;
    int reqCountY = 0;
    int reqCountZ = 0;
    
    /**
     * @brief Allocate MPI communication buffers based on current dimensions and neighbors.
     */
    void allocateMPIBuffers();

    /**
     * @brief Pack ghost zone data for MPI communication.
     * @param buffer Output buffer.
     * @param dim 0 => X, 1 => Y, 2 => Z.
     * @param sendPlus Whether to send to the plus side.
     */
    void packGhosts(std::vector<double> &buffer, int dim, bool sendPlus) const;

    /**
     * @brief Unpack received ghost zone data.
     * @param buffer Input buffer.
     * @param dim 0 => X, 1 => Y, 2 => Z.
     * @param plusSide Whether data is from the plus side.
     */
    void unpackGhosts(const std::vector<double> &buffer, int dim, bool plusSide);

    /**
     * @brief Initialize boundary conditions from ParameterSystem
     * @param params The parameter system containing boundary condition settings
     */
    void initBCsFromParameters(const class ParameterSystem &params);

    /**
     * @brief Get neighbor index in X dimension based on boundary conditions
     * @param i Current X index
     * @param leftSide true => i-1, false => i+1
     * @return Neighbor index accounting for boundary conditions
     */
    int getNeighborIndexX(int i, bool leftSide) const;

    /**
     * @brief Get neighbor index in Y dimension based on boundary conditions
     * @param j Current Y index
     * @param lowerSide true => j-1, false => j+1
     * @return Neighbor index accounting for boundary conditions
     */
    int getNeighborIndexY(int j, bool lowerSide) const;

    /**
     * @brief Get neighbor index in Z dimension based on boundary conditions
     * @param k Current Z index
     * @param lowerSide true => k-1, false => k+1
     * @return Neighbor index accounting for boundary conditions
     */
    int getNeighborIndexZ(int k, bool lowerSide) const;
};

}  // namespace agoge