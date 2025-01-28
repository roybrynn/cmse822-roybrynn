#pragma once

#include "Config.hpp"
#include "ParameterSystem.hpp"

namespace agoge {

/**
 * @class BoundaryManager
 * @brief A central place that owns boundary condition info for each face.
 *        Provides neighbor index logic for x,y,z directions.
 */
class BoundaryManager {
   public:
    /**
     * @brief Initialize BCs from ParameterSystem, e.g. "bc_xmin" => "outflow",
     * etc.
     */
    static void initBCsFromParameters(const ParameterSystem &params);

    /**
     * @brief Get neighbor index in X dimension for i +/- 1, using stored BC
     * info
     * @param i current index
     * @param Nx number of cells in x
     * @param leftSide true => i-1, false => i+1
     * @return neighbor index
     */
    static int getNeighborIndexX(int i, int Nx, bool leftSide);

    /**
     * @brief Similarly for Y dimension
     */
    static int getNeighborIndexY(int j, int Ny, bool lowerSide);

    /**
     * @brief Similarly for Z dimension
     */
    static int getNeighborIndexZ(int k, int Nz, bool lowerSide);

   private:
    // store BC for each face
    static config::BoundaryCondition bc_xmin;
    static config::BoundaryCondition bc_xmax;
    static config::BoundaryCondition bc_ymin;
    static config::BoundaryCondition bc_ymax;
    static config::BoundaryCondition bc_zmin;
    static config::BoundaryCondition bc_zmax;
};

}  // namespace agoge