#include "BoundaryManager.hpp"

#include <iostream>

namespace agoge {

// Static members
config::BoundaryCondition BoundaryManager::bc_xmin =
    config::BoundaryCondition::PERIODIC;
config::BoundaryCondition BoundaryManager::bc_xmax =
    config::BoundaryCondition::PERIODIC;
config::BoundaryCondition BoundaryManager::bc_ymin =
    config::BoundaryCondition::PERIODIC;
config::BoundaryCondition BoundaryManager::bc_ymax =
    config::BoundaryCondition::PERIODIC;
config::BoundaryCondition BoundaryManager::bc_zmin =
    config::BoundaryCondition::PERIODIC;
config::BoundaryCondition BoundaryManager::bc_zmax =
    config::BoundaryCondition::PERIODIC;

void BoundaryManager::initBCsFromParameters(const ParameterSystem &params) {
    bc_xmin = params.getBoundaryCondition("bc_xmin");
    bc_xmax = params.getBoundaryCondition("bc_xmax");
    bc_ymin = params.getBoundaryCondition("bc_ymin");
    bc_ymax = params.getBoundaryCondition("bc_ymax");
    bc_zmin = params.getBoundaryCondition("bc_zmin");
    bc_zmax = params.getBoundaryCondition("bc_zmax");

    std::cout << "[BoundaryManager] Initialized BCs:\n"
              << "  X: (" << int(bc_xmin) << "," << int(bc_xmax) << ")\n"
              << "  Y: (" << int(bc_ymin) << "," << int(bc_ymax) << ")\n"
              << "  Z: (" << int(bc_zmin) << "," << int(bc_zmax) << ")\n";
}

int BoundaryManager::getNeighborIndexX(int i, int Nx, bool leftSide) {
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

int BoundaryManager::getNeighborIndexY(int j, int Ny, bool lowerSide) {
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

int BoundaryManager::getNeighborIndexZ(int k, int Nz, bool lowerSide) {
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
