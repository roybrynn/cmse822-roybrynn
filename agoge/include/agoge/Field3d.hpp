// include/agoge/Field3d.hpp
#pragma once

#include <vector>

/**
 * @file Field3d.hpp
 * @brief Declaration of the Field3D structure for the Agoge solver.
 *
 * This structure stores the primary 3D fields (density, momentum components,
 * energy, and gravitational potential) for each cell in the computational
 * domain.
 */

namespace agoge {

/**
 * @struct Field3D
 * @brief Data structure storing the primary 3D fields for the Agoge solver.
 */
class Field3D {
   public:
    Field3D(int nx, int ny, int nz, double dx, double dy, double dz);

    // Indexing method
    int index(int i, int j, int k) const { return i + Nx * (j + Ny * k); }

    // Members
    /**
     * @brief Number of cells in the {x,y,z}-direction.
     */
    int Nx, Ny, Nz;
    /**
     * @brief Cell size in the x,y,z-direction.
     */
    double dx, dy, dz;
    /**
     * @brief Density/momenta/energy fields: rho[i + Nx*(j + Ny*k)].
     */
    std::vector<double> rho;
    std::vector<double> rhou;
    std::vector<double> rhov;
    std::vector<double> rhow;
    std::vector<double> E;
    std::vector<double> phi;

};

}  // namespace agoge
