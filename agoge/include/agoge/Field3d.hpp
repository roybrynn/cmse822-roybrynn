#pragma once

#include <vector>

/**
 * @file Field3d.hpp
 * @brief Declaration of the Field3D structure for the Agoge solver.
 *
 * This structure stores the primary 3D fields (density, momentum components,
 * energy, and gravitational potential) for each cell in the computational domain.
 */

namespace agoge {

/**
 * @struct Field3D
 * @brief Data structure storing the primary 3D fields for the Agoge solver.
 */
struct Field3D
{
    /**
     * @brief Number of cells in the x-direction.
     */
    int Nx;
    /**
     * @brief Number of cells in the y-direction.
     */
    int Ny;
    /**
     * @brief Number of cells in the z-direction.
     */
    int Nz;

    /**
     * @brief Cell size in the x-direction.
     */
    double dx;
    /**
     * @brief Cell size in the y-direction.
     */
    double dy;
    /**
     * @brief Cell size in the z-direction.
     */
    double dz;

    /**
     * @brief Density field: rho[i + Nx*(j + Ny*k)].
     */
    std::vector<double> rho;
    /**
     * @brief Momentum in the x-direction: rhou[i + Nx*(j + Ny*k)].
     */
    std::vector<double> rhou;
    /**
     * @brief Momentum in the y-direction: rhov[i + Nx*(j + Ny*k)].
     */
    std::vector<double> rhov;
    /**
     * @brief Momentum in the z-direction: rhow[i + Nx*(j + Ny*k)].
     */
    std::vector<double> rhow;
    /**
     * @brief Energy field: E[i + Nx*(j + Ny*k)].
     */
    std::vector<double> E;
    /**
     * @brief Gravitational potential: phi[i + Nx*(j + Ny*k)].
     */
    std::vector<double> phi;

    /**
     * @brief Constructor that allocates the SoA vectors based on the grid size.
     *
     * @param nx Number of cells in x
     * @param ny Number of cells in y
     * @param nz Number of cells in z
     * @param dx_ Cell size in x
     * @param dy_ Cell size in y
     * @param dz_ Cell size in z
     */
    Field3D(int nx, int ny, int nz, double dx_, double dy_, double dz_);

    /**
     * @brief Returns the 1D index in the data arrays for cell (i, j, k).
     *
     * @param i x-index
     * @param j y-index
     * @param k z-index
     * @return The flattened 1D index
     */
    inline int index(int i, int j, int k) const
    {
        return i + Nx * (j + Ny * k);
    }
};

} // namespace agoge
