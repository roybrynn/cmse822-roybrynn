#pragma once

#include "Field3d.hpp"
#include <string>

/**
 * @file HDF5_IO.hpp
 * @brief Declarations for HDF5 read/write routines in the Agoge project.
 *
 * Extra metadata (global_bbox, global_dimensions, mpi_decomposition, 
 * subdomain, local_dimensions, cell_size, and simulation_time) 
 * are stored in the /grid attribute group.
 */

namespace agoge {
namespace io {

/**
 * @brief Write the fields (`rho`, `rhou`, `rhov`, `rhow`, `E`, `phi`) from a
 * Field3D object to an HDF5 file.
 *
 * The data is written as 3D datasets of shape (Nz, Ny, Nx) in row-major
 * order, meaning index (k, j, i) is contiguous in the innermost dimension i.
 *
 * @param[in] Q            The Field3D data structure to write.
 * @param[in] filename     The output HDF5 file name.
 * @param[in] current_time The current simulation time (default: 0.0).
 */
void writeFieldHDF5(const Field3D &Q, const std::string &filename, double current_time = 0.0);

/**
 * @brief (Optional) Read the fields (`rho`, `rhou`, `rhov`, `rhow`, `E`, `phi`)
 * from an HDF5 file into a Field3D object.
 *
 * This function expects the same shape (Nz, Ny, Nx) that was used when writing.
 * If domain-size attributes are present, it can optionally check them.
 *
 * @param[out] Q         The Field3D data structure to populate.
 * @param[in]  filename  The input HDF5 file name.
 */
void readFieldHDF5(Field3D &Q, const std::string &filename);

/**
 * @brief Performs field I/O using HDF5, automatically naming the output
 *        directory and file using zero-padded epoch and rank.
 *
 * @param Q            The Field3D data to write.
 * @param problemName  The problem name.
 * @param rank         The MPI process rank.
 * @param current_time The current simulation time (default: 0.0).
 * @param outputDir    The base output directory (default: "output").
 */
void performFieldIO(const Field3D &Q, const std::string &problemName, int rank, 
                    double current_time = 0.0, const std::string &outputDir = "output");

} // namespace io
} // namespace agoge
