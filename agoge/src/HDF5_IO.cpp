/**
 * @file HDF5_IO.cpp
 * @brief Handles HDF5 input/output routines for the Agoge application.
 *
 * This file implements functions to write Field3D data and associated metadata
 * to HDF5 files.
 */

#include <H5Cpp.h>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <mpi.h>  // New include

#include <iostream>
#include <vector>

#include "Field3d.hpp"

namespace agoge {
namespace io {

// New: global IO epoch counter and total ranks
static int g_ioEpochCounter = 0;
static int g_totalRanks = 1;

// Modified init helper: Remove directories matching the output naming scheme.
void initOutputDirectory(const std::string &pattern) {
    // Remove matching directories and their contents.
    std::string removeCmd = "rm -rf " + pattern;
    system(removeCmd.c_str());
}

/**
 * @brief Writes the Field3D data to an HDF5 file.
 *
 * The function creates datasets for primary variables and stores both local 
 * and global grid metadata.
 *
 * @param Q The Field3D instance containing simulation data.
 * @param filename The output HDF5 file path.
 * @param current_time The current simulation time.
 */
void writeFieldHDF5(const Field3D &Q, const std::string &filename, double current_time = 0.0) {
    using namespace H5;
    try {
        // Create (or truncate) the file
        H5File file(filename.c_str(), H5F_ACC_TRUNC);

        // Write the field datasets (rho, etc.)
        auto writeDataset = [&](const std::string &name,
                                const std::vector<double> &data) {
            // Dimensions for the dataset: (Nz, Ny, Nx)
            hsize_t dims[3] = {static_cast<hsize_t>(Q.Nz),
                               static_cast<hsize_t>(Q.Ny),
                               static_cast<hsize_t>(Q.Nx)};
            DataSpace dataspace(3, dims);
            DataSet dataset =
                file.createDataSet(name, PredType::NATIVE_DOUBLE, dataspace);

            // Rearrange data from Field3D (stored as 1D using interiorIndex)
            // into a buffer with (k, j, i) ordering
            std::vector<double> buffer(Q.Nx * Q.Ny * Q.Nz);
            for (int k = 0; k < Q.Nz; ++k) {
                for (int j = 0; j < Q.Ny; ++j) {
                    for (int i = 0; i < Q.Nx; ++i) {
                        size_t idxQ = Q.interiorIndex(i, j, k);
                        size_t idxBuf = static_cast<size_t>(k) * (Q.Ny * Q.Nx) +
                                        static_cast<size_t>(j) * Q.Nx +
                                        static_cast<size_t>(i);
                        buffer[idxBuf] = data[idxQ];
                    }
                }
            }
            dataset.write(buffer.data(), PredType::NATIVE_DOUBLE);
        };

        // Write each field dataset.
        writeDataset("rho", Q.rho);
        writeDataset("rhou", Q.rhou);
        writeDataset("rhov", Q.rhov);
        writeDataset("rhow", Q.rhow);
        writeDataset("E", Q.E);
        writeDataset("phi", Q.phi);

        // Create a group for grid metadata and coordinate arrays.
        Group grid = file.createGroup("/grid");

        // Write interior bounding box as an attribute.
        hsize_t bbox_dims[1] = {6};
        DataSpace bbox_space(1, bbox_dims);
        const auto &bbox = Q.getBoundingBox();
        double bbox_data[6] = {bbox.xmin, bbox.xmax, bbox.ymin,
                                bbox.ymax, bbox.zmin, bbox.zmax};
        Attribute attr_bbox = grid.createAttribute("bounding_box", PredType::NATIVE_DOUBLE, bbox_space);
        attr_bbox.write(PredType::NATIVE_DOUBLE, bbox_data);

        // Existing: Store domain dimensions (interior) and cell sizes.
        hsize_t dims_attr[1] = {3};
        DataSpace dims_space(1, dims_attr);
        int domain_dims[3] = {Q.Nx, Q.Ny, Q.Nz};
        Attribute attr_dims = grid.createAttribute("domain_dimensions", PredType::NATIVE_INT, dims_space);
        attr_dims.write(PredType::NATIVE_INT, domain_dims);

        double cell_sizes[3] = {Q.dx, Q.dy, Q.dz};
        Attribute attr_size = grid.createAttribute("cell_size", PredType::NATIVE_DOUBLE, dims_space);
        attr_size.write(PredType::NATIVE_DOUBLE, cell_sizes);

        // --- New metadata attributes ---

        // Global bounding box (from Q.global_bbox)
        hsize_t gbox_dims[1] = {6};
        DataSpace gbox_space(1, gbox_dims);
        const auto &gb = Q.global_bbox;
        double gbox_data[6] = {gb.xmin, gb.xmax, gb.ymin, gb.ymax, gb.zmin, gb.zmax};
        Attribute attr_global_bbox = grid.createAttribute("global_bbox", PredType::NATIVE_DOUBLE, gbox_space);
        attr_global_bbox.write(PredType::NATIVE_DOUBLE, gbox_data);

        // Global domain dimensions (global interior sizes)
        int global_dims[3] = {Q.global_Nx, Q.global_Ny, Q.global_Nz};
        Attribute attr_global_dims = grid.createAttribute("global_dimensions", PredType::NATIVE_INT, dims_space);
        attr_global_dims.write(PredType::NATIVE_INT, global_dims);

        // MPI decomposition: total subdomains in x, y, z (Px, Py, Pz)
        int mpi_info[3] = {Q.Px, Q.Py, Q.Pz};
        Attribute attr_mpi = grid.createAttribute("mpi_decomposition", PredType::NATIVE_INT, dims_space);
        attr_mpi.write(PredType::NATIVE_INT, mpi_info);

        // Write Px, Py, and Pz as separate attributes.
        int px = Q.Px, py = Q.Py, pz = Q.Pz;
        DataSpace scalar_space = DataSpace(H5S_SCALAR);
        Attribute attr_px = grid.createAttribute("Px", PredType::NATIVE_INT, scalar_space);
        attr_px.write(PredType::NATIVE_INT, &px);
        Attribute attr_py = grid.createAttribute("Py", PredType::NATIVE_INT, scalar_space);
        attr_py.write(PredType::NATIVE_INT, &py);
        Attribute attr_pz = grid.createAttribute("Pz", PredType::NATIVE_INT, scalar_space);
        attr_pz.write(PredType::NATIVE_INT, &pz);

        // Subdomain indices for this rank: [subdomain_x, subdomain_y, subdomain_z]
        int subdomain[3] = {Q.subdomain_x, Q.subdomain_y, Q.subdomain_z};
        Attribute attr_subdomain = grid.createAttribute("subdomain", PredType::NATIVE_INT, dims_space);
        attr_subdomain.write(PredType::NATIVE_INT, subdomain);

        // Local dimensions (interior sizes on this rank)
        int local_interior[3] = {Q.Nx, Q.Ny, Q.Nz};
        Attribute attr_local = grid.createAttribute("local_dimensions", PredType::NATIVE_INT, dims_space);
        attr_local.write(PredType::NATIVE_INT, local_interior);

        // Add current simulation time as an attribute
        Attribute attr_time = grid.createAttribute("simulation_time", PredType::NATIVE_DOUBLE, scalar_space);
        attr_time.write(PredType::NATIVE_DOUBLE, &current_time);

    } catch (const H5::Exception &err) {
        std::cerr << "HDF5 write error: " << err.getDetailMsg() << std::endl;
        throw;
    }
}

/**
 * @brief Performs field I/O using HDF5, automatically naming the output
 *        directory and file using zero-padded epoch and rank.
 *
 * @param Q The Field3D data to write.
 * @param problemName The problem name.
 * @param rank The MPI process rank.
 * @param current_time The current simulation time.
 * @param outputDir The base output directory.
 */
void performFieldIO(const Field3D &Q, const std::string &problemName, int rank, 
                    double current_time, const std::string &outputDir) {
    MPI_Comm_size(MPI_COMM_WORLD, &g_totalRanks);

    // Build directory name e.g., "output_dir/Sedov_004ranks_0000"
    std::ostringstream dirStream;
    dirStream << outputDir << "/" << problemName << "_" << g_totalRanks << "ranks_" 
              << std::setw(4) << std::setfill('0') << g_ioEpochCounter;
    std::string dirName = dirStream.str();

    // On each output run, let only rank 0 remove directories matching the
    // naming scheme.
    if (rank == 0) {
        // Build the specific directory pattern for this epoch
        std::ostringstream patternStream;
        patternStream << outputDir << "/" << problemName << "_" << g_totalRanks
                      << "ranks_" << std::setw(4) << std::setfill('0')
                      << g_ioEpochCounter;
        std::string pattern = patternStream.str();

        // Remove any existing directory for this specific epoch
        initOutputDirectory(pattern);
    }
    
    // Synchronize all ranks so that removal is complete.
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Create directory (POSIX system call)
    std::string mkdirCmd = "mkdir -p " + dirName;
    system(mkdirCmd.c_str());
    
    // Build file name using output directory
    std::ostringstream fileStream;
    fileStream << dirName << "/rank_" << std::setw(4) << std::setfill('0') << rank << ".h5";
    std::string fileName = fileStream.str();
    
    if (rank == 0) std::cout << "Writing HDF5 files: " << dirName << " at time = " << current_time << "\n";
    writeFieldHDF5(Q, fileName, current_time);
    g_ioEpochCounter++;
}

}  // namespace io
}  // namespace agoge