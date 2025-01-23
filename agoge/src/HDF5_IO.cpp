// HDF5_IO.cpp
#include "HDF5_IO.hpp"

#include <H5Cpp.h>

#include <iostream>
#include <vector>

#include "Field3d.hpp"

namespace agoge {
namespace io {

void writeFieldHDF5(const Field3D &Q, const std::string &filename) {
    using namespace H5;

    // Create / open the file (truncate if exists)
    H5File file(filename.c_str(), H5F_ACC_TRUNC);

    // Dimensions for the 3D datasets: (Nz, Ny, Nx)
    hsize_t dims[3] = {static_cast<hsize_t>(Q.Nz), static_cast<hsize_t>(Q.Ny),
                       static_cast<hsize_t>(Q.Nx)};

    // Create a data space for the 3D shape
    DataSpace dataspace(3, dims);

    // Helper lambda to write a 3D dataset from Q.* vectors
    auto writeDataset = [&](const std::string &name,
                            const std::vector<double> &data) {
        // Create the dataset
        DataSet dataset =
            file.createDataSet(name, PredType::NATIVE_DOUBLE, dataspace);

        // Rearrange data from (i + Nx*(j + Ny*k)) to (k, j, i)
        std::vector<double> buffer(data.size());
        for (int k = 0; k < Q.Nz; ++k) {
            for (int j = 0; j < Q.Ny; ++j) {
                for (int i = 0; i < Q.Nx; ++i) {
                    size_t idxQ = i + Q.Nx * (j + Q.Ny * k);
                    size_t idxBuf = static_cast<size_t>(k) * (Q.Ny * Q.Nx) +
                                    static_cast<size_t>(j) * Q.Nx +
                                    static_cast<size_t>(i);
                    buffer[idxBuf] = data[idxQ];
                }
            }
        }

        // Write to the dataset
        dataset.write(buffer.data(), PredType::NATIVE_DOUBLE);
    };

    // Write each field
    writeDataset("rho", Q.rho);
    writeDataset("rhou", Q.rhou);
    writeDataset("rhov", Q.rhov);
    writeDataset("rhow", Q.rhow);
    writeDataset("E", Q.E);
    writeDataset("phi", Q.phi);

    // Store grid information as attributes at the root
    {
        // Define a group for grid metadata
        Group grid = file.createGroup("/grid");

        // Domain dimensions
        hsize_t adims_dims[1] = {3};
        DataSpace attr_space_dims(1, adims_dims);
        Attribute attr_dims = grid.createAttribute(
            "domain_dimensions", PredType::NATIVE_INT, attr_space_dims);
        int domain_dims[3] = {Q.Nx, Q.Ny, Q.Nz};
        attr_dims.write(PredType::NATIVE_INT, domain_dims);

        // Cell sizes
        hsize_t adims_size[1] = {3};
        DataSpace attr_space_size(1, adims_size);
        Attribute attr_size = grid.createAttribute(
            "cell_size", PredType::NATIVE_DOUBLE, attr_space_size);
        double cell_sizes[3] = {Q.dx, Q.dy, Q.dz};
        attr_size.write(PredType::NATIVE_DOUBLE, cell_sizes);

        // Bounding box (assuming origin at (0,0,0))
        hsize_t adims_bbox[2] = {3, 2};  // 3 dimensions, each with min and max
        DataSpace attr_space_bbox(2, adims_bbox);
        Attribute attr_bbox = grid.createAttribute(
            "bounding_box", PredType::NATIVE_DOUBLE, attr_space_bbox);
        double bbox[3][2] = {
            {0.0, Q.Nx * Q.dx}, {0.0, Q.Ny * Q.dy}, {0.0, Q.Nz * Q.dz}};
        attr_bbox.write(PredType::NATIVE_DOUBLE, bbox);
    }

    std::cout << "HDF5 file written to: " << filename << std::endl;
}

void readFieldHDF5(Field3D &Q, const std::string &filename) {
    using namespace H5;

    // Open the file in read-only mode
    H5File file(filename.c_str(), H5F_ACC_RDONLY);

    // Read grid metadata
    Group grid = file.openGroup("/grid");

    // Read domain dimensions
    Attribute attr_dims = grid.openAttribute("domain_dimensions");
    int domain_dims[3];
    attr_dims.read(PredType::NATIVE_INT, domain_dims);
    Q.Nx = domain_dims[0];
    Q.Ny = domain_dims[1];
    Q.Nz = domain_dims[2];

    // Read cell sizes
    Attribute attr_size = grid.openAttribute("cell_size");
    double cell_sizes[3];
    attr_size.read(PredType::NATIVE_DOUBLE, cell_sizes);
    Q.dx = cell_sizes[0];
    Q.dy = cell_sizes[1];
    Q.dz = cell_sizes[2];

    // Read bounding box (optional, for verification)
    Attribute attr_bbox = grid.openAttribute("bounding_box");
    double bbox[3][2];
    attr_bbox.read(PredType::NATIVE_DOUBLE, bbox);
    // You can use bbox if needed

    // Dimensions for the 3D datasets: (Nz, Ny, Nx)
    hsize_t dims[3];
    dims[0] = static_cast<hsize_t>(Q.Nz);
    dims[1] = static_cast<hsize_t>(Q.Ny);
    dims[2] = static_cast<hsize_t>(Q.Nx);

    // Create a data space for the 3D shape
    DataSpace dataspace(3, dims);

    // Helper lambda to read a 3D dataset into Q.* vectors
    auto readDataset = [&](const std::string &name, std::vector<double> &data) {
        // Open the dataset
        DataSet dataset = file.openDataSet(name);
        // Read the data into a buffer
        std::vector<double> buffer(Q.Nx * Q.Ny * Q.Nz);
        dataset.read(buffer.data(), PredType::NATIVE_DOUBLE);

        // Rearrange data from (k, j, i) to (i + Nx*(j + Ny*k))
        for (int k = 0; k < Q.Nz; ++k) {
            for (int j = 0; j < Q.Ny; ++j) {
                for (int i = 0; i < Q.Nx; ++i) {
                    size_t idxBuf = static_cast<size_t>(k) * (Q.Ny * Q.Nx) +
                                    static_cast<size_t>(j) * Q.Nx +
                                    static_cast<size_t>(i);
                    int idxQ = i + Q.Nx * (j + Q.Ny * k);
                    data[idxQ] = buffer[idxBuf];
                }
            }
        }
    };

    // Read each field
    readDataset("rho", Q.rho);
    readDataset("rhou", Q.rhou);
    readDataset("rhov", Q.rhov);
    readDataset("rhow", Q.rhow);
    readDataset("E", Q.E);
    readDataset("phi", Q.phi);

    std::cout << "HDF5 file read from: " << filename << std::endl;
}

}  // namespace io
}  // namespace agoge
