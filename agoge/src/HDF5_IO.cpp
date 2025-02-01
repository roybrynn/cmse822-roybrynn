#include "HDF5_IO.hpp"

#include <H5Cpp.h>

#include <iostream>
#include <vector>

#include "Field3d.hpp"

namespace agoge {
namespace io {

void writeFieldHDF5(const Field3D &Q, const std::string &filename) {
    using namespace H5;
    try {
        // Create (or truncate) the file
        H5File file(filename.c_str(), H5F_ACC_TRUNC);

        // Prepare a lambda to rearrange and write a 3D dataset.
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

        // Create a group for grid metadata and spatial coordinates.
        Group grid = file.createGroup("/grid");

        // Store domain dimensions as an attribute.
        hsize_t adims_dims[1] = {3};
        DataSpace attr_space_dims(1, adims_dims);
        int domain_dims[3] = {Q.Nx, Q.Ny, Q.Nz};
        Attribute attr_dims = grid.createAttribute(
            "domain_dimensions", PredType::NATIVE_INT, attr_space_dims);
        attr_dims.write(PredType::NATIVE_INT, domain_dims);

        // Store cell sizes as an attribute.
        hsize_t adims_size[1] = {3};
        DataSpace attr_space_size(1, adims_size);
        double cell_sizes[3] = {Q.dx, Q.dy, Q.dz};
        Attribute attr_size = grid.createAttribute(
            "cell_size", PredType::NATIVE_DOUBLE, attr_space_size);
        attr_size.write(PredType::NATIVE_DOUBLE, cell_sizes);

        // Store bounding box as an attribute.
        hsize_t bbox_dims[1] = {6};
        DataSpace bbox_space(1, bbox_dims);
        const BoundingBox &bbox = Q.getBoundingBox();
        double bbox_data[6] = {bbox.xmin, bbox.xmax, bbox.ymin,
                               bbox.ymax, bbox.zmin, bbox.zmax};
        Attribute attr_bbox = grid.createAttribute(
            "bounding_box", PredType::NATIVE_DOUBLE, bbox_space);
        attr_bbox.write(PredType::NATIVE_DOUBLE, bbox_data);

        // Write coordinate arrays computed from the built-in methods.
        // X coordinates
        std::vector<double> x_coords(Q.Nx);
        for (int i = 0; i < Q.Nx; ++i) {
            // Using the convention: center coordinate = bbox.xmin + (i+0.5)*dx
            x_coords[i] = bbox.xmin + (i + 0.5) * Q.dx;
        }
        hsize_t xdim[1] = {static_cast<hsize_t>(Q.Nx)};
        DataSpace x_space(1, xdim);
        DataSet x_dataset =
            grid.createDataSet("x", PredType::NATIVE_DOUBLE, x_space);
        x_dataset.write(x_coords.data(), PredType::NATIVE_DOUBLE);

        // Y coordinates
        std::vector<double> y_coords(Q.Ny);
        for (int j = 0; j < Q.Ny; ++j) {
            y_coords[j] = bbox.ymin + (j + 0.5) * Q.dy;
        }
        hsize_t ydim[1] = {static_cast<hsize_t>(Q.Ny)};
        DataSpace y_space(1, ydim);
        DataSet y_dataset =
            grid.createDataSet("y", PredType::NATIVE_DOUBLE, y_space);
        y_dataset.write(y_coords.data(), PredType::NATIVE_DOUBLE);

        // Z coordinates
        std::vector<double> z_coords(Q.Nz);
        for (int k = 0; k < Q.Nz; ++k) {
            z_coords[k] = bbox.zmin + (k + 0.5) * Q.dz;
        }
        hsize_t zdim[1] = {static_cast<hsize_t>(Q.Nz)};
        DataSpace z_space(1, zdim);
        DataSet z_dataset =
            grid.createDataSet("z", PredType::NATIVE_DOUBLE, z_space);
        z_dataset.write(z_coords.data(), PredType::NATIVE_DOUBLE);

        std::cout << "HDF5 file written to: " << filename << std::endl;
    } catch (const H5::Exception &err) {
        std::cerr << "HDF5 write error: " << err.getDetailMsg() << std::endl;
        throw;
    }
}

void readFieldHDF5(Field3D &Q, const std::string &filename) {
    using namespace H5;
    try {
        H5File file(filename.c_str(), H5F_ACC_RDONLY);

        // Read grid metadata, if desired.
        Group grid = file.openGroup("/grid");
        // (Reading of attributes and coordinate datasets can be performed here)
        // For brevity, detailed error checking and data rearrangement for the
        // fields is omitted.
        std::cout << "HDF5 file read from: " << filename << std::endl;
    } catch (const H5::Exception &err) {
        std::cerr << "HDF5 read error: " << err.getDetailMsg() << std::endl;
        throw;
    }
}

}  // namespace io
}  // namespace agoge