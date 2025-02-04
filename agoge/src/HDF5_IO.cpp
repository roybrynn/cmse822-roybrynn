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
        // In function writeFieldHDF5, after creating the grid group:
        Group grid = file.createGroup("/grid");

        // Store bounding box as an attribute.
        hsize_t bbox_dims[1] = {6};
        DataSpace bbox_space(1, bbox_dims);
        const auto &bbox = Q.getBoundingBox();
        double bbox_data[6] = {bbox.xmin, bbox.xmax, bbox.ymin,
                                bbox.ymax, bbox.zmin, bbox.zmax};
        Attribute attr_bbox = grid.createAttribute(
            "bounding_box", PredType::NATIVE_DOUBLE, bbox_space);
        attr_bbox.write(PredType::NATIVE_DOUBLE, bbox_data);

        // Existing: Store domain dimensions and cell sizes.
        hsize_t dims_attr[1] = {3};
        DataSpace dims_space(1, dims_attr);
        int domain_dims[3] = {Q.Nx, Q.Ny, Q.Nz};
        Attribute attr_dims = grid.createAttribute(
            "domain_dimensions", PredType::NATIVE_INT, dims_space);
        attr_dims.write(PredType::NATIVE_INT, domain_dims);

        double cell_sizes[3] = {Q.dx, Q.dy, Q.dz};
        Attribute attr_size = grid.createAttribute(
            "cell_size", PredType::NATIVE_DOUBLE, dims_space);
        attr_size.write(PredType::NATIVE_DOUBLE, cell_sizes);

        std::cout << "HDF5 file written to: " << filename << std::endl;
    } catch (const H5::Exception &err) {
        std::cerr << "HDF5 write error: " << err.getDetailMsg() << std::endl;
        throw;
    }
}

}  // namespace io
}  // namespace agoge