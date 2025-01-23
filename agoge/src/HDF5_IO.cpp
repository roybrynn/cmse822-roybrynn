#include "HDF5_IO.hpp"
#include "Field3d.hpp"

#include <H5Cpp.h>
#include <iostream>
#include <vector>

namespace agoge {
namespace io {

void writeFieldHDF5(const Field3D &Q, const std::string &filename)
{
    using namespace H5;

    // Create / open the file (truncate if exists)
    H5File file(filename.c_str(), H5F_ACC_TRUNC);

    // Dimensions for the 3D datasets: (Nz, Ny, Nx)
    // According to the problem statement, shape is (Nz, Ny, Nx)
    hsize_t dims[3];
    dims[0] = static_cast<hsize_t>(Q.Nz);
    dims[1] = static_cast<hsize_t>(Q.Ny);
    dims[2] = static_cast<hsize_t>(Q.Nx);

    // Create a data space for the 3D shape
    DataSpace dataspace(3, dims);

    // Helper lambda to write a 3D dataset from Q.* vectors
    auto writeDataset = [&](const std::string &name,
                            const std::vector<double> &data)
    {
        // Create the dataset
        DataSet dataset = file.createDataSet(
            name, PredType::NATIVE_DOUBLE, dataspace);

        // We need to transform the data from (i + Nx*(j + Ny*k)) indexing
        // to the shape (k, j, i) = (Nz, Ny, Nx) in row-major order.
        // We'll allocate a temporary buffer with that reordering.

        std::vector<double> buffer(data.size());
        // For each (k, j, i), copy data from the 1D index in Q
        // to the (k*Ny*Nx + j*Nx + i) offset in buffer.

        size_t idxBuf = 0;
        for (int k = 0; k < Q.Nz; ++k) {
            for (int j = 0; j < Q.Ny; ++j) {
                for (int i = 0; i < Q.Nx; ++i) {
                    // 1D index in Q.* is i + Nx*(j + Ny*k)
                    int idxQ = i + Q.Nx * (j + Q.Ny * k);
                    idxBuf = static_cast<size_t>(k) * (Q.Ny * Q.Nx)
                           + static_cast<size_t>(j) * Q.Nx
                           + static_cast<size_t>(i);
                    buffer[idxBuf] = data[idxQ];
                }
            }
        }

        // Write to the dataset
        dataset.write(buffer.data(), PredType::NATIVE_DOUBLE);
    };

    // Write each field
    writeDataset("rho",  Q.rho );
    writeDataset("rhou", Q.rhou);
    writeDataset("rhov", Q.rhov);
    writeDataset("rhow", Q.rhow);
    writeDataset("E",    Q.E   );
    writeDataset("phi",  Q.phi );

    // Optionally set attributes for domain sizes
    {
        // We can store Nx, Ny, Nz, dx, dy, dz as attributes on the root group.
        // Or you could store them on each dataset. We'll store on root for convenience.
        Group root = file.openGroup("/");
        // 1D attribute with 3 elements for Nx, Ny, Nz
        {
            hsize_t adims[1] = {3};
            DataSpace aspace(1, adims);
            int domain_dims[3] = {Q.Nx, Q.Ny, Q.Nz};
            Attribute attrDD = root.createAttribute(
                "domain_dimensions", PredType::NATIVE_INT, aspace);
            attrDD.write(PredType::NATIVE_INT, domain_dims);
        }
        // Another 1D attribute for dx, dy, dz
        {
            hsize_t adims[1] = {3};
            DataSpace aspace(1, adims);
            double domain_dx[3] = {Q.dx, Q.dy, Q.dz};
            Attribute attrDX = root.createAttribute(
                "cell_size", PredType::NATIVE_DOUBLE, aspace);
            attrDX.write(PredType::NATIVE_DOUBLE, domain_dx);
        }
    }

    std::cout << "HDF5 file written to: " << filename << std::endl;
}

void readFieldHDF5(Field3D &Q, const std::string &filename)
{
    using namespace H5;

    // Open the file in read-only mode
    H5File file(filename.c_str(), H5F_ACC_RDONLY);

    // We assume the shape is stored as (Nz, Ny, Nx). Let's read one dataset
    // to infer the shape, or assume Q.Nx, Q.Ny, Q.Nz are already set externally.

    // Helper lambda for reading one 3D dataset
    auto readDataset = [&](const std::string &name,
                           std::vector<double> &data)
    {
        DataSet dataset = file.openDataSet(name);
        DataSpace filespace = dataset.getSpace();

        // Check that the dataset shape matches (Nz, Ny, Nx)
        const int rank = filespace.getSimpleExtentNdims();
        if (rank != 3) {
            throw std::runtime_error("Dataset rank != 3 in " + name);
        }
        hsize_t dims[3];
        filespace.getSimpleExtentDims(dims, nullptr);
        if (static_cast<int>(dims[0]) != Q.Nz ||
            static_cast<int>(dims[1]) != Q.Ny ||
            static_cast<int>(dims[2]) != Q.Nx) {
            throw std::runtime_error("Dataset shape mismatch in " + name);
        }

        std::vector<double> buffer(Q.Nx * Q.Ny * Q.Nz);
        dataset.read(buffer.data(), PredType::NATIVE_DOUBLE);

        // Now map buffer (k, j, i) -> data[i + Nx*(j + Ny*k)]
        for (int k = 0; k < Q.Nz; ++k) {
            for (int j = 0; j < Q.Ny; ++j) {
                for (int i = 0; i < Q.Nx; ++i) {
                    size_t idxBuf = static_cast<size_t>(k) * (Q.Ny * Q.Nx)
                                  + static_cast<size_t>(j) * Q.Nx
                                  + static_cast<size_t>(i);
                    int idxQ = i + Q.Nx * (j + Q.Ny * k);
                    data[idxQ] = buffer[idxBuf];
                }
            }
        }
    };

    // Read each field
    readDataset("rho",  Q.rho );
    readDataset("rhou", Q.rhou);
    readDataset("rhov", Q.rhov);
    readDataset("rhow", Q.rhow);
    readDataset("E",    Q.E   );
    readDataset("phi",  Q.phi );

    // Optionally read attributes for domain sizes
    {
        Group root = file.openGroup("/");
        // domain_dimensions
        if (root.attrExists("domain_dimensions")) {
            Attribute attrDD = root.openAttribute("domain_dimensions");
            int domain_dims[3];
            attrDD.read(PredType::NATIVE_INT, domain_dims);
            if (domain_dims[0] != Q.Nx ||
                domain_dims[1] != Q.Ny ||
                domain_dims[2] != Q.Nz) {
                std::cerr << "Warning: Mismatch between Q.Nx,Ny,Nz and HDF5 domain_dimensions.\n";
            }
        }
        // cell_size
        if (root.attrExists("cell_size")) {
            Attribute attrDX = root.openAttribute("cell_size");
            double domain_dx[3];
            attrDX.read(PredType::NATIVE_DOUBLE, domain_dx);
            // Check or store them if desired
            if (domain_dx[0] != Q.dx ||
                domain_dx[1] != Q.dy ||
                domain_dx[2] != Q.dz) {
                std::cerr << "Warning: Mismatch between Q.dx,dy,dz and HDF5 cell_size.\n";
            }
        }
    }

    std::cout << "HDF5 file read from: " << filename << std::endl;
}

} // namespace io
} // namespace agoge
