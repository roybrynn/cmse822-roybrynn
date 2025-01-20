void writeFieldHDF5(const Field3D& Q, const std::string& filename)
{
    // Create HDF5 file
    H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);

    // Create a group that will store geometry info
    H5::Group rootGroup = file.createGroup("/Grid000000");

    // domain_dimensions attribute
    hsize_t dimsA[1] = {3};
    H5::DataSpace attrSpace(1, dimsA);
    int domain_dims[3] = {Q.Nz_, Q.Ny_, Q.Nx_};
    H5::Attribute attrDD = rootGroup.createAttribute("domain_dimensions", 
                             H5::PredType::NATIVE_INT, attrSpace);
    attrDD.write(H5::PredType::NATIVE_INT, domain_dims);

    // left_edge, right_edge
    double left[3]  = {0.0, 0.0, 0.0};
    double right[3] = {Q.dz_ * Q.Nz_, Q.dy_ * Q.Ny_, Q.dx_ * Q.Nx_};
    // or reorder if you interpret Nx -> x-dim, etc.

    H5::Attribute attrLE = rootGroup.createAttribute("left_edge", 
                             H5::PredType::NATIVE_DOUBLE, attrSpace);
    attrLE.write(H5::PredType::NATIVE_DOUBLE, left);

    H5::Attribute attrRE = rootGroup.createAttribute("right_edge", 
                             H5::PredType::NATIVE_DOUBLE, attrSpace);
    attrRE.write(H5::PredType::NATIVE_DOUBLE, right);

    // Create a 3D dataset under /Grid000000 for each field
    hsize_t dims[3] = {
        static_cast<hsize_t>(Q.Nz_),
        static_cast<hsize_t>(Q.Ny_),
        static_cast<hsize_t>(Q.Nx_)
    };
    H5::DataSpace dataspace(3, dims);

    // Helper to write a dataset in /Grid000000
    auto writeDataset = [&](const std::string &dsName, const std::vector<double>& data)
    {
        H5::DataSet dataset = rootGroup.createDataSet(
            dsName, 
            H5::PredType::NATIVE_DOUBLE, 
            dataspace
        );
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    };

    writeDataset("density", Q.rho);
    writeDataset("rhou",    Q.rhou);
    writeDataset("rhov",    Q.rhov);
    writeDataset("rhow",    Q.rhow);
    writeDataset("energy",  Q.E);
}

//------------------------------------------------------------
// A helper function to write a 3D field to HDF5
//------------------------------------------------------------
void writeFieldHDF5(const Field3D& Q, const std::string& filename)
{
    // Create HDF5 file (overwrite if exists)
    H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);

    // We will store each variable as a 3D dataset
    hsize_t dims[3] = { static_cast<hsize_t>(Q.Nz_), 
                        static_cast<hsize_t>(Q.Ny_), 
                        static_cast<hsize_t>(Q.Nx_) };

    // Create the data space for 3D
    H5::DataSpace dataspace(3, dims);

    // Helper lambda to write one dataset
    auto writeDataset = [&](const std::string &name, const std::vector<double> &data)
    {
        // Create dataset
        H5::DataSet dataset = file.createDataSet(
            name, 
            H5::PredType::NATIVE_DOUBLE, 
            dataspace);

        // We have data in [i + Nx*(j + Ny*k)] order, 
        // but HDF5 is going to interpret in (k, j, i) order, 
        // so we can just write it directly if we're consistent.
        // If you prefer to reorder, do it here. For now, we just write directly.
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    };

    writeDataset("rho",  Q.rho);
    writeDataset("rhou", Q.rhou);
    writeDataset("rhov", Q.rhov);
    writeDataset("rhow", Q.rhow);
    writeDataset("E",    Q.E);
}

