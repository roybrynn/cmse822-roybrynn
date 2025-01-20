//------------------------------------------------------------
// A simple Structure-of-Arrays container
//------------------------------------------------------------
struct Field3D
{
    // Each variable is a separate array of size Nx*Ny*Nz
    std::vector<double> rho;   // density
    std::vector<double> rhou;  // momentum in x
    std::vector<double> rhov;  // momentum in y
    std::vector<double> rhow;  // momentum in z
    std::vector<double> E;     // total energy

    // Grid dimensions
    int Nx_, Ny_, Nz_;
    double dx_, dy_, dz_;

    // Constructor allocates memory
    Field3D(int Nx, int Ny, int Nz, double dx, double dy, double dz)
        : Nx_(Nx), Ny_(Ny), Nz_(Nz),
          dx_(dx), dy_(dy), dz_(dz)
    {
        rho .resize(Nx_ * Ny_ * Nz_);
        rhou.resize(Nx_ * Ny_ * Nz_);
        rhov.resize(Nx_ * Ny_ * Nz_);
        rhow.resize(Nx_ * Ny_ * Nz_);
        E   .resize(Nx_ * Ny_ * Nz_);
    }

    // Inline helper for indexing
    inline int index(int i, int j, int k) const
    {
        return i + Nx_ * (j + Ny_ * k);
    }
};

