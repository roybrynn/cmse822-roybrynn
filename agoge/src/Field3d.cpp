#include "Field3D.hpp"

namespace agoge {

Field3D::Field3D(int nx, int ny, int nz, double dx_, double dy_, double dz_)
    : Nx(nx), Ny(ny), Nz(nz),
      dx(dx_), dy(dy_), dz(dz_)
{
    // Compute total number of cells
    const size_t N = static_cast<size_t>(Nx) * Ny * Nz;

    // Allocate each field
    rho .resize(N);
    rhou.resize(N);
    rhov.resize(N);
    rhow.resize(N);
    E   .resize(N);
    phi .resize(N);
}

} // namespace agoge
