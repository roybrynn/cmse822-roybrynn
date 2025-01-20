#pragma once

#include <cstddef>

namespace agoge {

// Holds info about domain size, BCs, etc.
struct Domain {
    int Nx, Ny, Nz;
    double Lx, Ly, Lz;
    // possible boundary condition enums, e.g. BCType bcXmin, bcXmax, ...

    Domain(int nx, int ny, int nz, double lx, double ly, double lz)
      : Nx(nx), Ny(ny), Nz(nz), Lx(lx), Ly(ly), Lz(lz)
    {}
};

} // namespace agoge
