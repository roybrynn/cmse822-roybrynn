#pragma once
#include <vector>
#include "Domain.hpp"

namespace agoge {

struct Field3D {
    // Arrays for rho, rhou, rhov, rhow, E, potential phi, etc.
    std::vector<double> rho, rhou, rhov, rhow, E, phi;

    // Store domain metadata (or keep separate)
    Domain domain;

    Field3D(const Domain &dom)
      : domain(dom)
    {
        size_t N = static_cast<size_t>(dom.Nx*dom.Ny*dom.Nz);
        rho .resize(N);
        rhou.resize(N);
        rhov.resize(N);
        rhow.resize(N);
        E   .resize(N);
        phi .resize(N, 0.0); // for gravity, if needed
    }

    // Helper to index a 3D coordinate
    inline int index(int i, int j, int k) const {
        return i + domain.Nx*(j + domain.Ny*k);
    }
};

} // namespace agoge
