// src/Field3d.cpp
#include "Field3d.hpp"

namespace agoge {

Field3D::Field3D(int nx, int ny, int nz, double dx, double dy, double dz,
                 agoge::ParameterSystem &params)
    : Nx(nx),
      Ny(ny),
      Nz(nz),
      dx(dx),
      dy(dy),
      dz(dz),
      rho(nx * ny * nz, 0.0),
      rhou(nx * ny * nz, 0.0),
      rhov(nx * ny * nz, 0.0),
      rhow(nx * ny * nz, 0.0),
      E(nx * ny * nz, 0.0),
      phi(nx * ny * nz, 0.0),
      params_(params) {}

}  // namespace agoge
