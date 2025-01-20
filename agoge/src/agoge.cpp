#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// HDF5 C++ API
#include "H5Cpp.h"

static const H5::H5FileAccPropList fapl = H5::FileAccPropList::DEFAULT;
static const H5::DataSetCreatPropList ds_creatplist = H5::DataSetCreatPropList::DEFAULT;
static const H5::DataSpace memspace = H5::DataSpace();

//------------------------------------------------------------
// PHYSICAL CONSTANTS AND PARAMETERS
//------------------------------------------------------------
static const double gamma_gas = 1.4;
static const double CFL       = 0.5;  // Courant-Friedrichs-Lewy number

// Grid sizes (small for demonstration)
static const int Nx = 64;
static const int Ny = 64;
static const int Nz = 64;

// Domain length in each direction
static const double Lx = 1.0;
static const double Ly = 1.0;
static const double Lz = 1.0;

// Derived spacing
static const double dx = Lx / Nx;
static const double dy = Ly / Ny;
static const double dz = Lz / Nz;

