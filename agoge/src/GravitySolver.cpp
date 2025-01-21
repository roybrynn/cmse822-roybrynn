#include "GravitySolver.hpp"
#include "Field3D.hpp"

#include <fftw3.h>
#include <cmath>
#include <vector>
#include <iostream>

namespace agoge {
namespace gravity {

void solvePoissonFFT(Field3D &Q)
{
    // Retrieve grid sizes and cell sizes
    const int Nx = Q.Nx;
    const int Ny = Q.Ny;
    const int Nz = Q.Nz;

    const double dx = Q.dx;
    const double dy = Q.dy;
    const double dz = Q.dz;

    // We'll assume G is 1.0 for demonstration, or you could fetch from config.
    // For dimensionless or code units, G can be set accordingly.
    constexpr double G = 1.0;

    // For convenience, define total physical lengths
    double Lx = Nx * dx;
    double Ly = Ny * dy;
    double Lz = Nz * dz;

    // 1. Allocate a real buffer for rho in real space
    //    and a complex buffer for the transform.
    // For an R2C 3D transform, the size of the complex array is Nx*(Ny)*(Nz/2+1) 
    // in row-major (z fastest) layout if we call the usual fftw_plan_dft_r2c_3d(Nz, Ny, Nx...).
    // But let's keep it consistent with an (z, y, x) transform order.
    // For demonstration, we'll do the standard FFTW approach:
    // forwardPlan = fftw_plan_dft_r2c_3d(Nz, Ny, Nx, real_in, complex_out, FFTW_ESTIMATE)
    // That means 'real_in' is of size Nx*Ny*Nz, 'complex_out' is Nx*Ny*(Nz/2+1).

    std::vector<double> realBuf(Nx * Ny * Nz, 0.0);
    fftw_complex* cplxBuf = fftw_alloc_complex(static_cast<size_t>(Nx) * Ny * (Nz/2 + 1));

    // 2. Copy Q.rho into the real buffer
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const int idx = Q.index(i, j, k);
                realBuf[idx] = Q.rho[idx];
            }
        }
    }

    // 3. Create the forward and backward FFT plans
    fftw_plan forwardPlan = fftw_plan_dft_r2c_3d(
        Nz, Ny, Nx,
        realBuf.data(),
        cplxBuf,
        FFTW_ESTIMATE
    );

    fftw_plan backwardPlan = fftw_plan_dft_c2r_3d(
        Nz, Ny, Nx,
        cplxBuf,
        realBuf.data(),
        FFTW_ESTIMATE
    );

    // 4. Execute forward transform (real -> complex)
    fftw_execute(forwardPlan);

    // 5. Multiply by the factor -4 pi G / |k|^2 in wave space, handle k=0 separately
    const double TWO_PI = 2.0 * M_PI;
    // We interpret wavenumbers as (kx, ky, kz):
    //   kx in [-Nx/2, Nx/2-1], ky in [-Ny/2, Ny/2-1], kz in [-Nz/2, Nz/2-1]
    // but FFTW indexes them in a certain layout. We'll do it in a standard way for R2C:
    //   For each kz in [0..Nz/2], ky in [0..Ny-1], kx in [0..Nx-1], but physically we adjust "aliases".

    auto idxCplx = [&](int kz, int ky, int kx) {
        // In R2C 3D layout, the kx index runs 0..(Nx), but physically only 0..(Nx/2) unique if Nx is even.
        // The total slice is Nx*(Ny)*(Nz/2+1)? Actually: (Nx*(Ny)*(Nz/2+1)) or a similar layout.
        return (kz * (Ny * (Nx/2 + 1))) + (ky * (Nx/2 + 1)) + kx;
    };

    for (int kz = 0; kz <= Nz/2; ++kz) {
        // k_z in 'alias' form
        int kzAlias = (kz <= Nz/2) ? kz : (kz - Nz);
        double kzVal = (kzAlias) * (TWO_PI / Lz);

        for (int ky = 0; ky < Ny; ++ky) {
            // alias for ky
            int kyAlias = (ky <= Ny/2) ? ky : (ky - Ny);
            double kyVal = (kyAlias) * (TWO_PI / Ly);

            for (int kx = 0; kx <= Nx/2; ++kx) {
                // alias for kx
                int kxAlias = (kx <= Nx/2) ? kx : (kx - Nx);
                double kxVal = (kxAlias) * (TWO_PI / Lx);

                const int idxK = idxCplx(kz, ky, kx);
                double k2 = kxVal*kxVal + kyVal*kyVal + kzVal*kzVal;

                if (k2 == 0.0) {
                    // k=0 mode: set to zero to avoid division by zero
                    cplxBuf[idxK][0] = 0.0;
                    cplxBuf[idxK][1] = 0.0;
                } else {
                    double factor = -4.0 * M_PI * G / k2;
                    double re = cplxBuf[idxK][0];
                    double im = cplxBuf[idxK][1];
                    // multiply
                    cplxBuf[idxK][0] = factor * re;
                    cplxBuf[idxK][1] = factor * im;
                }
            }
        }
    }

    // 6. Inverse FFT (complex -> real)
    fftw_execute(backwardPlan);

    // 7. Normalize, since FFTW doesn't automatically scale the inverse
    double invSize = 1.0 / static_cast<double>(Nx * Ny * Nz);
    for (int idx = 0; idx < Nx*Ny*Nz; ++idx) {
        realBuf[idx] *= invSize;
    }

    // 8. Copy back to Q.phi
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const int idx = Q.index(i, j, k);
                Q.phi[idx] = realBuf[idx];
            }
        }
    }

    // 9. Clean up
    fftw_destroy_plan(forwardPlan);
    fftw_destroy_plan(backwardPlan);
    fftw_free(cplxBuf);
}

} // namespace gravity
} // namespace agoge
