#include "GravitySolver.hpp"
#include "Field3D.hpp"

#include <fftw3.h>
#include <cmath>
#include <vector>
#include <iostream>

// NEW
#include "PerformanceMonitor.hpp"

namespace agoge {
namespace gravity {

void solvePoissonFFT(Field3D &Q)
{
    // Start timing
    PerformanceMonitor::instance().startTimer("solvePoissonFFT");

    const int Nx = Q.Nx;
    const int Ny = Q.Ny;
    const int Nz = Q.Nz;
    constexpr double G = 1.0;

    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;
    double Lx = Nx * dx;
    double Ly = Ny * dy;
    double Lz = Nz * dz;

    std::vector<double> realBuf(Nx * Ny * Nz, 0.0);
    fftw_complex* cplxBuf = fftw_alloc_complex(static_cast<size_t>(Nx)*Ny*(Nz/2 + 1));

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = Q.index(i,j,k);
                realBuf[idx] = Q.rho[idx];
            }
        }
    }

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

    fftw_execute(forwardPlan);

    const double TWO_PI = 2.0 * M_PI;
    auto idxCplx = [&](int kz, int ky, int kx) {
        return kz*(Ny*(Nx/2+1)) + ky*(Nx/2+1) + kx;
    };

    for(int kz = 0; kz <= Nz/2; kz++){
        int kzAlias = (kz <= Nz/2) ? kz : (kz - Nz);
        double kzVal = kzAlias * (TWO_PI / Lz);

        for(int ky = 0; ky < Ny; ky++){
            int kyAlias = (ky <= Ny/2) ? ky : (ky - Ny);
            double kyVal = kyAlias * (TWO_PI / Ly);

            for(int kx = 0; kx <= Nx/2; kx++){
                int kxAlias = (kx <= Nx/2) ? kx : (kx - Nx);
                double kxVal = kxAlias * (TWO_PI / Lx);

                double k2 = kxVal*kxVal + kyVal*kyVal + kzVal*kzVal;
                int idxK = idxCplx(kz, ky, kx);

                if(k2 == 0.0) {
                    cplxBuf[idxK][0] = 0.0;
                    cplxBuf[idxK][1] = 0.0;
                } else {
                    double factor = -4.0 * M_PI * G / k2;
                    double re = cplxBuf[idxK][0];
                    double im = cplxBuf[idxK][1];
                    cplxBuf[idxK][0] = factor * re;
                    cplxBuf[idxK][1] = factor * im;
                }
            }
        }
    }

    fftw_execute(backwardPlan);

    double invSize = 1.0 / double(Nx*Ny*Nz);
    for(int idx = 0; idx < Nx*Ny*Nz; idx++){
        realBuf[idx] *= invSize;
    }

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = Q.index(i,j,k);
                Q.phi[idx] = realBuf[idx];
            }
        }
    }

    fftw_destroy_plan(forwardPlan);
    fftw_destroy_plan(backwardPlan);
    fftw_free(cplxBuf);

    // Stop timing
    PerformanceMonitor::instance().stopTimer("solvePoissonFFT");
}

} // namespace gravity
} // namespace agoge
