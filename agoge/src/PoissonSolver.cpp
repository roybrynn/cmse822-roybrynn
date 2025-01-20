#include <fftw3.h>
#include <cmath>

// We'll assume Nx, Ny, Nz, Lx, Ly, Lz, etc. are global or in a struct.
// Also assume we have Q.rho storing density in real space.

void solvePoissonFFT(const Field3D &Q, std::vector<double> &phi)
{
    // 1) Allocate arrays for FFTW. We'll do a real-to-complex forward transform,
    //    so we need a complex array big enough to hold Nx*Ny*(Nz/2+1) if we're
    //    using an in-place R2C transform. The details vary by library & dimension.

    static fftw_complex *rhoK = nullptr; 
    static fftw_plan forwardPlan, backwardPlan; 
    static bool firstCall = true;

    if(firstCall) {
        // Allocate memory
        rhoK = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny*(Nz/2+1));

        // Create FFT plans (in-place or out-of-place). We'll assume out-of-place for clarity.
        // input array: a double* of size Nx*Ny*Nz
        // output array: an fftw_complex* of size Nx*Ny*(Nz/2+1)

        forwardPlan  = fftw_plan_dft_r2c_3d(Nz, Ny, Nx, /* real_in = */ nullptr, 
                                            /* complex_out = */ rhoK,
                                            FFTW_ESTIMATE);
        backwardPlan = fftw_plan_dft_c2r_3d(Nz, Ny, Nx, /* complex_in = */ rhoK,
                                            /* real_out = */ nullptr,
                                            FFTW_ESTIMATE);

        firstCall = false;
    }

    // We need a contiguous real array for the density, let's create a temp:
    static std::vector<double> realBuf(Nx*Ny*Nz, 0.0);
    // Copy Q.rho into realBuf
    for(int k=0; k<Nz; k++){
      for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
          int idx = Q.index(i,j,k);
          realBuf[idx] = Q.rho[idx];
        }
      }
    }

    // 2) Forward FFT of rho -> rhoK
    // But we must point the plan's real input to realBuf.data(), so let's do:
    fftw_execute_dft_r2c(forwardPlan, realBuf.data(), rhoK);

    // 3) Multiply in k-space by the factor: -4 pi G / k^2  (handle k=0 separately)
    const double G = 1.0;  // or pass in as parameter
    const double TWO_PI_X = 2.0*M_PI / Lx;
    const double TWO_PI_Y = 2.0*M_PI / Ly;
    const double TWO_PI_Z = 2.0*M_PI / Lz;

    for(int kz=0; kz<Nz; kz++){
      // correct “folded” index for FFT might be: if kz>Nz/2 => (kz - Nz)
      int kzalias = (kz <= Nz/2) ? kz : (kz - Nz);
      double kzVal = kzalias * TWO_PI_Z;

      for(int jy=0; jy<Ny; jy++){
        int jyalias = (jy <= Ny/2) ? jy : (jy - Ny);
        double kyVal = jyalias * TWO_PI_Y;

        for(int ix=0; ix<= Nx/2; ix++){
          // note: in r2c, x only goes up to Nx/2
          int ixalias = (ix <= Nx/2) ? ix : (ix - Nx);
          double kxVal = ixalias * TWO_PI_X;

          int idxK = (kz*(Ny*(Nx/2+1))) + (jy*(Nx/2+1)) + ix; // typical layout
          double k2 = kxVal*kxVal + kyVal*kyVal + kzVal*kzVal;

          if(k2 == 0.0){
            // k=0 mode => set to 0 or handle separately
            rhoK[idxK][0] = 0.0;
            rhoK[idxK][1] = 0.0;
          } else {
            double factor = -4.0*M_PI*G / k2;
            // multiply real & imaginary parts
            double realPart = rhoK[idxK][0];
            double imagPart = rhoK[idxK][1];
            rhoK[idxK][0] = factor * realPart;
            rhoK[idxK][1] = factor * imagPart;
          }
        }
      }
    }

    // 4) Inverse FFT -> phi array
    // We'll write results into realBuf, then copy to phi.
    fftw_execute_dft_c2r(backwardPlan, rhoK, realBuf.data());

    // 5) Normalize since FFTW does not do it automatically
    double norm = 1.0/( (double)(Nx*Ny*Nz) );
    for(int idx=0; idx<Nx*Ny*Nz; idx++){
      realBuf[idx] *= norm;
    }

    // 6) Copy result into phi array
    for(int k=0; k<Nz; k++){
      for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
          int idx = Q.index(i,j,k);
          phi[idx] = realBuf[idx];
        }
      }
    }
}
