#include "GravitySolver.hpp"
#include "Field3d.hpp"
#include "Config.hpp"

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>

namespace agoge {
namespace gravity {

/* ------------------------------------------------------------------------
   1. Helper: Check if a dimension is a power of two
   ------------------------------------------------------------------------ */
static bool isPowerOfTwo(int n)
{
    return (n > 0) && ((n & (n - 1)) == 0);
}

/* ------------------------------------------------------------------------
   2. Naive DFT forward + inverse
   (O(N^6)) approach
   ------------------------------------------------------------------------ */
static void naiveForwardDFT3D(const std::vector<double> &realIn,
                              std::vector<std::complex<double>> &cplxOut,
                              int Nx, int Ny, int Nz)
{
    cplxOut.assign(Nx*Ny*Nz, {0.0, 0.0});
    const double TWO_PI = 2.0 * M_PI;

    for(int kx = 0; kx < Nx; kx++){
      for(int ky = 0; ky < Ny; ky++){
        for(int kz = 0; kz < Nz; kz++){
          int idxK = kx + Nx*(ky + Ny*kz);
          std::complex<double> sumVal(0.0, 0.0);

          for(int x=0; x<Nx; x++){
            double alphaX = (double)(kx*x)/(double)Nx;
            for(int y=0; y<Ny; y++){
              double alphaY = (double)(ky*y)/(double)Ny;
              for(int z=0; z<Nz; z++){
                double alphaZ = (double)(kz*z)/(double)Nz;
                double phase = -TWO_PI*(alphaX + alphaY + alphaZ);

                int idxXYZ = x + Nx*(y + Ny*z);
                double val = realIn[idxXYZ];
                sumVal += val*std::complex<double>(std::cos(phase), std::sin(phase));
              }
            }
          }
          cplxOut[idxK] = sumVal;
        }
      }
    }
}

static void naiveInverseDFT3D(const std::vector<std::complex<double>> &cplxIn,
                              std::vector<double> &realOut,
                              int Nx, int Ny, int Nz)
{
    realOut.assign(Nx*Ny*Nz, 0.0);
    const double TWO_PI = 2.0*M_PI;
    double invSize = 1.0/double(Nx*Ny*Nz);

    for(int x=0; x<Nx; x++){
      for(int y=0; y<Ny; y++){
        for(int z=0; z<Nz; z++){
          int idxXYZ = x + Nx*(y + Ny*z);
          std::complex<double> sumVal(0.0, 0.0);

          for(int kx=0; kx<Nx; kx++){
            double alphaX = (double)(kx*x)/(double)Nx;
            for(int ky=0; ky<Ny; ky++){
              double alphaY = (double)(ky*y)/(double)Ny;
              for(int kz=0; kz<Nz; kz++){
                double alphaZ = (double)(kz*z)/(double)Nz;
                double phase = +TWO_PI*(alphaX + alphaY + alphaZ);

                int idxK = kx + Nx*(ky + Ny*kz);
                sumVal += cplxIn[idxK]*std::complex<double>(std::cos(phase), std::sin(phase));
              }
            }
          }
          realOut[idxXYZ] = sumVal.real()*invSize;
        }
      }
    }
}

/* ------------------------------------------------------------------------
   3. A basic Cooley–Tukey 1D in-place FFT (radix-2).
      We use complex<double> in a vector, and do separate code
      for forward/inverse.
   ------------------------------------------------------------------------ */
static void fft1D(std::complex<double> *data, int N, bool inverse)
{
    // This assumes N is a power of two, in-place, decimation-in-time.
    // For brevity, we skip advanced optimizations. We'll do standard Cooley-Tukey.
    // direction = +1 => forward, -1 => inverse
    int direction = inverse ? -1 : +1;

    // 1) Bit-reverse shuffle
    int j = 0;
    for(int i=0; i < N-1; i++){
        if(i < j) {
            std::swap(data[i], data[j]);
        }
        int m = N >> 1;
        while(j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    // 2) The iterative Cooley–Tukey
    for(int step = 1; step < N; step <<= 1) {
        int step2 = step << 1; // size of sub-FFT
        double theta = (direction * M_PI)/ double(step); // for e^{-i theta}

        std::complex<double> wstep(std::cos(theta), std::sin(theta));
        for(int k=0; k<N; k += step2) {
            std::complex<double> w(1.0, 0.0);
            for(int m=0; m<step; m++){
                std::complex<double> t = w * data[k+m+step];
                std::complex<double> u = data[k+m];
                data[k+m]       = u + t;
                data[k+m+step]  = u - t;
                w = w * wstep;
            }
        }
    }

    // 3) Normalize if inverse
    if(inverse) {
        double invN = 1.0/double(N);
        for(int i=0; i<N; i++) {
            data[i] *= invN;
        }
    }
}

/* ------------------------------------------------------------------------
   4. 3D transforms using the 1D cooley–tukey FFT.

   We do:
     - transform in x for each (y,z),
     - transform in y for each (x,z),
     - transform in z for each (x,y).
   We'll store data in an Nx*Ny*Nz complex array.
   We'll handle only power-of-two Nx,Ny,Nz for simplicity.
   ------------------------------------------------------------------------ */
static void forwardFFT3D(std::vector<std::complex<double>> &cplxIn, int Nx, int Ny, int Nz)
{
    // step 1: transform in x
    // for each y,z => apply 1D in x
    for(int z=0; z<Nz; z++){
      for(int y=0; y<Ny; y++){
        // gather pointer
        std::complex<double>* row = &cplxIn[ (0) + Nx*(y + Ny*z) ];
        fft1D(row, Nx, false);
      }
    }

    // step 2: transform in y
    // for each x,z => apply 1D in y, but we need to gather the data in a temp row
    {
      std::vector<std::complex<double>> temp(std::max(Nx,Ny));
      for(int z=0; z<Nz; z++){
        for(int x=0; x<Nx; x++){
          // gather the y-dimension
          for(int y=0; y<Ny; y++){
            temp[y] = cplxIn[x + Nx*(y + Ny*z)];
          }
          // 1D fft on temp
          fft1D(temp.data(), Ny, false);
          // store back
          for(int y=0; y<Ny; y++){
            cplxIn[x + Nx*(y + Ny*z)] = temp[y];
          }
        }
      }
    }

    // step 3: transform in z
    {
      std::vector<std::complex<double>> temp(std::max(Nz, Ny));
      for(int y=0; y<Ny; y++){
        for(int x=0; x<Nx; x++){
          // gather the z-dimension
          for(int z=0; z<Nz; z++){
            temp[z] = cplxIn[x + Nx*(y + Ny*z)];
          }
          // 1D fft
          fft1D(temp.data(), Nz, false);
          // store back
          for(int z=0; z<Nz; z++){
            cplxIn[x + Nx*(y + Ny*z)] = temp[z];
          }
        }
      }
    }
}

static void inverseFFT3D(std::vector<std::complex<double>> &cplxIn, int Nx, int Ny, int Nz)
{
    // step 1: inverse in z
    {
      std::vector<std::complex<double>> temp(std::max(Nz, Ny));
      for(int y=0; y<Ny; y++){
        for(int x=0; x<Nx; x++){
          // gather
          for(int z=0; z<Nz; z++){
            temp[z] = cplxIn[x + Nx*(y + Ny*z)];
          }
          fft1D(temp.data(), Nz, true);
          for(int z=0; z<Nz; z++){
            cplxIn[x + Nx*(y + Ny*z)] = temp[z];
          }
        }
      }
    }

    // step 2: inverse in y
    {
      std::vector<std::complex<double>> temp(std::max(Ny, Nx));
      for(int z=0; z<Nz; z++){
        for(int x=0; x<Nx; x++){
          // gather
          for(int y=0; y<Ny; y++){
            temp[y] = cplxIn[x + Nx*(y + Ny*z)];
          }
          fft1D(temp.data(), Ny, true);
          for(int y=0; y<Ny; y++){
            cplxIn[x + Nx*(y + Ny*z)] = temp[y];
          }
        }
      }
    }

    // step 3: inverse in x
    for(int z=0; z<Nz; z++){
      for(int y=0; y<Ny; y++){
        std::complex<double>* row = &cplxIn[ 0 + Nx*(y + Ny*z) ];
        fft1D(row, Nx, true);
      }
    }
}


/* ------------------------------------------------------------------------
   5. The main solvePoisson(...) function that decides which method to use.
   ------------------------------------------------------------------------ */

void solvePoisson(Field3D &Q, GravityMethod method)
{
    const int Nx = Q.Nx;
    const int Ny = Q.Ny;
    const int Nz = Q.Nz;
    double dx = Q.dx;
    double dy = Q.dy;
    double dz = Q.dz;
    double G  = agoge::config::G;

    // Copy Q.rho to a real buffer
    std::vector<double> realBuf(Nx*Ny*Nz, 0.0);
    for(int k=0; k<Nz; k++){
      for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
          int idx = Q.index(i,j,k);
          realBuf[idx] = Q.rho[idx];
        }
      }
    }

    // We'll store freq domain in cplxRho
    std::vector<std::complex<double>> cplxRho(Nx*Ny*Nz);

    if(method == GravityMethod::NAIVE_DFT)
    {
        // std::cout << "[GravitySolver] Using Naive DFT, O(N^6)\n";
        // 1) forward transform
        naiveForwardDFT3D(realBuf, cplxRho, Nx, Ny, Nz);
    }
    else // method == COOLEY_TUKEY
    {
        // std::cout << "[GravitySolver] Using Cooley-Tukey FFT, O(N^3 log N)\n";
        // check if Nx,Ny,Nz are power-of-two
        if(!isPowerOfTwo(Nx) || !isPowerOfTwo(Ny) || !isPowerOfTwo(Nz)) {
            std::cerr << "[GravitySolver] WARNING: Nx,Ny,Nz are not all powers of two => Cooley-Tukey won't handle general sizes here.\n";
        }
        // fill cplxRho with realBuf as complex
        for(int idx=0; idx<Nx*Ny*Nz; idx++){
            cplxRho[idx] = std::complex<double>(realBuf[idx], 0.0);
        }
        // forward 3D FFT
        forwardFFT3D(cplxRho, Nx, Ny, Nz);
    }

    // multiply by -4 pi G / k^2 in freq space
    double Lx = Nx*dx;
    double Ly = Ny*dy;
    double Lz = Nz*dz;
    const double TWO_PI = 2.0*M_PI;

    auto idxCplx = [&](int kx, int ky, int kz){
      return kx + Nx*(ky + Ny*kz);
    };

    for(int kx=0; kx<Nx; kx++){
      double kxVal = (kx <= Nx/2) ? kx : (kx - Nx);
      double alphaX = (TWO_PI/Lx)*kxVal;

      for(int ky=0; ky<Ny; ky++){
        double kyVal = (ky <= Ny/2) ? ky : (ky - Ny);
        double alphaY = (TWO_PI/Ly)*kyVal;

        for(int kz=0; kz<Nz; kz++){
          double kzVal = (kz <= Nz/2) ? kz : (kz - Nz);
          double alphaZ = (TWO_PI/Lz)*kzVal;

          double k2 = alphaX*alphaX + alphaY*alphaY + alphaZ*alphaZ;
          int idxK = idxCplx(kx, ky, kz);

          if(std::abs(k2) < 1.e-14) {
            cplxRho[idxK] = std::complex<double>(0.0, 0.0);
          } else {
            double factor = -4.0*M_PI*G / k2;
            cplxRho[idxK] *= factor;
          }
        }
      }
    }

    // inverse transform
    if(method == GravityMethod::NAIVE_DFT)
    {
        std::vector<double> tmpReal(Nx*Ny*Nz);
        naiveInverseDFT3D(cplxRho, tmpReal, Nx, Ny, Nz);
        // store in Q.phi
        for(int idx=0; idx<Nx*Ny*Nz; idx++){
            Q.phi[idx] = tmpReal[idx];
        }
    }
    else // method == COOLEY_TUKEY
    {
        // inverse 3D FFT
        inverseFFT3D(cplxRho, Nx, Ny, Nz);
        // store real part
        for(int idx=0; idx<Nx*Ny*Nz; idx++){
            Q.phi[idx] = cplxRho[idx].real();
        }
    }
    // std::cout << "[GravitySolver] Poisson solve completed.\n";
}

} // namespace gravity
} // namespace agoge
