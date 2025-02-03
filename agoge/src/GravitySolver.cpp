#include "GravitySolver.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "Config.hpp"
#include "Field3d.hpp"

namespace agoge {
namespace gravity {

//------------------------------------------------------------------------------
// Helper: Check if a number is a power of two.
static bool isPowerOfTwo(int n) { return (n > 0) && ((n & (n - 1)) == 0); }

//------------------------------------------------------------------------------
// Naive forward DFT (3D) --- O(N⁶). For educational purposes only.
static void naiveForwardDFT3D(const std::vector<double>& realIn,
                              std::vector<std::complex<double>>& cplxOut,
                              int Nx, int Ny, int Nz) {
    cplxOut.assign(Nx * Ny * Nz, {0.0, 0.0});
    const double TWO_PI = 2.0 * M_PI;
    for (int kx = 0; kx < Nx; kx++) {
        for (int ky = 0; ky < Ny; ky++) {
            for (int kz = 0; kz < Nz; kz++) {
                int idxK = kx + Nx * (ky + Ny * kz);
                std::complex<double> sumVal(0.0, 0.0);
                for (int x = 0; x < Nx; x++) {
                    double alphaX =
                        static_cast<double>(kx * x) / static_cast<double>(Nx);
                    for (int y = 0; y < Ny; y++) {
                        double alphaY = static_cast<double>(ky * y) /
                                        static_cast<double>(Ny);
                        for (int z = 0; z < Nz; z++) {
                            double alphaZ = static_cast<double>(kz * z) /
                                            static_cast<double>(Nz);
                            double phase = -TWO_PI * (alphaX + alphaY + alphaZ);
                            int idxXYZ = x + Nx * (y + Ny * z);
                            double val = realIn[idxXYZ];
                            sumVal +=
                                val * std::complex<double>(std::cos(phase),
                                                           std::sin(phase));
                        }
                    }
                }
                cplxOut[idxK] = sumVal;
            }
        }
    }
}

//------------------------------------------------------------------------------
// Naive inverse DFT (3D) --- O(N⁶). For educational purposes only.
static void naiveInverseDFT3D(const std::vector<std::complex<double>>& cplxIn,
                              std::vector<double>& realOut, int Nx, int Ny,
                              int Nz) {
    realOut.assign(Nx * Ny * Nz, 0.0);
    const double TWO_PI = 2.0 * M_PI;
    double invSize = 1.0 / static_cast<double>(Nx * Ny * Nz);
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                int idxXYZ = x + Nx * (y + Ny * z);
                std::complex<double> sumVal(0.0, 0.0);
                for (int kx = 0; kx < Nx; kx++) {
                    double alphaX =
                        static_cast<double>(kx * x) / static_cast<double>(Nx);
                    for (int ky = 0; ky < Ny; ky++) {
                        double alphaY = static_cast<double>(ky * y) /
                                        static_cast<double>(Ny);
                        for (int kz = 0; kz < Nz; kz++) {
                            double alphaZ = static_cast<double>(kz * z) /
                                            static_cast<double>(Nz);
                            double phase = TWO_PI * (alphaX + alphaY + alphaZ);
                            int idxK = kx + Nx * (ky + Ny * kz);
                            sumVal += cplxIn[idxK] *
                                      std::complex<double>(std::cos(phase),
                                                           std::sin(phase));
                        }
                    }
                }
                realOut[idxXYZ] = sumVal.real() * invSize;
            }
        }
    }
}

//------------------------------------------------------------------------------
// 1D FFT using the Cooley–Tukey algorithm (in-place, radix-2).
static void fft1D(std::complex<double>* data, int N, bool inverse) {
    int direction = inverse ? -1 : +1;

    // Bit-reversal permutation.
    int j = 0;
    for (int i = 0; i < N - 1; i++) {
        if (i < j) {
            std::swap(data[i], data[j]);
        }
        int m = N >> 1;
        while (j >= m && m > 0) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    // Cooley–Tukey iterative FFT.
    for (int step = 1; step < N; step <<= 1) {
        int step2 = step << 1;
        double theta = (direction * M_PI) / static_cast<double>(step);
        std::complex<double> wstep(std::cos(theta), std::sin(theta));
        for (int k = 0; k < N; k += step2) {
            std::complex<double> w(1.0, 0.0);
            for (int m = 0; m < step; m++) {
                std::complex<double> t = w * data[k + m + step];
                std::complex<double> u = data[k + m];
                data[k + m] = u + t;
                data[k + m + step] = u - t;
                w *= wstep;
            }
        }
    }

    // Normalize if performing an inverse FFT.
    if (inverse) {
        double invN = 1.0 / static_cast<double>(N);
        for (int i = 0; i < N; i++) {
            data[i] *= invN;
        }
    }
}

//------------------------------------------------------------------------------
// Forward 3D FFT: Apply successive 1D FFTs along x, then y, then z.
static void forwardFFT3D(std::vector<std::complex<double>>& cplxIn, int Nx,
                         int Ny, int Nz) {
    // FFT along x-direction for each (y,z)
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            fft1D(&cplxIn[0] + (0 + Nx * (y + Ny * z)), Nx, false);
        }
    }
    // FFT along y-direction for each (x,z)
    {
        std::vector<std::complex<double>> temp(std::max(Nx, Ny));
        for (int z = 0; z < Nz; z++) {
            for (int x = 0; x < Nx; x++) {
                for (int y = 0; y < Ny; y++) {
                    temp[y] = cplxIn[x + Nx * (y + Ny * z)];
                }
                fft1D(temp.data(), Ny, false);
                for (int y = 0; y < Ny; y++) {
                    cplxIn[x + Nx * (y + Ny * z)] = temp[y];
                }
            }
        }
    }
    // FFT along z-direction for each (x,y)
    {
        std::vector<std::complex<double>> temp(std::max(Nz, Ny));
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                for (int z = 0; z < Nz; z++) {
                    temp[z] = cplxIn[x + Nx * (y + Ny * z)];
                }
                fft1D(temp.data(), Nz, false);
                for (int z = 0; z < Nz; z++) {
                    cplxIn[x + Nx * (y + Ny * z)] = temp[z];
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
// Inverse 3D FFT: Apply successive inverse 1D FFTs along z, then y, then x.
static void inverseFFT3D(std::vector<std::complex<double>>& cplxIn, int Nx,
                         int Ny, int Nz) {
    // Inverse FFT along z-direction for each (x,y)
    {
        std::vector<std::complex<double>> temp(std::max(Nz, Ny));
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                for (int z = 0; z < Nz; z++) {
                    temp[z] = cplxIn[x + Nx * (y + Ny * z)];
                }
                fft1D(temp.data(), Nz, true);
                for (int z = 0; z < Nz; z++) {
                    cplxIn[x + Nx * (y + Ny * z)] = temp[z];
                }
            }
        }
    }
    // Inverse FFT along y-direction for each (x,z)
    {
        std::vector<std::complex<double>> temp(std::max(Ny, Nx));
        for (int z = 0; z < Nz; z++) {
            for (int x = 0; x < Nx; x++) {
                for (int y = 0; y < Ny; y++) {
                    temp[y] = cplxIn[x + Nx * (y + Ny * z)];
                }
                fft1D(temp.data(), Ny, true);
                for (int y = 0; y < Ny; y++) {
                    cplxIn[x + Nx * (y + Ny * z)] = temp[y];
                }
            }
        }
    }
    // Inverse FFT along x-direction for each (y,z)
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            fft1D(&cplxIn[0] + (0 + Nx * (y + Ny * z)), Nx, true);
        }
    }
}

//------------------------------------------------------------------------------
// Main function to solve Poisson's equation using either the naive DFT or
// Cooley–Tukey FFT approach.
void solvePoisson(Field3D& Q, GravityMethod method) {
    // Get interior dimensions and cell sizes.
    const int Nx = Q.Nx;
    const int Ny = Q.Ny;
    const int Nz = Q.Nz;
    const double dx = Q.dx;
    const double dy = Q.dy;
    const double dz = Q.dz;
    const double G = agoge::config::G;

    // Copy the interior density (rho) into a contiguous real buffer.
    std::vector<double> realBuf(Nx * Ny * Nz, 0.0);
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = Q.interiorIndex(i, j, k);
                realBuf[i + Nx * (j + Ny * k)] = Q.rho[idx];
            }
        }
    }

    // Prepare the frequency-domain buffer.
    std::vector<std::complex<double>> cplxRho(Nx * Ny * Nz);
    if (method == GravityMethod::NAIVE_DFT) {
        // std::cout << "[GravitySolver] Using Naive DFT (O(N⁶)).\n";
        naiveForwardDFT3D(realBuf, cplxRho, Nx, Ny, Nz);
    } else {  // COOLEY_TUKEY
        // std::cout << "[GravitySolver] Using Cooley–Tukey FFT (O(N³ log N)).\n";
        if (!isPowerOfTwo(Nx) || !isPowerOfTwo(Ny) || !isPowerOfTwo(Nz)) {
            std::cerr << "[GravitySolver] ERROR: Grid dimensions must be "
                         "powers of two for the Cooley–Tukey FFT.\n";
            return;
        }
        // Initialize cplxRho with the real density data.
        for (int idx = 0; idx < Nx * Ny * Nz; idx++) {
            cplxRho[idx] = std::complex<double>(realBuf[idx], 0.0);
        }
        forwardFFT3D(cplxRho, Nx, Ny, Nz);
    }

    // Apply the Green's function in Fourier space:
    // For each mode (except the zero mode) multiply by -4πG / k².
    const double Lx = Nx * dx;
    const double Ly = Ny * dy;
    const double Lz = Nz * dz;
    const double TWO_PI = 2.0 * M_PI;
    auto idxCplx = [Nx, Ny](int kx, int ky, int kz) -> int {
        return kx + Nx * (ky + Ny * kz);
    };

    for (int kx = 0; kx < Nx; kx++) {
        double kxVal = (kx <= Nx / 2) ? kx : kx - Nx;
        double alphaX = (TWO_PI / Lx) * kxVal;
        for (int ky = 0; ky < Ny; ky++) {
            double kyVal = (ky <= Ny / 2) ? ky : ky - Ny;
            double alphaY = (TWO_PI / Ly) * kyVal;
            for (int kz = 0; kz < Nz; kz++) {
                double kzVal = (kz <= Nz / 2) ? kz : kz - Nz;
                double alphaZ = (TWO_PI / Lz) * kzVal;
                double k2 = alphaX * alphaX + alphaY * alphaY + alphaZ * alphaZ;
                int idxK = idxCplx(kx, ky, kz);
                if (std::abs(k2) < 1.e-14) {
                    cplxRho[idxK] = {0.0, 0.0};
                } else {
                    double factor = -4.0 * M_PI * G / k2;
                    cplxRho[idxK] *= factor;
                }
            }
        }
    }

    // Inverse transform back to real space.
    if (method == GravityMethod::NAIVE_DFT) {
        std::vector<double> tmpReal;
        naiveInverseDFT3D(cplxRho, tmpReal, Nx, Ny, Nz);
        // Store the result in Q.phi (interior cells only).
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int interiorIdx = Q.interiorIndex(i, j, k);
                    Q.phi[interiorIdx] = tmpReal[i + Nx * (j + Ny * k)];
                }
            }
        }
    } else {  // COOLEY_TUKEY
        inverseFFT3D(cplxRho, Nx, Ny, Nz);
        // Store the real part of the inverse transform in Q.phi (interior cells
        // only).
        for (int k = 0; k < Nz; k++) {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int interiorIdx = Q.interiorIndex(i, j, k);
                    Q.phi[interiorIdx] = cplxRho[i + Nx * (j + Ny * k)].real();
                }
            }
        }
    }
    // std::cout << "[GravitySolver] Poisson solve completed.\n";
}

}  // namespace gravity
}  // namespace agoge
