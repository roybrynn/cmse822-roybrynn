/**
 * @file GravitySolver.cpp
 * @brief Implementation of gravity solvers for the Agoge application.
 *
 * This file implements methods for solving the Poisson equation using
 * both naive DFT and the Cooley–Tukey FFT algorithm, along with necessary
 * helper routines.
 */

#include "GravitySolver.hpp"
#include "Config.hpp"
#include "Field3d.hpp"

#include <cmath>
#include <iostream>
#include <vector>
#include <complex>

namespace agoge
{
    namespace gravity
    {

        static void naiveForwardDFT3D(const std::vector<double> &realIn,
                                      std::vector<std::complex<double>> &cplxOut,
                                      int Nx, int Ny, int Nz)
        {
            cplxOut.assign(Nx * Ny * Nz, {0.0, 0.0});
            const double TWO_PI = 2.0 * M_PI;
            for (int kx = 0; kx < Nx; kx++)
            {
                for (int ky = 0; ky < Ny; ky++)
                {
                    for (int kz = 0; kz < Nz; kz++)
                    {
                        int idxK = kx + Nx * (ky + Ny * kz);
                        std::complex<double> sumVal(0.0, 0.0);
                        for (int x = 0; x < Nx; x++)
                        {
                            double alphaX =
                                static_cast<double>(kx * x) / static_cast<double>(Nx);
                            for (int y = 0; y < Ny; y++)
                            {
                                double alphaY = static_cast<double>(ky * y) /
                                                static_cast<double>(Ny);
                                for (int z = 0; z < Nz; z++)
                                {
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
        static void naiveInverseDFT3D(const std::vector<std::complex<double>> &cplxIn,
                                      std::vector<double> &realOut, int Nx, int Ny,
                                      int Nz)
        {
            realOut.assign(Nx * Ny * Nz, 0.0);
            const double TWO_PI = 2.0 * M_PI;
            double invSize = 1.0 / static_cast<double>(Nx * Ny * Nz);
            for (int x = 0; x < Nx; x++)
            {
                for (int y = 0; y < Ny; y++)
                {
                    for (int z = 0; z < Nz; z++)
                    {
                        int idxXYZ = x + Nx * (y + Ny * z);
                        std::complex<double> sumVal(0.0, 0.0);
                        for (int kx = 0; kx < Nx; kx++)
                        {
                            double alphaX =
                                static_cast<double>(kx * x) / static_cast<double>(Nx);
                            for (int ky = 0; ky < Ny; ky++)
                            {
                                double alphaY = static_cast<double>(ky * y) /
                                                static_cast<double>(Ny);
                                for (int kz = 0; kz < Nz; kz++)
                                {
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
        static void fft1D(std::complex<double> *data, int N, bool inverse)
        {
            int direction = inverse ? -1 : +1;

            // Bit-reversal permutation.
            int j = 0;
            for (int i = 0; i < N - 1; i++)
            {
                if (i < j)
                {
                    std::swap(data[i], data[j]);
                }
                int m = N >> 1;
                while (j >= m && m > 0)
                {
                    j -= m;
                    m >>= 1;
                }
                j += m;
            }

            // Cooley–Tukey iterative FFT.
            for (int step = 1; step < N; step <<= 1)
            {
                int step2 = step << 1;
                double theta = (direction * M_PI) / static_cast<double>(step);
                std::complex<double> wstep(std::cos(theta), std::sin(theta));
                for (int k = 0; k < N; k += step2)
                {
                    std::complex<double> w(1.0, 0.0);
                    for (int m = 0; m < step; m++)
                    {
                        std::complex<double> t = w * data[k + m + step];
                        std::complex<double> u = data[k + m];
                        data[k + m] = u + t;
                        data[k + m + step] = u - t;
                        w *= wstep;
                    }
                }
            }

            // Normalize if performing an inverse FFT.
            if (inverse)
            {
                double invN = 1.0 / static_cast<double>(N);
                for (int i = 0; i < N; i++)
                {
                    data[i] *= invN;
                }
            }
        }

        //------------------------------------------------------------------------------
        // Forward 3D FFT: Apply successive 1D FFTs along x, then y, then z.
        static void forwardFFT3D(std::vector<std::complex<double>> &cplxIn, int Nx,
                                 int Ny, int Nz)
        {
            // FFT along x-direction for each (y,z)
            for (int z = 0; z < Nz; z++)
            {
                for (int y = 0; y < Ny; y++)
                {
                    fft1D(&cplxIn[0] + (0 + Nx * (y + Ny * z)), Nx, false);
                }
            }
            // FFT along y-direction for each (x,z)
            {
                std::vector<std::complex<double>> temp(std::max(Nx, Ny));
                for (int z = 0; z < Nz; z++)
                {
                    for (int x = 0; x < Nx; x++)
                    {
                        for (int y = 0; y < Ny; y++)
                        {
                            temp[y] = cplxIn[x + Nx * (y + Ny * z)];
                        }
                        fft1D(temp.data(), Ny, false);
                        for (int y = 0; y < Ny; y++)
                        {
                            cplxIn[x + Nx * (y + Ny * z)] = temp[y];
                        }
                    }
                }
            }
            // FFT along z-direction for each (x,y)
            {
                std::vector<std::complex<double>> temp(std::max(Nz, Ny));
                for (int y = 0; y < Ny; y++)
                {
                    for (int x = 0; x < Nx; x++)
                    {
                        for (int z = 0; z < Nz; z++)
                        {
                            temp[z] = cplxIn[x + Nx * (y + Ny * z)];
                        }
                        fft1D(temp.data(), Nz, false);
                        for (int z = 0; z < Nz; z++)
                        {
                            cplxIn[x + Nx * (y + Ny * z)] = temp[z];
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        // Inverse 3D FFT: Apply successive inverse 1D FFTs along z, then y, then x.
        static void inverseFFT3D(std::vector<std::complex<double>> &cplxIn, int Nx,
                                 int Ny, int Nz)
        {
            // Inverse FFT along z-direction for each (x,y)
            {
                std::vector<std::complex<double>> temp(std::max(Nz, Ny));
                for (int y = 0; y < Ny; y++)
                {
                    for (int x = 0; x < Nx; x++)
                    {
                        for (int z = 0; z < Nz; z++)
                        {
                            temp[z] = cplxIn[x + Nx * (y + Ny * z)];
                        }
                        fft1D(temp.data(), Nz, true);
                        for (int z = 0; z < Nz; z++)
                        {
                            cplxIn[x + Nx * (y + Ny * z)] = temp[z];
                        }
                    }
                }
            }
            // Inverse FFT along y-direction for each (x,z)
            {
                std::vector<std::complex<double>> temp(std::max(Ny, Nx));
                for (int z = 0; z < Nz; z++)
                {
                    for (int x = 0; x < Nx; x++)
                    {
                        for (int y = 0; y < Ny; y++)
                        {
                            temp[y] = cplxIn[x + Nx * (y + Ny * z)];
                        }
                        fft1D(temp.data(), Ny, true);
                        for (int y = 0; y < Ny; y++)
                        {
                            cplxIn[x + Nx * (y + Ny * z)] = temp[y];
                        }
                    }
                }
            }
            // Inverse FFT along x-direction for each (y,z)
            for (int z = 0; z < Nz; z++)
            {
                for (int y = 0; y < Ny; y++)
                {
                    fft1D(&cplxIn[0] + (0 + Nx * (y + Ny * z)), Nx, true);
                }
            }
        }

        static inline bool isPowerOfTwo(int n)
        {
            return (n > 0) && ((n & (n - 1)) == 0);
        }

        /**
         * @brief Subtract the mean density from realBuf, implementing the "Jeans swindle".
         */
        static void subtractMeanDensity(std::vector<double> &rhoBuf)
        {
            double sum = 0.0;
            for (double val : rhoBuf)
            {
                sum += val;
            }
            double mean = sum / double(rhoBuf.size());
            for (double &val : rhoBuf)
            {
                val = val - mean; // subtract mean
            }
        }

        /**
         * @brief Solves the Poisson equation for gravitational potential.
         *
         * Depending on the chosen method, this function applies either a naive DFT
         * or a Cooley-Tukey FFT approach to compute potential fluctuations.
         *
         * @param Q Field3D object containing density and ghost data.
         * @param method The gravity solution method.
         */
        void solvePoisson(Field3D &Q, GravityMethod method)
        {
            // 1) gather interior size
            int Nx = Q.Nx;
            int Ny = Q.Ny;
            int Nz = Q.Nz;
            double dx = Q.dx;
            double dy = Q.dy;
            double dz = Q.dz;
            double G = agoge::config::G;

            // 2) copy interior density into realBuf
            std::vector<double> realBuf(Nx * Ny * Nz, 0.0);
            for (int k = 0; k < Nz; k++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int i = 0; i < Nx; i++)
                    {
                        // map to interior index
                        int iG = i + config::Ng;
                        int jG = j + config::Ng;
                        int kG = k + config::Ng;
                        int iField = Q.index(iG, jG, kG);
                        realBuf[i + Nx * (j + Ny * k)] = Q.rho[iField];
                    }
                }
            }

            // 3) Subtract mean => "Jeans swindle"
            // multiply by 4 pi G as well => so the transform is 4 pi G * (rho-mean)
            subtractMeanDensity(realBuf);
            for (double &val : realBuf)
            {
                val *= (4.0 * M_PI * G);
            }

            // 4) transform to frequency domain
            std::vector<std::complex<double>> cplxData(Nx * Ny * Nz);
            if (method == GravityMethod::NAIVE_DFT)
            {
                naiveForwardDFT3D(realBuf, cplxData, Nx, Ny, Nz);
            }
            else
            {
                // check power-of-two
                if (!isPowerOfTwo(Nx) || !isPowerOfTwo(Ny) || !isPowerOfTwo(Nz))
                {
                    std::cerr << "[GravitySolver] Nx,Ny,Nz must be powers of two for Cooley-Tukey.\n";
                    return;
                }
                // fill cplxData
                for (int idx = 0; idx < Nx * Ny * Nz; idx++)
                {
                    cplxData[idx] = {realBuf[idx], 0.0};
                }
                forwardFFT3D(cplxData, Nx, Ny, Nz);
            }

            // 5) multiply by -1/k^2 except k=0 => 0
            double Lx = Nx * dx;
            double Ly = Ny * dy;
            double Lz = Nz * dz;
            auto idx3D = [&](int kx, int ky, int kz)
            {
                return kx + Nx * (ky + Ny * kz);
            };
            const double TWO_PI = 2.0 * M_PI;
            for (int kx = 0; kx < Nx; kx++)
            {
                // wave number shift
                double kxVal = (kx <= Nx / 2) ? kx : (kx - Nx);
                double alphaX = (TWO_PI / Lx) * kxVal;
                for (int ky = 0; ky < Ny; ky++)
                {
                    double kyVal = (ky <= Ny / 2) ? ky : (ky - Ny);
                    double alphaY = (TWO_PI / Ly) * kyVal;
                    for (int kz = 0; kz < Nz; kz++)
                    {
                        double kzVal = (kz <= Nz / 2) ? kz : (kz - Nz);
                        double alphaZ = (TWO_PI / Lz) * kzVal;
                        double k2 = alphaX * alphaX + alphaY * alphaY + alphaZ * alphaZ;
                        int idxK = idx3D(kx, ky, kz);
                        if (std::abs(k2) < 1.e-100)
                        {
                            // set DC=0 => no infinite offset
                            cplxData[idxK] = {0.0, 0.0};
                        }
                        else
                        {
                            // multiply by -1/k^2
                            double factor = -1.0 / k2;
                            cplxData[idxK] *= factor;
                        }
                    }
                }
            }

            // 6) inverse transform
            if (method == GravityMethod::NAIVE_DFT)
            {
                std::vector<double> tmpReal;
                naiveInverseDFT3D(cplxData, tmpReal, Nx, Ny, Nz);
                // store in Q.phi
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            int iG = i + config::Ng;
                            int jG = j + config::Ng;
                            int kG = k + config::Ng;
                            int iField = Q.index(iG, jG, kG);
                            int idx = i + Nx * (j + Ny * k);
                            Q.phi[iField] = tmpReal[idx];
                        }
                    }
                }
            }
            else
            {
                inverseFFT3D(cplxData, Nx, Ny, Nz);
                // store real part
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            int iG = i + config::Ng;
                            int jG = j + config::Ng;
                            int kG = k + config::Ng;
                            int iField = Q.index(iG, jG, kG);

                            int idx = i + Nx * (j + Ny * k);
                            Q.phi[iField] = cplxData[idx].real();
                        }
                    }
                }
            }

            // finished
            // Q.phi is now the potential fluctuations with net zero mass in periodic box
        }

    } // namespace gravity
} // namespace agoge