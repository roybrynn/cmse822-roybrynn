#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "mm_utils.hpp"

// Constants
constexpr double TOLERANCE = 0.001;
constexpr int DEF_SIZE = 1000;
constexpr int MAX_ITERS = 100000;
constexpr double LARGE = 1000000.0;

int main(int argc, char **argv) {
    int Ndim = (argc == 2) ? std::atoi(argv[1]) : DEF_SIZE;
    std::cout << "Matrix dimension (Ndim) = " << Ndim << std::endl;

    // Use std::vector to allocate storage dynamically.
    std::vector<TYPE> A(Ndim * Ndim);
    std::vector<TYPE> b(Ndim);
    std::vector<TYPE> xnew(Ndim, 0.0);
    std::vector<TYPE> xold(Ndim, 0.0);

    // Generate diagonally dominant matrix A
    initDiagDomNearIdentityMatrix(Ndim, A.data());

    // Initialize b with random values (between 0.0 and 0.5)
    for (int i = 0; i < Ndim; ++i) {
        b[i] = static_cast<TYPE>(std::rand() % 51) / 100.0;
    }

    // Implementing Jacobi Iteration
    double start_time = omp_get_wtime();
    TYPE conv = LARGE;
    int iters = 0;

    while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
        ++iters;

        #pragma omp parallel for
        for (int i = 0; i < Ndim; ++i) {
            TYPE sum = TYPE{0.0};
            for (int j = 0; j < Ndim; ++j) {
                if (i != j)
                    sum += A[i*Ndim + j]*xold[j];
            }
            xnew[i] = (b[i] - sum) / A[i*Ndim +i];
        }

        // Compute convergence criterion (Euclidean norm of difference)
        conv = static_cast<TYPE>(0.0);
        #pragma omp parallel for reduction(+: conv)
        for (int i = 0; i < Ndim; ++i) {
            TYPE tmp = xnew[i] - xold[i];
            conv += tmp * tmp;
        }
        conv = static_cast<TYPE>(std::sqrt(conv));

        // Swap vectors for next iteration
        std::swap(xold, xnew);
    }

    double elapsed_time = omp_get_wtime() - start_time;
    std::cout << "Converged after " << iters << " iterations in "
                << elapsed_time << " seconds with final convergence = "
                << conv << std::endl;

    TYPE err = 0.0, chksum = 0.0;

    for (int i = 0; i < Ndim; ++i) {
        xold[i] = static_cast<TYPE>(0.0);
        for (int j = 0; j < Ndim; ++j)
            xold[i] += A[i * Ndim + j] * xnew[j];

        TYPE diff = xold[i] - b[i];
        chksum += xnew[i];
        err += diff * diff;
    }
    err = static_cast<TYPE>(std::sqrt(err));

    std::cout << "Solution verification: Error = " << err
                << ", Checksum = " << chksum << std::endl;

    if (err > TOLERANCE)
        std::cout << "WARNING: Solution error exceeds tolerance!" << std::endl;

    return 0;
}