#include <cmath>
#include <cstdlib>  // for std::atoi and std::rand
#include <iostream>
#include <random>
#include <vector>

using TYPE = double;

/// \brief Generates a random, strictly diagonally dominant matrix and scales it
/// to be near the identity matrix.
/// \param n Number of rows/columns (matrix is n x n)
/// \param A Pointer to the matrix data (row-major order)
void initDiagDomNearIdentityMatrix(int n, TYPE* A) {
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<TYPE> dist(0.0,
                                              0.022);  // (rand() % 23) / 1000.0

    for (int i = 0; i < n; ++i) {
        TYPE sum = static_cast<TYPE>(0.0);
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = dist(rng);
            sum += A[i * n + j];
        }
        A[i * n + i] += sum;
        // Scale the row so the final matrix is almost an identity matrix
        for (int j = 0; j < n; ++j) {
            A[i * n + j] /= sum;
        }
    }
}
