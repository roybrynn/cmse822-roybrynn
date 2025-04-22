#include <iostream>
#include <vector>
#include <cassert>

// Serial 2D Matrix Transpose (row-major storage)
void transpose2D(const std::vector<double>& A, std::vector<double>& B,
                 int rows, int cols) {
    assert(A.size() == static_cast<size_t>(rows * cols));
    B.resize(rows * cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // A is in row-major order: A[i + rows * j]
            // Transpose: B[j + cols * i]
            B[j + cols * i] = A[i + cols * j];
        }
    }
}

void transpose3D(const std::vector<double>& A, std::vector<double>& B,
    int nx, int ny, int nz) {
    assert(A.size() == static_cast<size_t>(nx*ny*nz));
    B.resize(nx*ny*nz);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k)
            // A is in row-major order: A[i + rows * j]
            // Transpose: B[j + cols * i]
            B[j + ny * (i + nx * k)] = A[ i + nx * (j + ny * k)];
        }
    }
}

int main() {
    // int rows = 3, cols = 4;
    // std::vector<double> A(rows * cols);
    // std::vector<double> B;
    
    // // Initialize A with some values
    // for (int i = 0; i < rows * cols; ++i)
    //     A[i] = i + 1;
    
    // transpose2D(A, B, rows, cols);
    
    // // Print the transposed matrix
    // std::cout << "Transposed 2D Matrix:" << std::endl;
    // for (int i = 0; i < cols; ++i) {
    //     for (int j = 0; j < rows; ++j) {
    //         std::cout << B[i + cols * j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // return 0;

    int nx = 3, ny = 4, nz = 3;
    std::vector<double> A(nx * ny * nz);
    std::vector<double> B;

    for (int i = 0; i < nx * ny * nz; ++i)
        A[i] = i + 1;

    for (int k = 0; k < nz; ++k){
        for (int j = 0; j < ny; ++j){
            for (int i = 0; i < nx; ++i) {
                std::cout << A[i + nx * (j + ny * k)] << " ";
            }
            std::cout << std::endl;

        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;
    transpose3D(A, B, nx, ny, nz);

    std::cout << "Transposed 3D Matrix:" << std::endl;
    std::swap(nx, ny);
    for (int k = 0; k < nz; ++k){
        for (int j = 0; j < ny; ++j){
            for (int i = 0; i < nx; ++i) {
                std::cout << B[i + nx * (j + ny * k)] << " ";
            }
        std::cout << std::endl;

        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;
}