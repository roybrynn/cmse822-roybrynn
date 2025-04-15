#include <iostream>
#include <vector>
#include <chrono>

constexpr int N = 100000;
constexpr double TOL = 1e-4;

int main() {
    std::vector<double> a(N), b(N), c(N), res(N);

    // Fill the arrays using simple for loops.
    for (int i = 0; i < N; ++i) {
        a[i] = static_cast<double>(i);
        b[i] = 2.0f * static_cast<double>(i);
        res[i] = a[i] + b[i];
    }

    // Add two vectors using a simple loop and measure time.
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        c[i] = a[i] + b[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Test results.
    int err = 0;    
    for (int i = 0; i < N; ++i) {
        double diff = c[i] - res[i];
        if (diff * diff > TOL) {
            ++err;
        }
    }

    std::cout << "Vectors added with " << err << " errors\n";
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    double flops = 2.0 * N / elapsed.count();
    std::cout << "FLOP rate: " << flops << " FLOP/s\n";
    return 0;
}
