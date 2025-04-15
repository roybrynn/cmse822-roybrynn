// MPI Parallel Pi Calculation Using Custom MPI Data Types

#include <mpi.h>

#include <cstddef>  // for offsetof
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

constexpr int DARTS = 10000;  // darts per round
constexpr int ROUNDS = 100;   // number of rounds

// Custom struct for results
template <typename T>
struct PiResult {
    T pi_value;
    int darts_thrown;
};

// Dartboard function to approximate Pi
double compute_pi(int darts) {
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    int score = 0;
    for (int i = 0; i < darts; ++i) {
        auto x = dist(gen);
        auto y = dist(gen);
        if (x * x + y * y <= 1.0) {
            ++score;
        }
    }
    return 4.0 * static_cast<double>(score) / darts;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create custom MPI datatype
    MPI_Datatype mpi_pi_result_type;
    int block_lengths[2] = {1, 1};
    MPI_Aint displacements[2] = {offsetof(PiResult<double>, pi_value),
                                 offsetof(PiResult<double>, darts_thrown)};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};

    MPI_Type_create_struct(2, block_lengths, displacements, types,
                           &mpi_pi_result_type);
    MPI_Type_commit(&mpi_pi_result_type);

    // Each process calculates Pi over a subset of rounds
    int local_rounds = ROUNDS / size;
    std::vector<PiResult<double>> local_results(local_rounds);

    for (int i = 0; i < local_rounds; ++i) {
        double pi_estimate = compute_pi(DARTS);
        local_results[i] = PiResult<double>{pi_estimate, DARTS};
    }

    // Gather results at root process
    std::vector<PiResult<double>> all_results;
    if (rank == 0) {
        all_results.resize(ROUNDS);
    }

    MPI_Gather(local_results.data(), local_rounds, mpi_pi_result_type,
               all_results.data(), local_rounds, mpi_pi_result_type, 0,
               MPI_COMM_WORLD);

    // Process 0 computes the overall average
    if (rank == 0) {
        double sum_pi =
            std::accumulate(all_results.begin(), all_results.end(), 0.0,
                            [](double acc, const PiResult<double>& res) {
                                return acc + res.pi_value;
                            });

        double ave_pi = sum_pi / ROUNDS;

        std::cout
            << "Parallel MPI Pi Calculation using custom MPI datatypes:\n";
        for (size_t i = 0; i < all_results.size(); ++i) {
            std::cout << "After " << (i + 1) * DARTS
                      << " throws, average value of Pi = "
                      << all_results[i].pi_value << "\n";
        }

        std::cout << "\nOverall average Pi after " << ROUNDS * DARTS
                  << " throws = " << ave_pi << '\n';
        std::cout << "Real value of Pi: 3.1415926535897\n";
    }

    MPI_Type_free(&mpi_pi_result_type);
    MPI_Finalize();
    return 0;
}
