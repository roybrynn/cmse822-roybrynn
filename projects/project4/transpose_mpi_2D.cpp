#include <mpi.h>
#include <vector>
#include <iostream>

// Assume matrix dimensions are divisible by np (sqrt(world_size))
void parallelTranspose2D(std::vector<double>& localBlock, int blockSize, MPI_Comm comm) {
    int commSize;
    MPI_Comm_size(comm, &commSize);
    
    // Each process sends its block as a contiguous vector
    int sendCount = blockSize * blockSize;
    std::vector<double> recvBlock(sendCount);
    
    // Use MPI_Alltoall when each process sends an equal amount to every other process
    // For simplicity, assume a 1D reordering of blocks. In a real application,
    // you must calculate sendcounts and displacements based on the 2D process grid.
    MPI_Alltoall(localBlock.data(), sendCount, MPI_DOUBLE,
                 recvBlock.data(), sendCount, MPI_DOUBLE, comm);
    
    // After the communication, you may need to rearrange data locally.
    // Here we simply overwrite the local block.
    localBlock.swap(recvBlock);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    // Assume each MPI process holds a block of the global matrix.
    // Define blockSize (e.g., 128)
    int blockSize = 128;
    std::vector<double> localBlock(blockSize * blockSize);
    
    // Initialize localBlock with some data (omitted for brevity)
    
    // Perform parallel transpose on MPI_COMM_WORLD (or a subcommunicator if desired)
    parallelTranspose2D(localBlock, blockSize, MPI_COMM_WORLD);
    
    MPI_Finalize();
    return 0;
}