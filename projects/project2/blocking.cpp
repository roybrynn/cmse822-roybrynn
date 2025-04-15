#include <mpi.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// Helper function to convert (x, y, z) into a 1D index for our arrays.
// Here, Nz_with_halo = local_nz + 2 (one ghost layer at each end).
inline int index3D(int x, int y, int z, int Nx, int Ny, int Nz_with_halo) {
    return (z * Ny * Nx) + (y * Nx) + x;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Global domain dimensions (for demonstration):
    // We'll decompose along z: Nx, Ny remain the same for each rank.
    const int Nx = 16;
    const int Ny = 16;
    const int Nz = 32;

    // Assume size divides Nz exactly for brevity:
    const int local_nz = Nz / size;

    // Each rank stores local_nz plus 2 ghost layers in z:
    const int Nz_with_halo = local_nz + 2;

    // Allocate arrays for current data and result data
    // in row-major storage:
    std::vector<double> data(Nx * Ny * Nz_with_halo, 0.0);
    std::vector<double> result(Nx * Ny * Nz_with_halo, 0.0);

    // Initialize interior cells with some artificial values
    for (int z = 1; z <= local_nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx = index3D(x, y, z, Nx, Ny, Nz_with_halo);
                data[idx] = 100.0 * rank + z;  // e.g., a simple formula
            }
        }
    }

    // Identify neighbor ranks in the z-direction:
    // up_rank = rank-1; down_rank = rank+1,
    // with MPI_PROC_NULL for out-of-bounds
    int up_rank = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    // For a single exchange, we need to send:
    //  - The plane at z=1  to up_rank
    //  - The plane at z=local_nz to down_rank
    //
    // We'll store them in buffers (send_up/down, recv_up/down)
    std::vector<double> send_up(Nx * Ny), recv_up(Nx * Ny);
    std::vector<double> send_down(Nx * Ny), recv_down(Nx * Ny);

    // Pack the "top" interior plane (z=1) into send_up
    {
        int zplane = 1;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx_local = index3D(x, y, zplane, Nx, Ny, Nz_with_halo);
                int idx_plane = y * Nx + x;
                send_up[idx_plane] = data[idx_local];
            }
        }
    }

    // Pack the "bottom" interior plane (z=local_nz) into send_down
    {
        int zplane = local_nz;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx_local = index3D(x, y, zplane, Nx, Ny, Nz_with_halo);
                int idx_plane = y * Nx + x;
                send_down[idx_plane] = data[idx_local];
            }
        }
    }

    // ===========================
    // EXCHANGE WITH UP NEIGHBOR
    // ===========================
    // 1) If we have an up_rank, we send z=1 to up_rank:
    if (up_rank != MPI_PROC_NULL) {
        MPI_Send(send_up.data(), Nx * Ny, MPI_DOUBLE, up_rank, 101,
                 MPI_COMM_WORLD);
    }
    // 2) If we have an up_rank, receive from up_rank into recv_up:
    if (up_rank != MPI_PROC_NULL) {
        MPI_Recv(recv_up.data(), Nx * Ny, MPI_DOUBLE, up_rank, 102,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // ===========================
    // EXCHANGE WITH DOWN NEIGHBOR
    // ===========================
    // 3) If we have a down_rank, we send z=local_nz to down_rank:
    if (down_rank != MPI_PROC_NULL) {
        MPI_Send(send_down.data(), Nx * Ny, MPI_DOUBLE, down_rank, 102,
                 MPI_COMM_WORLD);
    }
    // 4) If we have a down_rank, receive from down_rank into recv_down:
    if (down_rank != MPI_PROC_NULL) {
        MPI_Recv(recv_down.data(), Nx * Ny, MPI_DOUBLE, down_rank, 101,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Because we used blocking calls, at this point we have updated the
    // top and bottom "ghost" planes in recv_up and recv_down.
    // Next, we copy them into our local arrays at z=0 and z=local_nz+1.

    // Store ghost data from up_rank into z=0
    if (up_rank != MPI_PROC_NULL) {
        int zghost = 0;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx_local = index3D(x, y, zghost, Nx, Ny, Nz_with_halo);
                int idx_plane = y * Nx + x;
                data[idx_local] = recv_up[idx_plane];
            }
        }
    }
    // Store ghost data from down_rank into z=local_nz+1
    if (down_rank != MPI_PROC_NULL) {
        int zghost = local_nz + 1;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx_local = index3D(x, y, zghost, Nx, Ny, Nz_with_halo);
                int idx_plane = y * Nx + x;
                data[idx_local] = recv_down[idx_plane];
            }
        }
    }

    // ===========================
    // 3-POINT STENCIL CALCULATION
    // ===========================
    // We'll do a simple average over z: data[z-1], data[z], data[z+1].
    // First, update interior (z=2..local_nz-1) ignoring ghost
    for (int z = 2; z < local_nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idxZ = index3D(x, y, z, Nx, Ny, Nz_with_halo);
                int idxZm1 = index3D(x, y, z - 1, Nx, Ny, Nz_with_halo);
                int idxZp1 = index3D(x, y, z + 1, Nx, Ny, Nz_with_halo);
                result[idxZ] = (data[idxZm1] + data[idxZ] + data[idxZp1]) / 3.0;
            }
        }
    }

    // Then handle boundary planes z=1 and z=local_nz
    // (which depend on z=0 or z=local_nz+1 for ghost data)
    {
        int z = 1;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idxZ = index3D(x, y, z, Nx, Ny, Nz_with_halo);
                int idxZm1 = index3D(x, y, z - 1, Nx, Ny, Nz_with_halo);
                int idxZp1 = index3D(x, y, z + 1, Nx, Ny, Nz_with_halo);
                result[idxZ] = (data[idxZm1] + data[idxZ] + data[idxZp1]) / 3.0;
            }
        }
    }
    {
        int z = local_nz;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idxZ = index3D(x, y, z, Nx, Ny, Nz_with_halo);
                int idxZm1 = index3D(x, y, z - 1, Nx, Ny, Nz_with_halo);
                int idxZp1 = index3D(x, y, z + 1, Nx, Ny, Nz_with_halo);
                result[idxZ] = (data[idxZm1] + data[idxZ] + data[idxZp1]) / 3.0;
            }
        }
    }

    // (Optional) Print or collect some data for verification
    if (rank == 0) {
        int corner_idx = index3D(0, 0, 1, Nx, Ny, Nz_with_halo);
        std::cout << "Rank 0 corner result: " << result[corner_idx] << "\n";
    }

    MPI_Finalize();
    return 0;
}