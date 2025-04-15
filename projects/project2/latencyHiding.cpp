#include <mpi.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// Helper function to convert (x, y, z) into a 1D index.
// We store data with an extra boundary in z: local_nz + 2 ghost layers.
inline int index3D(int x, int y, int z, int Nx, int Ny, int Nz_with_halo) {
    // Nx, Ny are the original domain's x and y extents
    // Nz_with_halo = local_nz + 2 for the ghost layers in z
    return (z * Ny * Nx) + (y * Nx) + x;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Global 3D domain dimensions (for demonstration)
    // Decomposed along the z-axis. Nx, Ny remain the same on each rank.
    const int Nx = 16;
    const int Ny = 16;
    const int Nz = 32;

    // For simplicity, assume Nx, Ny, Nz are divisible by 1,1,size respectively
    // so that local_nz = Nz / size is integer.
    const int local_nz = Nz / size;

    // Each process will own z in [rank*local_nz, (rank+1)*local_nz - 1]
    // plus 1 "ghost" layer at the top, 1 "ghost" layer at the bottom
    // => total z-extent = local_nz + 2.
    const int Nz_with_halo = local_nz + 2;

    // Allocate arrays for "current" and "next" time-step or for a simple
    // average. For a real FD solver, you'd manage multiple arrays or apply
    // in-place.
    std::vector<double> data(Nx * Ny * Nz_with_halo, 0.0);
    std::vector<double> result(Nx * Ny * Nz_with_halo, 0.0);

    // Initialize local domain (e.g., some artificial values).
    // We ignore the ghost layers for now, just fill the interior.
    for (int z = 1; z <= local_nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx = index3D(x, y, z, Nx, Ny, Nz_with_halo);
                // Example initialization: combine local indices & rank
                data[idx] = 100.0 * rank + z;
            }
        }
    }

    // Determine neighbor ranks in the z direction.
    // In 1D decomposition, "up" neighbor is rank-1, "down" neighbor is rank+1
    // use MPI_PROC_NULL for boundaries
    int up_rank = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    // We'll exchange the planes z=1 (top interior) and z=local_nz (bottom
    // interior). Each plane has Nx*Ny elements.

    // Helper arrays for sending/receiving
    std::vector<double> send_up(Nx * Ny), recv_up(Nx * Ny);
    std::vector<double> send_down(Nx * Ny), recv_down(Nx * Ny);

    // Prepare for halo exchange: define MPI_Requests.
    // We'll do nonblocking sends/recvs for demonstration.
    MPI_Request reqs[4];

    // For demonstration, do a single halo exchange + stencil iteration
    // (In a real code, you'd do this in a time-stepping loop.)

    // 1) Pack data to send "up" (to rank-1).
    //    This is the plane z=1 in local storage (the top interior plane).
    {
        int zplane = 1;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int local_idx = index3D(x, y, zplane, Nx, Ny, Nz_with_halo);
                int plane_idx = y * Nx + x;
                send_up[plane_idx] = data[local_idx];
            }
        }
    }

    // 2) Pack data to send "down" (to rank+1).
    //    This is the plane z=local_nz in local storage (bottom interior plane).
    {
        int zplane = local_nz;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int local_idx = index3D(x, y, zplane, Nx, Ny, Nz_with_halo);
                int plane_idx = y * Nx + x;
                send_down[plane_idx] = data[local_idx];
            }
        }
    }

    // Post Irecvs first (standard best practice)
    MPI_Irecv(recv_up.data(), Nx * Ny, MPI_DOUBLE, up_rank, 101, MPI_COMM_WORLD,
              &reqs[0]);
    MPI_Irecv(recv_down.data(), Nx * Ny, MPI_DOUBLE, down_rank, 102,
              MPI_COMM_WORLD, &reqs[1]);

    // Then post Isends
    MPI_Isend(send_down.data(), Nx * Ny, MPI_DOUBLE, down_rank, 101,
              MPI_COMM_WORLD, &reqs[2]);
    MPI_Isend(send_up.data(), Nx * Ny, MPI_DOUBLE, up_rank, 102, MPI_COMM_WORLD,
              &reqs[3]);

    // Now we can do some local interior work while messages are in-flight.
    // For example, compute a 3-point stencil for z=2..(local_nz-1).
    // We'll average only in z-direction for demonstration:
    for (int z = 2; z < local_nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                // Stencil:  (data[z-1] + data[z] + data[z+1]) / 3
                // ignoring x,y dimension for simplicity
                int idxZ = index3D(x, y, z, Nx, Ny, Nz_with_halo);
                int idxZm1 = index3D(x, y, z - 1, Nx, Ny, Nz_with_halo);
                int idxZp1 = index3D(x, y, z + 1, Nx, Ny, Nz_with_halo);
                result[idxZ] = (data[idxZm1] + data[idxZ] + data[idxZp1]) / 3.0;
            }
        }
    }

    // Wait for the communications to complete, so we can update boundary planes
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // Unpack the newly received ghost data into z=0 (from up_rank) and
    // z=local_nz+1 (from down_rank) (In a real code, you'd store them in these
    // ghost layers.) Then we can finish the boundary planes in the stencil
    // calculation.
    if (up_rank != MPI_PROC_NULL) {
        int zghost = 0;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int local_idx = index3D(x, y, zghost, Nx, Ny, Nz_with_halo);
                int plane_idx = y * Nx + x;
                data[local_idx] = recv_up[plane_idx];
            }
        }
    }
    if (down_rank != MPI_PROC_NULL) {
        int zghost = local_nz + 1;
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int local_idx = index3D(x, y, zghost, Nx, Ny, Nz_with_halo);
                int plane_idx = y * Nx + x;
                data[local_idx] = recv_down[plane_idx];
            }
        }
    }

    // Now finalize the boundary planes z=1 and z=local_nz with the 3-point
    // stencil, including the newly received ghost cells at z=0 or z=local_nz+1.
    // For minimal demonstration, let's handle those edges:
    // Update z=1 plane (which depends on z=0)
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
    // Update z=local_nz plane (which depends on z=local_nz+1)
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

    // Optional: print or gather some stats to verify correctness
    // For demonstration, let's just print rank 0's corner values
    if (rank == 0) {
        int corner_idx = index3D(0, 0, 1, Nx, Ny, Nz_with_halo);
        std::cout << "Rank " << rank
                  << ": result corner = " << result[corner_idx] << "\n";
    }

    MPI_Finalize();
    return 0;
}