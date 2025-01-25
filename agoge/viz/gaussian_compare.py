#!/usr/bin/env python3

"""
gaussian_compare.py
Compare the final state of a Gaussian pulse advected in the x-direction
to the analytic solution. Compute L1 and L2 errors in density.

Now uses functions from agoge_viz.py to read HDF5 results and optionally plot comparison.
"""

import sys
import math
import numpy as np
import h5py

# Added import from agoge_viz
# (assuming agoge_viz.py is in the same directory or PYTHONPATH)
import agoge_viz as av


def analytic_gaussian(nx, ny, nz, dx, dy, dz, amp, sigma, t, u0, base_rho=1.0):
    """
    Return a 3D numpy array of the analytic density at time t.
    The pulse was centered at (0.5*nx*dx, 0.5*ny*dy, 0.5*nz*dz), and advected
    in +x direction with velocity u0. We use periodic boundaries, assisted
    by a helper approach if needed.
    """
    xC = 0.5 * nx * dx
    yC = 0.5 * ny * dy
    zC = 0.5 * nz * dz

    # shift center in x
    xCenterShifted = xC + u0 * t
    Lx = nx * dx
    Ly = ny * dy
    Lz = nz * dz

    rho_ana = np.zeros((nz, ny, nx), dtype=np.float64)

    for k in range(nz):
        zpos = (k + 0.5) * dz
        # wrap zpos around [0, Lz]
        while zpos < 0.0:   zpos += Lz
        while zpos >= Lz:  zpos -= Lz

        for j in range(ny):
            ypos = (j + 0.5) * dy
            # wrap ypos
            while ypos < 0.0:   ypos += Ly
            while ypos >= Ly:  ypos -= Ly

            for i in range(nx):
                xpos = (i + 0.5) * dx
                # wrap xpos
                while xpos < 0.0:   xpos += Lx
                while xpos >= Lx:  xpos -= Lx

                # measure distance in x from shifted center, with periodic wrap
                dx_ = xpos - xCenterShifted
                while dx_ >  0.5 * Lx: dx_ -= Lx
                while dx_ < -0.5 * Lx: dx_ += Lx

                # measure distance in y from yC, with wrap
                dy_ = ypos - yC
                while dy_ >  0.5 * Ly: dy_ -= Ly
                while dy_ < -0.5 * Ly: dy_ += Ly

                # measure distance in z from zC, with wrap
                dz_ = zpos - zC
                while dz_ >  0.5 * Lz: dz_ -= Lz
                while dz_ < -0.5 * Lz: dz_ += Lz

                r2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_
                bump = amp * math.exp(-r2 / (2 * sigma * sigma))
                rho_ana[k, j, i] = base_rho + bump

    return rho_ana


def main():
    """
    Compare the final simulation output (from HDF5) to the analytic solution,
    then print L1 and L2 errors. Optionally plot with functions from agoge_viz.
    Usage:
      ./gaussian_compare.py <filename.h5> <amp> <sigma> <u0> <time>
    """
    if len(sys.argv) < 6:
        print("Usage: gaussian_compare.py <hdf5_file> <amp> <sigma> <u0> <time>")
        sys.exit(1)

    filename = sys.argv[1]
    amp      = float(sys.argv[2])
    sigma    = float(sys.argv[3])
    u0       = float(sys.argv[4])
    t        = float(sys.argv[5])

    # Use agoge_viz to read the HDF5 data
    data = av.read_agoge_hdf5(filename)
    rho_final = data["rho"]
    Nx, Ny, Nz = data["domain_dimensions"]
    dx, dy, dz = data["cell_size"]

    # Compute analytic
    rho_analytic = analytic_gaussian(Nx, Ny, Nz, dx, dy, dz, amp, sigma, t, u0)

    # compute L1 and L2 errors
    diff = rho_final - rho_analytic
    L1 = np.sum(np.abs(diff)) / (Nx * Ny * Nz)
    L2 = math.sqrt(np.sum(diff**2) / (Nx * Ny * Nz))

    print(f"L1 error= {L1:.6e}, L2 error= {L2:.6e}")

    # Optional: create a line plot of final vs. analytic along x at y=Ny//2,z=Nz//2
    # (just as an example - user can adjust indices)
    mid_j = Ny // 2
    mid_k = Nz // 2
    # plot final
    av.plot_line(data, field="rho", axis="x", index_j=mid_j, index_k=mid_k)
    # plot analytic
    # reuse the same approach by creating a "data_analytic" dictionary to pass
    data_analytic = {
        "rho": rho_analytic,
        "domain_dimensions": data["domain_dimensions"],
        "cell_size": data["cell_size"],
        "bounding_box": data["bounding_box"]
    }
    av.plot_line(data_analytic, field="rho", axis="x", index_j=mid_j, index_k=mid_k)

    av.plt.show()

if __name__ == "__main__":
    main()