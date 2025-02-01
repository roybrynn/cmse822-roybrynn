#!/usr/bin/env python3
"""
agoge_viz.py

Visualize Field3D data from Agoge HDF5 files.
The file stores fields (rho, rhou, rhov, rhow, E, phi) as 3D datasets and
includes spatial coordinate arrays (x, y, z) and grid metadata in the '/grid' group.

Usage:
    python agoge_viz.py <filename.h5> [--field rho|rhou|rhov|rhow|E|phi] [--slice axis index]

Optional:
    --slice: Provide an axis (x, y, or z) and an index for slicing.
    If yt is installed, an additional projection plot is generated.
"""

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Optionally import yt for volumetric rendering.
try:
    import yt
except ImportError:
    yt = None

def load_agoge_data(filename):
    data = {}
    with h5py.File(filename, 'r') as f:
        # Read field datasets from the file's root group.
        for field in ["rho", "rhou", "rhov", "rhow", "E", "phi"]:
            data[field] = np.array(f[field])
        # Read spatial coordinates and metadata from the '/grid' group.
        grid = f["/grid"]
        for coord in ["x", "y", "z"]:
            if coord in grid:
                data[coord] = np.array(grid[coord])
            else:
                data[coord] = None
        # Optionally, load additional attributes such as domain dimensions.
        if "domain_dimensions" in grid.attrs:
            data["domain_dimensions"] = grid.attrs["domain_dimensions"]
        if "cell_size" in grid.attrs:
            data["cell_size"] = grid.attrs["cell_size"]
    return data

def plot_field(field_data, x, y, z, axis='z', index=None, field_name="Field"):
    # Field data is stored as (Nz, Ny, Nx) (row-major order).
    if axis == 'z':
        if index is None:
            index = field_data.shape[0] // 2
        slice_data = field_data[index, :, :]
        X, Y = np.meshgrid(x, y, indexing='xy')
        plt.figure()
        plt.contourf(X, Y, slice_data, cmap='viridis')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"{field_name} at z-slice index {index}")
        plt.colorbar()
    elif axis == 'y':
        if index is None:
            index = field_data.shape[1] // 2
        slice_data = field_data[:, index, :]
        X, Z = np.meshgrid(x, z, indexing='xy')
        plt.figure()
        plt.contourf(X, Z, slice_data, cmap='viridis')
        plt.xlabel("x")
        plt.ylabel("z")
        plt.title(f"{field_name} at y-slice index {index}")
        plt.colorbar()
    elif axis == 'x':
        if index is None:
            index = field_data.shape[2] // 2
        slice_data = field_data[:, :, index]
        Y, Z = np.meshgrid(y, z, indexing='xy')
        plt.figure()
        plt.contourf(Y, Z, slice_data, cmap='viridis')
        plt.xlabel("y")
        plt.ylabel("z")
        plt.title(f"{field_name} at x-slice index {index}")
        plt.colorbar()
    else:
        raise ValueError("Axis must be one of 'x', 'y', or 'z'")
    plt.show()

def plot_line(x, data_list, labels, title="", xlabel="x", ylabel="Field Value", styles=None):
    """
    Plot 1D line data where x is common for each dataset.

    Parameters:
        x (array-like): 1D coordinate array.
        data_list (list of array-like): List of 1D data arrays to plot.
        labels (list of str): Labels for each line.
        title (str): Plot title.
        xlabel (str): x-axis label.
        ylabel (str): y-axis label.
        styles (list of str, optional): Matplotlib style formats (e.g., 'r-', 'bo--').
    """
    plt.figure(figsize=(8, 5))
    for i, data in enumerate(data_list):
        if styles is not None and i < len(styles):
            plt.plot(x, data, styles[i], linewidth=2, label=labels[i])
        else:
            plt.plot(x, data, linewidth=2, label=labels[i])
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Visualize Agoge HDF5 field data.")
    parser.add_argument("filename", help="Path to the HDF5 file")
    parser.add_argument("--field", default="rho",
                        choices=["rho", "rhou", "rhov", "rhow", "E", "phi"],
                        help="Select which field to visualize")
    parser.add_argument("--slice", nargs=2, metavar=("AXIS", "INDEX"),
                        help="Specify a slice: axis (x, y, or z) and index")
    args = parser.parse_args()

    data = load_agoge_data(args.filename)
    field_data = data[args.field]

    # Use coordinate arrays from file.
    if data["x"] is None or data["y"] is None or data["z"] is None:
        raise RuntimeError("Coordinate arrays not found in HDF5 file.")

    # Determine slicing parameters.
    if args.slice:
        axis = args.slice[0].lower()
        index = int(args.slice[1])
    else:
        axis = 'z'
        index = field_data.shape[0] // 2

    # Plot the selected field slice.
    plot_field(field_data, data["x"], data["y"], data["z"],
               axis=axis, index=index, field_name=args.field)

    # Advanced visualization using yt, if available.
    if yt is not None:
        try:
            ds = yt.load(args.filename)
            p = yt.ProjectionPlot(ds, 'z', args.field)
            proj_filename = f"projection_{args.field}.png"
            p.save(proj_filename)
            print(f"Projection plot saved to {proj_filename}")
        except Exception as e:
            print("yt encountered an error:", e)
    else:
        print("yt is not installed; skipping advanced visualization.")

if __name__ == "__main__":
    main()