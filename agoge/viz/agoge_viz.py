"""
agoge_viz.py
============

A simple Python module for reading and plotting Agoge simulation outputs.

Requires:
  - h5py (for HDF5 IO)
  - matplotlib (for plotting)
  - numpy (for array manipulation)

Example usage in a Jupyter notebook:
------------------------------------
    import agoge_viz as av

    # Read data from "agoge_final.h5"
    data = av.read_agoge_hdf5("agoge_final.h5")

    # data is a dict with fields: ["rho", "rhou", "rhov", "rhow", "E", "phi"]
    # plus metadata like data["domain_dimensions"], data["cell_size"], etc.

    # Make a 1D line plot of density through the domain along X at j=0, k=0
    av.plot_line(data, field="rho", axis="x", index_j=0, index_k=0)

    # Make a 2D slice plot of density in the z=0 plane
    av.plot_slice(data, field="rho", plane="z", plane_index=0)
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

def read_agoge_hdf5(filename):
    """
    Read Agoge output from an HDF5 file.

    Returns a dictionary containing:
      - 'rho', 'rhou', 'rhov', 'rhow', 'E', 'phi': 3D numpy arrays
      - 'domain_dimensions': tuple (Nx, Ny, Nz)
      - 'cell_size': tuple (dx, dy, dz)
      - 'bounding_box': array shape (3,2), each row is [min, max] for x,y,z
    """
    data = {}
    with h5py.File(filename, "r") as f:
        # Access the '/grid' group
        grid = f.get("/grid")
        if grid is None:
            raise KeyError("Group '/grid' does not exist in the HDF5 file.")

        # Read attributes instead of datasets
        try:
            domain_dims = np.array(grid.attrs["domain_dimensions"])
            Nx, Ny, Nz = domain_dims
            data["domain_dimensions"] = (Nx, Ny, Nz)
        except KeyError:
            raise KeyError("Attribute 'domain_dimensions' not found in '/grid'.")

        try:
            cell_size = np.array(grid.attrs["cell_size"])
            dx, dy, dz = cell_size
            data["cell_size"] = (dx, dy, dz)
        except KeyError:
            raise KeyError("Attribute 'cell_size' not found in '/grid'.")

        try:
            bbox = np.array(grid.attrs["bounding_box"])  # shape (3,2)
            data["bounding_box"] = bbox
        except KeyError:
            raise KeyError("Attribute 'bounding_box' not found in '/grid'.")

        # Read fields: stored in shape (Nz, Ny, Nx)
        def read_field(name):
            arr = np.array(f[name])  # shape (Nz,Ny,Nx)
            return arr

        data["rho"]  = read_field("rho")
        data["rhou"] = read_field("rhou")
        data["rhov"] = read_field("rhov")
        data["rhow"] = read_field("rhow")
        data["E"]    = read_field("E")
        data["phi"]  = read_field("phi")

    return data



def plot_line(data, field="rho", axis="x", index_j=0, index_k=0, figsize=(6,4), ax=None):
    """
    Plot a 1D line of 'field' data along a given axis ('x','y','z'),
    holding the other two indices fixed.

    data: dictionary from read_agoge_hdf5
    field: which variable to plot (e.g. "rho", "phi", etc.)
    axis: one of "x","y","z"
    index_j, index_k: specify which slice or line. Meaning depends on axis:
        if axis='x', then we fix j=index_j, k=index_k
        if axis='y', then we fix i=index_j, k=index_k
        if axis='z', then we fix i=index_j, j=index_k
      (this is a bit hacky, but it's a simple flexible approach.)
    """
    arr = data[field]  # shape (Nz, Ny, Nx)
    Nx, Ny, Nz = data["domain_dimensions"]
    dx, dy, dz = data["cell_size"]
    bbox = data["bounding_box"]

    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()

    # xarr -> coordinates for plot
    if axis == "x":
        # gather line arr[k=j,k=index_k, i=0..Nx-1]
        line = arr[index_k, index_j, :]  # watch dimension ordering
        # dimension ordering: arr is (z,y,x) => arr[z, y, x]
        # so index_k => z, index_j => y, axis => x
        xvals = np.linspace(bbox[0,0], bbox[0,1], Nx)
        label_axis = "x"
    elif axis == "y":
        # line = arr[z=index_k, y=0..Ny-1, x=index_j]
        # not super intuitive, but let's do it carefully
        line = arr[index_k, :, index_j]
        xvals = np.linspace(bbox[1,0], bbox[1,1], Ny)
        label_axis = "y"
    elif axis == "z":
        # line = arr[z=0..Nz-1, y=index_k, x=index_j]
        line = arr[:, index_k, index_j]
        xvals = np.linspace(bbox[2,0], bbox[2,1], Nz)
        label_axis = "z"
    else:
        raise ValueError(f"Unknown axis '{axis}'")

    ax.plot(xvals, line, '-x', label=field)
    ax.set_xlabel(f"{label_axis}-coordinate")
    ax.set_ylabel(field)
    ax.set_title(f"Line plot of {field} along axis={axis}")
    plt.legend()
    plt.grid(True)
    


def plot_slice(data, field="rho", plane="z", plane_index=0, figsize=(6,5),
               cmap="viridis"):
    """
    Plot a 2D slice of 'field' from data. The slice is plane='x','y','z' at plane_index
    plane='z' => arr[z=plane_index, :, :]
    plane='y' => arr[:, y=plane_index, :]
    plane='x' => arr[:, :, x=plane_index]
    """
    arr = data[field]  # shape (Nz,Ny,Nx)
    Nx, Ny, Nz = data["domain_dimensions"]
    dx, dy, dz = data["cell_size"]
    bbox = data["bounding_box"]

    if plane == "z":
        if plane_index < 0 or plane_index >= Nz:
            raise ValueError(f"plane_index out of range for z: {plane_index}")
        slice2d = arr[plane_index,:,:]  # shape (Ny,Nx)
        extent = [bbox[0,0], bbox[0,1], bbox[1,0], bbox[1,1]]  # x from 0..Nx, y from 0..Ny
        xlab = "x"
        ylab = "y"
    elif plane == "y":
        if plane_index < 0 or plane_index >= Ny:
            raise ValueError(f"plane_index out of range for y: {plane_index}")
        # shape => (Nz, Nx)
        # We reorder so y is fixed => slice2d[z, x]
        # But let's do shape => (z,x)
        slice2d = arr[:,plane_index,:]  # shape (Nz, Nx)
        extent = [bbox[0,0], bbox[0,1], bbox[2,0], bbox[2,1]]  # x=0..Nx, z=0..Nz
        xlab = "x"
        ylab = "z"
    elif plane == "x":
        if plane_index < 0 or plane_index >= Nx:
            raise ValueError(f"plane_index out of range for x: {plane_index}")
        # shape => (Nz, Ny)
        # x is fixed => slice2d[z,y]
        slice2d = arr[:,:,plane_index]  # shape (Ny=?), we want shape (Nz, Ny?)
        # Actually shape => (Nz, Ny)
        # We'll interpret the horizontal as y and vertical as z
        extent = [bbox[1,0], bbox[1,1], bbox[2,0], bbox[2,1]]  # y=0..Ny, z=0..Nz
        xlab = "y"
        ylab = "z"
    else:
        raise ValueError(f"Unknown plane: {plane} (expected x,y,z)")

    # For matplotlib's imshow, we typically do imshow(slice2d, extent, origin='lower') 
    # but note slice2d shape => (height, width)
    # We'll assume shape => (Ny, Nx) or something similar => we must be careful about orientation.

    plt.figure(figsize=figsize)
    plt.imshow(slice2d, cmap=cmap, extent=extent, origin='lower', aspect='auto')
    cbar = plt.colorbar()
    cbar.set_label(field)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(f"Slice of {field} on plane={plane}, index={plane_index}")
    plt.show()


# Extensions:
#   - You can add specialized plot functions for derived quantities (e.g., velocity magnitude)
#   - Add 3D rendering or isosurface with e.g. Mayavi, etc.
#   - For chunked reading on large data, you could add lazy/partial read logic.
