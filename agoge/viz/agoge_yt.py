# agoge_yt.py

import h5py
import numpy as np
import yt
from IPython.display import display

def load_agoge_hdf5(filename):
    """
    Load Agoge HDF5 data and convert it into a yt dataset.

    Parameters:
        filename (str): Path to the HDF5 file.

    Returns:
        yt.Dataset: A yt dataset object.
    """
    # Open the HDF5 file
    with h5py.File(filename, 'r') as f:
        # Access the '/grid' group
        if '/grid' not in f:
            raise KeyError("The HDF5 file does not contain a '/grid' group.")

        grid = f['/grid']
        
        # Read attributes from the '/grid' group
        try:
            domain_dims = grid.attrs['domain_dimensions']  # [Nx, Ny, Nz]
            cell_size = grid.attrs['cell_size']            # [dx, dy, dz]
            bbox = grid.attrs['bounding_box']              # [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
        except KeyError as e:
            raise KeyError(f"Missing attribute {e} in '/grid' group.")

        Nx, Ny, Nz = domain_dims
        dx, dy, dz = cell_size
        xmin, xmax = bbox[0]
        ymin, ymax = bbox[1]
        zmin, zmax = bbox[2]

        # Read datasets
        required_fields = ['rho', 'rhou', 'rhov', 'rhow', 'E', 'phi']
        data = {}
        for field in required_fields:
            if field not in f:
                raise KeyError(f"Dataset '{field}' not found in the HDF5 file.")
            data[field] = f[field][:].transpose((2, 1, 0))  # (Nz, Ny, Nx) -> (Nx, Ny, Nz)

        # Compute velocity components
        with np.errstate(divide='ignore', invalid='ignore'):
            velocity_x = np.where(data['rho'] != 0, data['rhou'] / data['rho'], 0)
            velocity_y = np.where(data['rho'] != 0, data['rhov'] / data['rho'], 0)
            velocity_z = np.where(data['rho'] != 0, data['rhow'] / data['rho'], 0)

        data['velocity_x'] = velocity_x
        data['velocity_y'] = velocity_y
        data['velocity_z'] = velocity_z

    # Define the grid dimensions and physical dimensions
    domain = np.array([Nx, Ny, Nz], dtype=np.int32)
    bbox_np = np.array([
        [xmin, xmax],
        [ymin, ymax],
        [zmin, zmax]
    ])

    # Load data into yt's uniform grid
    ds = yt.load_uniform_grid(
        data,
        domain,
        length_unit="code_length",  # Replace with actual units if known, e.g., "cm", "m"
        bbox=bbox_np,
        nprocs=1,
        geometry="cartesian"
    )

    return ds

def visualize_agoge_data(filename, output_prefix="output"):
    """
    Load Agoge HDF5 data and create visualizations using yt.

    Parameters:
        filename (str): Path to the HDF5 file.
        output_prefix (str): Prefix for the output image files.
    """
    ds = load_agoge_hdf5(filename)
    print(ds)

    # Slice Plot for Density
    rho_slice = yt.SlicePlot(ds, 'z', 'rho')
    rho_slice.set_xlabel('x')
    rho_slice.set_ylabel('y')
    rho_slice.save(f"{output_prefix}_rho_slice_z.png")

    # Projection Plot for Energy
    E_projection = yt.ProjectionPlot(ds, 'z', 'E')
    E_projection.set_xlabel('x')
    E_projection.set_ylabel('y')
    E_projection.save(f"{output_prefix}_E_projection_z.png")

    # Vector Slice Plot for Velocity
    velocity_slice = yt.SlicePlot(ds, 'z', ['velocity_x', 'velocity_y', 'velocity_z'])
    velocity_slice.annotate_velocity()
    velocity_slice.set_xlabel('x')
    velocity_slice.set_ylabel('y')
    velocity_slice.save(f"{output_prefix}_velocity_slice_z.png")

    print("Visualizations saved successfully.")
