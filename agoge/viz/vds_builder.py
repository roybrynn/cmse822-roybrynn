"""
@file vds_builder.py
@brief Builds an HDF5 Virtual Data Set (VDS) that aggregates output files from MPI ranks.
Uses h5py's VDS API to automatically create VDS files for all fields.
"""

import os
import h5py
import numpy as np
import glob

def get_metadata(file_path):
    with h5py.File(file_path, 'r') as f:
        grid = f["/grid"]
        local_dims = np.array(grid.attrs["local_dimensions"], dtype=int)
        subdomain = np.array(grid.attrs["subdomain"], dtype=int)
    return local_dims, subdomain

def get_available_fields(file_path):
    """
    Detect all available field datasets in the first rank file.
    
    Parameters:
        file_path (str): Path to a representative HDF5 file
        
    Returns:
        list: Names of available field datasets
    """
    available_fields = []
    with h5py.File(file_path, 'r') as f:
        # Get all keys from the root of the file
        for key in f.keys():
            # Skip the "grid" group or any other groups - we only want datasets
            if not isinstance(f[key], h5py.Dataset):
                continue
            available_fields.append(key)
    
    return available_fields

def build_vds_for_field(input_dir, output_file, field, global_dims, global_bbox, cell_size, rank_files):
    """
    Build a VDS file for a specific field.
    
    Parameters:
        input_dir (str): Directory containing rank files
        output_file (str): Output VDS file path
        field (str): Field name to map
        global_dims (ndarray): Global domain dimensions
        global_bbox (ndarray): Global domain bounding box
        cell_size (ndarray): Cell size information
        rank_files (list): List of rank file paths
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Create virtual layout for shape (Nz, Ny, Nx)
        layout = h5py.VirtualLayout(shape=(global_dims[2], global_dims[1], global_dims[0]), dtype='float64')
        
        # Iterate over rank files and map them into the layout
        for file in rank_files:
            local_dims, subdomain = get_metadata(file)
            offset = (subdomain[0] * local_dims[0],
                     subdomain[1] * local_dims[1],
                     subdomain[2] * local_dims[2])
            
            # Check if this field exists in this file
            with h5py.File(file, 'r') as f:
                if field not in f:
                    print(f"Warning: Field '{field}' not found in {file}. Skipping.")
                    continue
            
            # Create virtual source
            vsource = h5py.VirtualSource(file, field, shape=(local_dims[2], local_dims[1], local_dims[0]))
            layout[offset[2]:offset[2]+local_dims[2],
                 offset[1]:offset[1]+local_dims[1],
                 offset[0]:offset[0]+local_dims[0]] = vsource
        
        # Create the file if it doesn't exist or open it if it does
        mode = 'a' if os.path.exists(output_file) else 'w'
        with h5py.File(output_file, mode, libver='latest') as f:
            # Create virtual dataset
            f.create_virtual_dataset(field, layout)
            
            # Add grid metadata if it doesn't exist yet
            if "/grid" not in f:
                grid = f.create_group("grid")
                grid.attrs.create("global_bbox", global_bbox)
                grid.attrs.create("global_dimensions", global_dims)
                grid.attrs.create("cell_size", cell_size)
                
                # Calculate and store coordinate arrays
                xmin, xmax, ymin, ymax, zmin, zmax = global_bbox
                dx, dy, dz = cell_size
                
                # Create coordinate datasets for cell centers
                x_coords = np.array([xmin + (i + 0.5) * dx for i in range(global_dims[0])])
                y_coords = np.array([ymin + (j + 0.5) * dy for j in range(global_dims[1])])
                z_coords = np.array([zmin + (k + 0.5) * dz for k in range(global_dims[2])])
                
                grid.create_dataset("x_coords", data=x_coords)
                grid.create_dataset("y_coords", data=y_coords)
                grid.create_dataset("z_coords", data=z_coords)
                
                # Store coordinate edges for visualization
                x_edges = np.zeros(len(x_coords) + 1)
                y_edges = np.zeros(len(y_coords) + 1)
                z_edges = np.zeros(len(z_coords) + 1)
                
                # Calculate edge coordinates
                x_edges[:-1] = x_coords - 0.5 * dx
                x_edges[-1] = x_coords[-1] + 0.5 * dx
                
                y_edges[:-1] = y_coords - 0.5 * dy
                y_edges[-1] = y_coords[-1] + 0.5 * dy
                
                z_edges[:-1] = z_coords - 0.5 * dz
                z_edges[-1] = z_coords[-1] + 0.5 * dz
                
                grid.create_dataset("x_edges", data=x_edges)
                grid.create_dataset("y_edges", data=y_edges)
                grid.create_dataset("z_edges", data=z_edges)
        
        return True
    except Exception as e:
        print(f"Error creating VDS for field '{field}': {str(e)}")
        return False

def build_vds(input_dir, output_vds_file, field=None):
    """
    Build a VDS file for all available fields or a specific field.
    
    Parameters:
        input_dir (str): Directory containing rank files
        output_vds_file (str): Output VDS file path
        field (str, optional): Field name to map. If None, all available fields are mapped.
    """
    # Find rank files
    file_pattern = os.path.join(input_dir, "rank_*.h5")
    rank_files = glob.glob(file_pattern)
    if not rank_files:
        raise RuntimeError("No rank HDF5 files found in {}".format(input_dir))
    
    # Read global metadata from the first file
    with h5py.File(rank_files[0], 'r') as f:
        grid = f["/grid"]
        global_dims = np.array(grid.attrs["global_dimensions"], dtype=int)
        global_bbox = np.array(grid.attrs["global_bbox"], dtype=float)
        cell_size = np.array(grid.attrs["cell_size"], dtype=float)
    
    # Determine which fields to process
    if field is not None:
        # Process only the specified field
        fields = [field]
    else:
        # Auto-detect available fields from the first rank file
        fields = get_available_fields(rank_files[0])
        if not fields:
            raise RuntimeError("No field datasets found in {}".format(rank_files[0]))
    
    # Build VDS for each field
    successful_fields = []
    for field_name in fields:
        print(f"Building VDS for field: {field_name}")
        success = build_vds_for_field(
            input_dir, 
            output_vds_file, 
            field_name,
            global_dims,
            global_bbox,
            cell_size,
            rank_files
        )
        if success:
            successful_fields.append(field_name)
    
    if successful_fields:
        print(f"VDS created at: {output_vds_file} with fields: {', '.join(successful_fields)}")
    else:
        print(f"Failed to create VDS for any fields")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Build HDF5 Virtual Data Set from MPI rank files.")
    parser.add_argument("--input-dir", required=True, help="Directory containing output rank files.")
    parser.add_argument("--output", default="global_vds.h5", help="Output VDS file name")
    parser.add_argument("--field", default=None, help="Field name to map (default: all available fields)")
    args = parser.parse_args()
    build_vds(args.input_dir, args.output, args.field)
