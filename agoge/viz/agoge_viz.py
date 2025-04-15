"""
@file agoge_viz.py
@brief Object-oriented visualization and analysis framework for Agoge outputs.

This module provides classes for loading, analyzing, and visualizing Agoge simulation data.
It supports both individual snapshots and time series analysis.
"""

import glob
import h5py
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse
import importlib.util
import sys
from pathlib import Path


class AgogeOutput:
    """
    Class representing a single Agoge output snapshot.
    
    This class loads data from an Agoge output directory (either directly from
    rank files or through a VDS file) and provides methods for visualization
    and analysis.
    """
    
    def __init__(self, directory, build_vds=True):
        """
        Initialize by loading an Agoge output directory.
        
        Parameters:
            directory (str): Path to the Agoge output directory
            build_vds (bool): Whether to build a VDS file if it doesn't exist
        """
        self.directory = directory
        self.fields = {}
        self.coords = {}
        self.meta = {}
        self._loaded_fields = set()
        self._vds_file = os.path.join(directory, "global_vds.h5")
        self._vds_builder = None
        
        # Try to load or create the VDS file
        if build_vds:
            self._ensure_vds_file()
    
    def _ensure_vds_file(self):
        """Ensure that a VDS file exists for this output directory."""
        # Check if VDS file already exists
        if not os.path.exists(self._vds_file):
            # Check for rank files
            rank_files = glob.glob(os.path.join(self.directory, "rank_*.h5"))
            if not rank_files:
                raise RuntimeError(f"No rank files found in {self.directory}")
            
            # Import VDS builder dynamically
            if self._vds_builder is None:
                vds_builder_path = os.path.join(os.path.dirname(__file__), "vds_builder.py")
                if os.path.exists(vds_builder_path):
                    # Import the module dynamically
                    spec = importlib.util.spec_from_file_location("vds_builder", vds_builder_path)
                    self._vds_builder = importlib.util.module_from_spec(spec)
                    sys.modules["vds_builder"] = self._vds_builder
                    spec.loader.exec_module(self._vds_builder)
                else:
                    raise ImportError(f"VDS builder module not found at {vds_builder_path}")
            
            # Build the VDS file (for all fields)
            print(f"Creating VDS file for {self.directory}...")
            self._vds_builder.build_vds(self.directory, self._vds_file)
            print(f"VDS file created at {self._vds_file}")
    
    def _load_metadata(self):
        """Load metadata from VDS file."""
        if not self.meta:
            with h5py.File(self._vds_file, 'r') as f:
                grid = f["/grid"]
                
                # Load attributes
                for attr_name in grid.attrs:
                    self.meta[attr_name] = grid.attrs[attr_name]
                
                # Get simulation time if available, default to 0.0 if not
                self.simulation_time = grid.attrs.get("simulation_time", 0.0)
                
                # Load coordinate arrays
                if "x_coords" in grid:
                    self.coords["x"] = np.array(grid["x_coords"])
                if "y_coords" in grid:
                    self.coords["y"] = np.array(grid["y_coords"])
                if "z_coords" in grid:
                    self.coords["z"] = np.array(grid["z_coords"])
                if "x_edges" in grid:
                    self.coords["x_edges"] = np.array(grid["x_edges"])
                if "y_edges" in grid:
                    self.coords["y_edges"] = np.array(grid["y_edges"])
                if "z_edges" in grid:
                    self.coords["z_edges"] = np.array(grid["z_edges"])
                
                # If edge coordinates are missing, compute them
                if "x" in self.coords and "x_edges" not in self.coords:
                    x = self.coords["x"]
                    dx = x[1] - x[0] if len(x) > 1 else 1.0
                    x_edges = np.zeros(len(x) + 1)
                    x_edges[:-1] = x - 0.5 * dx
                    x_edges[-1] = x[-1] + 0.5 * dx
                    self.coords["x_edges"] = x_edges
                
                if "y" in self.coords and "y_edges" not in self.coords:
                    y = self.coords["y"]
                    dy = y[1] - y[0] if len(y) > 1 else 1.0
                    y_edges = np.zeros(len(y) + 1)
                    y_edges[:-1] = y - 0.5 * dy
                    y_edges[-1] = y[-1] + 0.5 * dy
                    self.coords["y_edges"] = y_edges
                
                if "z" in self.coords and "z_edges" not in self.coords:
                    z = self.coords["z"]
                    dz = z[1] - z[0] if len(z) > 1 else 1.0
                    z_edges = np.zeros(len(z) + 1)
                    z_edges[:-1] = z - 0.5 * dz
                    z_edges[-1] = z[-1] + 0.5 * dz
                    self.coords["z_edges"] = z_edges
                
                # Also store available field names
                self.available_fields = [key for key in f.keys() if key != "grid"]
    
    def get_simulation_time(self):
        """Return the simulation time of this output."""
        if not self.meta:
            self._load_metadata()
        return self.simulation_time

    def list_fields(self):
        """Return a list of available fields in this output, including derived fields."""
        if not hasattr(self, 'available_fields'):
            self._load_metadata()
        
        # Return both available fields and supported derived fields
        derived_fields = self.list_available_derived_fields()
        return self.available_fields + derived_fields
    
    def load_field(self, field_name):
        """
        Load a specific field from the VDS file. If the field is not directly available,
        try to compute it as a derived field.
        
        Parameters:
            field_name (str): Name of the field to load
            
        Returns:
            ndarray: The loaded or computed field data
        """
        if field_name in self.fields:
            return self.fields[field_name]
        
        # Ensure metadata is loaded
        if not self.meta:
            self._load_metadata()
        
        # Try to load the field from the VDS file
        try:
            with h5py.File(self._vds_file, 'r') as f:
                if field_name in f:
                    self.fields[field_name] = np.array(f[field_name])
                    self._loaded_fields.add(field_name)
                    return self.fields[field_name]
                
            # If field not found in file, try to compute as a derived field
            return self.compute_derived_field(field_name)
        except (KeyError, ValueError):
            # Try to compute as a derived field
            try:
                return self.compute_derived_field(field_name)
            except ValueError:
                raise KeyError(f"Field '{field_name}' not found in {self._vds_file} and is not a known derived field")
    
    def unload_field(self, field_name):
        """
        Unload a field from memory to save space.
        
        Parameters:
            field_name (str): Name of the field to unload
        """
        if field_name in self.fields:
            del self.fields[field_name]
            self._loaded_fields.remove(field_name)
    
    def unload_all_fields(self):
        """Unload all fields from memory to save space."""
        self.fields.clear()
        self._loaded_fields.clear()
    
    def plot_slice(self, field_name, axis='z', index=None, vmin=None, vmax=None, 
                   title=None, figsize=(8, 6), cmap="viridis"):
        """
        Plot a slice of a field along the specified axis.
        
        Parameters:
            field_name (str): Name of the field to plot
            axis (str): Axis along which to slice ('x', 'y', or 'z')
            index (int, optional): Index of the slice. If None, middle of the domain is used.
            vmin, vmax (float, optional): Range for colormap
            title (str, optional): Title for the plot
            figsize (tuple): Figure size (width, height) in inches
            cmap (str): Colormap name
            
        Returns:
            tuple: (fig, ax) matplotlib figure and axis objects
        """
        # Load the field if not already loaded
        if field_name not in self.fields:
            self.load_field(field_name)
        
        # Get field and coordinate data
        field = self.fields[field_name]
        x = self.coords.get("x", None)
        y = self.coords.get("y", None)
        z = self.coords.get("z", None)
        x_edges = self.coords.get("x_edges", None)
        y_edges = self.coords.get("y_edges", None)
        z_edges = self.coords.get("z_edges", None)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set default title if none provided
        if title is None:
            title = f"{field_name} (t = {self.get_simulation_time():.4f})"
        
        if axis.lower() == 'z':
            if index is None:
                index = field.shape[0] // 2
            
            slice_data = field[index, :, :]
            
            # Use provided edge coordinates if available, otherwise compute them
            if x_edges is None or y_edges is None:
                # Create mesh grid for pcolormesh - needs cell edges, not centers
                x_edges_local = np.zeros(len(x) + 1)
                y_edges_local = np.zeros(len(y) + 1)
                
                # Compute cell edges assuming uniform grid spacing
                dx = x[1] - x[0] if len(x) > 1 else 1.0
                dy = y[1] - y[0] if len(y) > 1 else 1.0
                
                # Extend coordinates to edges
                x_edges_local[:-1] = x - 0.5 * dx
                x_edges_local[-1] = x[-1] + 0.5 * dx
                y_edges_local[:-1] = y - 0.5 * dy
                y_edges_local[-1] = y[-1] + 0.5 * dy
            else:
                x_edges_local = x_edges
                y_edges_local = y_edges
            
            mesh = ax.pcolormesh(x_edges_local, y_edges_local, slice_data, 
                                cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
            
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_title(f"{title} at z = {z[index]:.4f}")
            fig.colorbar(mesh, ax=ax, label=field_name)
            ax.set_aspect('equal')
            
        elif axis.lower() == 'y':
            # Similar implementation for y-axis slicing
            if index is None:
                index = field.shape[1] // 2
            
            slice_data = field[:, index, :]
            
            if x_edges is None or z_edges is None:
                x_edges_local = np.zeros(len(x) + 1)
                z_edges_local = np.zeros(len(z) + 1)
                
                dx = x[1] - x[0] if len(x) > 1 else 1.0
                dz = z[1] - z[0] if len(z) > 1 else 1.0
                
                x_edges_local[:-1] = x - 0.5 * dx
                x_edges_local[-1] = x[-1] + 0.5 * dx
                z_edges_local[:-1] = z - 0.5 * dz
                z_edges_local[-1] = z[-1] + 0.5 * dz
            else:
                x_edges_local = x_edges
                z_edges_local = z_edges
            
            mesh = ax.pcolormesh(x_edges_local, z_edges_local, slice_data,
                               cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
            
            ax.set_xlabel("x")
            ax.set_ylabel("z")
            ax.set_title(f"{title} at y = {y[index]:.4f}")
            fig.colorbar(mesh, ax=ax, label=field_name)
            ax.set_aspect('equal')
        
        elif axis.lower() == 'x':
            # Similar implementation for x-axis slicing
            if index is None:
                index = field.shape[2] // 2
            
            slice_data = field[:, :, index]
            
            if y_edges is None or z_edges is None:
                y_edges_local = np.zeros(len(y) + 1)
                z_edges_local = np.zeros(len(z) + 1)
                
                dy = y[1] - y[0] if len(y) > 1 else 1.0
                dz = z[1] - z[0] if len(z) > 1 else 1.0
                
                y_edges_local[:-1] = y - 0.5 * dy
                y_edges_local[-1] = y[-1] + 0.5 * dy
                z_edges_local[:-1] = z - 0.5 * dz
                z_edges_local[-1] = z[-1] + 0.5 * dz
            else:
                y_edges_local = y_edges
                z_edges_local = z_edges
            
            mesh = ax.pcolormesh(y_edges_local, z_edges_local, slice_data,
                               cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
            
            ax.set_xlabel("y")
            ax.set_ylabel("z")
            ax.set_title(f"{title} at x = {x[index]:.4f}")
            fig.colorbar(mesh, ax=ax, label=field_name)
            ax.set_aspect('equal')
            
        else:
            raise ValueError("Axis must be one of 'x', 'y', or 'z'")
        
        plt.tight_layout()
        return fig, ax
    
    def plot_line(self, field_name, axis='x', index_y=None, index_z=None, 
                  figsize=(8, 6), color='blue', linewidth=2, marker=None):
        """
        Plot a 1D line through the field along a coordinate axis.
        
        Parameters:
            field_name (str): Name of the field to plot
            axis (str): Axis along which to plot ('x', 'y', or 'z')
            index_y, index_z (int): Indices for other dimensions
            figsize (tuple): Figure size (width, height) in inches
            color (str): Line color
            linewidth (float): Line width
            marker (str): Marker style
        
        Returns:
            tuple: (fig, ax) matplotlib figure and axis objects
        """
        # Load the field if not already loaded
        if field_name not in self.fields:
            self.load_field(field_name)
        
        # Get field and coordinate data
        field = self.fields[field_name]
        x = self.coords.get("x", None)
        y = self.coords.get("y", None)
        z = self.coords.get("z", None)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        if axis.lower() == 'x':
            if index_y is None:
                index_y = field.shape[1] // 2
            if index_z is None:
                index_z = field.shape[0] // 2
            
            line_data = field[index_z, index_y, :]
            ax.plot(x, line_data, color=color, linewidth=linewidth, marker=marker)
            ax.set_xlabel("x")
            ax.set_title(f"{field_name} along x-axis at y={y[index_y]:.4f}, z={z[index_z]:.4f}")
            
        elif axis.lower() == 'y':
            if index_z is None:
                index_z = field.shape[0] // 2
            if index_x is None:
                index_x = field.shape[2] // 2
                
            line_data = field[index_z, :, index_x]
            ax.plot(y, line_data, color=color, linewidth=linewidth, marker=marker)
            ax.set_xlabel("y")
            ax.set_title(f"{field_name} along y-axis at x={x[index_x]:.4f}, z={z[index_z]:.4f}")
            
        elif axis.lower() == 'z':
            if index_y is None:
                index_y = field.shape[1] // 2
            if index_x is None:
                index_x = field.shape[2] // 2
                
            line_data = field[:, index_y, index_x]
            ax.plot(z, line_data, color=color, linewidth=linewidth, marker=marker)
            ax.set_xlabel("z")
            ax.set_title(f"{field_name} along z-axis at x={x[index_x]:.4f}, y={y[index_y]:.4f}")
            
        else:
            raise ValueError("Axis must be one of 'x', 'y', or 'z'")
            
        ax.set_ylabel(field_name)
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        return fig, ax
    
    def compare_with(self, other_output, field_name, axis='z', index=None, 
                    log_scale=False, figsize=(15, 5)):
        """
        Compare this output with another output for a specific field.
        
        Parameters:
            other_output (AgogeOutput): Another AgogeOutput to compare with
            field_name (str): Field name to compare
            axis (str): Axis for slicing ('x', 'y', 'z')
            index (int): Index for slicing
            log_scale (bool): Use logarithmic scale for error plot
            figsize (tuple): Figure size (width, height) in inches
            
        Returns:
            tuple: (fig, axes) matplotlib figure and axes objects
        """
        # Load fields if needed
        if field_name not in self.fields:
            self.load_field(field_name)
        if field_name not in other_output.fields:
            other_output.load_field(field_name)
        
        field1 = self.fields[field_name]
        field2 = other_output.fields[field_name]
        
        # Check dimensions
        if field1.shape != field2.shape:
            print(f"Warning: Fields have different shapes: {field1.shape} vs {field2.shape}")
            # Find the minimum dimensions
            min_shape = tuple(min(d1, d2) for d1, d2 in zip(field1.shape, field2.shape))
            print(f"Truncating to common shape: {min_shape}")
            # Truncate to match minimum shape
            field1 = field1[:min_shape[0], :min_shape[1], :min_shape[2]]
            field2 = field2[:min_shape[0], :min_shape[1], :min_shape[2]]
        
        # Compute absolute error
        abs_error = np.abs(field1 - field2)
        
        # Compute relative error with epsilon to avoid division by zero
        epsilon = 1e-12
        rel_error = abs_error / (np.maximum(np.abs(field1), np.abs(field2)) + epsilon)
        
        # Create a three-panel figure
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        # Plot field1
        if axis.lower() == 'z':
            if index is None:
                index = field1.shape[0] // 2
                
            # First field
            slice1 = field1[index, :, :]
            mesh1 = axes[0].pcolormesh(
                self.coords.get("x_edges", None), 
                self.coords.get("y_edges", None), 
                slice1, 
                shading='flat'
            )
            axes[0].set_title(f"{field_name} (First)")
            axes[0].set_aspect('equal')
            fig.colorbar(mesh1, ax=axes[0])
            
            # Second field
            slice2 = field2[index, :, :]
            mesh2 = axes[1].pcolormesh(
                other_output.coords.get("x_edges", None), 
                other_output.coords.get("y_edges", None), 
                slice2, 
                shading='flat'
            )
            axes[1].set_title(f"{field_name} (Second)")
            axes[1].set_aspect('equal')
            fig.colorbar(mesh2, ax=axes[1])
            
            # Error
            if log_scale:
                error_to_plot = np.log10(abs_error[index, :, :] + epsilon)
                error_title = f"log10(|Δ{field_name}|)"
            else:
                error_to_plot = abs_error[index, :, :]
                error_title = f"|Δ{field_name}|"
                
            mesh_err = axes[2].pcolormesh(
                self.coords.get("x_edges", None), 
                self.coords.get("y_edges", None), 
                error_to_plot, 
                shading='flat'
            )
            axes[2].set_title(error_title)
            axes[2].set_aspect('equal')
            fig.colorbar(mesh_err, ax=axes[2])
            
        # Add similar cases for 'y' and 'x' axes (omitted for brevity)
            
        plt.tight_layout()
        
        # Print statistics
        print(f"Comparison of '{field_name}':")
        print(f"  Max absolute error:  {np.max(abs_error):.6e}")
        print(f"  Mean absolute error: {np.mean(abs_error):.6e}")
        print(f"  Max relative error:  {np.max(rel_error):.6e}")
        print(f"  Mean relative error: {np.mean(rel_error):.6e}")
        
        return fig, axes

    def compute_derived_field(self, field_name):
        """
        Compute a derived field based on the field name.
        
        Parameters:
            field_name (str): Name of the derived field to compute
            
        Returns:
            ndarray: The computed derived field
            
        Supported derived fields:
            - u: x-velocity (rhou/rho)
            - v: y-velocity (rhov/rho)
            - w: z-velocity (rhow/rho)
            - vel_magnitude: Velocity magnitude sqrt(u^2 + v^2 + w^2)
        """
        # Check if it's already computed
        if field_name in self.fields:
            return self.fields[field_name]
            
        # Handle velocity components
        if field_name == 'u':
            # Load necessary fields if not already loaded
            if 'rhou' not in self.fields:
                self.load_field('rhou')
            if 'rho' not in self.fields:
                self.load_field('rho')
            
            # Avoid division by zero
            epsilon = 1e-12
            u = self.fields['rhou'] / np.maximum(self.fields['rho'], epsilon)
            self.fields[field_name] = u
            self._loaded_fields.add(field_name)
            return u
            
        elif field_name == 'v':
            if 'rhov' not in self.fields:
                self.load_field('rhov')
            if 'rho' not in self.fields:
                self.load_field('rho')
            
            epsilon = 1e-12
            v = self.fields['rhov'] / np.maximum(self.fields['rho'], epsilon)
            self.fields[field_name] = v
            self._loaded_fields.add(field_name)
            return v
            
        elif field_name == 'w':
            if 'rhow' not in self.fields:
                self.load_field('rhow')
            if 'rho' not in self.fields:
                self.load_field('rho')
            
            epsilon = 1e-12
            w = self.fields['rhow'] / np.maximum(self.fields['rho'], epsilon)
            self.fields[field_name] = w
            self._loaded_fields.add(field_name)
            return w
            
        elif field_name == 'vel_magnitude':
            # Compute velocity magnitude
            u = self.compute_derived_field('u')
            v = self.compute_derived_field('v')
            w = self.compute_derived_field('w')
            
            vel_mag = np.sqrt(u**2 + v**2 + w**2)
            self.fields[field_name] = vel_mag
            self._loaded_fields.add(field_name)
            return vel_mag
        
        # Add more derived fields as needed here
        
        else:
            raise ValueError(f"Unknown derived field: {field_name}")
    
    def list_available_derived_fields(self):
        """
        Return a list of supported derived fields.
        
        Returns:
            list: Names of supported derived fields
        """
        return ['u', 'v', 'w', 'vel_magnitude']


class AgogeTimeSeries:
    """
    Class for managing and visualizing a time series of Agoge outputs.
    """
    
    def __init__(self, base_directory, pattern="*ranks_*"):
        """
        Initialize by scanning for Agoge output directories in the base directory.
        
        Parameters:
            base_directory (str): Base directory containing Agoge output directories
            pattern (str): Glob pattern to match output directories
        """
        self.base_directory = base_directory
        self.pattern = pattern
        self.outputs = []
        self.timestamps = []
        self._scan_directories()
    
    def _scan_directories(self):
        """Scan for output directories matching the pattern."""
        # Create full pattern
        full_pattern = os.path.join(self.base_directory, self.pattern)
        
        # Find matching directories
        matching_dirs = sorted(glob.glob(full_pattern))
        
        if not matching_dirs:
            print(f"Warning: No directories matching '{full_pattern}' found.")
            return
            
        print(f"Found {len(matching_dirs)} output directories.")
        
        # Extract timestamps from directory names and load simulation times
        self.simulation_times = []
        for dir_path in matching_dirs:
            # Try to extract a timestamp from the directory name
            try:
                dir_name = os.path.basename(dir_path)
                timestamp = int(dir_name.split('_')[-1])
                self.timestamps.append(timestamp)
                self.outputs.append(dir_path)
                
                # Try to load the simulation time from the output
                try:
                    # Create temporary AgogeOutput to get simulation time
                    temp_output = AgogeOutput(dir_path, build_vds=False)
                    sim_time = temp_output.get_simulation_time()
                    self.simulation_times.append(sim_time)
                except:
                    # If we can't load the simulation time, use the timestamp as a fallback
                    self.simulation_times.append(float(timestamp))
                    
            except (IndexError, ValueError):
                # If timestamp extraction fails, use the index as a substitute
                idx = len(self.timestamps)
                self.timestamps.append(idx)
                self.outputs.append(dir_path)
                self.simulation_times.append(float(idx))
    
    def load_output(self, index):
        """
        Load an output at the specified index.
        
        Parameters:
            index (int): Index of the output to load
            
        Returns:
            AgogeOutput: The loaded output
        """
        if index < 0 or index >= len(self.outputs):
            raise IndexError(f"Index {index} out of range for {len(self.outputs)} outputs")
            
        return AgogeOutput(self.outputs[index])
    
    def animate_field(self, field_name, axis='z', index=None, vmin=None, vmax=None,
                      figsize=(8, 6), cmap="viridis", interval=200, use_sim_time=True):
        """
        Create an animation of a field across all time steps.
        
        Parameters:
            field_name (str): Name of the field to animate
            axis (str): Axis along which to slice ('x', 'y', or 'z')
            index (int): Index for slicing
            vmin, vmax (float): Range for colormap
            figsize (tuple): Figure size (width, height) in inches
            cmap (str): Colormap name
            interval (int): Animation interval in milliseconds
            use_sim_time (bool): Whether to use simulation time (True) or output index (False) in title
            
        Returns:
            matplotlib.animation.FuncAnimation: Animation object
        """
        # Load first output to get metadata and coordinate system
        first_output = self.load_output(0)
        first_output.load_field(field_name)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        if axis.lower() == 'z':
            if index is None:
                index = first_output.fields[field_name].shape[0] // 2
                
            slice_data = first_output.fields[field_name][index, :, :]
            x_edges = first_output.coords.get("x_edges")
            y_edges = first_output.coords.get("y_edges")
            
            # Initial plot
            mesh = ax.pcolormesh(x_edges, y_edges, slice_data, cmap=cmap,
                               vmin=vmin, vmax=vmax, shading='flat')
            
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            title = ax.set_title(f"{field_name} at z = {first_output.coords['z'][index]:.4f}, t = {self.timestamps[0]}")
            fig.colorbar(mesh, ax=ax, label=field_name)
            ax.set_aspect('equal')
            
            # Define update function for animation
            def update(frame):
                output = self.load_output(frame)
                output.load_field(field_name)
                slice_data = output.fields[field_name][index, :, :]
                mesh.set_array(slice_data.ravel())
                
                if use_sim_time:
                    time_val = self.simulation_times[frame]
                else:
                    time_val = self.timestamps[frame]
                    
                title.set_text(f"{field_name} at z = {output.coords['z'][index]:.4f}, t = {time_val:.4f}")
                output.unload_all_fields()  # Free memory
                return mesh,
                
        # Add similar cases for 'y' and 'x' axes (omitted for brevity)
        
        # Create animation
        anim = FuncAnimation(fig, update, frames=len(self.outputs),
                          interval=interval, blit=True)
        
        plt.tight_layout()
        return anim
        
    def plot_time_evolution(self, field_name, stat='mean', location=None, figsize=(10, 6), 
                           use_sim_time=True):
        """
        Plot the time evolution of a field statistic or at a specific point location.
        
        Parameters:
            field_name (str): Name of the field
            stat (str): Statistic to plot ('max', 'min', 'mean', 'std', 'sum')
                        or 'point' to sample at a specific location
            location (tuple): (x_idx, y_idx, z_idx) for point sampling when stat='point'
            figsize (tuple): Figure size (width, height) in inches
            use_sim_time (bool): Whether to use simulation time (True) or output index (False) for x-axis
            
        Returns:
            tuple: (fig, ax) matplotlib figure and axis objects
        """
        values = []
        
        for i in range(len(self.outputs)):
            output = self.load_output(i)
            output.load_field(field_name)
            field = output.fields[field_name]
            
            if stat == 'max':
                values.append(np.max(field))
            elif stat == 'min':
                values.append(np.min(field))
            elif stat == 'mean':
                values.append(np.mean(field))
            elif stat == 'std':
                values.append(np.std(field))
            elif stat == 'sum':
                values.append(np.sum(field))
            elif stat == 'point':
                if location is None:
                    # Default to center of domain
                    z_idx = field.shape[0] // 2
                    y_idx = field.shape[1] // 2
                    x_idx = field.shape[2] // 2
                else:
                    z_idx, y_idx, x_idx = location
                values.append(field[z_idx, y_idx, x_idx])
            else:
                raise ValueError(f"Unknown statistic: {stat}")
            
            output.unload_all_fields()  # Free memory
        
        fig, ax = plt.subplots(figsize=figsize)
        
        x_values = self.simulation_times if use_sim_time else self.timestamps
        ax.plot(x_values, values, marker='o', linestyle='-', linewidth=2)
        ax.set_xlabel("Simulation Time" if use_sim_time else "Output Index", fontsize=12)
        
        if stat == 'point':
            if location:
                ax.set_ylabel(f"{field_name} at ({location[2]},{location[1]},{location[0]})", fontsize=12)
                ax.set_title(f"Time Evolution of {field_name} at Point ({location[2]},{location[1]},{location[0]})")
            else:
                ax.set_ylabel(f"{field_name} at center", fontsize=12)
                ax.set_title(f"Time Evolution of {field_name} at Domain Center")
        else:
            ax.set_ylabel(f"{stat.capitalize()} of {field_name}", fontsize=12)
            ax.set_title(f"Time Evolution of {stat.capitalize()} of {field_name}")
            
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        return fig, ax
        
    def compute_convergence(self, field_name, reference_index=-1, norm='L2'):
        """
        Compute the convergence history of a field relative to a reference solution.
        
        Parameters:
            field_name (str): Name of the field to analyze
            reference_index (int): Index of the reference solution (default: last output)
            norm (str): Norm to use for comparison ('L1', 'L2', 'Linf')
            
        Returns:
            tuple: (errors, convergence_rates)
        """
        if len(self.outputs) < 2:
            raise ValueError("Need at least two outputs to compute convergence")
            
        # Load the reference solution
        reference_output = self.load_output(reference_index)
        reference_output.load_field(field_name)
        reference_field = reference_output.fields[field_name]
        
        errors = []
        for i in range(len(self.outputs)):
            if i == reference_index:
                continue
            output = self.load_output(i)
            output.load_field(field_name)
            field = output.fields[field_name]
            
            if norm == 'L1':
                error = np.sum(np.abs(field - reference_field))
            elif norm == 'L2':
                error = np.sqrt(np.sum((field - reference_field)**2))
            elif norm == 'Linf':
                error = np.max(np.abs(field - reference_field))
            else:
                raise ValueError(f"Unknown norm: {norm}")
            
            errors.append(error)
            output.unload_all_fields()  # Free memory
        
        # Compute convergence rates
        convergence_rates = []
        for i in range(1, len(errors)):
            rate = np.log(errors[i] / errors[i-1]) / np.log(self.timestamps[i] / self.timestamps[i-1])
            convergence_rates.append(rate)
        
        return errors, convergence_rates


class AgogeViz:
    """
    Class for visualizing and analyzing Agoge outputs.
    """
    
    def __init__(self, outputs, compare_module=None):
        """
        Initialize AgogeViz with a list of outputs and an optional compare module.
        """
        self.outputs = outputs
        self.compare_module = compare_module
    
    def plot_slice(self, field_name, slice_axis='z', slice_index=0):
        """
        Simple slice plot along a given axis (x, y, or z).
        """
        # Example: load data for the specified field, select slice, and display
        pass

    def plot_line(self, field_name, axis='x', index_tuple=(0, 0)):
        """
        Plot a 1D line of the data across one axis at fixed other coordinates.
        """
        pass

    def compare_outputs(self, idx1, idx2, field_name):
        """
        Bitwise comparison of two outputs using the compare module if provided.
        """
        if self.compare_module:
            pass
        else:
            raise RuntimeError("No compare module supplied.")


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct and visualize global Agoge field data."
    )
    parser.add_argument("directory", nargs="?", help="Directory containing HDF5 files (or global_vds.h5)")
    parser.add_argument("--compare", nargs=2, metavar=("DIR1", "DIR2"), 
                        help="Compare global fields from two directories")
    parser.add_argument("--field", default="rho", help="Field name to visualize or compare")
    parser.add_argument("--axis", default="z", choices=["x", "y", "z"], help="Slicing axis")
    parser.add_argument("--index", type=int, default=None, help="Slice index")
    args = parser.parse_args()

    # New: If compare flag is provided, perform comparison and exit.
    if args.compare:
        output1 = AgogeOutput(args.compare[0])
        output2 = AgogeOutput(args.compare[1])
        output1.compare_with(output2, field_name=args.field, axis=args.axis, index=args.index)
        return

    if not args.directory:
        print("Error: Directory not specified.")
        return

    output = AgogeOutput(args.directory)
    output.plot_slice(field_name=args.field, axis=args.axis, index=args.index)

if __name__ == "__main__":
    main()