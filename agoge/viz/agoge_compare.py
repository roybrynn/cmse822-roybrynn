#!/usr/bin/env python3
"""
Compare two Agoge simulation outputs to quantify differences between them.
This program leverages the AgogeOutput class from agoge_viz module for data loading 
and visualization.
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add parent directory to path to ensure we can import agoge_viz
script_dir = Path(__file__).resolve().parent
if script_dir not in sys.path:
    sys.path.append(str(script_dir))

from agoge_viz import AgogeOutput

def compare_fields(dir1, dir2, field="rho", visualize=False, output_file=None, 
                  axis='z', index=None, log_scale=False):
    """
    Compare the specified field between two Agoge output directories.
    
    Parameters:
        dir1 (str): Path to first Agoge output directory
        dir2 (str): Path to second Agoge output directory
        field (str): Field name to compare (default: "rho")
        visualize (bool): Whether to generate visualization plots
        output_file (str): Path to save comparison results (optional)
        axis (str): Slicing axis for visualization ('x', 'y', 'z')
        index (int): Slice index for visualization (default: middle of domain)
        log_scale (bool): Use logarithmic scale for error plots
        
    Returns:
        tuple: (absolute_error, relative_error) arrays
    """
    print(f"Loading data from {dir1} and {dir2}...")
    
    # Create AgogeOutput objects for both directories
    output1 = AgogeOutput(dir1)
    output2 = AgogeOutput(dir2)
    
    # Load the specified field from both outputs
    field1 = output1.load_field(field)
    field2 = output2.load_field(field)
    
    # Check if the fields have the same shape
    if field1.shape != field2.shape:
        print(f"Warning: Fields have different shapes: {field1.shape} vs {field2.shape}")
        # Find the minimum dimensions for comparison
        min_shape = tuple(min(d1, d2) for d1, d2 in zip(field1.shape, field2.shape))
        print(f"Truncating to common shape: {min_shape}")
        # Truncate arrays to match minimum shape
        field1 = field1[:min_shape[0], :min_shape[1], :min_shape[2]]
        field2 = field2[:min_shape[0], :min_shape[1], :min_shape[2]]
    
    # Compute absolute error
    abs_error = np.abs(field1 - field2)
    
    # Compute relative error with epsilon to avoid division by zero
    epsilon = 1e-12
    rel_error = abs_error / (np.maximum(np.abs(field1), np.abs(field2)) + epsilon)
    
    # Print statistics
    print(f"Comparison of '{field}' between:")
    print(f"  - {dir1}")
    print(f"  - {dir2}")
    print(f"Shape: {field1.shape}")
    print(f"Statistics:")
    print(f"  Max absolute error:  {np.max(abs_error):.6e}")
    print(f"  Mean absolute error: {np.mean(abs_error):.6e}")
    print(f"  Max relative error:  {np.max(rel_error):.6e}")
    print(f"  Mean relative error: {np.mean(rel_error):.6e}")
    
    # Save comparison results to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(f"Comparison of '{field}' between {dir1} and {dir2}\n")
            f.write(f"Shape: {field1.shape}\n")
            f.write(f"Max absolute error:  {np.max(abs_error):.6e}\n")
            f.write(f"Mean absolute error: {np.mean(abs_error):.6e}\n")
            f.write(f"Max relative error:  {np.max(rel_error):.6e}\n")
            f.write(f"Mean relative error: {np.mean(rel_error):.6e}\n")
        print(f"Results saved to {output_file}")
    
    # Visualize if requested
    if visualize:
        # Use the built-in comparison method from AgogeOutput
        fig, axes = output1.compare_with(output2, field, axis=axis, index=index, log_scale=log_scale)
        plt.show()
    
    # Free memory by unloading fields
    output1.unload_all_fields()
    output2.unload_all_fields()
    
    return abs_error, rel_error

def compare_multiple_fields(dir1, dir2, fields=None, visualize=False, output_dir=None, **kwargs):
    """
    Compare multiple fields between two Agoge output directories.
    
    Parameters:
        dir1 (str): Path to first Agoge output directory
        dir2 (str): Path to second Agoge output directory
        fields (list): List of field names to compare (default: ["rho", "E", "phi"])
        visualize (bool): Whether to generate visualization plots
        output_dir (str): Directory to save comparison results
        **kwargs: Additional arguments passed to compare_fields
    """
    if fields is None:
        # Create AgogeOutput temporarily to get available fields
        temp_output = AgogeOutput(dir1)
        available_fields = temp_output.list_fields()
        # Default to standard hydrodynamic fields if available
        standard_fields = ["rho", "E", "phi"]
        fields = [f for f in standard_fields if f in available_fields]
        if not fields:
            # Use whatever fields are available
            fields = available_fields
        temp_output.unload_all_fields()
        
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    results = {}
    for field in fields:
        print(f"\nComparing field: {field}")
        
        output_file = None
        if output_dir:
            output_file = os.path.join(output_dir, f"{field}_comparison.txt")
            
        abs_error, rel_error = compare_fields(
            dir1, dir2, field=field, 
            visualize=visualize, 
            output_file=output_file,
            **kwargs
        )
        
        results[field] = {
            'abs_error': abs_error,
            'rel_error': rel_error,
            'max_abs': np.max(abs_error),
            'mean_abs': np.mean(abs_error),
            'max_rel': np.max(rel_error),
            'mean_rel': np.mean(rel_error)
        }
        
    # Generate summary if multiple fields
    if len(fields) > 1 and output_dir:
        summary_file = os.path.join(output_dir, "summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Summary comparison between {dir1} and {dir2}\n")
            f.write(f"{'Field':<10} {'Max Abs Error':<20} {'Mean Abs Error':<20} {'Max Rel Error':<20} {'Mean Rel Error':<20}\n")
            f.write("-" * 90 + "\n")
            
            for field in fields:
                res = results[field]
                f.write(f"{field:<10} {res['max_abs']:<20.6e} {res['mean_abs']:<20.6e} {res['max_rel']:<20.6e} {res['mean_rel']:<20.6e}\n")
                
        print(f"\nSummary saved to {summary_file}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Compare two Agoge output directories to quantify differences."
    )
    parser.add_argument("dir1", help="First Agoge output directory")
    parser.add_argument("dir2", help="Second Agoge output directory")
    parser.add_argument("--fields", nargs='+', 
                      help="Fields to compare (default: all available fields)")
    parser.add_argument("--visualize", action="store_true", 
                      help="Generate visualizations of the fields and their differences")
    parser.add_argument("--output-dir", 
                      help="Directory to save comparison results")
    parser.add_argument("--axis", default="z", choices=["x", "y", "z"], 
                      help="Slicing axis for visualization")
    parser.add_argument("--index", type=int, 
                      help="Slice index for visualization (default: middle of domain)")
    parser.add_argument("--log-scale", action="store_true",
                      help="Use logarithmic scale for error plots")
    parser.add_argument("--threshold", type=float, default=1e-12,
                      help="Threshold for considering errors significant (default: 1e-12)")
    
    args = parser.parse_args()
    
    # Compare the specified fields
    results = compare_multiple_fields(
        args.dir1, args.dir2,
        fields=args.fields,
        visualize=args.visualize,
        output_dir=args.output_dir,
        axis=args.axis,
        index=args.index,
        log_scale=args.log_scale
    )
    
    # Determine the overall success/failure based on the maximum error
    max_error = max(res['max_abs'] for res in results.values())
    if max_error <= args.threshold:
        print(f"\n✅ Tests PASSED: Maximum error ({max_error:.6e}) is below threshold ({args.threshold:.6e})")
        sys.exit(0)
    else:
        print(f"\n❌ Tests FAILED: Maximum error ({max_error:.6e}) exceeds threshold ({args.threshold:.6e})")
        sys.exit(1)

if __name__ == "__main__":
    main()
