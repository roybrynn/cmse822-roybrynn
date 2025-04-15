#!/usr/bin/env python3
"""
Generate performance plots from agoge scaling study data.
This script creates various plots to visualize strong and weak scaling 
performance, parallel efficiency, and communication overhead.

Usage:
    python plot_performance.py [options]

Options:
    --input-dir <dir>       Directory containing CSV data files (default: '.')
    --output-dir <dir>      Directory for output plots (default: same as input-dir)
    --prefix <prefix>       Prefix for plot filenames (default: '')
    --comparative           Generate comparative plots (requires both blocking/nonblocking data)
    --format <format>       Plot file format (png, pdf, svg) (default: png)
    --no-display            Don't display plots, only save them
    --help                  Show this help message
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# Set high-quality plot defaults
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

def load_data(input_dir, prefix=""):
    """Load CSV files containing scaling data."""
    data = {}
    prefix = f"{prefix}_" if prefix else ""
    
    # Try to load strong scaling data
    strong_file = os.path.join(input_dir, f"{prefix}strong_scaling.csv")
    if os.path.exists(strong_file):
        data['strong'] = pd.read_csv(strong_file)
    
    # Try to load weak scaling data
    weak_file = os.path.join(input_dir, f"{prefix}weak_scaling.csv")
    if os.path.exists(weak_file):
        data['weak'] = pd.read_csv(weak_file)
    
    # Try to load efficiency data
    eff_file = os.path.join(input_dir, f"{prefix}scaling_efficiency.csv")
    if os.path.exists(eff_file):
        data['efficiency'] = pd.read_csv(eff_file)
    
    # Try to load combined data
    combined_file = os.path.join(input_dir, f"{prefix}combined_scaling.csv")
    if os.path.exists(combined_file):
        data['combined'] = pd.read_csv(combined_file)
    
    return data

def plot_strong_scaling(data, output_dir, prefix="", fmt='png', show=True):
    """Generate strong scaling performance plots."""
    if 'strong' not in data:
        print("No strong scaling data found.")
        return
    
    df = data['strong']
    output_prefix = f"{prefix}_" if prefix else ""
    
    # Group by resolution
    resolutions = df['GridSize'].unique()
    
    # Plot 1: Raw Performance (Zone Updates per Second)
    plt.figure()
    for res in resolutions:
        res_data = df[df['GridSize'] == res]
        plt.plot(res_data['Processors'], res_data['ZoneUpdates(M)'], 
                 marker='o', linestyle='-', label=f"{res}³")
    
    plt.xscale('log', base=2)
    plt.yscale('log', base=10)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.xlabel("Number of Processors")
    plt.ylabel("Zone Updates per Second (M)")
    plt.title("Strong Scaling: Raw Performance")
    plt.legend(title="Grid Size")
    
    # Add ideal scaling line from first point
    if len(resolutions) > 0:
        first_res = resolutions[0]
        first_data = df[df['GridSize'] == first_res]
        if not first_data.empty:
            p0 = first_data['Processors'].min()
            u0 = first_data[first_data['Processors'] == p0]['ZoneUpdates(M)'].values[0]
            x = np.array([p0, first_data['Processors'].max()])
            y = u0 * x / p0
            plt.plot(x, y, 'k--', alpha=0.5, label="Ideal Scaling")
            plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_prefix}strong_raw_performance.{fmt}"))
    if show:
        plt.show()
    else:
        plt.close()
    
    # Plot 2: Normalized Performance (Zone Updates per core-second)
    plt.figure()
    for res in resolutions:
        res_data = df[df['GridSize'] == res]
        plt.plot(res_data['Processors'], 
                 res_data['ZoneUpdates(M)'] / res_data['Processors'], 
                 marker='o', linestyle='-', label=f"{res}³")
    
    plt.xscale('log', base=2)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.xlabel("Number of Processors")
    plt.ylabel("Zone Updates per Core-Second (M)")
    plt.title("Strong Scaling: Normalized Performance")
    plt.legend(title="Grid Size")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_prefix}strong_normalized_performance.{fmt}"))
    if show:
        plt.show()
    else:
        plt.close()
    
    # Plot 3: Parallel Efficiency
    plt.figure()
    for res in resolutions:
        res_data = df[df['GridSize'] == res].sort_values('Processors')
        if len(res_data) > 0:
            # Use first point as reference
            base_proc = res_data['Processors'].iloc[0]
            base_perf = res_data['ZoneUpdates(M)'].iloc[0] * base_proc
            
            # Calculate efficiency
            efficiency = [(perf * proc) / base_perf 
                          for perf, proc in zip(res_data['ZoneUpdates(M)'], res_data['Processors'])]
            
            plt.plot(res_data['Processors'], efficiency, 
                     marker='o', linestyle='-', label=f"{res}³")
    
    plt.xscale('log', base=2)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    plt.xlabel("Number of Processors")
    plt.ylabel("Parallel Efficiency")
    plt.title("Strong Scaling: Parallel Efficiency")
    plt.legend(title="Grid Size")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_prefix}strong_parallel_efficiency.{fmt}"))
    if show:
        plt.show()
    else:
        plt.close()
    
    # Plot 4: Communication Overhead (if available)
    if 'CommunicationTime' in df.columns and 'TimeStepSeconds' in df.columns:
        plt.figure()
        for res in resolutions:
            res_data = df[df['GridSize'] == res]
            comm_overhead = res_data['CommunicationTime'].astype(float) / res_data['TimeStepSeconds'].astype(float)
            plt.plot(res_data['Processors'], comm_overhead, 
                    marker='o', linestyle='-', label=f"{res}³")
        
        plt.xscale('log', base=2)
        plt.grid(True, which="both", ls="--", alpha=0.7)
        plt.xlabel("Number of Processors")
        plt.ylabel("Communication Time / Total Time")
        plt.title("Strong Scaling: Communication Overhead")
        plt.legend(title="Grid Size")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{output_prefix}strong_communication_overhead.{fmt}"))
        if show:
            plt.show()
        else:
            plt.close()

def plot_weak_scaling(data, output_dir, prefix="", fmt='png', show=True):
    """Generate weak scaling performance plots."""
    if 'weak' not in data:
        print("No weak scaling data found.")
        return
    
    df = data['weak']
    output_prefix = f"{prefix}_" if prefix else ""
    
    # Group by baseline resolution
    baselines = df['Baseline'].unique()
    
    # Plot 1: Raw Performance (Zone Updates per Second)
    plt.figure()
    for baseline in baselines:
        baseline_data = df[df['Baseline'] == baseline]
        plt.plot(baseline_data['Processors'], baseline_data['ZoneUpdates(M)'], 
                 marker='o', linestyle='-', label=f"{baseline}³ per process")
    
    plt.xscale('log', base=2)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.xlabel("Number of Processors")
    plt.ylabel("Zone Updates per Second (M)")
    plt.title("Weak Scaling: Raw Performance")
    plt.legend(title="Zones per Process")
    
    # Add ideal scaling line from first point
    if len(baselines) > 0:
        first_baseline = baselines[0]
        first_data = df[df['Baseline'] == first_baseline]
        if not first_data.empty:
            p0 = first_data['Processors'].min()
            u0 = first_data[first_data['Processors'] == p0]['ZoneUpdates(M)'].values[0]
            x = np.array([p0, first_data['Processors'].max()])
            y = np.full_like(x, u0)  # Constant horizontal line for ideal case
            plt.plot(x, y, 'k--', alpha=0.5, label="Ideal Scaling")
            plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_prefix}weak_raw_performance.{fmt}"))
    if show:
        plt.show()
    else:
        plt.close()
    
    # Plot 2: Normalized Performance (Zone Updates per processor)
    plt.figure()
    for baseline in baselines:
        baseline_data = df[df['Baseline'] == baseline]
        plt.plot(baseline_data['Processors'], 
                 baseline_data['ZoneUpdates(M)'] / baseline_data['Processors'], 
                 marker='o', linestyle='-', label=f"{baseline}³ per process")
    
    plt.xscale('log', base=2)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.xlabel("Number of Processors")
    plt.ylabel("Zone Updates per Core-Second (M)")
    plt.title("Weak Scaling: Normalized Performance")
    plt.legend(title="Zones per Process")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_prefix}weak_normalized_performance.{fmt}"))
    if show:
        plt.show()
    else:
        plt.close()
    
    # Plot 3: Parallel Efficiency
    plt.figure()
    for baseline in baselines:
        baseline_data = df[df['Baseline'] == baseline].sort_values('Processors')
        if len(baseline_data) > 0:
            # Use first point as reference
            base_perf = baseline_data['ZoneUpdates(M)'].iloc[0]
            
            # Calculate efficiency
            efficiency = [perf / base_perf for perf in baseline_data['ZoneUpdates(M)']]
            
            plt.plot(baseline_data['Processors'], efficiency, 
                     marker='o', linestyle='-', label=f"{baseline}³ per process")
    
    plt.xscale('log', base=2)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    plt.xlabel("Number of Processors")
    plt.ylabel("Parallel Efficiency")
    plt.title("Weak Scaling: Parallel Efficiency")
    plt.legend(title="Zones per Process")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{output_prefix}weak_parallel_efficiency.{fmt}"))
    if show:
        plt.show()
    else:
        plt.close()
    
    # Plot 4: Communication Overhead (if available)
    if 'CommunicationTime' in df.columns and 'TimeStepSeconds' in df.columns:
        plt.figure()
        for baseline in baselines:
            baseline_data = df[df['Baseline'] == baseline]
            comm_overhead = baseline_data['CommunicationTime'].astype(float) / baseline_data['TimeStepSeconds'].astype(float)
            plt.plot(baseline_data['Processors'], comm_overhead, 
                    marker='o', linestyle='-', label=f"{baseline}³ per process")
        
        plt.xscale('log', base=2)
        plt.grid(True, which="both", ls="--", alpha=0.7)
        plt.xlabel("Number of Processors")
        plt.ylabel("Communication Time / Total Time")
        plt.title("Weak Scaling: Communication Overhead")
        plt.legend(title="Zones per Process")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{output_prefix}weak_communication_overhead.{fmt}"))
        if show:
            plt.show()
        else:
            plt.close()

def plot_comparative(blocking_dir, nonblocking_dir, output_dir, fmt='png', show=True):
    """Generate comparative plots between blocking and non-blocking communication."""
    # Load both datasets
    blocking_data = load_data(blocking_dir, prefix="blocking")
    nonblocking_data = load_data(nonblocking_dir, prefix="nonblocking")
    
    # Ensure we have data to compare
    if not blocking_data or not nonblocking_data:
        print("Insufficient data for comparative plots.")
        return
    
    # Compare strong scaling if available
    if 'strong' in blocking_data and 'strong' in nonblocking_data:
        # Extract common resolutions
        blocking_res = set(blocking_data['strong']['GridSize'].unique())
        nonblocking_res = set(nonblocking_data['strong']['GridSize'].unique())
        common_res = blocking_res.intersection(nonblocking_res)
        
        # Plot comparative performance
        if common_res:
            # Plot speedup: nonblocking over blocking
            plt.figure(figsize=(10, 6))
            
            for res in sorted(common_res):
                b_data = blocking_data['strong'][blocking_data['strong']['GridSize'] == res]
                nb_data = nonblocking_data['strong'][nonblocking_data['strong']['GridSize'] == res]
                
                # Merge on processor count
                merged = pd.merge(b_data, nb_data, on='Processors', suffixes=('_blocking', '_nonblocking'))
                
                # Calculate speedup
                speedup = merged['ZoneUpdates(M)_nonblocking'] / merged['ZoneUpdates(M)_blocking']
                
                plt.plot(merged['Processors'], speedup, marker='o', linestyle='-', 
                         label=f"Grid {res}³")
            
            plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
            plt.xscale('log', base=2)
            plt.grid(True, which="both", ls="--", alpha=0.7)
            plt.xlabel("Number of Processors")
            plt.ylabel("Speedup (Non-blocking / Blocking)")
            plt.title("Strong Scaling: Non-blocking vs. Blocking Speedup")
            plt.legend(title="Grid Size")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"comparative_strong_speedup.{fmt}"))
            if show:
                plt.show()
            else:
                plt.close()
            
            # Plot communication overhead comparison
            plt.figure(figsize=(10, 6))
            
            for res in sorted(common_res):
                b_data = blocking_data['strong'][blocking_data['strong']['GridSize'] == res]
                nb_data = nonblocking_data['strong'][nonblocking_data['strong']['GridSize'] == res]
                
                # Calculate communication overhead
                b_overhead = b_data['CommunicationTime'].astype(float) / b_data['TimeStepSeconds'].astype(float)
                nb_overhead = nb_data['CommunicationTime'].astype(float) / nb_data['TimeStepSeconds'].astype(float)
                
                # Find matching processor counts
                common_procs = set(b_data['Processors']).intersection(set(nb_data['Processors']))
                
                # Filter and sort
                b_filtered = pd.DataFrame({'Processors': b_data['Processors'], 'Overhead': b_overhead})
                b_filtered = b_filtered[b_filtered['Processors'].isin(common_procs)].sort_values('Processors')
                
                nb_filtered = pd.DataFrame({'Processors': nb_data['Processors'], 'Overhead': nb_overhead})
                nb_filtered = nb_filtered[nb_filtered['Processors'].isin(common_procs)].sort_values('Processors')
                
                # Plot both lines for this resolution
                plt.plot(b_filtered['Processors'], b_filtered['Overhead'], marker='o', linestyle='--', 
                         label=f"Blocking, {res}³")
                plt.plot(nb_filtered['Processors'], nb_filtered['Overhead'], marker='s', linestyle='-', 
                         label=f"Non-blocking, {res}³")
            
            plt.xscale('log', base=2)
            plt.grid(True, which="both", ls="--", alpha=0.7)
            plt.xlabel("Number of Processors")
            plt.ylabel("Communication Overhead")
            plt.title("Strong Scaling: Communication Overhead Comparison")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"comparative_strong_comm_overhead.{fmt}"))
            if show:
                plt.show()
            else:
                plt.close()
    
    # Compare weak scaling if available
    if 'weak' in blocking_data and 'weak' in nonblocking_data:
        # Extract common baselines
        blocking_baselines = set(blocking_data['weak']['Baseline'].unique())
        nonblocking_baselines = set(nonblocking_data['weak']['Baseline'].unique())
        common_baselines = blocking_baselines.intersection(nonblocking_baselines)
        
        # Plot comparative performance
        if common_baselines:
            # Plot speedup: nonblocking over blocking
            plt.figure(figsize=(10, 6))
            
            for baseline in sorted(common_baselines):
                b_data = blocking_data['weak'][blocking_data['weak']['Baseline'] == baseline]
                nb_data = nonblocking_data['weak'][nonblocking_data['weak']['Baseline'] == baseline]
                
                # Merge on processor count
                merged = pd.merge(b_data, nb_data, on='Processors', suffixes=('_blocking', '_nonblocking'))
                
                # Calculate speedup
                speedup = merged['ZoneUpdates(M)_nonblocking'] / merged['ZoneUpdates(M)_blocking']
                
                plt.plot(merged['Processors'], speedup, marker='o', linestyle='-', 
                         label=f"Baseline {baseline}³")
            
            plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
            plt.xscale('log', base=2)
            plt.grid(True, which="both", ls="--", alpha=0.7)
            plt.xlabel("Number of Processors")
            plt.ylabel("Speedup (Non-blocking / Blocking)")
            plt.title("Weak Scaling: Non-blocking vs. Blocking Speedup")
            plt.legend(title="Zones per Process")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"comparative_weak_speedup.{fmt}"))
            if show:
                plt.show()
            else:
                plt.close()
            
            # Plot communication overhead comparison
            plt.figure(figsize=(10, 6))
            
            for baseline in sorted(common_baselines):
                b_data = blocking_data['weak'][blocking_data['weak']['Baseline'] == baseline]
                nb_data = nonblocking_data['weak'][nonblocking_data['weak']['Baseline'] == baseline]
                
                # Calculate communication overhead
                b_overhead = b_data['CommunicationTime'].astype(float) / b_data['TimeStepSeconds'].astype(float)
                nb_overhead = nb_data['CommunicationTime'].astype(float) / nb_data['TimeStepSeconds'].astype(float)
                
                # Find matching processor counts
                common_procs = set(b_data['Processors']).intersection(set(nb_data['Processors']))
                
                # Filter and sort
                b_filtered = pd.DataFrame({'Processors': b_data['Processors'], 'Overhead': b_overhead})
                b_filtered = b_filtered[b_filtered['Processors'].isin(common_procs)].sort_values('Processors')
                
                nb_filtered = pd.DataFrame({'Processors': nb_data['Processors'], 'Overhead': nb_overhead})
                nb_filtered = nb_filtered[nb_filtered['Processors'].isin(common_procs)].sort_values('Processors')
                
                # Plot both lines for this baseline
                plt.plot(b_filtered['Processors'], b_filtered['Overhead'], marker='o', linestyle='--', 
                         label=f"Blocking, {baseline}³")
                plt.plot(nb_filtered['Processors'], nb_filtered['Overhead'], marker='s', linestyle='-', 
                         label=f"Non-blocking, {baseline}³")
            
            plt.xscale('log', base=2)
            plt.grid(True, which="both", ls="--", alpha=0.7)
            plt.xlabel("Number of Processors")
            plt.ylabel("Communication Overhead")
            plt.title("Weak Scaling: Communication Overhead Comparison")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"comparative_weak_comm_overhead.{fmt}"))
            if show:
                plt.show()
            else:
                plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Generate performance plots from agoge scaling study data."
    )
    parser.add_argument("--input-dir", default=".", 
                      help="Directory containing CSV data files (default: '.')")
    parser.add_argument("--output-dir", default=None, 
                      help="Directory for output plots (default: same as input-dir)")
    parser.add_argument("--prefix", default="", 
                      help="Prefix for plot filenames (default: '')")
    parser.add_argument("--comparative", action="store_true", 
                      help="Generate comparative plots (requires blocking/nonblocking data)")
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"], 
                      help="Plot file format (default: png)")
    parser.add_argument("--no-display", action="store_true", 
                      help="Don't display plots, only save them")
    
    args = parser.parse_args()
    
    # Use input_dir as output_dir if not specified
    output_dir = args.output_dir if args.output_dir else args.input_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Display settings
    show_plots = not args.no_display
    
    if args.comparative:
        # Look for blocking and nonblocking subdirectories
        blocking_dir = os.path.join(args.input_dir, "blocking")
        nonblocking_dir = os.path.join(args.input_dir, "nonblocking")
        
        # Check if the directories exist
        if os.path.isdir(blocking_dir) and os.path.isdir(nonblocking_dir):
            plot_comparative(blocking_dir, nonblocking_dir, output_dir, args.format, show_plots)
        else:
            print(f"Cannot find blocking and nonblocking directories in {args.input_dir}")
            print("Expected structure: {input_dir}/blocking and {input_dir}/nonblocking")
    else:
        # Load data from the specified directory
        data = load_data(args.input_dir, args.prefix)
        
        if not data:
            print(f"No scaling data found in {args.input_dir}")
            return
        
        # Generate plots
        if 'strong' in data:
            plot_strong_scaling(data, output_dir, args.prefix, args.format, show_plots)
        
        if 'weak' in data:
            plot_weak_scaling(data, output_dir, args.prefix, args.format, show_plots)
        
        print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    main()