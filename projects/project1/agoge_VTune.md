# Profiling agoge's GravityCollapse Problem on Intel18 Nodes using Intel VTune

This tutorial guides you through profiling the **agoge** application's GravityCollapse problem on the Intel18 nodes of the MSU HPCC using Intel® VTune™ Profiler. You will profile two key components separately:
- **EulerSolver**
- **GravitySolver**

After profiling, you will extract performance metrics to compute the arithmetic intensity (AI) and place your data on a roofline model plot for the Intel18 nodes.

---

## Part I. Preparation

### 1. Accessing the Intel18 Nodes

Log in to the `dev-intel18` node using SSH or via OnDemand. For example, from your terminal, log in to the HPCC gateway node then:
```bash
ssh your_username@dev-intel18
```

2. Load the Required Modules

Before building or profiling, load the modules for HDF5 and VTun. Doing this will automatically load the necessary compilers and software stack. On dev-intel18, use:
 
```bash
module load HDF5
module load VTune
```

1. Build the agoge Code

Ensure that your code is built with both optimization and debug symbols so VTune can map performance data to your source lines. Check the `CXXFLAGS` in the Makefile to verify, then build agoge on intel18:

```bash
agoge> make clean; make 
```

This should produce the `agoge_run` binary. 

## Part II. Profiling with VTune CLI

We will use VTune’s HPC Performance Characterization analysis to capture floating-point operation counts (GFLOPS) and memory bandwidth data. 

### 1. Profiling the EulerSolver

a. Run VTune Collection for EulerSolver

To run with _just_ the finite-difference EulerSolver, set the folloing in GravityCollapse.yaml:

```yaml
use_gravity: false
do_euler_update: true
```

Then, to run a profiling session, invoke the following:

```bash
vtune -collect hpc-performance -result-dir vtune_euler_result -- ./agoge_run problems/GravityCollapse.yaml 
```

 • Explanation:
 • -collect hpc-performance: Collects HPC performance data (FLOP counts, memory traffic, etc.).
 • -result-dir vtune_euler_result: Directory where results will be saved.
 

b. Generate a Summary Report

After collection, generate a summary report:

```bash
vtune -report summary -result-dir vtune_euler_result
```

This report will show overall performance metrics such as DP GFLOPS and average DRAM bandwidth. You should get something like this:

```vtune: Executing actions 75 % Generating a report                              Elapsed Time: 17.590s
    SP GFLOPS: 0.000
    DP GFLOPS: 1.820
    x87 GFLOPS: 0.000
    CPI Rate: 0.999
    Average CPU Frequency: 3.218 GHz
    Total Thread Count: 1
Effective CPU Utilization: 2.4%
 | The metric value is low, which may signal a poor logical CPU cores
 | utilization caused by load imbalance, threading runtime overhead, contended
 | synchronization, or thread/process underutilization. Explore sub-metrics to
 | estimate the efficiency of MPI and OpenMP parallelism or run the Locks and
 | Waits analysis to identify parallel bottlenecks for other parallel runtimes.
 |
    Average Effective CPU Utilization: 0.942 out of 40
Memory Bound: 8.7% of Pipeline Slots
    Cache Bound: 18.8% of Clockticks
    DRAM Bound: 13.2% of Clockticks
        DRAM Bandwidth Bound: 0.0% of Elapsed Time
    NUMA: % of Remote Accesses: 0.0%

    Bandwidth Utilization
    Bandwidth Domain             Platform Maximum  Observed Maximum  Average  % of Elapsed Time with High BW Utilization(%)
    ---------------------------  ----------------  ----------------  -------  ---------------------------------------------
    DRAM, GB/sec                 218                         35.900   10.159                                           0.0%
    DRAM Single-Package, GB/sec  109                         32.900    9.111                                           0.0%
Vectorization: 1.8% of Packed FP Operations
 | A significant fraction of floating point arithmetic instructions are scalar.
 | This indicates that the code was not fully vectorized. Use Intel Advisor to
 | see possible reasons why the code was not vectorized.
 |
    Instruction Mix
        SP FLOPs: 0.0% of uOps
            Packed: 0.0% from SP FP
                128-bit: 0.0% from SP FP
                256-bit: 0.0% from SP FP
                512-bit: 0.0% from SP FP
            Scalar: 0.0% from SP FP
        DP FLOPs: 27.5% of uOps
            Packed: 1.8% from DP FP
                128-bit: 0.0% from DP FP
                256-bit: 1.8% from DP FP
                512-bit: 0.0% from DP FP
            Scalar: 98.2% from DP FP
             | A significant fraction of floating point arithmetic instructions
             | are scalar. This indicates that the code was not fully
             | vectorized. Use Intel Advisor to see possible reasons why the
             | code was not vectorized.
             |
        x87 FLOPs: 0.0% of uOps
        Non-FP: 72.5% of uOps
    FP Arith/Mem Rd Instr. Ratio: 1.373
    FP Arith/Mem Wr Instr. Ratio: 3.974
Collection and Platform Info
    Application Command Line: ./agoge_run "problems/GravityCollapse.yaml"
    Operating System: 5.15.0-126-generic DISTRIB_ID=Ubuntu DISTRIB_RELEASE=22.04 DISTRIB_CODENAME=jammy DISTRIB_DESCRIPTION="Ubuntu 22.04.5 LTS"
    Computer Name: dev-intel18
    Result Size: 30.0 MB
    Collection start time: 20:08:17 04/02/2025 UTC
    Collection stop time: 20:08:35 04/02/2025 UTC
    Collector Type: Driverless Perf per-process sampling
    CPU
        Name: Intel(R) Xeon(R) Processor code named Skylake
        Frequency: 2.394 GHz
        Logical CPU Count: 40
        Max DRAM Single-Package Bandwidth: 109.000 GB/s
        LLC size: 28.8 MB
```

Repeat this profiling for different problem sizes: $32^3$, $64^3$, $128^3$, $256^3$, and $512^3$. Note that you should decrease the `sound_crossings` parameter by a factor of about 16 for each increased resolution in order to keep the wallclock runtime roughly constant.

c. Extract Key Metrics

From the report, note:
 • DP GFLOPS: The double-precision performance (e.g., “DP GFLOPS: 4.215”).
 • Average DRAM Bandwidth: The average memory traffic (e.g., “DRAM, GB/sec: Average 2.706”).

These metrics will be used to compute the arithmetic intensity (AI).

### 2. Profiling the GravitySolver

Repeat the procedure for GravitySolver. In the yaml file, set:

```yaml
use_gravity: true
do_euler_update: false
```

Then run the profiling for the same problem sizes as above: $32^3$, $64^3$, $128^3$, $256^3$, and $512^3$. E.g.,

```bash
vtune -collect hpc-performance -result-dir vtune_gravity_result -- ./agoge_run problems/GravityCollapse.yaml 
```

Then, generate its summary report:

```bash
vtune -report summary -result-dir vtune_gravity_result
```

Record the key metrics (DP GFLOPS and DRAM bandwidth) from this report.

## Part III. Calculating and Plotting Arithmetic Intensity

### 1. Compute the Arithmetic Intensity (AI)

Arithmetic intensity is defined as:
$$
\text{AI} = \frac{\text{FLOPs/s}}{\text{Memory Traffic (Bytes/s)}}
$$

For each solver:
 • Example for EulerSolver:
 • DP GFLOPS ≈ 4.215 × 10⁹ FLOPs/s.
 • DRAM Bandwidth ≈ 2.706 × 10⁹ Bytes/s.
Then:
$$
\text{AI}_{\text{Euler}} \approx \frac{4.215 \times 10^9}{2.706 \times 10^9} \approx 1.56 \text{ FLOPs/Byte}
$$
 • Repeat for GravitySolver:
Use its respective metrics to compute AI.

### 2. Plotting on a Roofline Model

Even though VTune CLI does not directly produce a roofline plot, you can use the data gathered from the Empirical_Roofline_Tool used in Part 1 of Project 1. Here is an example python script for plotting your own roofline model using the data from ERT.

```Python 
import numpy as np
import matplotlib.pyplot as plt

def plot_roofline(l1_bw=191.6, l2_bw=86.7, dram_bw=41.8, peak_perf=39.3):
    """
    Plots an Empirical Roofline Model with customizable L1, L2, and DRAM bandwidths, 
    as well as peak performance.
    
    Parameters:
        l1_bw (float): L1 cache bandwidth in GB/s.
        l2_bw (float): L2 cache bandwidth in GB/s.
        dram_bw (float): DRAM bandwidth in GB/s.
        peak_perf (float): Peak floating point performance in GFLOPs/sec.
    """
    # Define FLOP/Byte operational intensity range
    intensity = np.logspace(-2, 2, 100)  # Log scale from 0.01 to 100

    # Compute roofline limits
    l1_limit = np.minimum(intensity * l1_bw, peak_perf)
    l2_limit = np.minimum(intensity * l2_bw, peak_perf)
    dram_limit = np.minimum(intensity * dram_bw, peak_perf)

    # Plot setup
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("FLOPs / Byte")
    ax.set_ylabel("GFLOPs / sec")
    ax.set_title("Empirical Roofline Model")

    # Plot the roofline bands
    ax.plot(intensity, l1_limit, label=f"L1 - {l1_bw} GB/s", linestyle='-', linewidth=2, color='r')
    ax.plot(intensity, l2_limit, label=f"L2 - {l2_bw} GB/s", linestyle='-', linewidth=2, color='g')
    ax.plot(intensity, dram_limit, label=f"DRAM - {dram_bw} GB/s", linestyle='-', linewidth=2, color='b')

    # Peak Performance Line
    ax.axhline(y=peak_perf, color='black', linestyle='--', linewidth=2, label=f"Peak Perf: {peak_perf} GFLOPs/s")

    ax.scatter(0.2990692864529472, 1.446, s=100, label="Example Point")
    
    ax.scatter(1.9717563989408649, 4.468, s=100, label="Example Point")

    # Grid and legend
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend()
    
    # Show plot
    plt.show()

    return ax 

# Example usage with default values extracted from PDF
ax = plot_roofline()
```

Adjust peak_bandwidth and compute_bound based on the actual hardware specifications of the Intel18 nodes.
Place points on this plot corresponding to your measured AI and peak FLOP rate for the EulerSolver and GravitySolver at the various resolutions. Comment on your results w.r.t. the roofline model. Does this match your expectations? Why or why not?
