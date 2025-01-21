# Comprehensive Guide to Profiling Matrix Multiplication with Intel VTune on MSU ICER HPCC

Welcome to this detailed tutorial on using Intel VTune Profiler from the command line to analyze the performance of the `matrix_multiply` application from the oneAPI samples on the MSU ICER High-Performance Computing Cluster (HPCC). This guide leverages the [oneAPI `matrix_multiply` sample](https://github.com/oneapi-src/oneAPI-samples/tree/master/Tools/VTuneProfiler) to provide a practical, hands-on experience in identifying and resolving performance bottlenecks using VTune. Additionally, this updated tutorial incorporates essential troubleshooting steps and detailed explanations of VTune's output to enhance your profiling workflow.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Accessing MSU ICER HPCC](#accessing-msu-icer-hpcc)
3. [Cloning the `matrix_multiply` Sample](#cloning-the-matrix_multiply-sample)
4. [Setting Up the Environment](#setting-up-the-environment)
5. [Compiling the `matrix_multiply` Application](#compiling-the-matrix_multiply-application)
6. [Running VTune Analysis from the Command Line](#running-vtune-analysis-from-the-command-line)
7. [Understanding and Resolving Common VTune Errors](#understanding-and-resolving-common-vtune-errors)
8. [Collecting and Viewing Results](#collecting-and-viewing-results)
9. [Analyzing Performance Metrics](#analyzing-performance-metrics)
10. [Advanced VTune Features](#advanced-vtune-features)
11. [Best Practices and Tips](#best-practices-and-tips)
12. [Conclusion](#conclusion)

---

## Prerequisites

Before you begin, ensure you have the following:

- **MSU ICER HPCC Access**: A valid user account with SSH access to the HPCC.
- **Intel VTune Profiler**: Installed on the HPCC. If not, contact the system administrator.
- **Git**: Installed for cloning repositories.
- **Basic Knowledge**: Familiarity with command-line operations and compiling code on HPC systems.
- **Access to `dev-intel18` Node**: Students should log in to the `dev-intel18` node for performing this tutorial.

## Accessing MSU ICER HPCC

### 1. **Establish a Secure Connection**

Log in to the `dev-intel18` node of the HPCC using SSH. Replace `username` with your actual username.

```bash
ssh username@icer.hpc.msu.edu
```

### 2. **Navigate to Your Working Directory**

Create and move to a directory where you want to clone the sample code.

```bash
mkdir -p ~/vtune_profiling
cd ~/vtune_profiling
```

## Cloning the `matrix_multiply` Sample

### 1. **Clone the Repository**

Use `git` to clone the oneAPI `matrix_multiply` sample repository.

```bash
git clone https://github.com/oneapi-src/oneAPI-samples.git
```

### 2. **Navigate to the VTune Profiler Tools Sample**

Change directory to the `VTuneProfiler` sample.

```bash
cd oneAPI-samples/Tools/VTuneProfiler/matrix_multiply
```

### 3. **Explore the Sample**

List the contents to familiarize yourself with the project structure.

```bash
ls
```

You should see files like `matrix_multiply.cpp`, `Makefile`, and other related scripts.

## Setting Up the Environment

### 1. **Load Necessary Modules**

HPCC environments typically use module systems such as `Environment Modules` or `Lmod`. Load the Intel oneAPI compiler and VTune modules. Replace module names if different on your system.

```bash
module load intel/oneapi/2025.0
module load vtune/2025.0
```

### 2. **Verify Module Loading**

Ensure that the modules are loaded correctly.

```bash
module list
```

You should see entries for Intel oneAPI and VTune.

### 3. **Set Environment Variables** (if necessary)

Although loading modules usually sets the necessary environment variables, verify that `PATH` and `LD_LIBRARY_PATH` include VTune binaries and libraries.

```bash
echo $PATH
echo $LD_LIBRARY_PATH
```

If needed, manually export them:

```bash
export PATH=/opt/intel/oneapi/vtune/2025.0/bin64:$PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/vtune/2025.0/lib64:$LD_LIBRARY_PATH
```

## Compiling the `matrix_multiply` Application

Proper compilation is crucial for accurate profiling. The application should include debugging symbols and appropriate optimization flags.

### 1. **Review the Makefile**

Open the `Makefile` to ensure it includes necessary flags. Typically, the `Makefile` in the sample should already be configured, but verify the following:

```makefile
CXX = icpx
CXXFLAGS = -g -O2 -std=c++17
```

- `-g`: Includes debugging information.
- `-O2`: Enables a moderate level of optimization.
- `-std=c++17`: Specifies the C++ standard.

### 2. **Modify Compilation Flags for Compatibility**

Given the error encountered regarding CPU instruction sets, adjust the compilation flags to target a more generic architecture that doesn't assume Intel-specific instructions.

1. **Edit the Makefile**

   Open the `Makefile` in a text editor.

   ```bash
   nano Makefile
   ```

2. **Adjust `CXXFLAGS`**

   Replace `-march=native` and `-mtune=native` with flags that target a baseline architecture supporting the required instruction sets.

   **Recommended Flags:**

   ```makefile
   CXXFLAGS = -g -O2 -march=core2 -mtune=generic -std=c++17
   ```

   - `-march=core2`: Targets the Intel Core 2 architecture, which supports SSE, SSE2, and SSE3.
   - `-mtune=generic`: Optimizes the code for a generic CPU, improving portability.

   Alternatively, you can specify an AMD architecture or use lower optimization flags if needed.

3. **Save and Exit**

   If using `nano`, press `CTRL + O` to save and `CTRL + X` to exit.

### 3. **Clean Previous Builds**

Remove any existing compiled files to ensure a fresh build.

```bash
make clean
```

### 4. **Compile the Application**

Run `make` to build the `matrix_multiply` executable.

```bash
make
```

Upon successful compilation, an executable named `matrix` should be created.

### 5. **Verify the Executable**

Confirm that the executable has debugging symbols.

```bash
file matrix
```

**Expected Output:**

```
matrix: ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked, with debug_info, not stripped
```

## Running VTune Analysis from the Command Line

Intel VTune Profiler can be operated using the `vtune` command-line interface (CLI). To bypass `ptrace` restrictions, we'll configure VTune to use hardware-based sampling.

### 1. **Basic Hotspots Analysis with Hardware Sampling**

Identify which functions consume the most CPU time using hardware-based sampling.

```bash
vtune -collect hotspots -knob sampling-mode=hw -result-dir hotspots_hw_analysis ./matrix
```

- **Parameters:**
  - `-collect hotspots`: Initiates a hotspots analysis, focusing on CPU time consumption.
  - `-knob sampling-mode=hw`: Sets the sampling mode to hardware-based, reducing reliance on `ptrace`.
  - `-result-dir hotspots_hw_analysis`: Specifies the directory to store the analysis results.
  - `./matrix`: Runs your application.

**Sample Output:**

```
vtune: Warning: Access to /proc/kallsyms file is limited. Consider changing /proc/sys/kernel/kptr_restrict to 0 to enable resolution of OS kernel and kernel module symbols.

vtune: Warning: To profile kernel modules during the session, make sure they are available in the /lib/modules/kernel_version/ location.

vtune: Collection started. To stop the collection, either press CTRL-C or enter from another console window: vtune -r /mnt/ufs18/home-005/scouch/cmse822/oneAPI-samples/Tools/VTuneProfiler/matrix_multiply_c/hotspots_hw_analysis -command stop.

Addr of buf1 = 0x1458faeb7010
Offs of buf1 = 0x1458faeb7180
Addr of buf2 = 0x1458f8eb6010
Offs of buf2 = 0x1458f8eb61c0
Addr of buf3 = 0x1458f6eb5010
Offs of buf3 = 0x1458f6eb5100
Addr of buf4 = 0x1458f4eb4010
Offs of buf4 = 0x1458f4eb4140
Threads #: 16 Pthreads
Matrix size: 2048
Using multiply kernel: multiply1
Freq = 2.400000 GHz
Execution time = 3.457 seconds
MFLOPS: 4970.220 mflops
vtune: Collection stopped.
vtune: Using result path `/mnt/ufs18/home-005/scouch/cmse822/oneAPI-samples/Tools/VTuneProfiler/matrix_multiply_c/hotspots_hw_analysis'
vtune: Executing actions 75 % Generating a report                              Elapsed Time: 3.543s
    CPU Time: 41.528s
        Effective Time: 41.528s
        Spin Time: 0s
        Overhead Time: 0s
    Instructions Retired: 54,864,000,000
    Microarchitecture Usage: 21.6% of Pipeline Slots
     | You code efficiency on this platform is too low.
     |
     | Possible cause: memory stalls, instruction starvation, branch
     | misprediction or long latency instructions.
     |
     | Next steps: Run Microarchitecture Exploration analysis to identify the
     | cause of the low microarchitecture usage efficiency.
     |
        CPI Rate: 2.510
         | The CPI may be too high. This could be caused by issues such as
         | memory stalls, instruction starvation, branch misprediction or long
         | latency instructions. Explore the other hardware-related metrics to
         | identify what is causing high CPI.
         |
    Total Thread Count: 17
    Paused Time: 0s

Top Hotspots
Function                    Module     CPU Time  % of CPU Time(%)
--------------------------  ---------  --------  ----------------
multiply1                   matrix      40.576s             97.7%
[Outside any known module]  [Unknown]    0.942s              2.3%
init_arr                    matrix       0.010s              0.0%
Effective CPU Utilization: 29.3%
 | The metric value is low, which may signal a poor logical CPU cores
 | utilization caused by load imbalance, threading runtime overhead, contended
 | synchronization, or thread/process underutilization. Explore sub-metrics to
 | estimate the efficiency of MPI and OpenMP parallelism or run the Locks and
 | Waits analysis to identify parallel bottlenecks for other parallel runtimes.
 |
    Average Effective CPU Utilization: 11.722 out of 40
Collection and Platform Info
    Application Command Line: ./matrix
    Operating System: 5.15.0-126-generic DISTRIB_ID=Ubuntu DISTRIB_RELEASE=22.04 DISTRIB_CODENAME=jammy DISTRIB_DESCRIPTION="Ubuntu 22.04.5 LTS"
    Computer Name: dev-intel18
    Result Size: 9.7 MB
    Collection start time: 14:32:11 21/01/2025 UTC
    Collection stop time: 14:32:15 21/01/2025 UTC
    Collector Type: Driverless Perf per-process sampling
    CPU
        Name: Intel(R) Xeon(R) Processor code named Skylake
        Frequency: 2.394 GHz
        Logical CPU Count: 40
        LLC size: 28.8 MB
        Cache Allocation Technology
            Level 2 capability: not detected
            Level 3 capability: available

If you want to skip descriptions of detected performance issues in the report,
enter: vtune -report summary -report-knob show-issues=false -r <my_result_dir>.
Alternatively, you may view the report in the csv format: vtune -report
<report_name> -format=csv.
vtune: Executing actions 100 % done
```

### 2. **Concurrency Analysis**

Analyze threading and parallelism within the application.

```bash
vtune -collect concurrency -result-dir concurrency_analysis ./matrix
```

### 3. **Memory Access Analysis**

Investigate memory-related bottlenecks.

```bash
vtune -collect memory-access -result-dir memory_analysis ./matrix
```

### 4. **Specifying Sampling Options**

Customize sampling intervals to suit your profiling needs.

```bash
vtune -collect hotspots -knob sampling-interval=10 -result-dir custom_hotspots ./matrix
```

### 5. **Running in Batch Mode**

Suitable for automated scripts and scheduling.

```bash
vtune -collect hotspots -result-dir batch_analysis -report-output batch_report.txt ./matrix
```

## Understanding and Resolving Common VTune Errors

### **Error: `Cannot start data collection because the scope of ptrace system call is limited.`**

**Error Message:**

```
vtune: Error: Cannot start data collection because the scope of ptrace system call is limited. To enable profiling, please set /proc/sys/kernel/yama/ptrace_scope to 0. To make this change permanent, set kernel.yama.ptrace_scope to 0 in /etc/sysctl.d/10-ptrace.conf and reboot the machine.
```

**Meaning:**

VTune relies on the `ptrace` system call to attach to and monitor your application during profiling. The current system setting restricts the use of `ptrace`, preventing VTune from initiating data collection.

**Resolution:**

To prevent this error and continue profiling without modifying `ptrace_scope`, configure VTune to use hardware-based sampling, which minimizes reliance on `ptrace`. This was implemented in the command used above:

```bash
vtune -collect hotspots -knob sampling-mode=hw -result-dir hotspots_hw_analysis ./matrix
```

**Additional Recommendations:**

- **Run on `dev-intel18` Node:** Ensure you are profiling on the `dev-intel18` node where necessary permissions and configurations are in place.
- **Contact System Administrators:** If further issues arise, consult with system administrators for assistance with system settings.

## Collecting and Viewing Results

After profiling, you can generate and interpret reports to understand the performance characteristics of your application.

### 1. **Listing Available Results**

```bash
vtune -list-results
```

This command displays all available analysis results in the current directory.

### 2. **Generating a Hotspots Report**

```bash
vtune -report hotspots -result-dir hotspots_hw_analysis -format text > hotspots_report.txt
```

- **Parameters:**
  - `-report hotspots`: Specifies the type of report.
  - `-result-dir hotspots_hw_analysis`: Directory containing the analysis results.
  - `-format text`: Outputs the report in plain text format.
  - `> hotspots_report.txt`: Redirects the output to a file.

### 3. **Viewing the Report**

Use `less` or any text editor to view the generated report.

```bash
less hotspots_report.txt
```

### 4. **Transferring Results to Your Local Machine**

For detailed analysis, you might prefer viewing reports on your local machine.

```bash
scp username@icer.hpc.msu.edu:/home/username/vtune_profiling/oneAPI-samples/Tools/VTuneProfiler/matrix_multiply_c/hotspots_report.txt ~/local_directory/
```

Replace `username`, `icer.hpc.msu.edu`, and `~/local_directory/` with your actual username, cluster address, and desired local directory.

## Analyzing Performance Metrics

Understanding the output from VTune Profiler is crucial for identifying performance bottlenecks and optimizing your application. Here's a detailed explanation of a sample VTune output:

### **Sample VTune Output:**

```
CPU Time: 41.528s
    Effective Time: 41.528s
    Spin Time: 0s
    Overhead Time: 0s
Instructions Retired: 54,864,000,000
Microarchitecture Usage: 21.6% of Pipeline Slots
 | Your code efficiency on this platform is too low.
 |
 | Possible cause: memory stalls, instruction starvation, branch
 | misprediction or long latency instructions.
 |
 | Next steps: Run Microarchitecture Exploration analysis to identify the
 | cause of the low microarchitecture usage efficiency.
 |
    CPI Rate: 2.510
     | The CPI may be too high. This could be caused by issues such as
     | memory stalls, instruction starvation, branch misprediction or long
     | latency instructions. Explore the other hardware-related metrics to
     | identify what is causing high CPI.
     |
Total Thread Count: 17
Paused Time: 0s

Top Hotspots
Function                    Module     CPU Time  % of CPU Time(%)
--------------------------  ---------  --------  ----------------
multiply1                   matrix      40.576s             97.7%
[Outside any known module]  [Unknown]    0.942s              2.3%
init_arr                    matrix       0.010s              0.0%
Effective CPU Utilization: 29.3%
 | The metric value is low, which may signal a poor logical CPU cores
 | utilization caused by load imbalance, threading runtime overhead, contended
 | synchronization, or thread/process underutilization. Explore sub-metrics to
 | estimate the efficiency of MPI and OpenMP parallelism or run the Locks and
 | Waits analysis to identify parallel bottlenecks for other parallel runtimes.
 |
    Average Effective CPU Utilization: 11.722 out of 40
```

### **Detailed Breakdown:**

#### 1. **CPU Time Metrics**

- **CPU Time:** `41.528s`
  - **Effective Time:** `41.528s`
    - **Interpretation:** The CPU was actively executing instructions for the entire duration of the profiling session.
  - **Spin Time:** `0s`
    - **Interpretation:** No time was spent in busy-waiting loops, indicating efficient thread synchronization.
  - **Overhead Time:** `0s`
    - **Interpretation:** Minimal profiling overhead, ensuring that profiling data accurately reflects application performance.

- **Instructions Retired:** `54,864,000,000`
  - **Interpretation:** Total number of instructions executed by the CPU during profiling. A high number signifies intensive computational activity.

#### 2. **Microarchitecture Usage**

- **Microarchitecture Usage:** `21.6% of Pipeline Slots`
  - **Interpretation:** Indicates the efficiency of CPU pipeline utilization. Low usage (21.6%) suggests significant idle pipeline slots, meaning the CPU isn't being fully utilized.
  - **Possible Causes:**
    - Memory stalls
    - Instruction starvation
    - Branch misprediction
    - Long latency instructions

- **CPI Rate (Cycles Per Instruction):** `2.510`
  - **Interpretation:** Average number of CPU cycles taken per instruction. A higher CPI indicates inefficiency.
  - **Possible Causes:**
    - Memory stalls
    - Instruction starvation
    - Branch misprediction
    - Long latency instructions

#### 3. **Threading Metrics**

- **Total Thread Count:** `17`
  - **Interpretation:** The application spawned 17 threads, slightly exceeding the specified 16 pthreads, possibly including VTune's internal threads.

- **Effective CPU Utilization:** `29.3%`
  - **Interpretation:** Only 29.3% of the available logical CPU cores were effectively utilized, indicating underutilization of CPU resources.
  - **Average Effective CPU Utilization:** `11.722 out of 40`
    - **Interpretation:** On average, only approximately 11.7 out of 40 logical CPUs were actively utilized.

- **Possible Causes of Low CPU Utilization:**
  - Load imbalance
  - Threading runtime overhead
  - Contended synchronization
  - Thread/process underutilization

#### 4. **Top Hotspots**

```
Top Hotspots
Function                    Module     CPU Time  % of CPU Time(%)
--------------------------  ---------  --------  ----------------
multiply1                   matrix      40.576s             97.7%
[Outside any known module]  [Unknown]    0.942s              2.3%
init_arr                    matrix       0.010s              0.0%
```

- **multiply1 Function:**
  - **CPU Time:** `40.576s`
  - **% of CPU Time:** `97.7%`
  - **Interpretation:** Dominates CPU usage, accounting for nearly all execution time. Primary focus for performance optimization.

- **[Outside any known module]:**
  - **CPU Time:** `0.942s`
  - **% of CPU Time:** `2.3%`
  - **Interpretation:** Represents CPU time spent outside recognized modules, such as system calls or library functions.

- **init_arr Function:**
  - **CPU Time:** `0.010s`
  - **% of CPU Time:** `0.0%`
  - **Interpretation:** Minimal CPU usage, likely involved in initializing data structures.

## Recommendations and Next Steps

Based on the VTune Profiler output, here are actionable steps to optimize your `matrix_multiply` application:

### 1. **Optimize the `multiply1` Function**

Given that `multiply1` accounts for **97.7%** of CPU time, it's the primary target for optimization.

#### **Potential Optimization Strategies:**

- **Algorithmic Improvements:**
  - **Blocking (Tiling):** Break down large matrix operations into smaller blocks that fit into the cache, reducing cache misses.
  - **Strassen's Algorithm:** Implement advanced algorithms that reduce the computational complexity of matrix multiplication.

- **Parallelization Enhancements:**
  - **Thread Affinity:** Bind threads to specific CPU cores to reduce context switching and improve cache utilization.
  - **Load Balancing:** Ensure that work is evenly distributed among all threads to prevent some threads from being idle.

- **Vectorization:**
  - **SIMD Instructions:** Utilize Single Instruction, Multiple Data (SIMD) instructions to perform parallel operations on multiple data points.
  - **Compiler Flags:** Use optimization flags that enable auto-vectorization (e.g., `-O3`, `-xHost`).

- **Memory Access Optimization:**
  - **Data Layout:** Arrange data structures to maximize cache line utilization and minimize cache misses.
  - **Prefetching:** Implement prefetching to load data into the cache before it's needed, reducing memory latency.

### 2. **Address Low Microarchitecture Usage and High CPI**

- **Run Microarchitecture Exploration Analysis:**
  - **Purpose:** Gain deeper insights into pipeline stalls, cache misses, branch mispredictions, and other microarchitectural inefficiencies.
  - **Command:**
    ```bash
    vtune -collect microarchitecture-exploration -result-dir micro_exploration ./matrix
    ```

- **Optimize Based on Findings:**
  - **Memory Stalls:** Improve data locality and cache usage.
  - **Instruction Starvation:** Increase instruction-level parallelism.
  - **Branch Mispredictions:** Simplify branching logic or use branch prediction-friendly patterns.
  - **Long Latency Instructions:** Optimize or replace instructions that take multiple cycles to execute.

### 3. **Enhance Effective CPU Utilization**

- **Increase Parallelism:**
  - **More Threads:** Experiment with increasing the number of threads beyond 16, up to the number of logical CPUs (40), to better utilize available resources.

- **Reduce Synchronization Overhead:**
  - **Lock-Free Programming:** Implement lock-free data structures or algorithms to minimize waiting times.
  - **Batch Operations:** Reduce the frequency of synchronization points by batching operations.

- **Load Balancing:**
  - **Dynamic Scheduling:** Use dynamic scheduling techniques to allocate work based on thread performance and workload.
  - **Thread Pools:** Implement thread pools to manage and reuse threads efficiently.

### 4. **Investigate Unknown Modules Activity**

- **Identify `[Outside any known module]`:**
  - **Potential Areas:** System calls, dynamically loaded libraries, or background processes.
  - **Action:** Use additional profiling tools or VTune's advanced features to map these activities to specific functions or modules.

### 5. **Address Profiling Warnings**

- **Kernel Symbol Access:**
  - **Action:** If kernel-level profiling is essential, adjust system settings as previously discussed or consult with system administrators to enable necessary access.

### 6. **Iterative Profiling and Optimization**

- **Cycle of Improvement:**
  - **Profile:** Collect detailed performance data.
  - **Analyze:** Identify bottlenecks and inefficiencies.
  - **Optimize:** Implement code optimizations based on analysis.
  - **Re-profile:** Measure the impact of optimizations and identify new areas for improvement.

## Advanced VTune Features

Intel VTune offers advanced capabilities for in-depth performance analysis. Here are some features you can leverage with the `matrix_multiply` sample.

### 1. **Custom Analysis Types**

Create tailored analyses based on specific metrics relevant to your application.

```bash
vtune -collect general-exploration -knob sampling-mode=hw -result-dir custom_analysis ./matrix
```

- **Parameters:**
  - `-collect general-exploration`: A flexible analysis type for custom metrics.
  - `-knob sampling-mode=hw`: Ensures hardware-based sampling.
  - `-result-dir custom_analysis`: Specifies the directory to store results.

### 2. **Automated Reports**

Schedule VTune analyses and generate reports automatically, integrating them into scripts or CI/CD pipelines.

```bash
vtune -collect hotspots -knob sampling-mode=hw -result-dir auto_analysis -report-output auto_report.txt ./matrix
```

- This command runs the analysis and immediately generates a report, saving it to `auto_report.txt`.

### 3. **Integration with Build Systems**

Integrate VTune profiling into automated build and test pipelines using scripts. For example, create a shell script `profile.sh`:

```bash
#!/bin/bash

# Load modules
module load intel/oneapi/2025.0
module load vtune/2025.0

# Compile the application
make clean
make

# Run VTune analysis
vtune -collect hotspots -knob sampling-mode=hw -result-dir pipeline_analysis -report-output pipeline_report.txt ./matrix

# Notify completion
echo "Profiling complete. Report saved to pipeline_report.txt"
```

- **Make the script executable:**

  ```bash
  chmod +x profile.sh
  ```

- **Execute the script:**

  ```bash
  ./profile.sh
  ```

### 4. **Remote Analysis**

Profile applications running on different nodes if your HPCC setup supports remote profiling.

```bash
vtune -collect hotspots -knob sampling-mode=hw -result-dir remote_analysis -target remote_node ./matrix
```

- Replace `remote_node` with the actual target node's identifier.

## Best Practices and Tips

1. **Minimize Noise:**
   - Ensure that no other heavy processes are running on the node during profiling to obtain accurate results.

2. **Use Representative Workloads:**
   - Profile using real-world matrix sizes and scenarios to capture genuine performance characteristics.

3. **Iterative Profiling:**
   - Profile, optimize, and profile again to measure the impact of your changes.

4. **Leverage Multiple Analysis Types:**
   - Combine different analysis types (e.g., hotspots, concurrency, memory) for a comprehensive performance overview.

5. **Document Findings:**
   - Keep detailed notes of identified bottlenecks and the optimizations applied for future reference.

6. **Understand the Application:**
   - A deep understanding of the `matrix_multiply` application's structure and algorithms will aid in effective profiling and optimization.

7. **Utilize VTune Documentation:**
   - Refer to the [Intel VTune Profiler Documentation](https://www.intel.com/content/www/us/en/docs/vtune-profiler/) for advanced usage and troubleshooting.

8. **Run on `dev-intel18` Node:**
   - Ensure all profiling activities are conducted on the `dev-intel18` node to align with system configurations and permissions.

## Conclusion

Profiling and optimizing applications are critical steps in harnessing the full potential of high-performance computing systems like MSU ICER HPCC. Intel VTune Profiler, with its robust command-line interface, offers powerful tools to identify and resolve performance bottlenecks effectively. By following this updated tutorial and applying the techniques to the `matrix_multiply` sample, you can gain valuable insights into your application's performance and implement meaningful optimizations.

**Key Takeaways:**

- **Addressed Common Errors:** Configured VTune to use hardware-based sampling to avoid `ptrace` restrictions.
- **Detailed Output Analysis:** Provided comprehensive explanations of VTune's profiling output to facilitate informed optimizations.
- **Focused Optimization:** Emphasized optimizing the `multiply1` function as the primary hotspot for performance gains.
- **Advanced Techniques:** Introduced advanced VTune features and best practices to enhance profiling workflows.

**Remember:** Performance tuning is an iterative process. Continuously profile, analyze, and refine your code to achieve optimal efficiency.

Happy profiling!