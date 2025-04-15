# Project 6: Hybrid Parallelism in Agoge

## Overview

In this project, you will extend the Agoge codebase to implement hybrid parallelism using OpenMP for both CPU and GPU (offload) execution. You will focus on the EulerSolver, with a bonus challenge to parallelize the GravitySolver. You will conduct thorough correctness testing and performance analysis, including strong and weak scaling studies, efficiency calculations, and CPU vs. GPU comparisons. All work should be performed and benchmarked on the MSU HPCC, including multi-node and multi-GPU runs.

## 1. OpenMP Parallelism in EulerSolver

### 1.1. CPU Parallelism
- **Task:** Add OpenMP parallel regions to the core computational loops in `agoge/src/EulerSolver.cpp`.
- **Reference:** Focus on the main update and flux calculation routines. Use `#pragma omp parallel for` to parallelize outer loops over grid cells.
- **Hints:**
  - Ensure all shared and private variables are correctly specified.
  - Use OpenMP reduction for summing operations (e.g., conserved quantities).
  - Test with different thread counts (set via `OMP_NUM_THREADS`).

### 1.2. GPU Offload (OpenMP Target)
- **Task:** Implement OpenMP target offload for the same routines, enabling execution on GPUs.
- **Reference:** Use `#pragma omp target teams distribute parallel for` in `EulerSolver.cpp`.
- **Hints:**
  - Ensure data is mapped to/from the device as needed.
  - Test on HPCC GPU nodes (e.g., `amd20` or `intel18` with GPUs).
  - Use OpenMP device selection and offload environment variables.

### 1.3. Code Organization
- **Task:** Use preprocessor macros or CMake options to enable/disable OpenMP and GPU offload code paths.
- **Reference:** Update `agoge/src/EulerSolver.cpp` and `agoge/include/agoge/EulerSolver.hpp` as needed.
- **Hints:**
  - Guard OpenMP code with `#ifdef _OPENMP`.
  - Provide a serial fallback for environments without OpenMP.


## 2. Correctness Testing

- **Task:** Ensure all parallel versions produce correct results.
- **Reference:** Use existing test problems in `agoge/problems/` (e.g., `SodShockTube`).
- **Hints:**
  - Compare output fields (density, momentum, energy) between serial, OpenMP CPU, and GPU runs.
  - Use the `agoge/viz/agoge_compare.py` tool to check serial and parallel output. 

## 3. Performance Analysis

### 3.1. Strong and Weak Scaling Studies
- **Task:** Measure and analyze performance for both strong and weak scaling.
- **Reference:**
  - Use the performance report printed by Agoge (timing of major routines).
  - Run on different node types and with varying numbers of threads/GPUs.
- **Hints:**
  - For strong scaling: fix the problem size, vary the number of threads/GPUs.
  - For weak scaling: increase the problem size proportionally with resources.
  - Record wall-clock times and compute zone-updates per second.

### 3.2. Parallel Efficiency
- **Task:** Compute and report parallel efficiency for all configurations.
- **Formula:**
  - Efficiency = (Speedup) / (Number of threads/GPUs)
  - Speedup = (Serial time) / (Parallel time)
- **Hints:**
  - Present results in tables and plots.
  - Discuss any bottlenecks or scaling limitations observed.

### 3.3. CPU vs. GPU Comparison
- **Task:** Compare performance between CPU and GPU runs.
- **Reference:**
  - Use the same problem and grid size for both runs.
  - Analyze and discuss the results.

### 3.4. Multi-node and Multi-GPU Runs
- **Task:** Run Agoge on multiple nodes and/or with multiple GPUs per node.
- **Reference:**
  - Use MPI for distributed memory parallelism (see `agoge/src/main.cpp` and `agoge/include/agoge/Field3d.hpp`).
  - Combine MPI with OpenMP for hybrid runs.
- **Hints:**
  - Document your job submission scripts and environment settings.
  - Report any issues or special considerations for multi-node/multi-GPU execution.



## 4. Bonus: OpenMP in GravitySolver

- **Task:** Add OpenMP parallelism to `agoge/src/GravitySolver.cpp`.
- **Reference:**
  - Focus on FFT and Poisson solve routines.
  - Use OpenMP to parallelize over independent grid slices or dimensions.
- **Hints:**
  - Ensure thread safety of any third-party libraries (e.g., FFTW).
  - Test for correctness and performance as above.



## 5. Documentation and Submission

- **Task:** Document all code changes with Doxygen comments and inline explanations.
- **Reference:**
  - Follow the coding standards in `docs/codingStandards.md`.
  - Update or create a `REPORT.md` summarizing your methodology, results, and conclusions.
- **Hints:**
  - Include plots/tables of scaling and efficiency results.
  - Clearly state which Agoge code files you modified.
  - Submit your code and report via GitHub as instructed.



## 6. Additional Notes

- Use the MSU HPCC for all runs. Document the node types and modules used.
- If you encounter issues with OpenMP or GPU offload, consult the Agoge Slack/Discussion board or course staff.
- For visualization, use scripts in `agoge/viz/` or your own tools.
- For extra credit, optimize and document OpenMP parallelism in the GravitySolver.



**Good luck! May your parallel code be correct and your speedups be superlinear!**