Here's the **updated tutorial** for setting up, building, and running the **Empirical Roofline Toolkit (ERT) on HPCC at Michigan State University**. It incorporates all the fixes, optimizations, and best practices we’ve developed throughout our discussion.

---

# **Empirical Roofline Toolkit (ERT) Setup and Execution on HPCC at MSU**

## **Introduction**
The **Empirical Roofline Toolkit (ERT)** is used to measure and analyze computational performance using the **Roofline model**. This tutorial walks you through:
1. **Installing and configuring ERT** on **HPCC at MSU**.
2. **Building and running ERT** across different **node types** (`amd20`, `intel16`, `intel18`, etc.).
3. **Optimizing parameters** for different architectures.
4. **Converting `.ps` output to `.pdf`** for easier viewing.

---

## **1. Load Required Modules**
Since ERT requires **a modern compiler**, **HDF5**, **Python 3**, and **GNUplot**, first load the necessary modules.

Run the following in your HPCC terminal:
```bash
module purge
module load GCC/12.3.0
module load OpenMPI/4.1.5-GCC-12.3.0
module load HDF5/1.14.0-gompi-2023a
module load git/2.41.0-GCCcore-12.3.0-nodocs
module load Vim/9.1.0004-GCCcore-12.3.0
module load CMake/3.27.6-GCCcore-12.3.0
module load Python/3.10.12-GCCcore-12.3.0-bare
module load gnuplot
```
Ensure that you are using **GCC 12.3.0**:
```bash
gcc --version
```
If it still shows an **older version**, manually adjust your paths:
```bash
export PATH=/opt/software/GCCcore-12.3.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/software/GCCcore-12.3.0/lib64:$LD_LIBRARY_PATH
```

---

## **2. Clone and Set Up ERT**
Navigate to your working directory and clone the ERT repository:
```bash
cd ~/cmse822/
git clone https://bitbucket.org/berkeleylab/cs-roofline-toolkit.git
cd cs-roofline-toolkit/Empirical_Roofline_Tool-1.1.0
```

---

## **3. Configure ERT for HPCC**
Create a configuration file (`config.hpcc.msu.edu.01`) optimized for **serial execution** on HPCC.

```bash
# === ERT Configuration File for MSU HPCC ===

# Set where results should be stored (relative or absolute path)
ERT_RESULTS ERT_Results

# Theoretical peak values for reference
ERT_SPEC_GBYTES_DRAM 200
ERT_SPEC_GFLOPS 5000

# === Compilation Settings ===
ERT_CC g++
ERT_CFLAGS -O3 -march=native -std=c++11
ERT_LD g++
ERT_LDFLAGS -lstdc++
ERT_LDLIBS -lstdc++

# === ERT Microkernel Settings ===
ERT_DRIVER driver1
ERT_KERNEL kernel1
ERT_FLOPS 8
ERT_ALIGN 64

# === Parallelism Settings (Set to Serial Mode) ===
ERT_MPI False
ERT_OPENMP False
ERT_GPU False

# === Execution Settings ===
ERT_RUN ./ERT_CODE

# === Experiment Settings ===
ERT_NUM_EXPERIMENTS 3
ERT_MEMORY_MAX 1073741824
ERT_WORKING_SET_MIN 8192
ERT_TRIALS_MIN 10

# === Output and Post-processing ===
ERT_GNUPLOT gnuplot
```

---

## **4. Build ERT**
Now, build ERT with:
```bash
./ert --build ~/cmse822/cmse822-codex-private/projects/project1/config.hpcc.msu.edu.01
```
If you get an error about **undefined references to `std::ios_base::Init`**, ensure that `ERT_LDLIBS -lstdc++` is set in the config file.

---

## **5. Submit a SLURM Job to Run ERT**
Since we are benchmarking across different node types, submit jobs via **SLURM**.

### **SLURM Script (`run_ert.slurm`)**
```bash
#!/bin/bash
#SBATCH --job-name=ERT_Benchmark
#SBATCH --output=ERT_output_%j.txt
#SBATCH --error=ERT_error_%j.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=amd20  # Change this to test different architectures

module purge
module load GCC/12.3.0
module load gnuplot
module load Python/3.10.12-GCCcore-12.3.0-bare

cd $SLURM_SUBMIT_DIR
./ert --run ~/cmse822/cmse822-codex-private/projects/project1/config.hpcc.msu.edu.01
```

### **Submit the job for different architectures**
```bash
sbatch --constraint=amd20 run_ert.slurm
sbatch --constraint=intel16 run_ert.slurm
sbatch --constraint=intel18 run_ert.slurm
```

---

## **6. Analyze Results**
Once the job completes, check the results directory:
```bash
ls -lh ERT_Results/Run.*/
```

### **View the Roofline Plot**
ERT generates a **PostScript (`.ps`) file** for the Roofline plot:
```bash
display ERT_Results/Run.<N>/roofline.ps
```

---

## **7. Convert `.ps` to `.pdf`**
If you want a **PDF instead of PostScript**, convert it using:
```bash
ps2pdf ERT_Results/Run.<N>/roofline.ps ERT_Results/Run.<N>/roofline.pdf
```
To convert **all ERT-generated `.ps` files**:
```bash
find ERT_Results/ -name "*.ps" -exec ps2pdf {} {}.pdf \;
```

---

## **8. Tuning for Different Architectures**
| **Architecture** | **Optimized Configurations** |
|-----------------|-----------------------------|
| **AMD Nodes (`amd20`)** | `ERT_CFLAGS -O3 -march=znver2`, `ERT_ALIGN 64`, `ERT_FLOPS 16` |
| **Intel Skylake (`intel16`)** | `ERT_CFLAGS -O3 -march=skylake-avx512`, `ERT_ALIGN 128`, `ERT_FLOPS 16` |
| **Intel Cascade Lake (`intel18`)** | `ERT_CFLAGS -O3 -march=cascadelake`, `ERT_ALIGN 64`, `ERT_FLOPS 8` |
| **GPU Nodes (CUDA)** | `ERT_GPU True`, `ERT_CFLAGS -arch=sm_80` |

---

## **9. Troubleshooting**
### **Compiler Issues**
If `gcc --version` shows an older version after loading `GCC/12.3.0`, manually set it:
```bash
export CC=g++
export CXX=g++
```

### **Linker Errors (`std::ios_base::Init` missing)**
Ensure that `ERT_LDLIBS -lstdc++` is set in the config file.

### **Missing CUDA**
If `module load CUDA/12.2.2` fails, check available versions:
```bash
module spider CUDA
```
Use the closest version instead.

