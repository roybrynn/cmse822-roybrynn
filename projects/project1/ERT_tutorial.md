# **Empirical Roofline Toolkit (ERT) Setup and Execution on HPCC at MSU**

## **Introduction**

The **Empirical Roofline Toolkit (ERT)** is used to measure and analyze computational performance using the **Roofline model**. This tutorial walks you through:
1. **Installing and configuring ERT** on **HPCC at MSU**.
2. **Building and running ERT** across different **node types** (`amd20`, `intel16`, `intel18`, etc.).
3. **Optimizing parameters** for different architectures.
4. **Converting `.ps` output to `.pdf`** for easier viewing.

---

## **1. Load Required Modules**

Since ERT requires **a modern compiler**, **HDF5**, **Python 3**, and **GNUplot**, first load the necessary modules. Simply loading the `GNUPlot` module should suffice as all dependent modules (such as GCC) will be loaded as well.

Run the following in your HPCC terminal:
```bash
module purge
module load gnuplot
```
Ensure that you are using **GCC 12.3.0**:
```bash
gcc --version
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

Start from the configuration file (`project1/config.hpcc.msu.edu.01`) optimized for **serial execution** on HPCC.

---

## **4. Build and Run ERT**

Now, build and run ERT all in one go with:
```bash
./ert ~/<path/to/your/repo>/projects/project1/config.hpcc.msu.edu.01
```

If you get an error about **undefined references to `std::ios_base::Init`**, ensure that `ERT_LDLIBS -lstdc++` is set in the config file.

If you get an error about not being able to produce a graph:
```
FLOP count 1...
  Building ERT core code...
  Running ERT core code...
  Processing results...
/bin/sh: 1: gnuplot: not found
  Failure...
Unable to produce a 'graph1' for ERT_Results/Run.012/FLOPS.001

--- Making ERT individual graphs failed ---
```

verify that you have loaded the gnuplot module (see Step 1).

---

## **5. Submit a SLURM Job to Run ERT**

Since we are benchmarking across different node types, I recommend you submit jobs via **SLURM**.
SLURM accepts various _constraints_ to specify the node type to run on. Gather data for various node types and compare and contrast your results. 
A starter SLURM script is now provided in the `project1` directory.

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
./ert ~/cmse822/cmse822-codex-private/projects/project1/config.hpcc.msu.edu.01
```

### **Submit the job for different architectures**

For example:

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

Here are some _suggested_ optimization flags to try and tailor the ERT to different architectures. TBH, I have not noticed a difference when using them...

| **Architecture** | **Optimized Configurations** |
|-----------------|-----------------------------|
| **AMD Nodes (`amd20`)** | `ERT_CFLAGS -O3 -march=znver2`, `ERT_ALIGN 64`, `ERT_FLOPS 16` |
| **Intel Skylake (`intel16`)** | `ERT_CFLAGS -O3 -march=skylake-avx512`, `ERT_ALIGN 128`, `ERT_FLOPS 16` |
| **Intel Cascade Lake (`intel18`)** | `ERT_CFLAGS -O3 -march=cascadelake`, `ERT_ALIGN 64`, `ERT_FLOPS 8` |
| **GPU Nodes (CUDA)** | `ERT_GPU True`, `ERT_CFLAGS -arch=sm_80` |

