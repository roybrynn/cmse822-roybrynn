# Agoge

**Agoge** is a pedagogical 3D **compressible Euler** solver written in C++. It was developed as a teaching and research tool for advanced computational physics and engineering courses (e.g., CMSE 822) to demonstrate:

- High-performance computing (HPC) concepts
- Numerical methods for hyperbolic PDEs
- Modern software engineering practices using C++ and CMake
- Techniques for modular code design and extensibility
- Integration of self-gravity via an FFT-based Poisson solver
- Potential for GPU offloading and other advanced parallelization strategies

The name “agoge”—inspired by the rigorous training system of ancient Sparta—reflects its role as a challenging training ground for mastering computational methods and software development in scientific computing.

---

## Key Features

1. **3D Compressible Euler Equations**
   - Second-order finite-difference (with minmod limiters) or higher-order schemes available.
   - Two-stage Runge-Kutta time stepping.
   - Supports periodic and simple boundary conditions.

2. **Self-Gravity (Optional)**
   - FFT-based Poisson solver to solve \(\nabla^2 \phi = 4\pi G\,\rho\)
   - Adds gravitational source terms (body force in momentum and work in energy equations) using spatially varying potential (\(\phi\)).

3. **Modular Code Structure**
   - Clear separation of data structures (`Field3D`), solver routines (`EulerSolver`), gravity routines (`GravitySolver`), and I/O (`HDF5_IO`).
   - Build system based on CMake for cross-platform compilation and easy extensibility.

4. **HDF5 I/O**
   - Write final solution fields (density, momentum components, energy, gravitational potential, etc.) to HDF5 files.
   - Output is compatible with popular visualization tools (e.g., yt, ParaView).

5. **Extensibility**
   - Designed to be a test bed for adding additional physics (e.g., viscosity, magnetic fields, multi-species flows) and numerical methods (e.g., higher-order or WENO schemes).   
   - Performance monitoring via a custom `PerformanceMonitor` module.

---

## Repository Layout

```
agoge/
├── include/
│   └── agoge/
│       ├── BoundaryManager.hpp
│       ├── Config.hpp
│       ├── Domain.hpp
│       ├── EulerSolver.hpp
│       ├── Field3d.hpp
│       ├── GravitySolver.hpp
│       ├── HDF5_IO.hpp
│       ├── ParameterSystem.hpp
│       ├── PerformanceMonitor.hpp
│       └── YamlParser.hpp
├── problems/
│   ├── GravityCollapse.cpp
│   ├── GravityCollapse.hpp
│   ├── GaussianPulse.cpp
│   ├── GaussianPulse.hpp
│   ├── Problem.hpp
│   ├── ProblemRegistry.cpp
│   └── ProblemRegistry.hpp
├── src/
│   ├── BoundaryManager.cpp
│   ├── Config.cpp
│   ├── EulerSolver.cpp
│   ├── Field3d.cpp
│   ├── GravitySolver.cpp
│   ├── HDF5_IO.cpp
│   ├── main.cpp
│   └── YamlParser.cpp
├── tests/
│   ├── CMakeLists.txt
│   ├── test_EulerSolver.cpp
│   └── test_Field3D.cpp
├── viz/
│   ├── agoge_yt.py
│   ├── agoge_viz.py
│   └── example.ipynb
├── projects/
│   └── project1/
│       └── README.md
├── .github/
│   └── workflows/
│       └── update_students.yml   # (For GitHub Classroom automation)
├── Makefile
├── README.md
├── design.md
└── instructions.md
```

- **`include/`**: Public header files for core components.
- **`src/`**: Implementation of the solver, including `main.cpp`.
- **`problems/`**: Problem definitions (e.g., SodShockTube, GravityCollapse) and a registry for selecting the active problem.
- **`tests/`**: Unit and integration tests.
- **`.github/workflows/`**: GitHub Actions workflow files for automated tasks (e.g., updating student repos in Classroom).
- **`viz/`**: Visualization scripts and example notebooks.

---

## Building Agoge on HPCC at MSU

Agoge is designed to run on high-performance computing clusters. At MSU’s HPCC, follow these steps:

1. **Access HPCC:**
   - Log in to HPCC using your preferred method (SSH, terminal, etc.).
   - Load the necessary modules (e.g., GCC, HDF5). For example:
     ```bash
     module purge
     module load HDF5
     ```

2. **Clone the Repository:**
   ```bash
   git clone https://github.com/your-username/agoge.git
   cd agoge
   ```

3. **Configure and Build:**
   From the top-level `agoge` directory, you may simply run `make` to build the code:
   
4. **Generate the Executable:**
   - The build process should produce an executable (e.g., `agoge_run`).

---

## Running Agoge on HPCC

Once built, you can run Agoge on HPCC either interactively or as a batch job.

### **Interactive Run:**

```bash
./agoge_run problems/SodShockTube.yaml
```

### **Batch Job Submission Example (SLURM):**

Create a job script (e.g., `job.sh`):

```bash
#!/bin/bash
#SBATCH --job-name=agoge_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:30:00
#SBATCH --partition=standard

# Load necessary modules
module purge
module load HDF5

# Navigate to the agoge directory
cd $SLURM_SUBMIT_DIR

# Run the executable with the desired problem configuration
./agoge_run ../problems/SodShockTube.yaml
```

Submit the job:

```bash
sbatch job.sh
```

---

## Running Self-Gravity & Advanced Features

- **Self-Gravity:**
  - Enable self-gravity simply by specifying a problem that uses the GravitySolver.
  - The solver uses an FFT-based Poisson solver (see `GravitySolver.cpp`) to compute the gravitational potential.  

- **GPU Offloading (Future Work):**
  - The code is structured to allow extensions such as GPU offloading via CUDA or OpenMP target directives.
  - This will require additional configuration and might involve separate compilation flags and modules on HPCC.

---

## Output & Visualization

- **Output Files:**
  - Agoge writes simulation data to HDF5 files (e.g., `agoge_final.h5`).
- **Visualization Tools:**
  - Use Python (with yt or matplotlib), ParaView, or similar HDF5-compatible tools.
  - See the example Jupyter notebook in the `viz` directory.

---

## Contributing

Contributions, bug reports, and feature requests are welcome. Please submit pull requests or open issues on GitHub. This project is used in CMSE 822 as a training ground for HPC and scientific software development.

---

## License

Agoge is released under the [MIT License](LICENSE). You are free to use, modify, and distribute Agoge for personal, educational, and commercial purposes. Please cite or link back to the repository if you use it in your research or teaching.

---

## Acknowledgments

- **CMSE 822 at MSU:** For inspiring the development of Agoge as an educational tool in parallel computing.
- **The HPC Community:** For contributions and best practices that have informed Agoge’s design.
- **Libraries and Tools:** Special thanks to the developers of HDF5, FFTW, and the GitHub Actions community.

---

## Contact

For questions or clarifications, please contact the course instructor or open an issue on the repository.
