# Agoge

**Agoge** is a pedagogical 3D **compressible Euler** solver written in C++. It is designed primarily for CMSE 822 (or similar advanced computational physics/engineering courses) to **demonstrate**:

- High-performance computing (HPC) concepts  
- Numerical methods for hyperbolic PDEs  
- Modular code design and build practices using CMake  
- Extensions to handle **self-gravity** and possible **GPU** offloading

The name “agoge” reflects the code’s purpose as a **training ground** for learning **computational methods** and **software engineering** best practices in scientific computing.

## Key Features

1. **3D Compressible Euler Equations**  
   - Second-order finite-difference (or optionally higher-order)  
   - Two-stage Runge-Kutta time stepping  
   - Periodic or simple boundary conditions  

2. **Self-Gravity (Optional)**  
   - FFT-based Poisson solver for \(\nabla^2 \phi = 4\pi G\,\rho\)  
   - Adds gravitational source terms to momentum & energy equations  

3. **Modular Code Structure**  
   - Separation of data structures (`Field3D`), solver routines (`EulerSolver`), and I/O (`HDF5_IO`)  
   - CMake-based build system for easy compilation and extensibility  

4. **HDF5 I/O**  
   - Example of writing out final states (density, momentum, energy, etc.) to HDF5  
   - Can be viewed with typical HPC visualization tools (yt, ParaView, etc.)  

5. **Room for Extensions**  
   - GPU offloading (e.g., with CUDA, AMReX, or OpenMP target)  
   - Additional physics (viscosity, MHD, multi-species, etc.)  
   - More sophisticated boundary conditions, limiters, Riemann solvers  

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
├── CMakeLists.txt
├── Makefile
├── README.md
├── design.md
└── instructions.md
```

- **`include/`**: Public headers declaring data structures and solver interfaces.  
- **`src/`**: Implementation files, including the main driver (`main.cpp`).  
- **`tests/`**: Unit and integration tests.  
- **`CMakeLists.txt`**: Entry point for the build system.

---

## Building Agoge

You’ll need a **C++17**-compatible compiler, **CMake 3.15+**, and optionally **HDF5** (if you want HDF5 I/O) and **FFTW** (if you want self-gravity).

1. **Clone** this repository:
   ```bash
   git clone https://github.com/username/agoge.git
   cd agoge
   ```
2. **Configure** and **build**:
   ```bash
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make
   ```
3. **Optional**: If using HDF5 or FFTW, ensure they are installed and visible to CMake, e.g.:
   ```bash
   cmake .. -DHDF5_ROOT=/path/to/hdf5 -DFFTW3_ROOT=/path/to/fftw
   make
   ```
4. **Tests** (if enabled):
   ```bash
   ctest --output-on-failure
   ```

This should produce an executable, for instance:  
```
build/agoge_main
```

---

## Running Agoge

By default, `agoge_main` will:

1. **Initialize** a 3D domain (e.g., \(64^3\) or as specified in `main.cpp` or a config file).  
2. **Set** initial conditions (a uniform field with small perturbations, or a user-defined profile).  
3. **Time-step** the compressible Euler equations with or without self-gravity.  
4. **Write** final output to `agoge_final.h5` (if HDF5 is enabled).

Example usage:
```bash
./agoge_main
```
Depending on how the code is written, you might pass arguments like the grid size, time steps, or use a parameter file. If you’re using the **Poisson solver** for self-gravity, ensure it’s enabled during compilation (e.g., `-DAGOGE_USE_GRAVITY=ON` in CMake).

---

## Output & Visualization

- The solver writes fields (density, momentum components, energy, potential, etc.) to **HDF5** files.  
- You can visualize them with **yt**, **ParaView**, or any HDF5-compatible post-processing tool.  
- For instance, with Python and yt:
  ```python
  import yt
  ds = yt.load("agoge_final.h5")
  slc = yt.SlicePlot(ds, "z", "density")
  slc.save()
  ```

---

## Extending & Customizing

- **Higher-Order Fluxes**: Replace or augment the finite-difference flux with more advanced Riemann solvers, limiters, or WENO schemes.  
- **Boundary Conditions**: Add outflow, reflecting, or custom BC logic in `EulerSolver.cpp`.  
- **GPU Offload**: Use frameworks like **AMReX**, CUDA, or OpenMP target offload. Annotate the solver routines to run on devices, and manage data.  
- **Additional Physics**: Viscous terms, magnetic fields, multi-species, or advanced EOS can be added by extending the data structures and PDE solvers.

---

## Contributing

Contributions, bug reports, and feature requests are welcome! Please submit a pull request or open an issue. This project is used in an educational setting (e.g., **CMSE 822**), so it’s a great test bed for learning HPC and scientific software engineering best practices.

---

## License

[MIT License](LICENSE) — Agoge is free for personal, educational, and commercial use. If you use or modify Agoge for your own research or teaching, a citation or link back to this repository is greatly appreciated.

---

## Acknowledgments

- **CMSE 822** at [Your Institution] for providing the impetus to develop a pedagogical compressible flow code.  
- [FFTW](http://fftw.org/) and [HDF5](https://portal.hdfgroup.org/) for core HPC libraries.  
- Everyone in the HPC and astrophysics communities whose open-source contributions inspired Agoge’s structure.
