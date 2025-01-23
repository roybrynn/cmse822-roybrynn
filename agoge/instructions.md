# Agoge QuickStart

## 1. Directory Layout (Minimal)

```
agoge/
├── CMakeLists.txt
├── src/
│   ├── CMakeLists.txt
│   ├── main.cpp
│   ├── PerformanceMonitor.hpp
│   ├── PerformanceMonitor.cpp
│   ├── Field3D.hpp / Field3D.cpp
│   ├── EulerSolver.hpp / EulerSolver.cpp
│   ├── GravitySolver.hpp / GravitySolver.cpp
│   ├── HDF5_IO.hpp / HDF5_IO.cpp
│   └── Config.hpp
├── problems/
│   ├── CMakeLists.txt (optional subdir build file)
│   ├── Problem.hpp
│   ├── ProblemFactory.hpp / ProblemFactory.cpp
│   ├── SodShockTube.hpp / SodShockTube.cpp
│   ├── GravityCollapse.hpp / GravityCollapse.cpp
│   └── ... (more problems)
```

1. **`src/`**: Core solver (Euler, Gravity, HDF5 I/O, `main.cpp`) plus performance monitor.  
2. **`problems/`**: Abstract `Problem` interface, factory, and problem-specific setups.

---

## 2. Building with CMake

1. **Install** libraries: HDF5, FFTW (if required).  
2. **Clone** or place `agoge/` on your system.  
3. **Configure & Build**:
   ```bash
   cd agoge
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make
   ```
4. This should produce an **executable**, e.g., `agoge_run` (depending on your `src/CMakeLists.txt` setup).

---

## 3. Running Agoge

1. **Select a Problem** by name (e.g., “sod”, “collapse”).  
2. **Optional**: Provide grid dimensions and domain size.  
3. **Example**:
   ```bash
   ./agoge_run sod 128 1 1 1.0 1.0 1.0
   ```
   This selects the “sod” problem, with `Nx=128`, `Ny=1`, `Nz=1`, domain lengths = 1.0.  
4. Agoge will:
   - Initialize the chosen problem’s fields.  
   - Optionally solve Poisson’s equation if the problem requests gravity.  
   - Evolve with an RK2 integrator and adaptive time step (CFL).  
   - Write **`agoge_final.h5`** on completion.  
   - Print a **performance report** (timing of major routines).

---

## 4. Inspecting Output

- **HDF5** file `agoge_final.h5` contains density, momentum, energy, etc.  
- Visualize with tools like **yt**, **ParaView**, or **h5dump**.  
- If you need to customize I/O or run multiple output dumps, adjust `HDF5_IO.cpp` or add calls in the main loop.
