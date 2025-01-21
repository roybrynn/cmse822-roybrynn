Below is an **example layout** for a fully functional **Agoge** solver directory, using **CMake** to build and link all the source files (including FFTW and HDF5). You can clone this structure, copy in the exact source/header contents from our earlier steps, and build the code directly.

```
agoge/
├── CMakeLists.txt
├── src/
│   ├── CMakeLists.txt
│   ├── Config.hpp
│   ├── Field3D.hpp
│   ├── Field3D.cpp
│   ├── EulerSolver.hpp
│   ├── EulerSolver.cpp
│   ├── GravitySolver.hpp
│   ├── GravitySolver.cpp
│   ├── HDF5_IO.hpp
│   ├── HDF5_IO.cpp
│   └── main.cpp
└── README.md  (optional)
```

Below are the **CMake** files plus a reminder of the **source/header** files. Make any tweaks necessary for your environment (e.g., adjusting library names, paths, or version requirements).

---

## 1. Top-Level `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.15)
project(agoge LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optionally find FFTW and HDF5 libraries
# (You can remove these if you manually specify them later, or if your code doesn't need them)
find_package(FFTW3 REQUIRED)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)

# Add the 'src' subdirectory which contains our source code
add_subdirectory(src)

# (Optional) If you want tests or other directories, add_subdirectory(tests) here
```

**Notes**:
- If your system doesn’t have CMake modules for FFTW/HDF5, you might need to set environment variables or pass flags like `-DFFTW3_ROOT=/path/to/fftw -DHDF5_ROOT=/path/to/hdf5` when running `cmake`.
- Alternatively, you could comment out the `find_package` lines if you don’t require them or are manually specifying library paths.

---

## 2. `src/CMakeLists.txt`

```cmake
# Go into the 'src' directory for building the core library and the executable
file(GLOB AGOGE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    # Alternatively list them explicitly:
    # Field3D.cpp
    # EulerSolver.cpp
    # GravitySolver.cpp
    # HDF5_IO.cpp
    # (main.cpp is built separately for the executable below)
)

# Remove 'main.cpp' from the library if it is included in the wildcard
list(REMOVE_ITEM AGOGE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

# Create a library from the solver source files
add_library(agoge_lib ${AGOGE_SOURCES})

# Add link dependencies
# We'll assume find_package(FFTW3 REQUIRED) created FFTW3::fftw3 target
# and find_package(HDF5 ...) created HDF5::HDF5, HDF5::HDF5_CPP, etc.
target_link_libraries(agoge_lib
    PUBLIC
        FFTW3::fftw3        # or FFTW3::fftw3f if using float version
        HDF5::HDF5
        HDF5::HDF5_CPP
        HDF5::HDF5_HL
)

# Make sure the library can see the headers in this directory
target_include_directories(agoge_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create the main executable
add_executable(agoge_run main.cpp)
target_link_libraries(agoge_run PRIVATE agoge_lib)

# Optionally set RPATH or other linking details here if needed.
```

**Notes**:
- The code looks for all `.cpp` files in `src/` and puts them in `agoge_lib`, except `main.cpp` which is used for the main executable.
- We assume both FFTW and HDF5 library targets are named in a standard way (e.g., `FFTW3::fftw3`, `HDF5::HDF5_CPP`). Depending on your system, you might need to tweak these names or add additional libraries.

---

## 3. Source and Header Files

Here is a final recap of *all* the source/header files that should appear in your `src/` directory. The contents match what we generated previously. (You may copy-paste exactly if you wish.)

1. **`Config.hpp`**  
   - Declares global constants (e.g., `gamma_gas`, `G`, `use_gravity`).  

2. **`Field3D.hpp`** / **`Field3D.cpp`**  
   - Declares and defines a SoA data structure with `rho`, `rhou`, `rhov`, `rhow`, `E`, `phi`, along with `(Nx, Ny, Nz, dx, dy, dz)`.  

3. **`EulerSolver.hpp`** / **`EulerSolver.cpp`**  
   - Declares and defines `computeL(...)` for \(-\nabla \cdot F\) plus gravitational source terms if needed.  
   - Provides a 2-stage `runRK2(...)` integrator.  

4. **`GravitySolver.hpp`** / **`GravitySolver.cpp`**  
   - Declares and defines `solvePoissonFFT(...)`, which uses **FFTW** to solve \(\nabla^2 \phi = 4 \pi G \rho\) under periodic BCs.  

5. **`HDF5_IO.hpp`** / **`HDF5_IO.cpp`**  
   - Declares and defines `writeFieldHDF5(...)` (and `readFieldHDF5(...)` if desired) using the **HDF5 C++** API.  
   - Writes/reads the 3D fields `(rho, rhou, rhov, rhow, E, phi)` in `(Nz, Ny, Nx)` shape.  

6. **`main.cpp`**  
   - A simple driver that sets up `Field3D`, initializes, runs an RK2 loop (with optional `solvePoissonFFT`), and writes out `agoge_final.h5`.  

**Example**:
```
src/
 ├── Config.hpp
 ├── Field3D.hpp
 ├── Field3D.cpp
 ├── EulerSolver.hpp
 ├── EulerSolver.cpp
 ├── GravitySolver.hpp
 ├── GravitySolver.cpp
 ├── HDF5_IO.hpp
 ├── HDF5_IO.cpp
 └── main.cpp
```

**Implementation**: The exact code for each file has been given in previous steps (the ChatGPT-generated solutions). You can place them in these filenames and build them with the provided CMake setup.

---

## 4. Building and Running

1. **Install** or **load** modules for FFTW3 and HDF5 with C++ support on your system.  
2. **Clone or copy** the `agoge` directory onto your machine.  
3. **Configure** and **build**:
   ```bash
   cd agoge
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make
   ```
4. **Run**:
   ```bash
   ./src/agoge_run
   ```
   or the executable may appear as simply `./agoge_run` in the `build/src/` directory (depending on CMake version and generator).

You can pass optional command-line arguments if your `main.cpp` expects them (like `Nx Ny Nz Lx Ly Lz`). If everything is set up correctly, the code will run, output progress statements, and write an **HDF5** file named **`agoge_final.h5`**.

---

### Summary

With these **two** CMakeLists files and the **ten** source/header files, you have a **complete** HPC-ready, modular **Agoge** codebase that demonstrates:

- **Compressible Euler** equations with or without gravity.  
- **FFT-based Poisson solver** for self-gravity.  
- **HDF5** output for final (or intermediate) states.  

Simply copy the provided content into the corresponding files, adjust any library paths or `find_package` directives for your environment, and run the steps above to build and execute your **Agoge** solver. Happy computing!