# Compile-Time Constants Plan for Agoge

## Overview

The plan is to convert the runtime parameters `Nx`, `Ny`, `Nz`, and `Nghost` to compile-time constants using preprocessor macros and modern C++ techniques. This will allow for better performance through compiler optimizations, fixed-size array optimizations, and potentially improved vectorization, while maintaining build-time flexibility.

## Augment Header File

Add to the #file:agoge/include/agoge/Config.hpp new compile-time domain configuration:

```cpp
// Config.hpp would contain:
// Default grid dimensions if not provided at compile time
#ifndef AGOGE_NX
#define AGOGE_NX 64
#endif

#ifndef AGOGE_NY
#define AGOGE_NY 64
#endif

#ifndef AGOGE_NZ
#define AGOGE_NZ 64
#endif

#ifndef AGOGE_NGHOST
#define AGOGE_NGHOST 1
#endif

namespace agoge {
namespace config {
    // C++ constants that wrap the preprocessor definitions for type safety
    constexpr int DEFAULT_NX = AGOGE_NX;
    constexpr int DEFAULT_NY = AGOGE_NY;
    constexpr int DEFAULT_NZ = AGOGE_NZ;
    constexpr int DEFAULT_NGHOST = AGOGE_NGHOST;


    // existing code....
}
}
```

## Files to Modify

### Core Class Modifications

1. **include/agoge/Field3d.hpp**:
   - Convert `Field3D` to template class with grid dimensions as template parameters
   - Replace dynamic `std::vector` with `std::array` where appropriate
   - Update indexing functions to leverage compile-time constants

2. **src/Field3d.cpp**:
   - Implement template specializations
   - Optimize boundary condition loops with compile-time dimensions
   - Optimize MPI buffer allocation with fixed sizes

### Dependent Files

3. **include/agoge/Config.hpp**:
   - Add references to the new compile-time constants
   - Maintain backward compatibility with runtime configuration

4. **src/EulerSolver.cpp and EulerSolver.hpp**:
   - Update to work with templated Field3D class
   - Optimize loops for fixed-size arrays

5. **src/GravitySolver.cpp and GravitySolver.hpp**:
   - Update Poisson solver to use template parameters
   - Optimize FFT and DFT operations for fixed grid sizes

6. **src/main.cpp**:
   - Create appropriate Field3D template instantiation
   - Maintain compatibility with runtime parameters via command line
   - Add compile-time configuration selection

7. **Problem implementations**:
   - Update all problem classes to work with templated Field3D

## Implementation Considerations

1. **Performance Optimizations**:
   - Fixed-size arrays allow better compiler loop unrolling and vectorization
   - Indexing operations can be optimized at compile time
   - Array boundary checks can be eliminated

2. **Backward Compatibility**:
   - Provide type aliases for common grid sizes
   - Consider keeping a non-templated version for runtime flexibility

3. **Memory Management**:
   - Large grids may need special handling to avoid stack overflows
   - Consider specialized allocators for large compile-time arrays

4. **MPI Communication**:
   - Use fixed-size buffers for MPI operations
   - Optimize packing/unpacking operations with compile-time loop bounds

This refactoring will significantly change the core data structures but could yield substantial performance improvements for typical problem sizes in the agoge application.

<!-- copilot::ignore:start -->
## Prompts 

Collecting workspace information# Refactoring Prompts for Compile-Time Constants Implementation

## 1. Field3d.hpp Refactoring Prompt

```
Refactor the Field3D class in include/agoge/Field3d.hpp to use compile-time constants for grid dimensions. Convert it to a template class with parameters for NX, NY, NZ, and NGHOST. Replace std::vector with std::array for performance optimization.

Specifically:
1. Add include for the new Config.hpp with compile-time constants
2. Convert Field3D to template class: template<int NX=DEFAULT_NX, int NY=DEFAULT_NY, int NZ=DEFAULT_NZ, int NGHOST=DEFAULT_NGHOST>
3. Replace std::vector<double> fields with std::array<double, (NX+2*NGHOST)*(NY+2*NGHOST)*(NZ+2*NGHOST)>
4. Optimize indexing functions to leverage compile-time constants
5. Update constructor to properly initialize the fixed-size arrays
6. Add static constexpr members for total sizes including ghost zones

Keep compatibility with the rest of the codebase and ensure the refactored class maintains all functionality of the original class. Show me only the refactored Field3d.hpp file.
```

## 2. Field3d.cpp Refactoring Prompt

```
Refactor the Field3d.cpp implementation file to work with the new templated Field3D class. The template version needs special handling since implementation must be available in the header or explicitly instantiated.

Please:
1. Move implementation details from Field3d.cpp into Field3d.hpp inside the template class
2. Create explicit template instantiations for common grid sizes (64³, 128³, 256³)
3. Optimize boundary condition loops with compile-time loop bounds using template parameters
4. Update the MPI buffer allocation function to use fixed-size arrays when appropriate
5. Optimize applyBCs() method with compile-time constants for loop bounds
6. Update the packing/unpacking routines to use compile-time dimensions

Ensure all dependencies and includes are properly maintained. Show me the changes needed for Field3d.cpp to work with the new templated version.
```

## 3. Config.hpp Refactoring Prompt

```
Modify the include/agoge/Config.hpp file to add compile-time constants for grid dimensions and maintain backward compatibility with runtime configuration. 

Please:
1. Add the preprocessor definitions for AGOGE_NX, AGOGE_NY, AGOGE_NZ, and AGOGE_NGHOST with default values of 64, 64, 64, and 1
2. Add constexpr constants in the agoge::config namespace that wrap these preprocessor definitions
3. Keep all existing functionality in the file
4. Add utility templates or functions to help with compile-time grid dimension calculations
5. Add comments explaining the purpose of these constants and how to override them

Make sure the changes integrate smoothly with the existing code and maintain backward compatibility. Show me the updated Config.hpp file.
```

## 4. EulerSolver Refactoring Prompt

```
Update the EulerSolver.hpp and EulerSolver.cpp files to work with the new templated Field3D class. Optimize loops and calculations using compile-time constants.

Please:
1. Update function signatures to accept templated Field3D instances
2. Convert computational loops to use compile-time bounds where possible
3. Apply template specializations for common grid sizes to optimize critical routines
4. Update any EulerSolver-specific structures to work with compile-time dimensions
5. Ensure backward compatibility through appropriate type aliases or overloads
6. Add constexpr calculations where possible for flux computations

Focus on performance-critical sections like the flux calculations and conserved-to-primitive variable conversions. Show me both the updated EulerSolver.hpp and the key changes to EulerSolver.cpp.
```

## 5. GravitySolver Refactoring Prompt

```
Refactor the GravitySolver.hpp and GravitySolver.cpp files to leverage the templated Field3D class and optimize FFT and DFT operations with compile-time knowledge of grid dimensions.

Specifically:
1. Update function signatures to accept templated Field3D instances
2. Optimize FFT and DFT operations with compile-time grid sizes
3. Add specializations for power-of-two grid dimensions that can leverage faster FFT algorithms
4. Modify the Poisson solver to use template parameters for dimension-specific optimizations
5. Use compile-time constants to optimize memory allocation and access patterns
6. Update Green's function calculation to leverage compile-time optimizations

Focus on the performance-critical FFT operations and ensure that the optimizations don't reduce numerical accuracy. Show me the key changes needed for both GravitySolver.hpp and GravitySolver.cpp.
```

## 6. Main.cpp Refactoring Prompt

```
Update src/main.cpp to work with the new templated Field3D class, maintaining compatibility with runtime parameters while enabling compile-time optimizations.

Please implement:
1. Create appropriate Field3D template instantiation based on compile-time constants
2. Add logic to select between runtime and compile-time configurations based on input parameters
3. If runtime dimensions match compile-time constants, use the optimized template version
4. Otherwise, fall back to a runtime-specified version with dynamic arrays
5. Update all Field3D usages throughout main.cpp to work with the templated version
6. Add appropriate type deduction and forwarding to maintain interface compatibility

Ensure the changes are minimally invasive to the existing codebase and preserve all current functionality. Show me the key changes to src/main.cpp needed for integration with the templated Field3D.
```

## 7. Problem Classes Refactoring Prompt

```
Update the problem implementation classes to work with the templated Field3D class. Focus on maintaining compatibility while enabling compile-time optimizations.

For each problem class:
1. Update function signatures to accept templated Field3D instances
2. Modify initialization functions to properly setup templated fields
3. Update any problem-specific calculations to leverage compile-time dimensions
4. Ensure boundary condition application is compatible with the templated field
5. Add appropriate type deduction and forwarding to maintain interface compatibility
6. Test for compatibility with both compile-time and runtime field configurations

Focus on the Problem base class, ProblemRegistry, and key implementations like SodShockTube and GravityCollapse. Show me the key changes needed for the problem classes to work with the templated Field3D.
```

## 8. BoundaryManager Refactoring Prompt

```
Update the BoundaryManager.cpp file to handle templated Field3D classes and optimize boundary condition application using compile-time constants.

Please:
1. Update function signatures to accept templated Field3D instances
2. Convert boundary condition loops to use compile-time bounds where possible
3. Optimize MPI communication patterns with compile-time knowledge of field dimensions
4. Add specialized implementations for common boundary types (periodic, outflow, etc.)
5. Ensure boundary condition handling is compatible with both compile-time and runtime fields
6. Update any utility functions used in boundary management to work with templates

Focus on performance-critical sections like ghost cell updates and MPI exchanges. Show me the key changes needed for BoundaryManager.cpp to work with the templated Field3D.
```
<!-- copilot::ignore:end -->
