# CMSE/CSE 822 Project 4: GPU Computing

OMP offload of Ben's exercises.

## OpenMP Vector Addition 

- vector addition: in class

## OpenMP Jacobi Solver 

In this project, you

- Jacobi solver: on your own. NO starter code. Read S5.5. 

## **Implementing a Jacobi Iterative Solver**

### **Step 1: Understanding the Jacobi Method**
**Textbook Reference:**  
- See **Section 5.5.3** (Computational Form) in Eijkhout's "Intro to HPC".  
- **Equation (5.16)** illustrates the Jacobi iterative scheme explicitly:

\[
x_i^{(t+1)} = a_{ii}^{-1}\left(\sum_{j\ne i}a_{ij}x_j^{(t)} + b_i\right)
\]

This scheme is stationary, as every iteration follows the same update process without dependence on the iteration number (**Section 5.5.1**).

The following steps will walk you through implementing the Jacobi solver in serial. 

### **Step 2: Declaring Constants and Vectors**

Create a new file `jacobi_solver.cpp` and include headers and constants:

```cpp
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "mm_utils.hpp"

// Constants
constexpr double TOLERANCE = 0.001;
constexpr int DEF_SIZE = 1000;
constexpr int MAX_ITERS = 100000;
constexpr double LARGE = 1000000.0;
```

---

### **Step 3: Parsing Inputs and Initialization**

Handle matrix dimension from command line and initialize vectors:

```cpp
int main(int argc, char **argv) {
    int Ndim = (argc == 2) ? std::atoi(argv[1]) : DEF_SIZE;
    std::cout << "Matrix dimension (Ndim) = " << Ndim << std::endl;

    // Use std::vector to allocate storage dynamically.
    std::vector<TYPE> A(Ndim * Ndim);
    std::vector<TYPE> b(Ndim);
    std::vector<TYPE> xnew(Ndim, 0.0);
    std::vector<TYPE> xold(Ndim, 0.0);

    // Generate diagonally dominant matrix A
    initDiagDomNearIdentityMatrix(Ndim, A.data());

    // Initialize b with random values (between 0.0 and 0.5)
    for (int i = 0; i < Ndim; ++i) {
        b[i] = static_cast<TYPE>(std::rand() % 51) / 100.0;
    }
```

---

### **Step 4: Implementing the Jacobi Iteration**

- **Textbook Reference:** Follow the iterative scheme from **Section 5.5.3**, specifically equation (5.16).

```cpp
    double start_time = omp_get_wtime();
    TYPE conv = LARGE;
    int iters = 0;

    while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
        ++iters;

        // Compute new iteration
        for (int i = 0; i < Ndim; ++i) {
            xnew[i] = static_cast<TYPE>(0.0);
            for (int j = 0; j < Ndim; ++j) {
                if (i != j)
                    xnew[i] += <INSERT YOUR CODE HERE>;
            }
            xnew[i] = <INSERT YOUR CODE HERE>;
        }

        // Compute convergence criterion (Euclidean norm of difference)
        conv = static_cast<TYPE>(0.0);
        for (int i = 0; i < Ndim; ++i) {
            TYPE tmp = xnew[i] - xold[i];
            conv += tmp * tmp;
        }
        conv = static_cast<TYPE>(std::sqrt(conv));

        // Swap vectors for next iteration
        std::swap(xold, xnew);
    }

    double elapsed_time = omp_get_wtime() - start_time;
    std::cout << "Converged after " << iters << " iterations in "
              << elapsed_time << " seconds with final convergence = "
              << conv << std::endl;
```

- **Convergence Criterion:** (Textbook, **Section 5.5.7**)  
You test convergence by measuring the Euclidean norm of the difference between successive iterates. The iteration stops when the computed difference is below the tolerance.

### **Step 5: Verifying the Solution**

Multiply the solution vector by the original matrix and compare with vector `b`:

```cpp
    TYPE err = 0.0, chksum = 0.0;

    for (int i = 0; i < Ndim; ++i) {
        xold[i] = static_cast<TYPE>(0.0);
        for (int j = 0; j < Ndim; ++j)
            xold[i] += A[i * Ndim + j] * xnew[j];

        TYPE diff = xold[i] - b[i];
        chksum += xnew[i];
        err += diff * diff;
    }
    err = static_cast<TYPE>(std::sqrt(err));

    std::cout << "Solution verification: Error = " << err
              << ", Checksum = " << chksum << std::endl;

    if (err > TOLERANCE)
        std::cout << "WARNING: Solution error exceeds tolerance!" << std::endl;

    return 0;
}
```


### **Step 6: Compilation and Execution**

Compile using a modern C++20 compiler:

```sh
g++ -std=c++20 -fopenmp jacobi_solver.cpp -o jacobi_solver
```

Run with default or specified dimension:

```sh
./jacobi_solver
# or with specified dimension
./jacobi_solver 2500
```


## **Parallelization with OpenMP**

Starting from your serial Jacobi solver, we will now proceed with implementing parallelization via OpenMP for both the CPU and GPU. 

### **Step 1: Understanding the Jacobi Method (Review)**

The Jacobi solver iteratively solves a system of linear equations:

\[
A x = b
\]

This method splits the coefficient matrix \( A \) into three parts:

- Lower triangular part \( L \)
- Diagonal part \( D \)
- Upper triangular part \( U \)

The original equation is rearranged to isolate the unknown vector \( x \):

\[
D x = b - (L + U) x \quad\rightarrow\quad x_{\text{new}} = D^{-1}(b - (L + U)x_{\text{old}})
\]

The iterative Jacobi algorithm updates \( x \) repeatedly until convergence is reached, typically determined by a tolerance on the change in successive \( x \).

- **Advantages**: Simple implementation and easy correctness verification.
- **Disadvantages**: Convergence can be slow, and it requires \( A \) to be diagonally dominant to guarantee convergence.


### **Step 2: Starting with Your Serial Implementation**

Begin with your working serial Jacobi solver. The basic iterative structure is:

```cpp
while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
    iters++;

    for (int i = 0; i < Ndim; i++) {
        xnew[i] = 0.0;
        for (int j = 0; j < Ndim; j++) {
            if (i != j)
                xnew[i] += <YOUR CODE HERE>;
        }
        xnew[i] = <YOUR CODE HERE>;
    }

    conv = 0.0;
    for (int i = 0; i < Ndim; i++) {
        TYPE tmp = xnew[i] - xold[i];
        conv += tmp * tmp;
    }
    conv = sqrt(conv);

    std::swap(xold, xnew);
}
```

Verify correctness thoroughly before parallelizing by checking:

- Multiply the solution vector by \( A \) and ensure the result closely matches \( b \).

### **Step 3: Initial Parallelization using OpenMP**

To begin parallelization, use OpenMP’s GPU offloading directives. OpenMP’s `target` directive offloads code to a GPU device, and `loop` specifies loops for parallel execution:

### Implementation

Insert OpenMP pragmas around computational loops:

```cpp
while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
    iters++;

    #pragma omp target map(tofrom: xnew[0:Ndim], xold[0:Ndim]) \
                       map(to: A[0:Ndim*Ndim], b[0:Ndim])
    #pragma omp loop
    for (int i = 0; i < Ndim; i++) {
        xnew[i] = 0.0;
        for (int j = 0; j < Ndim; j++) {
            if (i != j)
                xnew[i] += <YOUR CODE HERE>;
        }
        xnew[i] = <YOUR CODE HERE>;
    }

    conv = 0.0;
    #pragma omp target map(to: xnew[0:Ndim], xold[0:Ndim]) \
                       map(tofrom: conv)
    #pragma omp loop reduction(+:conv)
    for (int i = 0; i < Ndim; i++) {
        TYPE tmp = xnew[i] - xold[i];
        conv += tmp * tmp;
    }
    conv = sqrt(conv);

    std::swap(xold, xnew);
}
```

#### Testing
- Validate parallel correctness by comparing outputs with your serial implementation.

#### Performance Analysis
- Record runtime.
- Note: This initial parallelization often leads to poor performance due to heavy data transfers between the host (CPU) and GPU each iteration.


### **Step 4: Improving Data Transfer using `target data` Regions**

Repeated copying of data between host and GPU is costly. To reduce data movement overhead, create a persistent data region on the GPU using OpenMP’s `target data` directive:

#### Implementation

Surround the iterative loop with a `target data` region:

```cpp
#pragma omp target data map(tofrom: xold[0:Ndim], xnew[0:Ndim]) \
                        map(to: A[0:Ndim*Ndim], b[0:Ndim])
while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
    iters++;

    #pragma omp target
    #pragma omp loop private(j)
    for (int i = 0; i < Ndim; i++) {
        xnew[i] = 0.0;
        for (int j = 0; j < Ndim; j++) {
            if (i != j)
                xnew[i] += <YOUR CODE HERE>;
        }
        xnew[i] = <YOUR CODE HERE>;
    }

    conv = 0.0;
    #pragma omp target map(tofrom: conv)
    #pragma omp loop reduction(+:conv)
    for (int i = 0; i < Ndim; i++) {
        TYPE tmp = xnew[i] - xold[i];
        conv += tmp * tmp;
    }
    conv = sqrt(conv);

    std::swap(xold, xnew);
}
```

#### Testing and Performance Analysis
- Verify correctness again.
- Measure runtime improvements. Expect a significant reduction due to fewer data transfers.


### **Step 5: Eliminating Branches (Branchless Implementation)**

GPU architectures perform poorly with branching. Conditional checks cause divergence among GPU threads, slowing execution significantly. Instead, convert branching conditions into arithmetic operations.

Replace:
```cpp
if (i != j)
    xnew[i] += A[i*Ndim + j] * xold[j];
```

with the branchless form:
```cpp
xnew[i] += A[i*Ndim + j] * xold[j] * static_cast<TYPE>(i != j);
```

#### Testing and Performance Analysis
- Test correctness carefully again.
- Measure performance improvements from reduced branching.

### **Step 6: Optimizing Memory Access Patterns (Coalescing)**

GPU memory access is fastest when threads access contiguous memory locations (coalesced access). Change your indexing patterns to ensure consecutive threads read consecutive memory addresses:

Change:
```cpp
xnew[i] += A[i*Ndim + j] * xold[j] * static_cast<TYPE>(i != j);
```

to coalesced access pattern:
```cpp
xnew[i] += A[j*Ndim + i] * xold[j] * static_cast<TYPE>(i != j);
```

This simple transposition drastically enhances performance by aligning memory access with GPU architecture.

#### Testing and Performance Analysis
- Validate correctness again thoroughly.
- Record runtime; expect further significant improvements.


### **Step 7: Final Validation and Comprehensive Performance Comparison**

Perform rigorous correctness validation:

- Compare GPU-parallelized solver results with your initial serial implementation across multiple test cases and matrix sizes.
- Benchmark serial versus each parallel optimization step clearly and quantitatively.


### **Step 8: Documenting Results and Analysis**

Summarize your optimizations and clearly state each step’s impact on runtime:

| Optimization Step             | Typical Performance Impact  |
|-------------------------------|-----------------------------|
| Initial parallelization       | Poor (due to heavy data movement) |
| Persistent data regions       | Large improvement (~5-7x speedup) |
| Branchless implementation     | Moderate improvement (~30%) |
| Coalesced memory accesses     | Significant improvement (~50%) |


## **Step 10: Submission and Presentation**

Prepare your final submission including:

- Well-commented code demonstrating each optimization step.
- Performance analysis with clear benchmarking results.
- Discussion of observed performance trends and conclusions.
