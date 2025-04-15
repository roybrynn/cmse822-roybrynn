# CMSE 822 Group Quest 3: Parallel Gravity and the FFT Fortress (Enhanced Code-Focused Draft)

## Setting

In the hidden depths of the computational labyrinth, you and your party confront the **FFT Dungeon**—its final guardians, the Fourier Phantoms, challenge you to parallelize the gravitational Poisson solver. This quest demands refactoring **`GravitySolver.cpp`** (and its companion `GravitySolver.hpp`) to invoke **MPI**-based parallel FFT. Achieve near-\( N \log N \) scaling, or the spectral gates will remain closed forever.

---

## Group Member Roles

To better navigate this refactor, assign clear roles:

- **Navigator**: Oversees to-do items, organizes sessions, ensures a smooth progress timeline.  
- **Lorekeeper**: Documents code changes, design rationale, and performance analyses.  
- **Bard**: Delivers a concise 10-minute presentation of your epic success, from domain decomposition to final benchmarks.

**Navigator**:  
**Lorekeeper**:  
**Bard**:

---

## Quest Objectives

Your overarching goal is to inject **MPI parallelism** into the existing **FFT-based Poisson solver** in `GravitySolver.cpp`, ensuring:

1. **Parallel 3D FFT**:  
   - Adapt `forwardFFT3D` and `inverseFFT3D` from single-threaded (serial) calls to a distributed-memory MPI approach.
   - Achieve near-\( N \log N \) time complexity as you scale the problem size.

2. **Efficient Domain Decomposition**:  
   - Implement a **slab** (or pencil) decomposition of the 3D grid.  
   - Use **MPI collectives** (e.g., `MPI_Alltoallv`) for transposing data between 1D FFT passes.

3. **Performance Demonstration**:  
   - Show **speedups** as you increase the number of ranks.  
   - Validate correctness with physical test problems, ensuring the final solutions match the serial reference.

---

## Setting the Scene: GravitySolver.cpp

In `GravitySolver.cpp` (see the **`agoge::gravity`** namespace), you have:

1. **Naive DFT Routines**  
   - `naiveForwardDFT3D` / `naiveInverseDFT3D` (O\(N^6\) complexity). Located around lines **30–150** in the snippet provided, these are purely educational but could also be parallelized if you want to compare performance.  

2. **Cooley–Tukey 1D FFT**  
   - The static function `fft1D(std::complex<double> *data, int N, bool inverse)` around lines **160–220**.  
   - Used by `forwardFFT3D` and `inverseFFT3D` to perform dimension-by-dimension transforms.

3. **Forward and Inverse 3D FFT**  
   - `forwardFFT3D` (~lines **260–320**)  
     - Performs **1D FFT** along x, then y, then z (or vice versa), by iterating through slices in memory.  
   - `inverseFFT3D` (~lines **340–400**)  
     - Reverses the transform in z, y, x.  
   - Both currently assume the entire array is local.

4. **Data Structures**  
   - `std::vector<std::complex<double>> cplxData(Nx * Ny * Nz)` for the frequency domain.  
   - `std::vector<double> realBuf(Nx * Ny * Nz)` for the real domain.  
   - Indices typically follow `idx = x + Nx * (y + Ny * z)`.

5. **`solvePoisson(...)`**  
   - Uses the 3D transforms to switch between real-space density (`realBuf`) and frequency-space potential (`cplxData`).  
   - Apply your parallel code inside or just after the forward/inverse FFT calls.  

**Hint**: You can keep the same function signatures but adapt them internally to handle parallel data.

---

## Starter Code Guidance and Hints (Specific to GravitySolver.cpp)

1. **Examine `forwardFFT3D`**  
   - This function loops over `(z, y, x)` in triple-nested loops to call `fft1D`. For a slab decomposition by `z`, each rank would own a chunk of `z` slices.  
   - After the first dimension transform, you’ll likely need to **redistribute** data so each rank has the correct slices for the next dimension. This is where `MPI_Alltoallv` can help.

2. **Examine `inverseFFT3D`**  
   - This function runs the inverse pass in the opposite dimension order (z, then y, then x).  
   - The same data redistribution logic from `forwardFFT3D` applies, but in reverse order.

3. **Check the Index Macros**  
   - The code uses `x + Nx * (y + Ny * z)` extensively. Ensure that after distribution, each rank modifies these index calculations to reflect local subdomain offsets (for example, `local_z_start = rank * (Nz / size)`).

4. **Slab vs. Pencil**  
   - A **slab** approach, focusing on distributing the `z` dimension, is simpler to implement in `forwardFFT3D` / `inverseFFT3D`. Start there.  
   - For advanced HPC: a **pencil** decomposition requires a 2D distribution (e.g., splitting x and y or y and z). This is more complex but can scale better for large \(p\).

---

## Plan Your Parallel Implementation

### Step 1: MPI Initialization and Data Setup

- **Insert MPI_Init** in your main driver (not necessarily in `GravitySolver.cpp`).  
- Inside `solvePoisson(...)`, figure out each rank’s local range in `z` (for slab decomposition):

```cpp
  int localNz = Nz / size;
  int zStart = rank * localNz;
  // Allocate local memory for cplxData, realBuf, etc.
```

- Hint: If using C++, you can define std::vector<std::complex<double>> localCplxData(localNx *localNy* localNz); or keep the full size and handle only the portion you own.

### Step 2: Parallelize the 1D FFT Calls

- For the first pass (say, along z):
- Each rank calls fft1D(...) on its sub-block of size localNz for each (x, y).
- No communication is needed for that step because each rank already holds its portion of z.

### Step 3: Transpose / Redistribute for Next Dimension (y)

- After finishing the z pass, data must be rearranged so each rank can do 1D FFT in y.
- For a slab approach, each rank will now need some part of all z slices in the appropriate layout for y.
- Use MPI_Alltoallv to share blocks among ranks. Create buffers that gather your local portion of the next dimension lines, then send them to the appropriate ranks.

### Step 4: 1D FFT Along y

- Once data is rearranged, call fft1D on each local chunk.
- Then transpose again for the x dimension.

### Step 5: 1D FFT Along x

- Perform local 1D FFT on the x dimension.
- The result is the fully transformed data in the frequency domain.

Step 6: Inverse FFT (Mirror Steps)

- For inverseFFT3D, invert the dimension order: x → transpose → y → transpose → z.
- Carefully maintain the same logic but in reverse.

## Integration Test (Trial of Strength): Dust Collapse

- A simplified check of gravitational collapse. The code should produce the same final distribution of potential and density in both serial and parallel modes.
- Ensure that when solvePoisson(...) is used, total energy remains nearly constant.
- Compare results from serial vs. parallel runs on small domain sizes (e.g., (16^3) or (32^3)).

Hint: Print sums (like (\sum \rho), (\sum \phi)) at each time step to ensure consistency across ranks.

## Submission Requirements

 1. Codebase

- A fully MPI-enabled GravitySolver.cpp demonstrating how you refactored forwardFFT3D / inverseFFT3D.
- Comments describing how each dimension’s data is distributed and how you used collectives.

 2. Quest Chronicle (Final Report)

- Document your domain decomposition approach (slab or pencil).
- Show performance scaling plots (run times vs. ( N \log N )) for different domain sizes (e.g., (64^3, 128^3, 256^3)) and various ranks.

 3. Presentation (The Bard’s Tale)

- A 10-minute talk with slides summarizing your parallelization journey: distribution plan, indexing modifications, performance data.

C++ Syntax Tip: When gathering sub-blocks for transposes, you can use:

```cpp
// Example of using a temporary buffer:
std::vector<std::complex<double>> sendBuf(numElementsToSend);
std::vector<int> sendCounts(size), sendDispls(size);
// fill sendBuf from localCplxData
MPI_Alltoallv(sendBuf.data(), sendCounts.data(), sendDispls.data(),
              yourMPIType, recvBuf.data(), recvCounts.data(),
              recvDispls.data(), yourMPIType, MPI_COMM_WORLD);
```

Make sure each rank carefully computes the correct offsets (sendDispls/recvDispls) before the call.

## Performance Challenge – “Conquer the Parallel FFT Fortress”

This is it—**the Final Boss Battle** against the **Parallel Performance Dragon**! Compete with your fellow adventurers to **achieve the best performance** for your parallel FFT-based Poisson solver under a standardized test:

- **4 nodes of `amd20` on HPCC**, with each node running multiple MPI ranks.

🏆 **Victory Conditions:**  
Maximize the **Figure of Merit (FoM)** using the legendary formula:  
\[
\text{FoM} = \frac{\text{(Total Grid Cells Transformed)} \times \text{(Time Steps)}}{\text{Wall-Clock Time}}
\]  

### Rules of Engagement

1. You may use **any parallel version** of your `GravitySolver.cpp` code to complete this trial.
2. Ensure correctness by validating the output of your solver against the serial reference implementation.
3. Optimize your implementation for **scalability** and **load balancing** across all MPI ranks.

### Grading

This challeng will be worth up to **10% extra credit toward your team's grade** for this project will be determined by your relative performance in this challenge:

- The team with the **highest FoM** will receive the full 10%.
- The team with the **lowest FoM** will receive 0%.
- All other teams will be scored proportionally based on their performance ranking.

### Tips for Victory

- Profile your code to identify bottlenecks in communication (e.g., `MPI_Alltoallv`) or computation (e.g., 1D FFT calls).
- Experiment with different domain decomposition strategies (e.g., slab vs. pencil) to improve scalability.
- Use optimized libraries like **FFTW-MPI** if allowed, or fine-tune your custom implementation.
- Test your solver on increasingly larger grids (e.g., \(128^3\), \(256^3\), \(512^3\)) to demonstrate scaling.

May your transforms be swift, your communication efficient, and your FoM legendary!

## Final Encouragement

Keep a watchful eye on indexing errors. Start small and gradually scale. Leverage `std::cout << "Rank: " << rank ...` prints to debug partial arrays. Once the basics work, measure performance on bigger domains to confirm near-\( N \log N \) scaling.

With this code-centered roadmap, you should bravely conquer the FFT Dungeon in `GravitySolver.cpp`, forging an MPI-based solver that’s both correct and fast. Good luck, adventurers—may your transforms be swift and your load balancing strong!
