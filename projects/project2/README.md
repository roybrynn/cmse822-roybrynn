# **Project 2: (Domain) Divide and Conquer**

## **Overview**  

In this project, you will embark on an epic quest to implement **MPI-based domain decomposition** in Agoge's **finite-difference (FD) hydrodynamic solver**, mastering the arcane arts of distributed-memory parallelism. You will wield **blocking and non-blocking point-to-point communication spells** to perform **halo exchanges**, channel the power of **latency hiding**, and decipher performance mysteries using **Intel VTune** and the **Empirical Roofline Toolkit (ERT)**.  

---

## **Project Breakdown**

### **Task 1: Implement MPI-Based Domain Decomposition**

Your first challenge is to **decompose** the computational grid into **MPI subdomains**, assigning each process dominion over its own **local portion** of the simulation domain.  
Each process must:  

1. Compute hydrodynamic updates on its local grid.  
2. Exchange **halo/ghost-cell** data with neighboring MPI ranks using **point-to-point communication**.  
3. Apply **boundary conditions** at process edges. 

In the [Appendix](#appendix-implementation-hints) of this document, I give you an example implementation plan for achieving this in Agoge. 
You may choose to follow this plan closely, or modify it as you see fit.    

🏰 **Trials of the Domain Decomposition Dungeon:**  

- Modify the **agoge FD solver** to support **distributed-memory parallelism** using MPI.  
- Implement **three communication strategies** for the **halo exchange**:
  1. **Blocking point-to-point MPI** (`MPI_Send`, `MPI_Recv`)  
  2. **Non-blocking MPI** (`MPI_Isend`, `MPI_Irecv` with `MPI_Waitall`)  
  3. **Latency hiding**: Overlap computation and communication  
- Verify correctness by running **small test cases** and **checking numerical consistency** with the serial solver.  

🔍 **Reference Scrolls (Exercises)**: 4.1 (basic message passing), 4.2 (scaling ping-pong communication), 4.13 & 4.14 (non-blocking MPI).  

---

### **Task 2: Performance Analysis with Empirical Roofline Toolkit (ERT)**

Venture forth to the **Oracle of Empirical Roofline**, where you will seek wisdom regarding the **memory bandwidth** and **floating-point performance** in parallel solver on **different node types** in **ICER HPCC**.  

📜 **Steps to complete:**  

1. **Run ERT** with MPI enabled on three different HPCC node types.  
2. Compare performance across:  
   - Single-node vs. multi-node execution.  
   - Different message sizes and decomposition strategies.  
3. Generate **roofline plots** and interpret where your solver lies on the **memory-bandwidth vs. compute-performance curve**.  

🔍 **Ancient Texts**: ERT GitHub repository ([link](https://github.com/ebugger/Empirical-Roofline-Toolkit)).  

---

### **Task 3: Profiling with Intel VTune**

Descend into the **Crypt of Intel VTune**, where hidden inefficiencies lurk in the darkness, waiting to be uncovered. Use this powerful relic to analyze:  

- **Communication overhead**: MPI synchronization, message transfer times.  
- **Computational load balance** across ranks.  
- **Memory access patterns**: cache misses, NUMA behavior.  

🔎 **Steps to complete:**  

1. **Run VTune CLI** on your **MPI** enable version of Agoge on the `intel18` nodes of HPCC (See tutorial in Project 1).  
2. Identify and **analyze performance bottlenecks**.  
3. Compare profiling results across different MPI implementations (blocking vs. non-blocking).  

---

### **Task 4: Scaling Studies & Sedov Explosion Test**

Take on the **Trials of the Sedov Explosion**, where your solver must scale with increasing power. Conduct **strong and weak scaling tests** to determine if your solver is fit for legendary battles.  

🔥 **Steps to complete:**  

1. **Strong scaling test**: Fix problem size, increase process count. Do this for just a single node but for three different node types.  
2. **Weak scaling test**: Increase problem size with process count. The problem size should increase proportional with process count.
3. Run tests on **multiple HPCC node types** (e.g., different CPUs). For weak scaling, ensure that you go beyond single-node scale.
4. Compare scaling results **(single-node vs. multi-node performance)**.  

📊 **Battle Records:**  

- Plots showing **runtime vs. number of processes**.  
- **Efficiency calculations** (speedup, parallel efficiency).  
- Discussion of **single-node vs. multi-node performance trends**.  

---

### **Task 5: Performance Challenge – “Maximum Performance” Optimization**

This is it—**the Final Boss Battle** against the **Parallel Performance Dragon**! Compete with your fellow adventurers to **achieve the best performance** under a standardized test:

- 4 nodes of `amd20` on HPCC. 

🏆 **Victory Conditions:**  
Maximize **zone-updates per second** using the legendary formula:  
\[
\text{FoM} = \frac{\text{(Total Grid Cells Updated)} \times \text{(Time Steps)}}{\text{Wall-Clock Time}}
\]  

You may use ANY parallel version of your agoge code to complete this trial and produce the large FoM possible. The final 10% of your team's grade for this project will be determined by your relative performance in this challenge. 
The team with the highest FoM will receive all 10%, while the slowest team will receive 0%, and the rest scored in order of their performance FoM. 

⚔️ **Optimization Strategies:**  

You may take _any_ approach you like to optimizing your code to produce the highest FoM, including optimizing the _serial_ portions of the code to achieve higher performance. Other optimization strategies might be:

- **Experiment with different MPI implementations** (blocking vs. non-blocking).  
- **Tune message sizes & domain decomposition strategy**.  
- **Overlap communication with computation (latency hiding)**.  
- **Utilize MPI process pinning and NUMA-aware execution**.  


📊 **Final Deliverables:**  

- **Optimized MPI solver code**.  
- **Performance analysis report** with a discussion of optimizations.  
- **Charts comparing different MPI approaches**.  

## Appendix: MPI Implementation Hints

Below is a **detailed implementation plan** for introducing domain decomposition and MPI parallelism into Agoge’s `EulerSolver` **without** heavily modifying the physics code. The plan focuses on augmenting **`Field3D`** to represent each sub-domain on an MPI rank (including ghost cells) and modifying **`applyBCs()`** to perform halo exchanges between ranks rather than local outflow/periodic copies only. We consider **three** approaches to ghost‐cell communication:

1. **Blocking MPI** calls (simple, direct).
2. **Non‐blocking MPI** calls (uses `MPI_Isend`/`MPI_Irecv`, then `MPI_Wait`).
3. **Non‐blocking / asynchronous** with **latency hiding** (overlap compute and comm).

---

### 1. Overview of the Required Changes

1. **Domain Decomposition** in `Field3D`:  
   - Each MPI rank will store a portion (sub‐block) of the global domain, plus ghost cells.  
   - `Field3D::Nx, Ny, Nz` become local interior sizes.  
   - The bounding box or `(xmin, xmax, …)` becomes local to that subdomain.  
   - `Field3D::applyBCs()` is extended to send/receive ghost‐zone data from neighbors in x, y, z directions.

2. **Minimal Changes** to the Euler/Gravity Solvers:  
   - They keep reading/writing ghost zones as though it is “serial.”  
   - The actual ghost data on each rank is updated via MPI in `Field3D::applyBCs()`.

3. **Main Program** (`main.cpp`) Gains MPI Initialization / Finalization:  
   - `MPI_Init(&argc, &argv)` at the beginning, `MPI_Finalize()` at the end.  
   - Possibly add a small **domain partition** routine that sets each rank’s sub‐domain sizes, bounding box, etc.  
   - Construct the local `Field3D` accordingly.

Below, we describe in detail how to implement each step, then how to do the communication in the three different ways.

---

#### 2. Step‐by‐Step Plan

#### 2.1. Introduce MPI Initialization and Partitioning in `main.cpp`

**(A)** **Add** an `#include <mpi.h>` at the top of `main.cpp`.  
**(B)** **Call** `MPI_Init` at the start and `MPI_Finalize` at the end.  

```cpp
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    // ... existing code ...

    MPI_Finalize();
    return 0;
}
```

**(C)** **Determine** how you want to partition the domain. E.g., a simple 1D decomposition in x:

```cpp
// E.g. read total Nx, rankCount from MPI
int rank, nprocs;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

// Suppose Nx is divisible by nprocs
int NxLocal = Nx / nprocs;
int xStart = rank * NxLocal;    // global start index
int xEnd   = xStart + NxLocal;  // global end index (non-inclusive if you want)
```

Similarly for `y` or `z` if doing 2D/3D decomposition. Or keep it simple (like 1D in x for demonstration).  

**(D)** **Construct** each rank’s local bounding box:
```cpp
double xRange = (xmax - xmin);
double localWidth = xRange / nprocs;
double localXmin = xmin + rank*localWidth;
double localXmax = localXmin + localWidth;

// Then build Field3D with NxLocal, local bounding box...
agoge::BoundingBox localBox = { localXmin, localXmax, ymin, ymax, zmin, zmax };
Field3D Qlocal(NxLocal, Ny, Nz, localBox, ghostCells);
```

**(E)** The rest of `main.cpp` proceeds, but it only deals with `Qlocal`—the PDE solver code works on that sub‐array.  

---

#### 2.2. Modifications to `Field3D.hpp`

**Goal**: Each rank’s `Field3D` describes its local subdomain. We add:

1. **Rank / neighbor** info: we might want `int leftRank, rightRank, upRank, downRank, etc.` or store a small array of neighbors.  
2. **A method** to do the MPI exchanges in `applyBCs()` or a separate method `_exchangeHalos()`.  

For example:

```cpp
class Field3D {
public:
    // existing members

    // rank neighbors
    int rankLeft, rankRight; // if 1D, or we have rankUp, rankDown, rankFront, etc.

    // each rank sub-domain offset in x,y,z
    int globalStartX;  // e.g. xStart
    int globalEndX;    // e.g. xEnd
    // possibly same for y,z if 2D/3D partition

    // new constructor that sets these from the partition logic
    Field3D(int nxLocal, int nyLocal, int nzLocal,
            const BoundingBox &bbox_in, int nghost,
            int leftR, int rightR, int gStartX, int gEndX, // etc
            ...);

    // updated applyBCs with MPI calls
    void applyBCsMPI();

private:
    // an internal function for MPI exchange
    void exchangeHalosBlocking();  // for the blocking approach
    void exchangeHalosNonBlocking();
    void exchangeHalosNonBlockingAsync();
};
```

**Your** existing `applyBCs()` that does outflow or periodic for the local ghost cells is still relevant for boundaries that do not cross MPI ranks. For domain boundaries that are purely physical, we do the old approach. For “internal” boundaries between ranks, we do MPI exchange.  

---

#### 2.3. “applyBCs()” or “applyBCsMPI()” Implementation

Inside `applyBCs()`, we do:

1. **Local** boundary fill for x‐min if `bc_xmin` is physical outflow or periodic that doesn’t cross ranks.  
2. **If** we have an MPI neighbor in x‐min direction, we pack the halo data and send it to rankLeft, then receive from rankLeft.  

**We’ll** show how to do that in the 3 communication modes below.

**IMPORTANT**: We are only doing PDE with ghost cells. The PDE code reads ghost data at each time step. So after we fill physical BC for real domain edges, we do **MPI** to fill interior subdomain boundaries.

---

### 3. Three Communication Approaches

We focus on a **1D** domain decomposition in x for clarity. The logic generalizes to 2D or 3D by exchanging up/down/front/back as well.

#### 3.1. **Case 1: Blocking MPI**  

**Algorithm** (for x direction):  
1. **Pack** your left ghost cells into a buffer, e.g. `sendBufLeft`, then call `MPI_Send` to rankLeft.  
2. **Receive** your left neighbor’s right ghost cells into `recvBufLeft` with `MPI_Recv`.  
3. Copy `recvBufLeft` into your local ghost region for x‐min.  
4. Do the same for x‐max side with rankRight.  
5. Possibly do `MPI_Barrier(MPI_COMM_WORLD)` if you want strict sync.

**Implementation** snippet:

```cpp
void Field3D::exchangeHalosBlocking()
{
    // Step A: pack the left halo of size (nyGhost*nzGhost) for each field
    if(rankLeft >= 0){ // means we have a left neighbor
       // pack e.g. for(var in [rho,rhou,rhov,rhow,E,phi])...
       // MPI_Send(sendBufLeft, count, MPI_DOUBLE, rankLeft, someTag, MPI_COMM_WORLD);

       // receive from rankLeft
       // MPI_Recv(recvBufLeft, count, MPI_DOUBLE, rankLeft, someTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

       // unpack recvBufLeft into local ghost zone
    }

    if(rankRight >=0){ 
       // do similarly
    }
}
```

**Pros**: Easy to code.  
**Cons**: Blocks until both sides are done sending.

---

#### 3.2. **Case 2: Non‐blocking MPI**  

We use `MPI_Isend`, `MPI_Irecv` and eventually `MPI_Waitall`:

```cpp
void Field3D::exchangeHalosNonBlocking()
{
    MPI_Request reqs[4]; // 2 sends, 2 recvs
    int nreq=0;

    if(rankLeft >=0){
       MPI_Isend(sendBufLeft, count, MPI_DOUBLE, rankLeft, 100, MPI_COMM_WORLD, &reqs[nreq++]);
       MPI_Irecv(recvBufLeft, count, MPI_DOUBLE, rankLeft, 101, MPI_COMM_WORLD, &reqs[nreq++]);
    }

    if(rankRight >=0){
       MPI_Isend(sendBufRight, count, MPI_DOUBLE, rankRight, 101, MPI_COMM_WORLD, &reqs[nreq++]);
       MPI_Irecv(recvBufRight, count, MPI_DOUBLE, rankRight, 100, MPI_COMM_WORLD, &reqs[nreq++]);
    }

    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    // then unpack
}
```

**Pros**: Doesn’t force each rank to block in order. Typically faster.  
**Cons**: We still do `MPI_Waitall()` *before* continuing the PDE code. So no overlap with compute.

---

#### 3.3. **Case 3: Non‐blocking/Asynchronous with Latency Hiding**  

Here, we try to overlap PDE compute with communication:

1. **Start** asynchronous sends/receives.  
2. **Immediately** do some local PDE compute steps that do NOT depend on halo data (e.g. interior updates).  
3. **Wait** for the communication to finish with `MPI_Waitall`.  
4. **Update** the PDE in the ghost regions once the data arrives.

**Implementation** outline in the PDE code or time step logic:

```cpp
// inside runRK2(...) or computeL(...):
// Step1: start halo exchange
field.exchangeHalosNonBlockingAsyncStart();

// Step2: compute fluxes in the interior region that doesn't need ghost data
// e.g. for i in [1..nx-2], ...
interiorComputation();

// Step3: wait for halo data
field.exchangeHalosNonBlockingAsyncFinish();

// Step4: compute fluxes near boundary that uses ghost data
boundaryComputation();
```

**Inside** the `Field3D`:

```cpp
void Field3D::exchangeHalosNonBlockingAsyncStart()
{
   // do MPI_Isend, MPI_Irecv with distinct requests
   // store them in e.g. `requestsLeftSend, requestsLeftRecv` members
}

void Field3D::exchangeHalosNonBlockingAsyncFinish()
{
   // call MPI_Waitall on the stored requests
   // then do the unpack
}
```

**Pros**: Potentially the highest performance, as you hide latencies behind interior computation.  
**Cons**: More complex code changes. The PDE solver needs a “two-phase” approach.

---

### 4. Minimal Changes to the Physics Code

You want to avoid rewriting Euler. The main PDE loops remain the same:

1. Where you previously called `applyBCs()` at the end of each stage, you now call something like:
   - `applyBCsMPIBlocking()`
   - or `applyBCsNonBlocking()` + `MPI_Wait...`
   - or the “async” variant in a 2‐phase approach.

2. The logic for PDE boundary conditions (like outflow or periodic at the global edges) still applies, but only if **rank** is at the domain boundary. If you’re an internal rank, you do an MPI exchange with your neighbor.

---

### 5. Summary of the Implementation Plan

1. **Add MPI**:
   1. `MPI_Init/MPI_Finalize` in `main.cpp`.  
   2. A small domain partition logic to compute local Nx, bounding box, plus neighbor ranks.

2. **Extend `Field3D`**:
   1. Add rank neighbor IDs (`rankLeft, rankRight`, etc.).  
   2. Add new methods for halo exchange:
      - `exchangeHalosBlocking()`
      - `exchangeHalosNonBlocking()`
      - `exchangeHalosNonBlockingAsyncStart()` + `exchangeHalosNonBlockingAsyncFinish()`.
   3. In `applyBCs()`, if we have a physical boundary, do local boundary logic. If we have an internal boundary, pack/unpack MPI buffers using one of the three patterns.

3. **Euler Solver**:
   - Keep calls to `Q.applyBCs()` or `Q.applyBCsMPI()`. Possibly rename or condition with a build flag, e.g., `#ifdef AGOGE_PARALLEL`.
   - For the “async” approach, reorganize the PDE time step or flux computation so you can do interior updates while waiting on halo data.

4. **Testing**:
   - Start with a single rank. That should reduce to the original code.  
   - Move to 2 ranks, test domain decomposition, check boundary.  
   - Scale to more ranks.
