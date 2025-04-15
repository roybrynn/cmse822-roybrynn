## I. Core Architectural Concepts (Ch. 1)

• **Von Neumann Architecture**  
  – *Concept*: Beyond noting that CPU, memory, and I/O are discrete components, this topic emphasizes the stored-program concept and the fetch/decode/execute cycle in a pipelined fashion [IntroHPC §1.1].  
  – *Video*: [Von Neumann Architecture – Neso Academy](https://www.youtube.com/watch?v=7mIY6n0TE8w)  
  *(This 11‐minute lecture covers the fundamentals of the stored-program concept and contrasts older architectures with modern RISC principles.)*

• **Pipelining Basics**  
  – *Concept*: Modern CPUs split operations into multiple stages (decode, execute, write-back) and use assembly-line processing. Discuss hazards (data, structural, control) and mitigation techniques like loop unrolling and branch prediction [IntroHPC §1.2.1].  
  – *Video*: [Pipelining in Computer Architecture – Neso Academy](https://www.youtube.com/watch?v=8FzvhmkB6mE)  
  *(A 10‐minute introduction to how pipelining works and what causes stalls in modern CPUs.)*

• **Cache/Memory Hierarchy**  
  – *Concept*: Explores the multiple levels of cache (L1, L2, L3) and main memory (DRAM), explaining how data locality and contiguous accesses affect performance. Also covers “capacity” and “conflict” misses [IntroHPC §1.3].  
  – *Video*: [Memory Hierarchy and Cache Basics – Computerphile](https://www.youtube.com/watch?v=O2L2Uv9pdDA)  
  *(In about 12 minutes, this video explains cache levels and their impact on performance.)*

• **Key Principle: Minimizing Data Movement**  
  – *Concept*: Emphasize that many HPC codes are memory-bound. The “memory wall” (the gap between CPU speed and DRAM bandwidth) means that reducing data movement is often more impactful than chasing higher CPU frequencies [IntroHPC §1.3.2].  
  – *Video*: [The Memory Wall – Data Movement Bottlenecks](https://www.youtube.com/watch?v=Qz9q3Fq2p0Q)  
  *(This concise 10‐minute talk discusses why and how data movement becomes the primary performance limiter.)*

• **Register Usage**  
  – *Concept*: Registers are the fastest storage available, yet they are very limited in number. Topics include loop-invariant code motion, register allocation, and the impact of register spills on performance [IntroHPC §1.3.4].  
  – *Video*: [CPU Registers Explained – Computerphile](https://www.youtube.com/watch?v=nAKjzTGG07g)  
  *(A brief 9‐minute explanation of what registers are and how compilers use them.)  
  – *Example*:  
    ```cpp
    #include <array>
    #include <iostream>
    
    constexpr size_t N = 1024;
    int main() {
        std::array<double, N> data{};
        for (size_t i = 0; i < N; ++i)
            data[i] = static_cast<double>(i);
        double sum = 0.0;
        for (size_t i = 0; i < N; ++i)
            sum += data[i];
        std::cout << "Sum = " << sum << '\n';
        return 0;
    }
    ```
    *(This short C++20 snippet illustrates local array usage that can benefit from keeping data in registers.)*

---

## II. Performance Metrics and Optimization Techniques (Ch. 6)

• **Peak Performance vs. Sustained Performance**  
  – *Concept*: Theoretical “peak FLOPs” derive from CPU specs, yet real applications typically achieve only a fraction of that due to stalls from memory latency or branch mispredictions [IntroHPC §6.1].  
  – *Video*: [Understanding CPU Performance Metrics](https://www.youtube.com/watch?v=oRbgrt1k6c8)  
  *(A 12‐minute lecture that breaks down how peak and sustained performance differ in practice.)*

• **Bandwidth Considerations**  
  – *Concept*: Memory bandwidth can be the throughput cap in HPC. Techniques like data blocking (tiling) and vector-friendly loops help to exploit caches and reduce memory transfers [IntroHPC §§6.2 & 6.7].  
  – *Video*: [Memory Bandwidth and Performance – A Short Guide](https://www.youtube.com/watch?v=2d5UQH6qK8s)  
  *(This 11‐minute video outlines how memory bandwidth influences real-world performance.)*

• **Latency Hiding**  
  – *Concept*: Out-of-order execution, hardware prefetching, and loop reordering/unrolling can hide memory latency and keep more instructions “in flight” [IntroHPC §§1.3.6 & 6.3.2].  
  – *Video*: [Hiding Latency in Modern CPUs](https://www.youtube.com/watch?v=3lC44jVnGzw)  
  *(In about 10 minutes, learn techniques used by modern processors to mask memory latency.)*

• **Compiler Optimizations**  
  – *Concept*: Flags like `-O3`, link-time optimizations, and auto-vectorization can dramatically boost performance—while cautioning that aggressive optimizations may affect floating-point consistency [IntroHPC §6.8].  
  – *Video*: [GCC Compiler Optimizations Explained](https://www.youtube.com/watch?v=28i6zRJpSgk)  
  *(A focused 13‐minute video that reviews key optimization flags and their impact on code performance.)*

• **Performance Measurement**  
  – *Concept*: Learn to measure execution time using C++’s `std::chrono` (or hardware counters) to compare naive versus optimized loops and diagnose bottlenecks [SciProg §24.8.2].  
  – *Video*: [C++ Chrono Library Tutorial](https://www.youtube.com/watch?v=Qf6w83VbJio)  
  *(A concise 10‐minute tutorial on using `std::chrono` for performance measurement.)*

• **Microbenchmarks**  
  – *Concept*: Small, focused tests (like summing an array or matrix-vector multiplication) help isolate whether performance limits come from CPU pipelines or cache bandwidth [IntroHPC §6.10].  
  – *Video*: [Writing Microbenchmarks in C++](https://www.youtube.com/watch?v=nG8vm8l7FoE)  
  *(This 12‐minute session walks through creating and interpreting microbenchmarks.)*

---

## III. Complexity Analysis (Ch. 15)

• **Big-O Notation**  
  – *Concept*: Understanding asymptotic scaling is essential—even when constants and data placement (cache effects) play a critical role in HPC [IntroHPC §15.1].  
  – *Video*: [Big-O Notation Explained – CS Dojo](https://www.youtube.com/watch?v=V6mKVRU1evU)  
  *(An 11‐minute explanation of algorithmic complexity.)*

• **Algorithmic Examples**  
  – *Concept*: Compare a blocked O(n²) matrix multiplication (optimized for cache reuse) with a naive O(n³) approach to show how tile size influences performance [IntroHPC §§7.4 & 6.9].  
  – *Video*: [Matrix Multiplication and Complexity](https://www.youtube.com/watch?v=J3dnhnEJ53k)  
  *(A 13‐minute breakdown of how algorithmic choices affect cache performance.)*

• **Amdahl’s Law (Preview)**  
  – *Concept*: Introduces how even a small serial fraction (f) limits overall speedup, foreshadowing the importance of minimizing serial bottlenecks in HPC [IntroHPC §2.2.3].  
  – *Video*: [Amdahl's Law Explained in 5 Minutes](https://www.youtube.com/watch?v=F2N9kBfj0vY)  
  *(A short, clear explanation of Amdahl’s Law and its implications.)*

---

## IV. Memory Hierarchies

• **Cache Levels**  
  – *Concept*: Understand the typical sizes and latencies (e.g., L1 ~32KB, L2 ~256KB, L3 ~10MB, with DRAM being 100× slower) and why sequential accesses reduce misses [IntroHPC §§1.3.5].  
  – *Video*: [Cache Memory & Its Levels – Explained](https://www.youtube.com/watch?v=LN6B19TRsFQ)  
  *(An 11‐minute video outlining cache levels and their performance impact.)*

• **Prefetching and Blocking**  
  – *Concept*: Blocking (tiling) transforms loops to exploit the limited number of cache lines and to hide memory latency via effective prefetching [IntroHPC §1.3.6].  
  – *Video*: [Cache Blocking and Loop Tiling](https://www.youtube.com/watch?v=SJdqF3gdX5s)  
  *(A 12‐minute explanation of these critical techniques.)*

• **NUMA (Non-Uniform Memory Access)**  
  – *Concept*: On multi-socket nodes, local versus remote memory access matters; ensuring local data placement can significantly boost performance [IntroHPC §1.5].  
  – *Video*: [Introduction to NUMA Architecture](https://www.youtube.com/watch?v=euSSKMtNI7Y)  
  *(A 13‐minute video that clarifies NUMA principles.)*

• **Real-World Demonstration: Cache Thrashing**  
  – *Concept*: Demonstrate how poor locality leads to cache thrashing and performance slowdowns [IntroHPC §6.5].  
  – *Video*: [Cache Thrashing Explained](https://www.youtube.com/watch?v=JOmfbxNqnc8)  
  *(An illustrative 10‐minute demonstration of cache thrashing.)*

---

## V. Hands-On Exercises and Their Purpose

1. **Exercises 1.1–1.5, 1.11, 1.14, 1.16, 1.17 [SciProg Ch. 1]**  
   – *Focus*: These exercises reinforce CPU design, register allocation, pipelining stalls, and loop transformations via small C++ programs.  
   – *Video*: [Efficient C++ Coding for HPC – Loop Optimization & Memory Management](https://www.youtube.com/watch?v=FwEGUAtd9jY)  
   *(A 14‐minute tutorial demonstrating how small optimizations can yield significant performance gains.)*

2. **Exercises 6.2–6.5 [SciProg Ch. 6]**  
   – *Focus*: These exercises highlight the interplay between memory bandwidth and peak FLOPs, emphasizing aligned allocations and cache-friendly data structures.  
   – *Video*: [Microbenchmarking & Performance Tuning in C++](https://www.youtube.com/watch?v=K1G4VTXJ15k)  
   *(A 13‐minute session on setting up and interpreting microbenchmarks.)*

---

## VI. Intel VTune Tutorial

• **Profiling Workflow & Bottleneck Identification**  
  – *Concept*: Compile code with debugging symbols (or moderate optimization), run VTune (or similar tools) to capture hotspots, and interpret memory traffic and CPI graphs to spot the “memory wall” [IntroHPC §6.4 & Tools Appendix].  
  – *Video*: [Getting Started with Intel VTune Profiler](https://www.youtube.com/watch?v=nPeXezfHnpw)  
  *(A 15‐minute introductory tutorial on using VTune for performance analysis.)*

• **Roofline Synergy**  
  – *Concept*: Combine hardware counters from VTune with the roofline model to check whether loops are compute- or memory-bound, guiding further optimizations [IntroHPC §6.8].  
  – *Video*: [The Roofline Model Explained](https://www.youtube.com/watch?v=oh2p5w5-1fI)  
  *(An accessible 12‐minute video that explains how to leverage the roofline model for performance tuning.)*

---

## VII. Context for the Group Project

• **Arithmetic Intensity**  
  – *Concept*: Demonstrate how to compute “FLOPs per byte” for kernels (e.g., `Y[j] += A[j][i]*B[i]`) to predict when a kernel will hit the memory bandwidth ceiling [IntroHPC §6.9].  
  – *Video*: [Understanding Arithmetic Intensity & the Roofline Model](https://www.youtube.com/watch?v=yBS15PKbbfM)  
  *(A concise 11‐minute explanation integrating arithmetic intensity with roofline analysis.)*

• **Stencil vs. Spectral Dwarfs**  
  – *Concept*: Compare stencil operations (with large data neighborhoods and poor cache reuse) against spectral (FFT-based) methods that have high numerical intensity but potential communication overheads [ParallelProg §9.1].  
  – *Video*: [Stencil Computations in HPC](https://www.youtube.com/watch?v=qCniC3zvGAY)  
  *(A focused 12‐minute lecture covering the pros and cons of each kernel type.)*

• **Roofline Model Usage & Communication Steps**  
  – *Concept*: Overlaying “peak CPU” and “peak memory” lines on performance graphs reveals inefficiencies and points to the next steps (e.g., further blocking or vectorization). Foreshadow the importance of efficient node-level kernels before tackling multi-node or GPU parallelism.  
  – *Videos*:  
    – [The Roofline Model Explained](https://www.youtube.com/watch?v=oh2p5w5-1fI) *(reused for clarity)*  
    – [MPI Communication Basics](https://www.youtube.com/watch?v=_RzN7X3W73Y)  
    *(A 10‐minute primer on basic communication principles in HPC.)*

---

## VIII. Overarching Skills and Transferable Knowledge

• **Reading HPC Literature**  
  – *Concept*: Develop the skill to efficiently read and digest scientific papers and technical documentation, vital for keeping up with HPC research and best practices [IntroHPC ch.6].  
  – *Video*: [How to Read a Scientific Paper Efficiently](https://www.youtube.com/watch?v=4k0_m8R_5iY)  
  *(A practical 10‐minute guide to navigating dense technical literature.)*

• **Software Instrumentation**  
  – *Concept*: Become comfortable using performance tools (VTune, perf, etc.) to measure instructions per cycle, cache miss rates, and other hardware counters.  
  – *Video*: [Introduction to Profiling Tools in HPC](https://www.youtube.com/watch?v=n5KypvZhQfs)  
  *(An 11‐minute overview of using software instrumentation effectively.)*

• **Data-Driven Tuning**  
  – *Concept*: Learn to iteratively refine loop structures, data layouts, and compiler flags based on measured performance metrics.  
  – *Video*: [Data-Driven Performance Optimization in C++](https://www.youtube.com/watch?v=X1R3wS9oMn8)  
  *(A 12‐minute session on making optimization decisions based on real data.)*

• **Collaboration & Communication**  
  – *Concept*: Group roles (Navigator, Lorekeeper, Bard) simulate real HPC team workflows and emphasize the importance of clear communication and effective version control practices [ParallelProg §2.1.2].  
  – *Video*: [Effective Team Collaboration in Software Projects](https://www.youtube.com/watch?v=f0A8ML0NFck)  
  *(A brief 10‐minute talk on best practices for team collaboration in technical projects.)*

---

## IX. Expected Outcomes

• **Mental Map of Single-Processor Architecture**  
  – A deep understanding of where registers, caches, and memory bottlenecks occur, equipping students with the tools to diagnose performance issues at the node level.

• **Proficiency in Performance Counters and Optimization Techniques**  
  – Skill in employing cache-blocking, loop fusion/fission, and vectorization to merge algorithmic design with hardware awareness.

• **Foundations for Parallel and Distributed HPC**  
  – Recognizing that efficient single-node kernels are the building blocks for scaling to multi-node and GPU environments.

---

## X. Conclusion

This module lays the groundwork by revealing how modern CPUs operate “under the hood” and demonstrating that minimizing data movement is key to high performance. The integration of complexity theory (Big-O, Amdahl’s Law) with practical hardware constraints paves the way for future exploration into parallel computing. With these carefully curated video lectures complementing in‐depth exercises and projects, students are well prepared to tackle both the theoretical and practical challenges of high-performance computing.

---

This detailed outline—with embedded video lectures—should serve as a comprehensive guide for graduate students in CMSE 822, ensuring they develop a robust, practical, and theoretical understanding of single-processor computing as the foundation for further HPC studies.