**Exercise 1.1**

**Prompt:** Compare the speed of a classical FPU and a pipelined one. Show that the result rate is dependent on the number of operations (‘n’): give a formula for $r(n)$, and for $r_\infty = \lim_{n \to \infty} r(n)$. What is the asymptotic improvement in $r$ over the non-pipelined case? Show that for $n = n_{1/2}$, $r(n) = r_\infty/2$. This is often used as the definition of $n_{1/2}$.

**Solution:**
1. For a non-pipelined FPU, the time for $n$ operations is:
   
   $$ t(n) = n \cdot \ell \cdot \tau $$
   
   where $\ell$ is the number of stages in the pipeline and $\tau$ is the clock cycle time. The result rate for non-pipelined is:

   $$ r_{\text{non-pipelined}} = \frac{1}{\ell \tau}. $$

2. For a pipelined FPU, the time for $n$ operations is:

   $$ t(n) = (\ell + n - 1) \cdot \tau. $$

   The result rate for pipelined execution is:

   $$ r(n) = \frac{n}{t(n)} = \frac{n}{(\ell + n - 1) \tau}. $$

3. Taking the limit as $n \to \infty$:

   $$ r_\infty = \lim_{n \to \infty} \frac{n}{(\ell + n - 1) \tau} = \frac{1}{\tau}. $$

   Thus, the asymptotic improvement is:

   $$ \frac{r_\infty}{r_{\text{non-pipelined}}} = \frac{1}{\tau} \cdot \ell \tau = \ell. $$

4. To find $n_{1/2}$, set $r(n_{1/2}) = \frac{r_\infty}{2}$:

   $$ \frac{n_{1/2}}{(\ell + n_{1/2} - 1) \tau} = \frac{1}{2\tau}. $$

   Solving for $n_{1/2}$:

   $$ n_{1/2} = \ell + 1. $$

---

**Exercise 1.2**

**Prompt:** Analyze the speedup and $n_{1/2}$ of linked triads.

**Solution:**
A linked triad operation involves performing multiple operations (e.g., addition and multiplication) where the output of one feeds directly into the next without going back to memory. The effective latency is reduced as stages of the pipeline overlap. Assuming each triad uses separate pipelines, the analysis follows a similar approach to Exercise 1.1:

1. Time for $n$ operations:
   
   $$ t(n) = (\ell + n - 1) \cdot \tau. $$

2. Speedup relative to a non-pipelined architecture:
   
   $$ \text{Speedup} = \frac{t_{\text{non-pipelined}}}{t_{\text{pipelined}}} = \frac{n \cdot \ell \cdot \tau}{(\ell + n - 1) \cdot \tau}. $$

3. Asymptotic behavior and $n_{1/2}$: same derivation as Exercise 1.1 applies, with improvements based on pipeline efficiency.

---

**Exercise 1.3**

**Prompt:** Analyze the speedup and $n_{1/2}$ of a processor with multiple pipelines operating in parallel. Suppose there are $p$ independent pipelines that execute the same instruction on a stream of operands.

**Solution:**
1. For $p$ pipelines, each can handle independent streams of data. The total throughput scales linearly with $p$, assuming no conflicts or bottlenecks.

2. Result rate for parallel pipelines:
   
   $$ r(n) = \frac{p \cdot n}{(\ell + n - 1) \cdot \tau}. $$

3. Asymptotically:

   $$ r_\infty = \frac{p}{\tau}. $$

4. At $n = n_{1/2}$:

   $$ r(n_{1/2}) = \frac{r_\infty}{2} \Rightarrow n_{1/2} = \frac{\ell + 1}{p}. $$

   Thus, increasing $p$ reduces $n_{1/2}$, improving performance.

---

**Exercise 1.4**

**Prompt:** Transform a dependent loop into one suitable for parallel execution (recursive doubling).

**Solution:**
1. Original loop:
   
   $$ x[i+1] = a[i] \cdot x[i] + b[i]. $$

2. Recursive doubling:

   Compute $x[i+2]$ directly from $x[i]$:

   $$ x[i+2] = a[i+1] \cdot (a[i] \cdot x[i] + b[i]) + b[i+1]. $$

3. Algorithm:
   - Step 1: Compute even-indexed terms (e.g., $x[i], x[i+2], x[i+4], \ldots$)
   - Step 2: Compute odd-indexed terms using interpolation.

4. Efficiency analysis:
   - Work complexity: $O(\log n)$ steps for doubling.
   - Communication cost: reduced compared to original loop.

---

**Exercise 1.5**

**Prompt:** Explain why L1 cache is smaller than L2 and L3 caches.

**Solution:**
1. Practical reason: L1 cache is closer to the processor core, requiring faster and smaller memory.
2. Theoretical reason: The speed of access decreases as the size increases due to physical constraints on memory density and signal propagation.

---

**Exercise 1.11**

**Prompt:** Analyze bank conflicts for recursive doubling and recursive halving algorithms.

**Solution:**
1. Recursive doubling:
   - Memory access pattern: stride doubles with each step.
   - Bank conflicts occur when stride aligns with bank size. Use prime-numbered banks to reduce conflicts.

2. Recursive halving:
   - Stride decreases with each step, reducing conflicts as computation proceeds.

3. Conclusion: Recursive halving is typically better for avoiding bank conflicts in both serial and parallel cases.

---

**Exercise 1.14**

**Prompt:** Compare memory-bound and compute-bound operations.

**Solution:**
1. Memory-bound: Limited by data transfer rates (e.g., vector addition).
2. Compute-bound: Limited by processor speed (e.g., matrix-matrix multiplication).
3. Optimization: Use tiling or blocking to improve memory locality and reduce cache misses.

---

**Exercise 1.16**

**Prompt:** Explore performance improvements using loop unrolling.

**Solution:**
1. Loop unrolling reduces overhead from loop control.
2. Example: Unroll a loop with stride 2:

   Original: $$ a[i] = b[i] + c[i], \quad \text{for } i = 0, n. $$

   Unrolled: $$ a[i] = b[i] + c[i], \quad a[i+1] = b[i+1] + c[i+1]. $$

3. Performance improves due to reduced branch instructions and better pipelining.

---

**Exercise 6.2**

**Prompt:** Analyze cache reuse in matrix-vector multiplication.

**Solution:**
Matrix-vector multiplication $y = Ax$:

1. Reuse opportunity: Vector $x$ is reused across rows of $A$.
2. Strategy: Partition $A$ into blocks that fit in cache.
3. Performance depends on matrix size and cache hierarchy.

---

**Exercise 6.3 and 6.4**

**Prompt:** Explore bandwidth and latency impacts on memory performance.

**Solution:**
1. Bandwidth: Affects sustained data transfer rates.
2. Latency: Influences initial data access delay.
3. Optimization: Overlap computation with data transfer using prefetching.

---

**Exercise 6.5**

**Prompt:** Evaluate performance using loop tiling for matrix operations.

**Solution:**
1. Tiling breaks large operations into cache-sized blocks.
2. Example: Multiply submatrices $C_{ij} = A_{ik} \cdot B_{kj}$.
3. Improves cache reuse and reduces memory traffic.
