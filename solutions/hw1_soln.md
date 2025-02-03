# **Solution Guide: Chapter 1 Exercises**

## **Exercise 1.1: Pipelining and Throughput**

> **Problem Restatement**  
> Compare the speeds (throughput) of a classical (non-pipelined) FPU vs. a pipelined FPU.
>
> 1. Show that the **result rate** $r(n)$ for the pipelined FPU depends on $n$.  
> 2. Compute the **asymptotic rate** $r_\infty$.  
> 3. Derive the **speedup** over the non-pipelined FPU.  
> 4. Find $n_{1/2}$, the number of operations needed to reach half of the asymptotic rate.

### **Solution**

1. **Non-pipelined**:  
   - Suppose each operation requires $\ell$ cycles. Then, for $n$ operations:
     $$T_{\mathrm{classic}}(n) = n (\ell\,\tau),$$
   where $\tau$ is the clock period.

2. **Pipelined**:  
   - A pipeline with $\ell$ stages finishes the first operation after $\ell$ cycles, but then **one result** (ideally) appears each cycle. A common simplified formula:
     $$T_{\mathrm{pipelined}}(n) = (\ell + n - 1)\,\tau.$$
   - **Result rate** $r(n)$:
     $$r(n) = \frac{n}{T_{\mathrm{pipelined}}(n)} = \frac{n}{(\ell + n -1)\,\tau}.$$

3. **Asymptotic rate** $r_{\infty}$:  
   $$r_{\infty} = \lim_{n\to\infty} r(n) 
       = \lim_{n\to\infty} \frac{n}{(\ell + n -1)\,\tau}
       = \frac{1}{\tau}.$$
   Thus for large $n$, the pipeline produces **1 result per cycle**.

4. **Non-pipelined vs. Pipelined**:  
   - The classical FPU’s rate is $\frac{1}{\ell\,\tau}$.  
   - The pipelined’s asymptotic rate is $\frac{1}{\tau}$.  
   - Speedup:
     $$\frac{r_\infty}{\frac{1}{\ell\,\tau}} = \ell.$$
   So the pipeline is up to $\ell$ times faster in the limit.

5. **Half the peak throughput** ($n_{1/2}$):  
   - We want $r(n_{1/2}) = \tfrac12 r_\infty = \tfrac12 \cdot \frac{1}{\tau}.$
   - Solve
     $$\frac{n_{1/2}}{(\ell + n_{1/2} -1)\,\tau} 
         = \frac{1}{2\,\tau}.$$
   - This yields $n_{1/2} = \ell -1$.  
   - So after $\ell -1$ operations, we are already at half of the pipeline’s max throughput.

---

## **Exercise 1.2: Linked Triads**

> **Problem Restatement**  
> We have an operation of the form  
> $$\texttt{a[i]} = \texttt{b[i]} + \texttt{c[i]} \times \texttt{d[i]},$$  
> which is effectively **one multiplication plus one addition** per iteration. On a non-pipelined FPU, one might naïvely pay a cost of $2\ell$ cycles per iteration. But modern FPUs *can fuse multiply-add* into one pipeline pass (or mostly overlap them) once fully loaded.

### **Solution (with Fused Multiply-Add Assumption)**

1. **Non-pipelined**:  
   - $2$ floating-point ops (FMA expanded) $\rightarrow$ each op $\ell$ cycles = $2\ell$ cycles per iteration.  
   - $n$ elements $\rightarrow T_{\mathrm{classic}}(n) = 2\,n\,(\ell\,\tau).$

2. **Pipelined (Fused)**:  
   - We assume the multiply-add can stream through a pipeline of depth $\ell$.  
   - Thus, after $\ell$ cycles of startup, *each iteration* can retire one FMA result per cycle.  
   - $$T_{\mathrm{pipelined}}(n) = (\ell + n - 1)\,\tau.$$

3. **Asymptotic Speedup**:  
   - Non-pipelined rate = $\frac{1}{2\,\ell\,\tau}$.  
   - Pipelined $\lim_{n\to\infty}$ rate = $\frac{1}{\tau}$.  
   - **Speedup**: 
     $$\frac{\frac{1}{\tau}}{\frac{1}{2\,\ell\,\tau}} = 2\,\ell.$$

4. **$n_{1/2}$ Computation**:  
   - Same style: solve 
     $$\frac{n_{1/2}}{(\ell + n_{1/2}-1)\,\tau} = \frac12 \cdot \frac{1}{\tau}.$$
   - We find again $n_{1/2} = \ell -1.$

Hence, with perfect FMA-pipelining, we can achieve up to $2\ell$ speedup for large $n$.

---

## **Exercise 1.3: Multiple Pipelines in Parallel**

> **Problem Restatement**  
> Suppose there are $p$ identical pipelines (say, a vector unit of width $p$), all doing the same operation in parallel. Each pipeline has depth $\ell$. Analyze the time, speedup, and half-peak point.

### **Solution**

1. **Non-pipelined**:  
   $$T_{\text{classic}}(n) = n (\ell\,\tau).$$

2. **$p$-way Pipelined**:  
   - Once started, we can produce $p$ results every cycle. Startup cost $\ell$ cycles.  
   $$T_{\mathrm{pipelined},p}(n) \,\approx\, (\ell + \frac{n}{p} - 1)\,\tau,$$
   assuming $\frac{n}{p}$ is integer for simplicity.

3. **Asymptotic Speedup**:  
   - Non-pipelined rate: $\frac{1}{\ell\,\tau}$.  
   - $p$-way pipeline rate: $\frac{p}{\tau}.$
   - Speedup $\approx \ell\,p.$

4. **$n_{1/2}$**:  
   - The pipeline array’s max throughput is $\frac{p}{\tau}$.  
   - We want half of that, i.e. $\frac{p}{2\tau}$.  
   - Setting 
     $$\frac{n_{1/2}}{(\ell + \tfrac{n_{1/2}}{p} -1)\,\tau} 
         = \frac{p}{2\,\tau}$$
     yields roughly $n_{1/2} = p(\ell -1)$ for integer-friendly approximations.  
   - So around $p(\ell-1)$ elements, we get half the pipeline array’s peak.

---

## **Exercise 1.4: Recursive Doubling of a Recurrence**

> **Problem Restatement**  
> We have a dependency
> ```cpp
> x[i+1] = a[i] * x[i] + b[i];
> ```
> which is fully sequential if done naively. We want to do **recursive doubling** to skip over intermediate elements and potentially achieve more parallelism.

### **Solution**

1. **Skip-1 Recurrence**:  
   From the original:
   $$x[i+1] = a[i]\,x[i] + b[i],$$
   we want $x[i+2]$ directly in terms of $x[i]$:
   $$x[i+2] 
       = a[i+1] \, x[i+1] + b[i+1]
       = a[i+1] \Bigl(a[i]\,x[i] + b[i]\Bigr) + b[i+1].$$
   So
   $$x[i+2] 
       = \bigl(a[i+1]\,a[i]\bigr)\,x[i]
         + \bigl(a[i+1]\,b[i] + b[i+1]\bigr).$$

2. **General Skip-$2^k$**:  
   Repeating that logic, we can derive $x[i + 4]$ in terms of $x[i]$, etc. In general,
   $$x[i + 2^k] 
       = \alpha_k(i)\,x[i] + \beta_k(i),$$
   for some precomputed coefficients $\alpha_k$, $\beta_k$.

3. **Filling in Missing (e.g. Odd) Indices**  
   - **Steps**:
     1. **Precompute** enough of these skip formulas, e.g. $\alpha_1(i)$, $\beta_1(i)$ so that `x[i+2]` can be found from `x[i]`.  
     2. **Compute even indices** in a separate pass (i.e. `i = 0,2,4,...`) in parallel, using the skip-2 formula.  
     3. **Fill in** `x[i+1], x[i+3], ...`: for that, you might use the original direct formula. Because now `x[i]` is final, no further dependence.

   Concretely, one might do:

   ```cpp
   // 1) Coeffs for skip-2
   for (i=0; i<n; i+=2) {
     alpha_2[i] = a[i+1]*a[i];
     beta_2[i]  = a[i+1]*b[i] + b[i+1];
   }
   // 2) Compute x[2], x[4], ...
   for (i=0; i<n; i+=2) {
     x[i+2] = alpha_2[i]*x[i] + beta_2[i];
   }
   // 3) Fill in x[1], x[3], ...
   for (i=0; i<n; i+=2) {
     x[i+1] = a[i]*x[i] + b[i];
   }
   ```

4. **Time Analysis**  
   - **Naive**: $T_0(n) \approx n \times (\text{latency of each iteration}).$ No concurrency.  
   - **Recursive Doubling**: We have some overhead of building skip formulas. But then we can do a “parallel chunk” of updates. Potentially an $O(\log n)$ step method if we keep doubling the skip.  

5. **Why Overhead May be Small**  
   - Once we break the tight data dependency, we can exploit pipeline concurrency or multi-core concurrency.  
   - The one-time overhead of computing skip-coefficients or filling in missing elements is small compared to the gain of parallelizing the main body for large $n$.

---

## **Exercise 1.5: Why L1 < L2 < L3**

> **Problem Restatement**  
> Explain both **practical** and **theoretical** reasons why L1 cache is smaller than L2, which is smaller than L3.

### **Solution**

1. **Practical (Cost/Power)**  
   - Faster caches (esp. L1) must use very power-hungry, expensive SRAM cells physically near the CPU core. Making them bigger would increase chip area and cost significantly.

2. **Theoretical (Latency/Access Time)**  
   - Accessing a larger memory structure is intrinsically slower: more lines, bigger tag directories, more wire delay, etc.  
   - L1 is kept minimal to get **1–2 cycle** access. L2 and L3 are bigger, but slower.

Hence each cache level is larger but slower.

---

## **Exercise 1.11: Recursive Doubling and Bank Conflicts**

> **Problem Restatement**  
> We have the summation:
> ```cpp
> for (s=2; s<2*n; s*=2)
>   for (i=0; i<n-s/2; i+=s)
>     x[i] += x[i+s/2];
> ```
> We want to see whether this leads to bank conflicts if memory is laid out in interleaved banks of size $2^k$. Also compare with *recursive halving*:
> ```cpp
> for (s=(n+1)/2; s>1; s/=2)
>   for (i=0; i<n; i++)
>     x[i] += x[i+s];
> ```
> Which might be better in a parallel scenario?

### **Solution**

- Suppose memory is interleaved into banks in chunks of $2^k$. That means:
  - An address $A$ belongs to bank $(A \bmod 2^k)$.  
  - In a parallel loop, multiple threads writing to `x[i]` for different `i` can conflict if those addresses share the same remainder mod $2^k$.

1. **Recursive Doubling Pattern**:
   - On iteration $\ell$, `s = 2^\ell`, so the offset is `s/2 = 2^{\ell-1}`.  
   - Then we do `x[i] += x[i + 2^{\ell-1}]` for `i` jumping by `s = 2^\ell`.  
   - If $2^{\ell-1}$ is a multiple of $2^k$, i.e. $\ell-1 \ge k$, then `x[i]` and `x[i + s/2]` may map to the **same** remainder mod $2^k$ for certain `i`.  
   - If multiple threads try to do that at once for many `i`, bank conflicts can occur systematically.

2. **Recursive Halving Pattern**:
   - `s` starts around `n/2` and is divided by 2 each iteration.  
   - The offset `s` might not be a neat power of 2 at each step (especially if `(n+1)/2` is not a perfect power of 2).  
   - This can reduce or at least randomize the alignment mod $2^k$, mitigating large-scale collisions.

3. **Which is Better?**  
   - The doubling approach can suffer more systematic collisions for large powers-of-two data sizes.  
   - The halving approach has a different offset pattern and can lead to fewer conflicts.  
   - In a **parallel** environment, bank conflicts matter more, so halving can be better.  

---

## **Exercise 1.14: Matrix-Matrix Product and Data Reuse**

> **Problem Restatement**  
> The matrix multiplication $C \leftarrow A \times B$ has high mathematical reuse, but naive code may not exploit it. Why?

### **Solution**

1. **Mathematical Reuse**  
   - Each element of $A$ is used in multiple dot products with columns of $B$. Each element of $B$ is used with all rows of $A$. Potentially $O(n)$ uses per element.

2. **Naive Code**  
   - For example:
     ```cpp
     for (i=0; i<n; i++)
       for (j=0; j<n; j++)
         for (k=0; k<n; k++)
           C[i][j] += A[i][k]*B[k][j];
     ```
   - The access pattern might do poor caching for `B[k][j]`, or reload the same row of `A` repeatedly if `k` changes fastest.

3. **What Determines Effective Reuse**  
   - **Loop ordering**: e.g., `(i, k, j)` might hold `A[i][k]` in register across the innermost loop in `j`.  
   - **Blocking/Tiling**: partition `A,B,C` into sub-blocks that fit in cache to maximize re-use in each block.

Hence naive triple loops may not realize the potential reuse without careful ordering or tiling.

---

## **Exercise 1.16: Summation with Tree-Based Reduction**

> **Problem Restatement**  
> Compare a **tree-based** summation code:
> ```cpp
> for (s=2; s<2*n; s*=2)
>   for (i=0; i<n-s/2; i+=s)
>     x[i] += x[i + s/2];
> sum = x[0];
> ```
> vs. a **standard** summation:
> ```cpp
> sum = 0;
> for (i=0; i<n; i++)
>   sum += x[i];
> ```
> in terms of **spatial and temporal locality**.

### **Solution**

1. **Standard Summation**  
   - **Spatial locality**: excellent, we read `x` in sequential order (stride=1). Each cache line is well-used.  
   - **Temporal locality**: minimal for `x`. We only read each element once. The accumulator `sum` can stay in a register, so no repeated memory writes for it.

2. **Tree-Based Summation**  
   - The stride changes each pass: `i += s`, with `s` doubling each iteration. This means we might skip large blocks in memory. So **spatial locality** can be poor. Many partial lines are loaded.  
   - Possibly some re-use of `x[i]` from pass to pass, but if $n$ is big, we likely evict it from cache before it’s reused. So **temporal locality** is also questionable.  
   - If done in parallel, it is a classic reduce pattern but can degrade locality unless carefully blocked.

Conclusion: the tree-based code is more about parallelism; it is typically **worse** for simple single-thread cache locality than the straightforward loop.

Below is a detailed solution guide in **Markdown** format for the listed exercises from Chapter 6 of _IntroToHPC_ by Victor Eijkhout. If there is any point requiring additional clarification for completeness or correctness, I will pose questions inline. Please let me know if you have any specific additional context or code snippets that differ from what was included above.

---

# **Solutions for Chapter 6 Exercises**

## Exercise 6.2

> **Exercise statement**  
> While the strategy just sketched will demonstrate the existence of cache sizes, it will not report the maximal bandwidth that the cache supports. What is the problem and how would you fix it?

### **Solution Explanation**

1. **Context of the Problem**  
   - The “strategy” referenced is typically creating a random-access or pseudo-random traversal of an array (e.g., using a permutation of indices) in order to identify cache size boundaries.  
   - Such random or indirect accesses defeat hardware prefetching and prevent contiguous reuse of cache lines.

2. **Why This Fails to Show *Maximum* Cache Bandwidth**  
   - **Maximum bandwidth** is usually observed with **contiguous** (or at least stride-1) memory access, where the memory interface and cache prefetchers can stream data at the highest rate.  
   - In the random-access test, many cache lines are only partially utilized and rarely get re-accessed. The result is that a large fraction of the time is spent on address translation, TLB misses, or suboptimal prefetch.

3. **How to Fix It**  
   - **Run a separate streaming test** that **fully utilizes** each cache line and benefits from prefetching to measure the **peak** or near-peak bandwidth.  
   - Concretely:  
     1. **Use a contiguous loop** over the array:  
        ```cpp
        for (int i = 0; i < N; i++) {
          sum += array[i];
        }
        ```  
     2. **Ensure the array is large enough** so that it exceeds L1 and L2 (so we measure main-memory bandwidth or L3 bandwidth).  
     3. **Consider measuring smaller subarrays** that fit in L1 or L2 to measure the maximum bandwidth specifically for those caches if desired.  
   - Thus, we may need **two different benchmarks**:  
     - One (random or strided) to pinpoint cache size thresholds.  
     - Another (streaming/contiguous) to measure peak cache or main-memory bandwidth.

> **Answer summary**:  
> The random-/strided-access method identifies cache boundaries but shows artificially low bandwidth. To measure true peak bandwidth, use a streaming (contiguous) kernel that efficiently exploits cache lines and hardware prefetching.

---

## Exercise 6.3

> **Exercise statement**  
> Give an example of a doubly-nested loop where the loops can be exchanged; give an example where this can not be done. If at all possible, use practical examples from this book.

### **Solution Explanation**

We consider two situations:

1. **Loops That *Can* Be Exchanged**  
   - Typically, if a nested loop has **no data dependencies** between iterations in the outer loop and the inner loop, the order can be swapped.  
   - **Example**: Summation of a 2D array’s elements:
     ```cpp
     double sum = 0.0;
     for (int i = 0; i < N; i++) {
       for (int j = 0; j < M; j++) {
         sum += A[i][j];
       }
     }
     ```
     - Because each `(i, j)` access is independent, we can exchange the loops:
     ```cpp
     double sum = 0.0;
     for (int j = 0; j < M; j++) {
       for (int i = 0; i < N; i++) {
         sum += A[i][j];
       }
     }
     ```

2. **Loops That *Cannot* Be Exchanged**  
   - If there is a **data dependency** where each outer iteration depends on results of the inner iteration (or vice versa), then the loop order is fixed for correctness.  
   - **Example**: A triangular iteration pattern where updating `A[i][j]` relies on a previously updated entry in the same row or column:
     ```cpp
     for (int i = 1; i < N; i++) {
       for (int j = 1; j < i; j++) {
         A[i][j] = A[i-1][j] + A[i][j-1];
       }
     }
     ```
     - Here, `A[i][j]` depends on `A[i-1][j]` (the previous row) and `A[i][j-1]` (the previous column).  
     - Exchanging the loops (swapping `i` and `j`) breaks the dependency chain and gives incorrect results.  
     - Thus, the loop order is **not** interchangeable.

> **Answer summary**:  
> - **Exchangeable**: Independent loops over a 2D array sum or similar.  
> - **Non-exchangeable**: Loops containing data dependencies across i–j iteration space (common in dynamic programming or triangular matrix updates).

---

## Exercise 6.4

> **Exercise statement**  
> If your loop describes the \((i,j)\) indices of a two-dimensional array, it is often best to let the \(i\)-index be in the inner loop for Fortran, and the \(j\)-index inner for C. Can you come up with at least two reasons why this is possibly better for performance?

### **Solution Explanation**

Recall:
- **Fortran** uses **column-major** layout: consecutive memory locations vary the **first index** (`i`).
- **C** uses **row-major** layout: consecutive memory locations vary the **last index** (`j`).

1. **Reason 1: Memory Contiguity and Prefetch**  
   - In Fortran, `A(i,j)` has `i` as the fastest-varying index. Looping over `i` in the innermost loop results in contiguous access, improving cache hits and prefetch efficiency.  
   - In C, `A[i][j]` has `j` as the fastest-varying index. Looping over `j` in the innermost loop yields contiguous memory access.

2. **Reason 2: Vectorization and Pipeline Efficiency**  
   - Modern compilers can **vectorize** inner loops that access contiguous memory. This includes auto-vectorization of short inner loops that step through consecutive addresses.  
   - If you place the “correct” index in the inner loop, the stride becomes 1, enabling more effective use of SIMD instructions.

> **Answer summary**:  
> - Cache line usage is optimized by aligning your innermost index with the storage order.  
> - Compilers more readily generate vectorized code with stride-1 loops.

---

## Exercise 6.5

> **Exercise statement**  
> For this reason, loop tiling is also known as *cache blocking*. The block size depends on how much data is accessed in the loop body; ideally you would try to make data reused in L1 cache, but it is also possible to block for L2 reuse (though that is typically less optimal).  
>
> **Analyze this example**. When is `x` brought into cache, when is it reused, and when is it flushed? What is the required cache size in this example? Rewrite this example using:
> ```cpp
> #define L1SIZE 65536
> ```

### **Solution Explanation**

Consider an example loop where we repeatedly operate on an array (or sub-block of an array) `x` in blocks of size `bs`. Pseudocode:

```cpp
double x[...]; // large array
int N = ...;   // total size

// Suppose we have some blocking approach:
for (int blockStart = 0; blockStart < N; blockStart += bs) {
  // This block is small enough to fit in L1
  // so we can reuse it effectively
  for (int i = blockStart; i < blockStart + bs; i++) {
    x[i] = 2.3 * x[i] + 1.2; 
  }
  // Possibly do repeated passes over the same block
  // ...
}
```

1. **When is `x` Brought into Cache?**  
   - On the **first access** within each block. The CPU loads the corresponding cache lines from main memory or a higher-level cache (L2/L3) into the L1 cache.

2. **When is it Reused?**  
   - **During the iteration(s) that still touch the same block**. As long as the block fits into L1 (and we have not moved on to a new block that evicts this data), further accesses stay in L1 without requiring reload from main memory.

3. **When is it Flushed (Evicted)?**  
   - Once we finish work on that block and move to the next block, new data loads eventually evict old data from L1 if the L1 size is insufficient to contain both blocks simultaneously.  
   - If the block size `bs` plus overhead for other arrays/data is near or larger than L1, the older block’s data is quickly flushed.

4. **Required Cache Size in This Example**  
   - The **ideal** is that each block `bs` (in bytes) should be **≤** `L1SIZE` minus some overhead for other variables, stack usage, etc.  
   - If `bs * sizeof(double)` is the primary occupant, then `bs * 8 ≤ L1SIZE`, if each element is a `double` (8 bytes).

5. **Rewrite with `#define L1SIZE 65536`**  
   - Suppose we know L1 size is 64 KiB. We can pick `bs` so that `bs * 8` fits well below 65536:
     ```cpp
     #define L1SIZE 65536
     #define DSIZE   8        // size of double in bytes
     // For safety, use half the L1 size or some fraction
     // to account for overhead, etc.
     #define FRACTION 0.5

     int bs = (int)( (FRACTION * L1SIZE) / DSIZE );

     // Then tile:
     for (int blockStart = 0; blockStart < N; blockStart += bs) {
       int blockEnd = min(blockStart + bs, N);
       for (int i = blockStart; i < blockEnd; i++) {
         x[i] = 2.3 * x[i] + 1.2; 
       }
     }
     ```
   - This ensures each block is small enough to remain in the L1 cache while being processed.

> **Answer summary**:  
> - `x` is loaded into L1 on the first access of each block, reused during that block’s iteration, and flushed upon moving to the next block that displaces it from L1.  
> - The block size should be chosen so that it (plus any overhead) fits in the L1 cache (64 KiB in the example).  
> - The rewrite with `L1SIZE` simply uses a compile-time definition to set a block size that respects the L1 capacity.

