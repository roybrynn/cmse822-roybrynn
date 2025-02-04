# CMSE 822 Group Quest 1: The Path to Peak Performance

## Setting

In this Quest, you will join the legendary Profilers' Guild, tasked with guiding the ancient code of Agoge through a labyrinth of memory-bound dungeons and compute-bound trials to achieve ultimate speed.

## Group Member Roles

First thing upon accepting this assignment, specify who will be the _Naviagator_, _Lorekeeper_, and the group _Bard_ for this project. The Navigator will be responsible for setting group meetings and ensuring that all aspects of the project are addressed equitably amongst the group. The Lorekeeper will be responsible for ensuring that the project report is prepared according to the requirements and that all prompts are explicitly and clearly addressed. The Bard will be responsible for giving a short (10 minute max) presentation in class on the group's work upon completion of the project.

**Navigator:** 

**Lorekeeper:** 

**Bard:**  

## Learning Objectives

- Develop collaborative skills using GitHub.
- Gain familiarity with using HPCC.
- Understand concepts of arithmetic intensity and memory hierarchies. 
- Understand how to interpret a roofline model and use it to estimate performance.
- Understand two basic computational kernel types (dwarfs): stencil operations and spectral methods
- Gain experience using Intel VTune for performance analysis.

Review the section in [HPSC](https://cmse822.github.io/assets/EijkhoutIntroToHPC2020.pdf) on computing arithmetic intensity for given compute kernels. 
Then, as a group, compute the arithmetic intensities of the following kernels in units of FLOPs/byte, assuming 8 bytes per float.

```C
  Y[j] += Y[j] + A[j][i] * B[i]
```

```C
  s += A[i] * A[i]
```

```C
  s += A[i] * B[i]
```

```C
  Y[i] = A[i] + C*B[i]
```

Included a table in your project report listing the arithmetic intensities for these kernels.

## Part 1: The Roofline Model 

In this part, you will explore the roofline model for analyzing the interplay between arithmetic intensity and memory bandwidth for architectures with complex memory hierarchies. 

1. Reference the materials on the Roofline Performance model at <https://crd.lbl.gov/divisions/amcr/computer-science-amcr/par/research/roofline/>. In particular, look through ["Roofline: An Insightful Visual Performance Model for Floating-Point Programs and Multicore Architectures"](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2008/EECS-2008-134.pdf) and the slides at <https://crd.lbl.gov/assets/pubs_presos/parlab08-roofline-talk.pdf>.
2. Clone the CS Roofline Toolkit, `git clone https://bitbucket.org/berkeleylab/cs-roofline-toolkit.git`. Modify one of the config files in `Empirical_Roofline_Tool-1.1.0/Config` as necessary for the machine you are using.
3. Run the ERT in serial mode on at least _3 different_ node types on HPCC. Report the peak performances and bandwidths (for all caches levels as well as DRAM). Where is the "ridge point" of the roofline for the various cases?
4. Consider the four FP kernels in "Roofline: An Insightful Visual Performance Model for Floating-Point Programs and Multicore Architectures" (see their Table 2). Assuming the high end of operational (i.e., "arithmetic") intensity, how would these kernels perform on the platforms you are testing? What optimization strategy would you recommend to increase performance of these kernels?
5. Address the same questions in (4) for the four kernels given in the Warm-up above. 

To your project write-up, add your plots of the roofline model for the systems you tested, and responses addressing the above questions. 

## Part 2: Enter the Agoge 

### Step 1: Clone your repo on HPCC 

If you haven't already, clone your cmse822 work repo on HPCC. You may need to make a personal access token (PAT) on GH for this:

#### Steps to make a PAT:

- **Generate a Personal Access Token (PAT):**
   - Go to your GitHub account and navigate to **Settings** > **Developer Settings** > **Personal Access Tokens** > **Tokens (classic)**.
   - Click **Generate new token** and select appropriate scopes, such as `repo` (for full repository access).
   - Copy the generated token (you will not see it again).

- **Use the PAT for Cloning:**
   When Git prompts for your password, use the **PAT** instead of your GitHub password.

   ```bash
   git clone https://github.com/cmse822/<your repo name>.git
   ```

   - **Username:** Enter(your GitHub username.
   - **Password:** Paste the PAT.

- **Alternative: Use SSH Authentication**
   If you want to avoid entering a token every time:
   - Set up an SSH key:
     ```bash
     ssh-keygen -t ed25519 -C "your_email@example.com"
     ```
   - Add the public key (`~/.ssh/id_ed25519.pub`) to your GitHub account under **Settings** > **SSH and GPG Keys**.
   - Use the SSH URL for cloning:
     ```bash
     git clone git@github.com:smcouch/cmse822-codex-private.git
     ```

- **Cache the Token (Optional):**
   If you don't want to enter your PAT repeatedly, enable Git credential caching:
   ```bash
   git config --global credential.helper cache
   ```

### Step 2: Build and run agoge 

Go in to the `agoge` directory of your repos and run `make`. If you have the proper modules loaded (default _should_ work), you should get the `agoge_run` binary. Try the Sod shock tube problem:
```bash
./agoge_run problems/SodShockTube.yaml 
```

### Step 3: Estimate big-O scalings

Now you will use the GravityCollapse problem to estimate the big-O scaling of the EulerSolver (a stencil operation) and the GravitySolver (an FFT spectral method). 

- run GravityCollapse with varied domain size (in factors of 2 up to as high as you are able to run)
- plot the zone updates per second for the two code units and estimate the computational complexity from the logarithmic slope of these lines. 

### Step 4: Performance profiling

Now, use Intel VTune to profile agoge, identifying hotspots and getting an empirical measure of the flop rates of the different code units. Consider your results in the context of the roofline models that you produced int he first part of this project.

## Deliverables

1. **Code Repository**: A GitHub repository containing all project code and documentation.
2. **Final Report**: A detailed report explaining the problem, methodology, results, and conclusions.
3. **Presentation**: A presentation summarizing the project and findings.

## Evaluation Criteria

- Code Quality: Readability, organization, and documentation.
- Methodology: Appropriateness and correctness of the computational methods used.
- Results: Accuracy and significance of the results obtained.
- Collaboration: Effective use of version control and teamwork.
- Presentation: Clarity and thoroughness of the final report and presentation.


## Resources

- [GitHub Guides](https://guides.github.com/)
- [Markdown Guide](https://www.markdownguide.org/)
- [Python Documentation](https://docs.python.org/3/)
