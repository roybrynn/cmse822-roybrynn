Project 2: Implement p2p MPI in FD Unit

Objectives:
    - Parallelize the FD solver using domain decomposition and blocking p2p MPI
    - Develop and run first complex MPI code
Dwarf: 
    - structured grid
MPI topic: 
    - point-to-point, halo exchange, blocking
    - running MPI programs 
Unit Tests:
    - Sedov but scaling, strong and weak
Challenge:
Most zone-updates per wall-clock second, but parallel now