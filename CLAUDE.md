# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Environment

- **Working Environment**: Windows terminal with Fortran compilation/execution in WSL
- **Important**: For compilation, debugging, or running code, manually operate in WSL environment
- **Compilation Rule**: Claude should only modify build files (Makefile, Compile) but let user handle actual compilation in WSL
- **Language Rule**: All code comments and documentation files (.md) must be written in English. Conversations with user can be in Chinese.
- **Git Policy**: Commit changes to src/ files immediately after modifications with brief update description

## Build Commands

- **Build executable**: `cd src && make` - Compiles the Fortran source files into `bosonDQMC.out`
- **Clean build artifacts**: `cd src && make clean` - Removes object files
- **Full build and deploy**: `./auto.sh` - Builds, cleans, copies executable to test directory, and submits SLURM job
- **Quick build from src**: `cd src && ./auto.sh` - Builds and copies executable to test directory

## Architecture Overview

This is a Bosonic Projector Quantum Monte Carlo (BPQMC) implementation written in Fortran 90 with MPI parallelization.

### Core Components

- **Main Program** (`main.f90`): MPI-based DQMC simulation driver that orchestrates warm-up, sweeps, and measurements
- **Model System** (`model.f90`): Manages lattice structure (kagome), operators (kinetic and Hubbard), auxiliary field configurations, and initial state wave functions
- **Initial State** (`initial_state.f90`): Manages trial wave functions for projector algorithm, including ground state calculation and energy gap computation
- **Local Sweeps** (`local_sweep.f90`): Implements local Monte Carlo updates and measurement collection
- **Lattice** (`lattice.f90`): Kagome lattice implementation with geometric operations
- **Fields** (`fields.f90`): Auxiliary field configuration management
- **Process Matrix** (`process_matrix.f90`): Defines imaginary-time propagators
- **Local U Operations** (`localU.f90`): Handles imaginary-time propagation and Monte Carlo updates
- **Matrix Multiplication** (`multiply.f90`): Handles imaginary-time propagation operations
- **Observables**: Equal-time (`obser_equal.f90`) and imaginary-time (`obser_tau.f90`) measurements
- **Stabilization** (`stabilization.f90`): Numerical stabilization for matrix operations

### Build System

- Uses a two-layer Makefile system: outer `Makefile` handles compiler detection and dependency setup, inner `Compile` handles actual compilation
- Automatically detects MPI Fortran wrappers (mpiifort, mpiifx, mpifort, mpif90, gfortran)
- Links against external libraries in `/home/*/Lib_90_new/` for modules, numerical routines, and random number generation
- Produces `bosonDQMC.out` executable

### Execution Environment

- **Test Directory**: Contains SLURM job script (`dqmc`) and configuration files
- **MPI Execution**: Run via `mpirun -np N ./bosonDQMC.out` where N is number of processes
- **Configuration**: Uses `confin.txt`, `paramC_sets.txt`, and `seeds.txt` for simulation parameters

### Key Features

- Supports both thermalization and measurement phases
- Implements local and global update schemes (global currently disabled)
- Fourier transform analysis for momentum-space observables
- Matrix stabilization for numerical stability
- Configurable imaginary-time measurements

## Development Context

The current program in the src/ directory implements finite-temperature determinant quantum Monte Carlo algorithm. We are gradually modifying it into a zero-temperature projector quantum Monte Carlo (PQMC) algorithm. Key changes made:

- **Algorithm Type**: Converting from finite-temperature DQMC to zero-temperature projector QMC
- **Ensemble**: Changed from grand canonical (fixed chemical potential μ) to canonical (fixed boson particle number Nbos)
- **Hamiltonian**: Removed chemical potential term from the kinetic Hamiltonian in non_interact.f90

Future modifications will implement the program according to the README.md specifications for the projector algorithm.

## Current PQMC Conversion Progress

### Stage 1: Data Structure Modifications [COMPLETED]
- **Modified `process_matrix.f90`**:
  - Updated `Propagator` type: UUR(Ndim,1), UUL(1,Ndim) for vector operations
  - Added Gbar(Ndim,Ndim) to store Ḡ = G - I 
  - Kept Gr for temporary compatibility during transition
  - Removed VUR, VUL, DUR, DUL arrays (no longer needed for PQMC)
  - Updated `WrapList` type: URlist(Ndim,1,0:Nst), ULlist(1,Ndim,0:Nst)
  - Modified initialization to use Init%PR and Init%PL trial wave functions

### Stage 2: Numerical Stabilization Algorithm [COMPLETED]
- **Modified `stabilization.f90`**:
  - Replaced matrix QR/UDV decomposition with LAPACK-optimized vector normalization
  - Implemented `stab_UR` and `stab_UL` using DZNRM2 and ZDSCAL for efficient normalization
  - Rewrote `stab_green` to calculate Gbar using ZGERC for outer product operations
  - Updated all Wrap functions (Wrap_pre, Wrap_L, Wrap_R) to use Gbar with numerical stability checks
  - Removed obsolete functions: QDRP_decompose, complex stab_green_big implementation
  - Added compatibility layer: Gr = Gbar + I where needed

### Stage 3: Propagation Algorithm Adaptation [COMPLETED]
- **Modified `operator_Hubbard.f90`**:
  - Updated `opU_mmult_R` and `opU_mmult_L` to use assumed-shape arrays `dimension(:,:)`
  - Automatic compatibility with both matrix (Ndim,Ndim) and vector (Ndim,1)/(1,Ndim) operations
  - Use `size()` function for dynamic array bounds, eliminating need for separate vector functions
- **Modified `non_interact.f90`**:
  - Updated `opT_mmult_R` and `opT_mmult_L` to use assumed-shape arrays `dimension(:,:)`
  - Dynamic temporary array allocation using `dimension(size(Mat,1), size(Mat,2))`
  - Unified interface handles all matrix and vector operations seamlessly
- **Zero modification required**: `multiply.f90` and `localU.f90` work without changes due to elegant Fortran overloading

### Next Steps:
- Stage 4: Gbar Integration - Gradually replace Gr usage with Gbar throughout the codebase
- Stage 5: Metropolis Algorithm Adaptation - Modify LocalU_metro according to readme.md PQMC specifications  
- Stage 6: Observable Calculation Adaptation - Rewrite obser_equal.f90 to use Gbar-based measurements
- Stage 7: Sweep Flow Modification - Update local_sweep.f90 to calculate observables only at Ltrot/2 (middle imaginary time)
- Stage 8: Gr Dependency Elimination - Remove all Gr dependencies, make Gbar the primary Green's function
- Stage 9: Final Integration and Documentation - Complete PQMC conversion with comprehensive testing

**Important**: Each stage must include testing/validation (manual in WSL) and documentation updates to CLAUDE.md, followed by git commits and GitHub synchronization