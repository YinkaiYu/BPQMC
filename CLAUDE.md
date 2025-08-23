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
- **Debug build with Intel compiler**: `cd src && make FC="mpiifort -fc=ifx" FFLAGS='-check all -traceback -c -I$(HOME)/Modules'` - Compiles with runtime checks and traceback for debugging
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

### Stage 4: Gbar Integration [COMPLETED]
- **Modified `multiply.f90`**:
  - Updated all propagation functions (`propU_L`, `propU_R`, `propT_L`, `propT_R`) to sync Gr operations with Gbar
  - Each function now applies the same operator multiplication to both Gr and Gbar matrices
- **Modified `localU.f90`**:
  - Updated propagation functions (`LocalU_prop_L`, `LocalU_prop_R`) to sync Gr operations with Gbar
  - LocalU_metro function kept unchanged (reserved for Stage 5)
- **Modified `globalK.f90`**:
  - Updated propagation functions (`GlobalK_prop_L`, `GlobalK_prop_R`) to sync Gr operations with Gbar
  - ratioK_fermion function kept unchanged (global updates currently disabled)
- **Modified `dynamics.f90`**:
  - Updated `Dyn_reset` to initialize PropGreen using `Prop%Gbar + ZKRON` instead of `Prop%Gr`
  - Maintains consistency with Gbar-based PQMC approach

### Stage 5: Metropolis Algorithm Adaptation [COMPLETED]
- **Modified `localU.f90`**:
  - Replaced DQMC determinant ratio with PQMC Pfaffian-based acceptance ratio `[r_b]^{N_b}`
  - Implemented `r_b = 1 + Δ_i/N_b * Gbar_{ii}` calculation using boson particle number Nbos
  - Updated Gbar using two-step PQMC formula: `Gbar' = (1+Δ)/r_b * Gbar` with BLAS operations
  - Modified function interface to use Gbar instead of Gr matrix
  - Updated LocalU_prop_L and LocalU_prop_R function calls to pass Prop%Gbar
  - Added CalcBasic module import for Nbos access
  - Corrected acceptance probability to `|ratio_exp * ratio_Pfa|^2` with conjugate term
- **Modified `readme.md`**:
  - Added mathematical definition: `r_b = P†B(2θ,τ)(1+Δ)B(τ,0)P / P†B(2θ,τ)B(τ,0)P = 1 + Δ_i/N_b * Gbar_{ii}`
  - Documented the relationship between matrix and scalar expressions for r_b

### Stage 6: Observable Calculation Adaptation [COMPLETED]
- **Modified `obser_equal.f90`**:
  - Replaced `Prop%Gr` with `Prop%Gbar + ZKRON` for Green's function calculations
  - Updated `Grupc = transpose(Prop%Gbar)` for PQMC-compatible measurement formulas
  - Simplified conjugate Green's functions using direct Gbar operations: `Grdoc = dconjg(transpose(Prop%Gbar))`
  - Maintains all physical observables (density, kinetic energy, correlations) with Gbar-based calculations
  - Mathematical equivalence: `G = Gbar + I` where `Gbar = G - I` stores the deviation from identity matrix

### Stage 7: Sweep Flow Modification [COMPLETED]
- **Modified `local_sweep.f90`**:
  - Updated `Local_sweep_L`: Observable calculations now occur only at `nt = Ltrot/2` (middle imaginary time)
  - Updated `Local_sweep_R`: Observable calculations now occur only at `nt = Ltrot/2` (middle imaginary time)
  - Reduced measurement frequency from `2*Ltrot*Nsweep` to `2*Nsweep` per bin for computational efficiency
  - PQMC algorithm requirement: Measurements at middle imaginary time provide optimal ground state projection
  - Maintained `Nobs` counter logic to properly track reduced measurement frequency

### Next Steps:
- Stage 8: Gr Dependency Elimination - Remove all Gr dependencies, make Gbar the primary Green's function
- Stage 9: Final Integration and Documentation - Complete PQMC conversion with comprehensive testing

**Important**: Each stage must include testing/validation (manual in WSL) and documentation updates to CLAUDE.md, followed by git commits and GitHub synchronization