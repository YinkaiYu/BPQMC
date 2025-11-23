# Repository Guidelines

## Environment & Workflow
- Work from Windows Terminal but perform compilation, debugging, and execution inside WSL.
- Modify build files when needed, yet leave full builds to the user unless explicitly requested.
- Keep cluster runs and submitted jobs under `test/`; use WSL for MPI launches.
- Documentation and code comments stay in English; user interactions can be Chinese.
- Commit changes to `src/` files immediately after edits with concise messages; avoid committing large outputs.

## Project Structure & Modules
- `src/`: Fortran 90 sources and build logic (`Makefile`, `Compile`).
- `test/`: Run area, SLURM script `dqmc`, configs (`confin.txt`, `paramC_sets.txt`, `seeds.txt`).
- `app/`: Archived/deployed binaries (optional).
- `data/`: Inputs or artifacts not created by builds.
- `auto.sh`: Chain build → copy binary → submit SLURM job.

## Build, Run, and Development
- Standard build: `cd src && make` → creates `src/BPQMC.out` (auto-detects MPI wrapper).
- Clean: `cd src && make clean`
- Debug (Intel ifx via mpiifort):
  `cd src && make FC="mpiifort -fc=ifx" FFLAGS='-check all -traceback -c -I$(HOME)/Modules'`
- Quick rebuild and deploy from `src/`: `cd src && ./auto.sh`
- Local run (1 rank): `cd test && mpirun -np 1 ./BPQMC.out`
- Submit to cluster: `cd test && sbatch dqmc`
- Helper: `./auto.sh` (build → clean → copy to `test/` → submit).
- Toolchain overrides: `make -C src HOME=/path/to/Lib_90_new`, `LDFLAGS='-lmkl'`, `SUFFIX='-heap-arrays -fopenmp'`.

## Build System Notes
- Two-layer Makefile setup: top-level `Makefile` selects compilers/dependencies, inner `Compile` performs the Fortran builds.
- Automatically detects MPI Fortran wrappers (`mpiifort`, `mpiifx`, `mpifort`, `mpif90`, `gfortran`) and links modules from `/home/*/Lib_90_new/`.
- Output executable is `BPQMC.out`; auxiliary scripts copy it into `test/` for runs.

## Execution Environment
- `test/` directory holds SLURM job script `dqmc` plus runtime inputs (`confin.txt`, `paramC_sets.txt`, `seeds.txt`).
- Run MPI jobs with `mpirun -np N ./BPQMC.out` (N processes) or via SLURM submission.
- Generated observables and logs should remain inside `test/` unless explicitly archived.

## Coding Style & Naming
- Fortran 90/95; 4-space indent; aim ≤ 100 columns.
- Filenames: lowercase_with_underscores (e.g., `process_matrix.f90`).
- Add `implicit none`; specify `intent(in|out|inout)`; prefer descriptive snake_case names.
- Keep docs/comments concise and in English.

## Testing Guidelines
- No formal unit tests; validate numerics by small runs in `test/`.
- Check `output.log`, `info.txt`, and observables for regressions; use fixed seeds in `test/seeds.txt`.
- Avoid committing large binaries; include minimal samples only when necessary.

## Commits & Pull Requests
- Commits: imperative subject (≤ 72 chars), optional scope (e.g., `lattice:`), concise body explaining why.
- PRs: describe rationale, parameter changes, performance notes; attach small before/after logs and link issues.

## Architecture Overview
- `main.f90`: MPI-enabled driver that manages warm-up, sweeps, measurements, and high-level control flow.
- `model.f90`: Defines lattice geometry, kinetic/Hubbard operators, auxiliary fields, and initializes trial wave functions.
- `initial_state.f90`: Builds projector trial states, evaluates gaps, and prepares left/right wave functions.
- `lattice.f90`: Encapsulates lattice construction and geometric utilities (kagome structure).
- `fields.f90`: Handles auxiliary field representation and updates.
- `process_matrix.f90`: Maintains propagator data (Gbar, wrap lists), imaginary-time propagation, and stabilization hooks.
- `localU.f90`: Implements local imaginary-time propagation and Metropolis updates.
- `local_sweep.f90`: Orchestrates left/right sweeps, measurement scheduling, and counter management.
- `multiply.f90`: Performs propagation matrix multiplications using stabilized routines.
- `obser_equal.f90`, `obser_tau.f90`: Collect equal-time and imaginary-time observables.
- `stabilization.f90`: Provides numerical stabilization for propagation matrices via LAPACK vector normalization.
- `globalK.f90` (currently disabled globally): Implements optional global updates using Gbar-based conversions.

## Key Features
- Supports projector QMC with distinct thermalization and measurement phases.
- Local and (optional) global update strategies; global flow currently off by default.
- Uses `Gbar = G - I` storage for Green’s function, reconstructing full `G` only on demand.
- Includes Fourier analysis utilities for momentum observables and matrix stabilization for numerical robustness.
- Configurable imaginary-time measurements with focus on mid-projection sampling.

## Development Context & PQMC Transition
- Code base transitioned from finite-temperature DQMC to zero-temperature PQMC with canonical ensemble (fixed boson number).
- Chemical potential terms were removed from kinetic Hamiltonian during the transition.
- Propagator structure now stores only `Gbar`, reducing memory footprint and aligning all modules (`process_matrix`, `stabilization`, `multiply`, `localU`, `local_sweep`, `globalK`) to operate directly on `Gbar`.
- Measurement routines (`obser_equal`, others) reconstruct `G = Gbar + I` as needed, keeping physical observables intact.
- Each major stage of the conversion should be documented here and validated manually in WSL after changes.

## Current Refactor: Rank-1 Propagation (Bosonic PQMC)

This section tracks the ongoing refactor to move from a Green’s-function-propagation algorithm to a rank-1 wavefunction propagation scheme, as now documented in `readme.md`.  The high-level goal is to:
- store and propagate only the left/right rank-1 wavefunctions at each time slice;
- reconstruct `G` (or `Gbar`) on demand for measurements;
- completely remove the need to maintain an `N×N` `Gbar` during propagation.

### Plan & Milestones

0. **Algorithm documentation (DONE)**  
   - Update `readme.md` to describe the rank-1 propagation scheme, the role of `UUL`, `UUR` and the overlap scalar, and how to reconstruct `G` / `Gbar` at measurement times.

1. **Wavefunction storage & normalization (DONE)**  
   - In `local_sweep.f90`, change the time-slice handling so that the rank-1 wavefunctions are normalised and stored at *every* imaginary-time slice (no more `Nwrap`-based stabilisation interval from `calc_basic.f90`).  
   - Use the `stabilization.f90` utilities and `process_matrix.f90` data structures to maintain full `URlist` and `ULlist` for later use.  
   - After this step, `Gbar` is still maintained as before, but wavefunction storage is sufficient for a rank-1-only algorithm.

2. **Decouple local updates from `Gbar` – `LocalU_metro` (DONE)**  
   - Extend `type :: Propagator` in `process_matrix.f90` to include a complex scalar `overlap` that tracks the inner product of `UUL` and `UUR`.  
   - Rewrite `LocalU_metro` in `localU.f90` to depend only on `UUL`, `UUR`, and `overlap` for computing local Metropolis ratios, while still updating `Gbar` in parallel for cross-checking.  
   - Ensure that when an auxiliary-field update is accepted, `UUR` and `overlap` are updated consistently with the new configuration.

3. **Adapt propagation routines to the new interface (DONE)**  
   - Update `LocalU_prop_L`, `LocalU_prop_R`, and any related routines in `localU.f90` to:  
     - retrieve `UUR` (or `UUL`) from `URlist`/`ULlist`;  
     - compute the overlap between `UUL` and `UUR`;  
     - loop over sites to call `LocalU_metro` purely in terms of wavefunctions;  
     - apply the imaginary-time propagation with `Op_U%mmult_L` / `Op_U%mmult_R` in the appropriate order.  
   - Keep `Gbar` propagation active for now so that behaviour can be compared to the original algorithm.

4. **Verify local updates without `Gbar` (TODO)**  
   - Search for and remove residual dependencies on `Gbar` in the local-update logic, beyond the parallel-maintained version used for checks.  
   - Run paired tests (same seeds in `test/`) comparing observables or acceptance statistics between:  
     - the original `Gbar`-based update path;  
     - the new rank-1-wavefunction-based path.  
   - Only proceed once the two paths agree to the desired numerical tolerance.

5. **Remove `Gbar` from equal-time measurements (TODO)**  
   - In `obser_equal.f90` (`Obs_equal_calc`), reconstruct `G`/`Gbar` using `UUL`, `UUR`, and `overlap` following the formulas in `readme.md` instead of using stored `Gbar`.  
   - Pay special attention to the `Ltrot == 0` branch in `local_sweep.f90` to ensure that measurements at zero imaginary-time length are still handled correctly.  
   - Perform parallel checks: compare reconstructed `Gbar` with the still-maintained `Gbar` to verify correctness.

6. **Clean up remaining `Gbar` dependencies (TODO)**  
   - Search the codebase for any other uses of `Gbar` (e.g. global updates, imaginary-time correlation functions).  
   - For features that are currently disabled or not actively maintained, introduce lightweight placeholders or compile-time guards so the code still compiles without requiring a fully propagated `Gbar`.

7. **Remove redundant `Gbar` maintenance (TODO)**  
   - Once all functional dependencies on `Gbar` are removed and verified, strip out the parallel `Gbar` updates from propagation and local-update routines.  
   - Verify that the resulting implementation uses only rank-1 wavefunctions and overlaps in the main Monte Carlo loop, and that performance scales as expected (matrix–vector only, no matrix–matrix updates).

8. **Documentation and cleanup (TODO)**  
   - Update `readme.md` and any in-code comments to reflect the final design (rank-1 propagation, overlap handling, storage strategy).  
   - Remove temporary notes from this section, keeping only a concise summary of the final architecture and any non-obvious design choices.  
   - If helpful, record a short “developer note” here on how to re-enable a `Gbar`-like path for debugging, without reintroducing it into the hot path.
