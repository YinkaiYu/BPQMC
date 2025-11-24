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
   - `local_sweep.f90` stores normalised left/right wavefunctions at wrap intervals (`Nwrap`) and always at `nt=Ltrot`.  
   - `WrapList` holds full `URlist/ULlist`; `Gbar` still propagated in parallel for checks.

2. **Decouple local updates from `Gbar` – `LocalU_metro` (DONE)**  
   - `Propagator` carries scalar `overlap`; `LocalU_metro` uses only `UUL/UUR/overlap` for ratios while still updating `Gbar` for comparison.  
   - Acceptance updates both `UUR` and `overlap` consistently.

3. **Propagate both wavefunctions explicitly (DONE)**  
   - `LocalU_prop_L`: compute overlap → Metropolis over sites → propagate `Gbar`/`UUL`/`UUR` (L: `mmult_L` +1, `mmult_R` -1).  
   - `LocalU_prop_R`: propagate `Gbar`/`UUL`/`UUR` (R: `mmult_R` +1, `mmult_L` -1) → compute overlap → Metropolis over sites.  
   - All other propagation routines simultaneously advance `UUL`/`UUR` (no per-U WrList helpers). Wrapping steps re-orthonormalise and log differences between propagated vectors and stored lists to monitor drift.

4. **Verify local updates without `Gbar` (DONE)**  
   - Wrapping now occurs every `Nwrap` slices (plus `nt=Ltrot`), with stored/propagated vector diffs logged; `Gbar` maintained in parallel for cross-checks.

5. **Remove `Gbar` from equal-time measurements (DONE)**  
   - `obser_equal.f90` now reconstructs `G`/`Gbar` directly from `UUL`, `UUR`, and the propagated `overlap` (rank-1 outer product); the `Ltrot == 0` path uses the same reconstruction.  
   - The measurement path no longer depends on the stored `Gbar`; reconstruction is now the sole source.

6. **Clean up remaining `Gbar` dependencies (DONE)**  
   - Imaginary-time (`dynamics.f90`) now seeds time-sliced Greens via reconstructed rank-1 `G` instead of stored `Gbar`.  
   - Global update path (`globalK.f90`) is stubbed as disabled for the rank-1 flow so compilation does not rely on `Gbar`.

7. **Remove redundant `Gbar` maintenance (DONE)**  
   - Local propagation and Metropolis updates no longer multiply/update `Gbar`; only rank-1 `UUL/UUR/overlap` are propagated (`localU.f90`, `multiply.f90`, `local_sweep.f90`, `stabilization.f90`).  
   - `Prop%Gbar` is left allocated for compatibility but is no longer touched in the hot path; Monte Carlo now runs purely on vector operations.

8. **Documentation and cleanup (TODO)**  
   - Update `readme.md` and any in-code comments to reflect the final design (rank-1 propagation, overlap handling, storage strategy).  
   - Remove temporary notes from this section, keeping only a concise summary of the final architecture and any non-obvious design choices.  
   - If helpful, record a short “developer note” here on how to re-enable a `Gbar`-like path for debugging, without reintroducing it into the hot path.
