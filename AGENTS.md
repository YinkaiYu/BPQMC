# Repository Guidelines

## Environment & Workflow
- Work from Windows Terminal but perform compilation, debugging, and execution inside WSL.
- Modify build files when needed, yet leave full builds to the user unless explicitly requested.
- Keep cluster runs and submitted jobs under `test/`; use WSL for MPI launches.
- Documentation and code comments stay in English; user interactions can be Chinese.
- Commit changes to `src/` files immediately after edits with concise messages; avoid committing large outputs.

## Project Structure & Modules
- `src/`: Fortran 90 sources and build logic (`Makefile`, `Compile`).
- `test/`: Active run area for SLURM/runtime inputs plus production/HPC HMC scripts.
- `test/README.md` documents the cleaned test layout. Top-level `test/` is now production/HPC-facing; archived benchmark/check tools live in subdirectories there.
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
- `paramC_sets.txt` may now start with a lattice header: `kagome` or `triangular`.
- `test/production_hmc.py` is the main driver for production tuning, staged HMC thermalization runs, summary collection, and report rendering.
- `test/render_hmc_report.py` converts production summary JSON into PNG plots and a Markdown summary.
- `test/hmc_report_template.ipynb` is the notebook entry point for interactive post-processing.
- `test/archive_small_benchmark/` contains the closed triangular `L=6` local-vs-HMC correctness campaign:
  `small_hmc_benchmark.py`, `render_small_hmc_stage_report.py`, `render_small_hmc_full_grid.py`, and the archived benchmark notebooks.
- `test/archive_manual_checks/` contains the older Fortran-side force/ratio/round-trip/check programs and their outputs.
- `test/runtime_outputs/root_run_20260319/` contains the last root-level manual run dump that used to clutter `test/`.
- The small-parameter correctness campaign is considered closed. Keep the final visible report at `data/triangular_hmc_small_benchmark/full_grid_progress_v3` and keep raw campaign data under `data/triangular_hmc_small_benchmark/archive_20260319`.

## Coding Style & Naming
- Fortran 90/95; 4-space indent; aim ≤ 100 columns.
- Filenames: lowercase_with_underscores (e.g., `process_matrix.f90`).
- Add `implicit none`; specify `intent(in|out|inout)`; prefer descriptive snake_case names.
- Keep docs/comments concise and in English.

## Testing Guidelines
- No formal unit tests; validate numerics by small runs in `test/`.
- Check `output.log`, `info.txt`, and observables for regressions; use fixed seeds in `test/seeds.txt`.
- Avoid committing large binaries; include minimal samples only when necessary.
- For HMC work, keep the test loop modular:
  - compile after each source-module change
  - run a small smoke test
  - run `test/production_hmc.py tune`
  - run `test/production_hmc.py stage` for long HMC-only thermalization traces
  - run `test/production_hmc.py collect` if a finished stage needs its summaries rebuilt from `runs/`
  - run `test/production_hmc.py report` or `test/render_hmc_report.py`
  - run `test/production_hmc.py benchmark` only when a legacy direct production comparison is still needed
- The triangular small-parameter correctness campaign is archived under `test/archive_small_benchmark/`.
  Use it only as a historical reference or if the final report must be regenerated.
- Use debug builds when a new lattice path crashes:
  - `cd src && make clean && make FFLAGS='-O0 -g -traceback -check all -fpe0 -c -I/home/yyk/Lib_90_new/Modules'`
- Production triangular targets in this branch are typically:
  - `L=12` for the full benchmark grid
  - `L=21` for spot checks
  - `beta=256`, `Dtau=0.001`, `U1=0`, `iniHam=5`, `iniTwist=1e-4`
- The active production ramp is slower and more granular than the final target:
  - it may start from `L=6`
  - it first gates on `squareOcc` and `IPR`
  - then it promotes gradually in `(Nbos, U2, beta, dtau, L)`
- Current workstation production examples already distinguish:
  - a healthy projector ladder at `L=6, Nbos=1e3, U2=1`, currently healthy through `beta=128`, `dtau=0.005`
  - a slow-drift rung at `L=6, Nbos=1e4, U2=1`
  - a stronger-coupling but stage-healthy rung at `L=6, Nbos=1e4, U2=10`
  - the current next rung is `L=6, Nbos=1e3, U2=1, beta=160, dtau=0.004`
- Do not stop after a single benchmark mismatch. First distinguish:
  - lattice-specific observable bug
  - HMC tuning/warm-up issue
  - stabilization or force inconsistency
- Do not trust a warm-start `confout.txt` blindly. `test/hmc_tools.py` validates its line count against `(Naux, Ndim, Ltrot)`. The archived small-benchmark driver also checks this before reusing a warm start.
- For strict small-benchmark health points, `doubleOcc` alone is not enough. Re-check thermal cut against at least `squareOcc`, `nearestOcc`, `IPR`, and `PF_Gamma`.

## Commits & Pull Requests
- Commits: imperative subject (≤ 72 chars), optional scope (e.g., `lattice:`), concise body explaining why.
- PRs: describe rationale, parameter changes, performance notes; attach small before/after logs and link issues.

## Architecture Overview
- `main.f90`: MPI-enabled driver that manages warm-up, sweeps, measurements, and high-level control flow.
- `model.f90`: Defines lattice geometry, kinetic/Hubbard operators, auxiliary fields, and initializes trial wave functions.
- `initial_state.f90`: Builds projector trial states, evaluates gaps, and prepares left/right rank-1 wave functions.
- `lattice.f90`: Encapsulates runtime-selectable lattice construction and geometric utilities for `kagome` and `triangular`.
- `fields.f90`: Handles auxiliary field representation and updates.
- `process_matrix.f90`: Stores propagated vectors (`UUR`, `UUL`), overlap scalar, wrap lists, and bookkeeping.
- `localU.f90`: Implements local imaginary-time propagation and Metropolis updates using only rank-1 data.
- `local_sweep.f90`: Orchestrates left/right sweeps, measurement scheduling, and counter management.
- `multiply.f90`: Propagates the rank-1 wavefunctions through kinetic and interaction operators.
- `obser_equal.f90`, `obser_tau.f90`: Collect equal-time and imaginary-time observables by reconstructing `G` on demand.
- `stabilization.f90`: Provides normalization/orthogonalization for the propagated vectors.
- `global_update.f90`: Active HMC implementation for the rank-1 flow.
- `globalK.f90`: Obsolete placeholder; not part of the active HMC path.

## Current Design (Rank-1 PQMC)
- Monte Carlo loop propagates only the left/right trial vectors and their overlap; no `N×N` Green’s matrix is stored.
- Local Metropolis ratios use the outer product of the current vectors; propagation is matrix–vector only.
- Equal-time and time-sliced Green’s functions are reconstructed on demand from `UUR`, `UUL`, and `overlap` when measuring.
- Wrapping/normalization keeps vectors well-conditioned and logs drift for diagnostics.
- Imaginary-time measurements and Fourier analysis operate on reconstructed observables while runtime outputs stay under `test/`.
- The code now chooses `Norb` and `Nbond` at runtime from `lattice_type` instead of compiling separate executables.
- The triangular path uses `Norb=1`, `Nbond=3`; the kagome path keeps `Norb=3`, `Nbond=2`.
- HMC uses the same propagation order as the local sweeps and keeps cumulative log norms for stable action evaluation.
- One full-size HMC force buffer has been removed to reduce peak memory usage.

## HMC and Triangular Workflow
- `test/production_hmc.py tune` is the main short-scan entry point for candidate HMC parameters.
- `test/production_hmc.py stage` runs HMC-only production stage jobs and stores sample-index traces.
- `test/production_hmc.py collect` rebuilds `production_stage*.json/csv` from completed `runs/`.
- `test/production_hmc.py report` renders Markdown + PNG summaries from `production_stage.json`, `production_tune.json`, or the legacy `production_benchmark.json`.
- `test/render_hmc_report.py` now understands three summary modes:
  - `tune`
  - `stage`
  - `benchmark`
- The stage workflow writes:
  - `production_stage.json`
  - `production_stage_cases.csv`
  - `production_stage_observables.csv`
  - `production_stage_runs.csv`
  - `production_stage_samples.csv`
  - `report/acceptance.png`
  - `report/ess_per_sec.png`
  - `report/tau_int.png`
  - `report/drift_ratios.png`
  - `report/trace_*.png`
  - `report/report.md`
- `test/hmc_report_template.ipynb` can be pointed at the same JSON for interactive visualization.
  It now supports `stage`, `tune`, and `benchmark` summaries and rewrites `report.md` image paths for notebook display.
- `test/archive_small_benchmark/small_hmc_benchmark.py` writes:
  - `small_tune.csv`, `small_tune.json`, `recommended_hmc.json`
  - `small_benchmark_cases.csv`, `small_benchmark_observables.csv`, `small_benchmark_samples.csv`
  - `small_benchmark_samples.csv` now keeps the key trace observables for all repeats, so repeated thermal-cut scans can be reconstructed without reopening each run manually
- `test/archive_small_benchmark/render_small_hmc_stage_report.py` reads one or more strict benchmark directories and writes:
  - `stage_cases.csv`
  - `stage_observables.csv`
  - `stage_samples.csv`
  - `stage_summary.json`
  - `overview.png`
  - `observable_*_n*.png`
  - `trace_*_rep*.png`
  - `tune_*.png`
  - `report.md`
- `test/archive_small_benchmark/render_small_hmc_full_grid.py` is a convenience wrapper around the stage renderer. It points at the repository's current best 15-point directory selection and defaults to `data/triangular_hmc_small_benchmark/full_grid_progress_v3`.
- The stage renderer depends on `pandas` and `matplotlib`; in this repository the recommended interpreter is `/home/yyk/conda/envs/notebook/bin/python`.
- Current staged small-benchmark reports live under `data/triangular_hmc_small_benchmark/stage_report_*`.
- The current staged checkpoint used during development is `data/triangular_hmc_small_benchmark/stage_report_v6a`.
- The current full-grid rendered visual report is `data/triangular_hmc_small_benchmark/full_grid_progress_v3`.
- The archived raw benchmark/tuning directories now live under `data/triangular_hmc_small_benchmark/archive_20260319`.
- `test/archive_small_benchmark/small_hmc_stage_report.ipynb` exposes both `NBOS_SELECT` and `TRACE_REPEAT` for interactive slicing of the rendered staged report.
- `HMC-REFERENCE.md` is the running handoff note for validated health points, unresolved slow modes,
  tuning heuristics, and production/HPC-oriented lessons.
- `HMC-REFERENCE.md` should always state the current healthy rung, next rung, and main blocker so the production ramp can resume after context loss.
- For current production bring-up, the primary decision is not local-vs-HMC agreement.
  The first question is whether `squareOcc` and `IPR` stabilize across sample index, seeds, and initial-state choices.
- A current tuning lesson already written into `HMC-REFERENCE.md`:
  on stronger-coupling rungs, shrinking `dt` can matter more than chasing the old acceptance target band.
- The renderer can also be used on unresolved benchmark directories to build progress reports with per-repeat traces.
  A current example is `data/triangular_hmc_small_benchmark/progress_report_n1000_u1e0_v2`.
  Another current example is `data/triangular_hmc_small_benchmark/progress_report_n100_u1e1_v3`.
- Future work should not treat the archived local-vs-HMC small benchmark as the main task.
  The benchmark proved small-parameter correctness well enough.
  The next priority is production-grade HMC thermalization and HPC-scale workflows, not further local-update comparison on difficult large-`Nbos`, large-`U2` points.

## Preferred Commit Granularity
- Commit source changes in small, reviewable chunks.
- Prefer separate commits for:
  - lattice/interface changes
  - HMC algorithm changes
  - observable or benchmark bug fixes
  - new test tooling
  - documentation updates
- Do not mix long-running output files with source commits.
