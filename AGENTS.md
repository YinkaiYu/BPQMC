# Repository Guidelines

## Environment & Workflow
- Work from Windows Terminal but perform compilation, debugging, and execution inside WSL.
- Modify build files when needed, yet leave full builds to the user unless explicitly requested.
- Keep cluster runs and submitted jobs under `test/`; use WSL for MPI launches.
- Documentation and code comments stay in English; user interactions can be Chinese.
- Commit changes to `src/` files immediately after edits with concise messages; avoid committing large outputs.

## Project Structure & Modules
- `src/`: Fortran 90 sources and build logic (`Makefile`, `Compile`).
- `test/`: Run area, SLURM script `dqmc`, configs (`confin.txt`, `paramC_sets.txt`, `seeds.txt`), HMC benchmark/tuning/report scripts.
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
- `test/production_hmc.py` is the main driver for production tuning and direct local-vs-HMC benchmarks.
- `test/render_hmc_report.py` converts benchmark JSON into PNG plots and a Markdown summary.
- `test/hmc_report_template.ipynb` is the notebook entry point for interactive post-processing.
- `test/small_hmc_benchmark.py` is the correctness-first triangular `L=6` workflow for local seeds, HMC tune scans, and strict local-vs-HMC benchmarks.
- `test/render_small_hmc_stage_report.py` aggregates passed small-benchmark directories into a staged PNG/Markdown report.
- `test/small_hmc_stage_report.ipynb` is the notebook entry point for those staged small-benchmark reports.

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
  - run `test/tune_hmc.py` or `test/production_hmc.py tune`
  - run direct `local` vs `HMC` comparison with `test/benchmark_hmc.py` or `test/production_hmc.py benchmark`
  - render figures with `test/render_hmc_report.py`
- For the triangular small-parameter correctness campaign, use:
  - `test/small_hmc_benchmark.py local-seed` to build validated `confout.txt` warm starts
  - `test/small_hmc_benchmark.py tune` for short `warm=0` full-field HMC scans
  - `test/small_hmc_benchmark.py benchmark` for long strict local-vs-HMC checks
  - `test/render_small_hmc_stage_report.py` to produce staged visual reports from the passing cases
- Use debug builds when a new lattice path crashes:
  - `cd src && make clean && make FFLAGS='-O0 -g -traceback -check all -fpe0 -c -I/home/yyk/Lib_90_new/Modules'`
- Production triangular targets in this branch are typically:
  - `L=12` for the full benchmark grid
  - `L=21` for spot checks
  - `beta=256`, `Dtau=0.001`, `U1=0`, `iniHam=5`, `iniTwist=1e-4`
- Do not stop after a single benchmark mismatch. First distinguish:
  - lattice-specific observable bug
  - HMC tuning/warm-up issue
  - stabilization or force inconsistency
- Do not trust a warm-start `confout.txt` blindly. `test/hmc_tools.py` now validates its line count against `(Naux, Ndim, Ltrot)` and `test/small_hmc_benchmark.py` will skip invalid warm starts.
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
- `test/tune_hmc.py` is still useful for small development scans on built-in parameter sets.
- `test/benchmark_hmc.py` is the fast regression benchmark for development-sized cases.
- `test/production_hmc.py tune` is the main production tuning entry point.
- `test/production_hmc.py benchmark` runs direct local-vs-HMC comparisons over a user-specified `(L, Nbos, U2)` grid and writes JSON/CSV summaries.
- `test/render_hmc_report.py` reads `production_benchmark.json` and writes:
  - `acceptance.png`
  - `ess_per_sec.png`
  - `tau_int.png`
  - `zscore_heatmap.png`
  - `report.md`
- `test/hmc_report_template.ipynb` can be pointed at the same JSON for interactive visualization.
- `test/small_hmc_benchmark.py` writes:
  - `small_tune.csv`, `small_tune.json`, `recommended_hmc.json`
  - `small_benchmark_cases.csv`, `small_benchmark_observables.csv`, `small_benchmark_samples.csv`
  - `small_benchmark_samples.csv` now keeps the key trace observables for all repeats, so repeated thermal-cut scans can be reconstructed without reopening each run manually
- `test/render_small_hmc_stage_report.py` reads one or more strict benchmark directories and writes:
  - `stage_cases.csv`
  - `stage_observables.csv`
  - `stage_samples.csv`
  - `stage_summary.json`
  - `overview.png`
  - `observable_*_n*.png`
  - `trace_*.png`
  - `tune_*.png`
  - `report.md`
- Current staged small-benchmark reports live under `data/triangular_hmc_small_benchmark/stage_report_*`.
- The current staged checkpoint used during development is `data/triangular_hmc_small_benchmark/stage_report_v6a`.
- `test/small_hmc_stage_report.ipynb` now exposes both `NBOS_SELECT` and `TRACE_REPEAT` for interactive slicing of the rendered staged report.
- `HMC-REFERENCE.md` is the running handoff note for validated health points, unresolved slow modes,
  tuning heuristics, and production/HPC-oriented lessons.

## Preferred Commit Granularity
- Commit source changes in small, reviewable chunks.
- Prefer separate commits for:
  - lattice/interface changes
  - HMC algorithm changes
  - observable or benchmark bug fixes
  - new test tooling
  - documentation updates
- Do not mix long-running output files with source commits.
