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
  It now also supports `--run-timeout-sec` on execute paths so a hung `mpirun` can be cut off
  once `info.txt` and the observables are already complete.
- `test/render_hmc_report.py` converts production summary JSON into PNG plots and a Markdown summary.
- `test/render_hmc_live_progress.py` converts raw in-flight stage traces into PNG plots and a Markdown summary before a repeat has fully finished.
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
  - run `test/production_hmc.py collect-tune` if a long tune finishes only partially but its completed `runs/` should already be summarized
  - run `test/production_hmc.py stage` for long HMC-only thermalization traces
  - run `test/production_hmc.py collect` if a finished or partially finished stage needs its summaries rebuilt from `runs/`
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
  - `test/production_hmc.py` can now expand multiple initial-state families with
    `--ini-type-values` and `--ini-ampl-values`
  - `test/production_hmc.py` can now also expand explicit seed-family lists with
    `--seed-base-values 50001,51001,52001`, which is the preferred way to probe
    `U2=1000` seed-family sensitivity without hand-copying many work roots;
    that explicit list is now used verbatim across the whole scan rather than shifted per grid point
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-uniform`
    to split the per-time-slice spatially uniform mode away from the residual `--hmc-mass`
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-shell1`
    to split the triangular lowest nonzero momentum shell away from both the residual
    `--hmc-mass` and the exact uniform mode
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-shell2`
    for the second triangular nonzero momentum shell
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-lowk`
    for a broader combined low-|k| split that groups the first three nonzero momentum shells
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-lowk-shells`
    to widen that grouped low-|k| subspace beyond the default three shells
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-midk`
    and `--hmc-mass-spatial-midk-shells`
    to split the next grouped momentum-shell band away from the lightest low-|k| block
  - `test/production_hmc.py` also supports `--hmc-mass-spatial-shell-map`
    to pass a comma-separated per-shell mass profile for the first few triangular
    nonzero momentum shells; when present, it supersedes the grouped `lowk` / `midk`
    shortcuts
- Current workstation production examples already distinguish:
  - a healthy projector ladder at `L=6, Nbos=1e3, U2=1`, currently healthy through `beta=160`, `dtau=0.004`
  - a partially collected but representative `beta=192`, `dtau=0.002` rung on the same baseline point,
    with `2/3` repeats currently complete and `stable_window`
  - the representative baseline ladder only moves `squareOcc` / `IPR` by about `7.6e-4` / `7.7e-4`
    across `beta=32 -> 192`, `dtau=0.01 -> 0.002`, so the current workstation priority has shifted
    from deeper projector scans to higher-`U2` scans
  - the old slow-drift rung at `L=6, Nbos=1e4, U2=1` is now mainly a historical scalar-mass reference
  - a stronger-coupling but stage-healthy rung at `L=6, Nbos=1e4, U2=10`
  - the current `U2`-ramp starting geometry is the stage-healthy preconditioned rung
    `L=6, Nbos=1e4, U2=1`, `mass=16`, `uniform_mass=1`, `nfrog=20`, `dt=0.008`
  - the current next rungs are:
    - `L=6, Nbos=1e4, U2=30`, `mass=16`, `uniform_mass=1`, `nfrog=12`, `dt=0.006`, now `healthy`
    - `L=6, Nbos=1e4, U2=100`, `mass=16`, `uniform_mass=1`, `nfrog=20`, `dt=0.001`, current long-stage `rep0` is `stable_window`
    - `L=6, Nbos=1e4, U2=300`, `mass=16`, `uniform_mass=1`, `nfrog=28`, `dt=0.0005`, stage launched after a conservative retune
    - `L=6, Nbos=1e4, U2=1000` is the current blocker:
      the old scalar/uniform-only `mass=16`, `uniform_mass=1`, `nfrog=20`, `dt=0.0004`
      reference stage is now formally `strong_drift`
      `shell1_mass=1` has also been closed as too aggressive
      shell1-only `shell1_mass=4` and `shell1_mass=8` completed stages are now also closed as `strong_drift`
      the first `shell2` family `shell1_mass=8`, `shell2_mass=4` is now also closed as `strong_drift`
      the lighter `shell2` follow-up `shell1_mass=4`, `shell2_mass=2` is no longer the main lead;
      all cheap candidates still kept `acceptance=1` with huge positive `DeltaH`
      the first combined low-|k| family `lowk_mass=4` is now also closed as `strong_drift`
      the completed three-shell broader-family ladder `lowk_mass=2 -> 1 -> 0.5`
      is also still `strong_drift`
      the widened grouped low-|k| basis remains the only serious active path, but the
      `uniform_mass=1` four-shell families `lowk_mass=0.5` and `lowk_mass=0.25`
      have both now completed as `strong_drift`
      the completed light-uniform four-shell lead
      `uniform_mass=0.5`, `lowk_mass=0.25`, `lowk_shells=4`, `nfrog=24`, `dt=0.00004`
      is now formally `strong_drift` at `2/3` repeats, with
      `squareOcc/IPR drift/span max≈0.772` and `repeat span≈0.414`
      widening that grouped low-|k| block to `lowk_shells=6` only made the tune more conservative
      the grouped-`midk` follow-ups are now also formally `strong_drift`
      (`midk_mass=1` and `midk_mass=0.5`, both with `midk_shells=2`)
      the current active lead is now the general shell-map preconditioner:
      - partial short-tune winner:
        `uniform_mass=0.5`,
        `shell_map=0.125,0.125,0.25,0.25,0.5,0.5`,
        `nfrog=20`, `dt=0.00004`
      - active retained-window follow-up: `512 / 256 / 32` with explicit seed families
        `50001,51001,52001`
  - the merged report entry is `data/triangular_hmc_production/overview/report.md`
  - that unified `report.md` now begins with `## Current Representative Stage Per Case`
    before the full historical stage table
  - the merged report now exposes per-case sections plus `bins`, `thermal_cut`, `warm`, and `samples/post`
    so short traces are easier to separate from genuinely deep thermalization runs
  - the merged report also renders per-case relative-change plots and a representative
    one-stage-per-`(beta, dtau)` convergence table, so `beta/dtau` trends can be inspected
    without mixing in old failed tuning attempts
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
- `test/production_hmc.py collect-tune` rebuilds `production_tune*.json/csv` from completed or partially completed tune `runs/`.
- `test/production_hmc.py stage` runs HMC-only production stage jobs and stores sample-index traces.
- `test/production_hmc.py collect` rebuilds `production_stage*.json/csv` from completed `runs/`.
- `test/production_hmc.py report` renders Markdown + PNG summaries from `production_stage.json`, `production_tune.json`, or the legacy `production_benchmark.json`.
- `test/production_hmc.py` can optionally pass `--hmc-mass-spatial-uniform` into the HMC kernel
  through the environment variable `BPQMC_HMC_MASS_SPATIAL_UNIFORM`.
- `test/production_hmc.py` can also pass `--hmc-mass-spatial-shell1`,
  `--hmc-mass-spatial-shell2`, `--hmc-mass-spatial-lowk`, and
  `--hmc-mass-spatial-midk` into the HMC kernel
  through `BPQMC_HMC_MASS_SPATIAL_SHELL1`, `BPQMC_HMC_MASS_SPATIAL_SHELL2`,
  `BPQMC_HMC_MASS_SPATIAL_LOWK`, and `BPQMC_HMC_MASS_SPATIAL_MIDK`.
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
- `test/render_hmc_live_progress.py` writes:
  - `live_progress/live_trace_*.png`
  - `live_progress/report.md`
  - `live_progress/summary.json`
  - its summary table now reports the worst per-run `drift/span` among the currently visible repeats,
    not a smoothed or longest-trace surrogate
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
- `PRODUCTION-RAMP.md` is the shorter persistent action checklist for the current production ramp.
- `HMC-REFERENCE.md` should always state the current healthy rung, next rung, and main blocker so the production ramp can resume after context loss.
- `PRODUCTION-RAMP.md` should be kept in sync whenever the active next rung or blocker changes.
- For production data review, prefer `test/render_hmc_production_overview.py` over opening many per-directory reports by hand.
- The unified `overview/report.md` now also contains a `## Live Progress` section for in-flight large-`U2` stages.
- That unified live section now also uses the worst currently visible repeat, so in-flight
  drift is not understated by averaging across seed families.
- In the merged trace plots, the dashed line is only the configured `thermal_cut`, and the overlaid lines are independent stage runs rather than one continued chain.
- The production stage gate should be interpreted on the worst retained repeat, not only on the mean drift across repeats.
- The production stage gate should also be read against the cross-repeat retained-window mismatch,
  which now appears in the overview as `repeat span`.
- When judging convergence, do not just check whether a deeper rung runs.
  Explicitly inspect whether `squareOcc` and `IPR` stop changing as `beta` increases and as `dtau` decreases.
- `test/production_convergence_scan.py` is the convenience wrapper for those fixed-`beta` and fixed-`dtau` scans.
  Its defaults are intentionally less conservative than the quick bring-up stages:
  `bins=1024`, `thermal_cut=512`, `warm=512`.
- `test/production_hmc.py stage` / `collect` now use the same `1024 / 512 / 512` defaults unless
  overridden explicitly.
- `test/dqmc_production` mirrors that split:
  - `PROD_COMMAND=tune` / `collect-tune` default to `64 / 32 / 32`
  - `PROD_COMMAND=stage` / `collect` default to `1024 / 512 / 512`
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
