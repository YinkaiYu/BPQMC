# Test Directory Layout

This directory is now organized around the future production/HPC workflow.
The short persistent plan for that workflow lives in `../PRODUCTION-RAMP.md`;
the more detailed handoff note lives in `../HMC-REFERENCE.md`.

## Active Top-Level Files

These are the files intended for ongoing use:

- `BPQMC.out`: current runtime binary copied into `test/`
- `auto.sh`: local build-copy helper
- `dqmc`: default SLURM job script
- `dqmc_production`: production-oriented SLURM job script
- `confin.txt`, `paramC_sets.txt`, `seeds.txt`: canonical runtime inputs
- `hmc_tools.py`: shared Python helpers
- `production_hmc.py`: main production tuning / stage / collect / report driver
- `render_hmc_report.py`: production report renderer for `tune`, `stage`, and legacy `benchmark` summaries
- `render_hmc_production_overview.py`: merged renderer for the whole `data/triangular_hmc_production/` tree
- `render_hmc_live_progress.py`: live renderer for in-flight stage runs whose current repeat has not yet completed
  - writes `live_progress/summary.json` in addition to Markdown/PNG so the unified overview
    can reuse the exact same `thermal_cut` and worst-repeat live drift summary
  - active stage runs also write `stage_config.json` in the work root, so the live renderer
    and the unified overview can recover the configured `thermal_cut` even before a full
    `production_stage.json` exists
- `production_convergence_scan.py`: convenience wrapper for fixed-`beta` and fixed-`dtau` stage scans
  - defaults to `bins=1024`, `thermal_cut=512`, `warm=512` for less conservative convergence checks
- `production_hmc.py stage` and `production_hmc.py collect` now use the same `1024 / 512 / 512`
  defaults unless the command line overrides them
- `hmc_report_template.ipynb`: production report notebook template for `stage`, `tune`, and `benchmark` summaries
- `analyze_hmc_monitor.py`, `analyze_hmc_trace.py`: HMC trajectory diagnostics that may still be useful in production debugging

## Archived Benchmark Tools

The closed small-parameter local-vs-HMC correctness campaign has been moved to:

- `archive_small_benchmark/`

This archive contains the old benchmark/tuning/render scripts and notebooks, including:

- `small_hmc_benchmark.py`
- `render_small_hmc_stage_report.py`
- `render_small_hmc_full_grid.py`
- `small_hmc_stage_report.ipynb`
- `benchmark_hmc.py`
- `tune_hmc.py`

These tools are still runnable, but they are no longer the main workflow.

## Archived Manual Checks

Older Fortran-side validation programs and their outputs now live in:

- `archive_manual_checks/`

This includes the finite-difference / round-trip / comparison helpers used during HMC bring-up.

## Archived Runtime Outputs

Old top-level run products were moved to:

- `runtime_outputs/root_run_20260319/`

This keeps `test/` clean for future production runs while preserving the last root-level manual run dump.

## Production Focus

The benchmark phase is considered closed.
Future work should prioritize:

- staged thermalization bring-up from easier local workstation points toward production physics
- production-grade HMC thermalization
- preconditioning and integrator improvements
- SLURM/HPC workflows
- large-`L`, large-`Nbos`, large-`U2` HMC runs

The active production workflow is now:

- `production_hmc.py tune`: short HMC grid scans for candidate parameters
  - `--run-timeout-sec` is available for tune/stage/benchmark execution paths when `mpirun`
    occasionally hangs after the observable files and `info.txt` are already complete
- `production_hmc.py collect-tune`: rebuild `production_tune*.json/csv` from finished or partially finished tune `runs/`
- `production_hmc.py stage`: long HMC-only runs with sample-index trace capture for `squareOcc`, `IPR`, `doubleOcc`, `nearestOcc`
- `production_hmc.py collect`: rebuild `production_stage*.json/csv` from finished or partially finished `runs/` without rerunning QMC
- `production_hmc.py report`: render a Markdown + PNG report from `production_stage.json`, `production_tune.json`, or `production_benchmark.json`
- `production_hmc.py` also supports `--hmc-mass-spatial-uniform`
  - `0` keeps the old scalar-mass HMC
  - positive values split the per-time-slice spatially uniform mode away from the residual `--hmc-mass`
- `production_hmc.py` also supports `--hmc-mass-spatial-lowk`
  and `--hmc-mass-spatial-lowk-shells`
  - these create a grouped low-`|k|` block spanning the first few nonzero triangular momentum shells
- `production_hmc.py` also supports `--hmc-mass-spatial-midk`
  and `--hmc-mass-spatial-midk-shells`
  - these assign the next grouped momentum-shell band to a separate intermediate mass
- `production_hmc.py` also supports `--hmc-mass-spatial-shell-map`
  - pass a comma-separated per-shell mass profile for the first few triangular nonzero
    momentum shells
  - this explicit shell-map takes precedence over grouped `lowk` / `midk` bands
- `production_hmc.py` also supports `--hmc-hybrid-local-sweeps`
  - use this when a static HMC mass profile still shows strong retained-window drift
  - it inserts a small number of local Metropolis sweeps between HMC proposals
- `production_hmc.py` also supports `--hmc-mass-spatial-shell1`
  - `0` keeps the scalar/uniform-only geometry
  - positive values split the triangular lowest nonzero momentum shell away from the residual
    and uniform modes

For thermalization checks across different initial states, `production_hmc.py` also supports:

- `--ini-type-values`
- `--ini-ampl-values`

These options expand multiple initial-state families into distinct case names under the same
`tune` or `stage` work root.

For a single merged entry point over all production directories, run:

- `/home/yyk/conda/envs/notebook/bin/python test/render_hmc_production_overview.py --root data/triangular_hmc_production --output-dir data/triangular_hmc_production/overview`

The resulting report lives at:

- `data/triangular_hmc_production/overview/report.md`
- that same `report.md` now starts with `## Current Representative Stage Per Case`
  before the full historical stage table
- the in-flight large-`U2` live links are now also folded into that same `report.md` under `## Live Progress`
- that live section now reports the worst per-run `drift/span` currently visible, so seed-family
  differences are not hidden by averaging
- the report now includes per-case sections and `samples/post` coverage, so late slow-mode
  drift is easier to diagnose from one entry point
- the dashed marker in each trace is only the configured `thermal_cut`, and the overlaid
  lines are independent stage runs rather than one continued chain
- the production stage gate is judged on the worst retained repeat as well as the mean drift
- the merged stage tables now also expose `repeats`, `max drift`, and `repeat span`

The default SLURM entry point is `dqmc_production`, which now dispatches by `PROD_COMMAND`:

- `PROD_COMMAND=tune`
- `PROD_COMMAND=collect-tune`
- `PROD_COMMAND=stage`
- `PROD_COMMAND=collect`
- `PROD_COMMAND=report`
- `PROD_COMMAND=benchmark` is kept only for legacy direct local-vs-HMC comparisons
- by default it uses `64 / 32 / 32` for `tune` / `collect-tune`
- by default it uses `1024 / 512 / 512` for `stage` / `collect`

When a stage repeat is still running and `collect` cannot yet summarize it, render a live report with:

- `MPLBACKEND=Agg /home/yyk/conda/envs/notebook/bin/python test/render_hmc_live_progress.py --work-root <stage-root>`

For the active workstation ramp, keep the current healthy rung, next rung, and main blocker
documented in `../HMC-REFERENCE.md`, and keep the actionable ladder in `../PRODUCTION-RAMP.md`,
before moving on to a harder stage.
At the moment, the baseline `L=6, Nbos=1e3, U2=1` projector ladder has already been pushed
through `beta=160, dtau=0.004`, and the representative `beta=192, dtau=0.002` stage is already
showing `stable_window` on `2/3` completed repeats. The current workstation priority is therefore
the higher-`U2` ramp, not still deeper projector scans on the same easy baseline point.
The current `U2`-ramp starting geometry is the stage-healthy spatial-uniform preconditioned rung
`L=6, Nbos=1e4, U2=1` with `mass=16`, `uniform_mass=1`, `20 x 0.008`.
The immediate next production rungs are:

- `L=6, Nbos=1e4, U2=30`, `mass=16`, `uniform_mass=1`, `12 x 0.006`, now stage-healthy
- `L=6, Nbos=1e4, U2=100`, `mass=16`, `uniform_mass=1`, `20 x 0.001`, current long-stage `rep0` is healthy
- `L=6, Nbos=1e4, U2=300`, `mass=16`, `uniform_mass=1`, `28 x 0.0005`, stage launched after the conservative retune
- `L=6, Nbos=1e4, U2=1000` is the current blocker:
  the old scalar/uniform-only `20 x 0.0004` reference stage is now `strong_drift`,
  shell1-only `mk=4` and `mk=8` completed stages are also now `strong_drift`,
  the first `shell2` family `mk1=8`, `mk2=4` is also now `strong_drift`,
  the lighter `shell2` follow-up `mk1=4`, `mk2=2` also never produced a selective tune window,
  the first combined low-|k| family `m_lowk=4` is also now `strong_drift`,
  the completed three-shell broader-family ladder `m_lowk=2 -> 1 -> 0.5` is also still `strong_drift`,
  grouped-`midk` follow-ups are now also formally `strong_drift`,
  and the current active replacement path is the explicit shell-map preconditioner
  `shell_map=0.125,0.125,0.25,0.25,0.5,0.5` plus `jitter=8`.
  The current retained-window lead is now `hybrid_local_sweeps=2`,
  `nfrog=24`, `dt=0.000035`; the completed `512/256/32` stage is still `strong_drift`,
  but the deeper `1024/512/512` rerun is now the best live candidate,
  with post-cut `squareOcc/IPR drift/span≈0.185`.
  The current short-tune lead is also `hybrid_local_sweeps=2`,
  `nfrog=16`, `dt=0.00005`, with partial
  `acceptance≈0.979`, `tau_int≈10.731`, `ESS/sec≈0.085`.
  `hybrid_local_sweeps=1` is still retained as a comparison branch.
