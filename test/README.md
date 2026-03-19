# Test Directory Layout

This directory is now organized around the future production/HPC workflow.

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
- `production_hmc.py collect-tune`: rebuild `production_tune*.json/csv` from finished or partially finished tune `runs/`
- `production_hmc.py stage`: long HMC-only runs with sample-index trace capture for `squareOcc`, `IPR`, `doubleOcc`, `nearestOcc`
- `production_hmc.py collect`: rebuild `production_stage*.json/csv` from finished `runs/` without rerunning QMC
- `production_hmc.py report`: render a Markdown + PNG report from `production_stage.json`, `production_tune.json`, or `production_benchmark.json`

For thermalization checks across different initial states, `production_hmc.py` also supports:

- `--ini-type-values`
- `--ini-ampl-values`

These options expand multiple initial-state families into distinct case names under the same
`tune` or `stage` work root.

For a single merged entry point over all production directories, run:

- `/home/yyk/conda/envs/notebook/bin/python test/render_hmc_production_overview.py --root data/triangular_hmc_production --output-dir data/triangular_hmc_production/overview`

The resulting report lives at:

- `data/triangular_hmc_production/overview/report.md`
- the report now includes per-case sections and `samples/post` coverage, so late slow-mode
  drift is easier to diagnose from one entry point
- the dashed marker in each trace is only the configured `thermal_cut`, and the overlaid
  lines are independent stage runs rather than one continued chain
- the production stage gate is judged on the worst retained repeat as well as the mean drift

The default SLURM entry point is `dqmc_production`, which now dispatches by `PROD_COMMAND`:

- `PROD_COMMAND=tune`
- `PROD_COMMAND=collect-tune`
- `PROD_COMMAND=stage`
- `PROD_COMMAND=collect`
- `PROD_COMMAND=report`
- `PROD_COMMAND=benchmark` is kept only for legacy direct local-vs-HMC comparisons
- by default it uses `64 / 32 / 32` for `tune` / `collect-tune`
- by default it uses `1024 / 512 / 512` for `stage` / `collect`

For the active workstation ramp, keep the current healthy rung, next rung, and main blocker
documented in `../HMC-REFERENCE.md` before moving on to a harder stage.
At the moment, the baseline `L=6, Nbos=1e3, U2=1` projector ladder has already been pushed
through `beta=160, dtau=0.004`, and the next rung is `beta=192, dtau=0.002`.
