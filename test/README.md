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
- `production_hmc.py stage`: long HMC-only runs with sample-index trace capture for `squareOcc`, `IPR`, `doubleOcc`, `nearestOcc`
- `production_hmc.py collect`: rebuild `production_stage*.json/csv` from finished `runs/` without rerunning QMC
- `production_hmc.py report`: render a Markdown + PNG report from `production_stage.json`, `production_tune.json`, or `production_benchmark.json`

The default SLURM entry point is `dqmc_production`, which now dispatches by `PROD_COMMAND`:

- `PROD_COMMAND=tune`
- `PROD_COMMAND=stage`
- `PROD_COMMAND=collect`
- `PROD_COMMAND=report`
- `PROD_COMMAND=benchmark` is kept only for legacy direct local-vs-HMC comparisons

For the active workstation ramp, keep the current healthy rung, next rung, and main blocker
documented in `../HMC-REFERENCE.md` before moving on to a harder stage.
At the moment, the baseline `L=6, Nbos=1e3, U2=1` projector ladder has already been pushed
through `beta=128, dtau=0.005`, and the next rung is `beta=160, dtau=0.004`.
