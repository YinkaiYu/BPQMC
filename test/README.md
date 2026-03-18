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
- `production_hmc.py`: main production tuning / benchmark driver
- `render_hmc_report.py`: production report renderer
- `hmc_report_template.ipynb`: production report notebook template
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

- production-grade HMC thermalization
- preconditioning and integrator improvements
- SLURM/HPC workflows
- large-`L`, large-`Nbos`, large-`U2` HMC runs
