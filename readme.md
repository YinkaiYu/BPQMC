# Bosonic Projector QMC

This repository implements a rank-1 bosonic projector QMC solver with two sampling modes:

- local Metropolis updates
- HMC global updates

The code now supports both `kagome` and `triangular` lattices at runtime. The same executable can switch lattice geometry through the first line of `paramC_sets.txt`.

## Quick Start

Build in WSL:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
make -C src
```

Run the default example in `test/`:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC/test
cp ../src/BPQMC.out .
mpirun -np 1 ./BPQMC.out
```

The cleaned `test/` layout is documented in `test/README.md`.
Top-level `test/` is now reserved for active production/HPC workflows; archived
small-benchmark tools live under `test/archive_small_benchmark/`.
The current active workflow is HMC production bring-up, not local-vs-HMC correctness benchmarking.
The persistent production checklist now lives in `PRODUCTION-RAMP.md`, while
`HMC-REFERENCE.md` remains the more detailed running handoff note.
The single merged entry point for accumulated production results is
`data/triangular_hmc_production/overview/report.md`.
That same `report.md` now also carries a `## Current Representative Stage Per Case`
section for the current best retained-window stage at each fixed `(L, Nbos, U2)`,
plus a `## Live Progress` section for in-flight large-`U2` stages.

Run a triangular-lattice HMC smoke test in a fresh temporary directory:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
mkdir -p /tmp/bpqmc_triangular_smoke
cp src/BPQMC.out /tmp/bpqmc_triangular_smoke/
cat > /tmp/bpqmc_triangular_smoke/confin.txt <<'EOF'
0
EOF
cat > /tmp/bpqmc_triangular_smoke/seeds.txt <<'EOF'
12345
22345
32345
42345
52345
62345
72345
82345
EOF
cat > /tmp/bpqmc_triangular_smoke/paramC_sets.txt <<'EOF'
triangular
1.0 0.0 1.0 6
2 2 8 1.6
2 2 8
2 120 4 1.0
.false. 60
.true. 20 1.0 1.0
.true. 10 0.012 0
2 0.1 0.0 0.0
5 0.0001 0.0
EOF
cd /tmp/bpqmc_triangular_smoke
mpirun -np 1 ./BPQMC.out
sed -n '1,40p' info.txt
```

## Runtime Input

`paramC_sets.txt` now accepts an explicit lattice header:

```text
lattice_type
RT          RU1         RU2         Nbos
Nlx         Nly         Ltrot       Beta
NlxTherm    NlyTherm    LtrotTherm
Nwrap       Nbin        Nsweep      shiftLoc
is_tau      Nthermal
is_warm     Nwarm       shiftWarm1  shiftWarm2
is_global   Nfrog       hmc_dt      NfrogJitter hmc_mass hmc_block_tau hmc_block_sites
iniType     iniAmpl     iniBias1    iniBias2
iniHam      iniTwist    imbalance
```

Supported `lattice_type` values:

- `kagome`
- `triangular`

Backward compatibility:

- if the first line is numeric, the code falls back to the old format and assumes `kagome`
- if the HMC line is absent, the code falls back to local updates

Meaning of the HMC line:

- `is_global = .false.` keeps local updates
- `is_global = .true.` enables HMC
- `Nfrog` is the nominal leapfrog step count
- `hmc_dt` is the leapfrog step size
- `NfrogJitter` randomizes trajectory length uniformly in `max(1, Nfrog-NfrogJitter) ... Nfrog+NfrogJitter`
- `hmc_mass` is the leapfrog mass
- `hmc_block_tau` and `hmc_block_sites` are optional block sizes; `0` means full-field HMC

The parser is backward compatible. The following HMC line formats are accepted:

```text
is_global Nfrog hmc_dt
is_global Nfrog hmc_dt NfrogJitter
is_global Nfrog hmc_dt NfrogJitter hmc_mass
is_global Nfrog hmc_dt NfrogJitter hmc_mass hmc_block_tau
is_global Nfrog hmc_dt NfrogJitter hmc_mass hmc_block_tau hmc_block_sites
```

For production triangular runs used in this branch, the usual choice is:

- `U1 = 0`
- `iniHam = 5`
- `iniTwist = 1e-4`
- `beta = 256`
- `Dtau = beta / Ltrot = 0.001`

## Lattice Support

The executable switches lattice geometry at runtime.

`kagome`

- `Norb = 3`
- `Nbond = 2`
- keeps the original kagome geometry and observable layout

`triangular`

- `Norb = 1`
- `Nbond = 3`
- uses the triangular primitive vectors and nearest-neighbor bonds from `feature/triangular-lattice`
- supports `iniHam = 4` and `iniHam = 5`

`iniHam = 5` on the triangular lattice applies a directional Peierls phase on bond direction 1:

```math
Z_1 = t \, e^{i 2\pi \,\text{iniTwist}}, \qquad Z_2 = Z_3 = t.
```

## Build Notes

The build uses a two-layer Makefile setup:

- `src/Makefile` selects compiler and external libraries
- `src/Compile` builds the Fortran objects and links `BPQMC.out`

Useful commands:

```bash
cd src && make
cd src && make clean
cd src && make FFLAGS='-O0 -g -traceback -check all -fpe0 -c -I/home/yyk/Lib_90_new/Modules'
```

## Sampling and Benchmarks

### Active Production Workflow

The active production/HPC driver is now `test/production_hmc.py`:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py --help
```

Key subcommands:

- `tune`: short HMC grid scans
- `collect-tune`: rebuild tune summaries from finished or partially finished `runs/`
- `stage`: long HMC-only runs with thermalization traces
- `collect`: rebuild stage summaries from finished or partially finished `runs/`
- `report`: render Markdown + PNG summaries
- `benchmark`: legacy direct local-vs-HMC comparison, kept only as an archived-style helper
- `render_hmc_live_progress.py`: render Markdown + PNG traces directly from an in-flight
  stage `runs/` tree before a repeat has finished

Optional production preconditioning knob:

- `--hmc-mass-spatial-uniform`
  - this keeps the input file format unchanged and drives an environment-backed HMC preconditioner
  - `hmc_mass` remains the residual-mode mass
  - `hmc-mass-spatial-uniform` sets a lighter or heavier mass for the per-time-slice spatially uniform mode
  - `0` disables the split and recovers the old scalar-mass HMC
- `--hmc-mass-spatial-shell1`
  - triangular-only experimental preconditioner for the lowest nonzero momentum shell
  - `hmc_mass` remains the residual-mode mass
  - `hmc-mass-spatial-uniform` still controls the exact spatially uniform mode
  - `hmc-mass-spatial-shell1` assigns a separate mass to the first real cosine/sine shell
  - `0` disables the shell split and recovers the scalar/uniform-only geometry

Example: short production tune scan on a conservative triangular ladder rung:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py tune \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 1000 \
  --u2-values 1 \
  --beta 32 \
  --dtau 0.01 \
  --bins 64 \
  --thermal-cut 32 \
  --warm 32 \
  --grid '8:0.02,12:0.015,16:0.01' \
  --hmc-jitter 2 \
  --hmc-mass 4 \
  --hmc-mass-spatial-uniform 0 \
  --repeats 2 \
  --work-root data/triangular_hmc_production/l6_n1e3_u1_tune
```

Example: long HMC-only stage run that focuses on `squareOcc` and `IPR` thermalization:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py stage \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 1000 \
  --u2-values 1 \
  --beta 32 \
  --dtau 0.01 \
  --bins 1024 \
  --thermal-cut 512 \
  --warm 512 \
  --repeats 3 \
  --hmc-nfrog 12 \
  --hmc-dt 0.015 \
  --hmc-jitter 2 \
  --hmc-mass 4 \
  --hmc-mass-spatial-uniform 0 \
  --stage-label l6_n1e3_u1_beta32_dtau1em2 \
  --work-root data/triangular_hmc_production/l6_n1e3_u1_stage
```

Rebuild the stage summary without rerunning QMC:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py collect \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 1000 \
  --u2-values 1 \
  --beta 32 \
  --dtau 0.01 \
  --bins 1024 \
  --thermal-cut 512 \
  --warm 512 \
  --repeats 3 \
  --hmc-nfrog 12 \
  --hmc-dt 0.015 \
  --hmc-jitter 2 \
  --hmc-mass 4 \
  --hmc-mass-spatial-uniform 0 \
  --stage-label l6_n1e3_u1_beta32_dtau1em2 \
  --work-root data/triangular_hmc_production/l6_n1e3_u1_stage
```

Rebuild a tune summary from partially completed tune runs:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py collect-tune \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 1000 \
  --u2-values 1 \
  --beta 192 \
  --dtau 0.002 \
  --bins 64 \
  --warm 0 \
  --repeats 2 \
  --grid '8:0.02,12:0.015,16:0.01' \
  --hmc-jitter 2 \
  --hmc-mass 4 \
  --hmc-mass-spatial-uniform 0 \
  --work-root data/triangular_hmc_production/l6_n1e3_u1_beta192_dtau2em3_tune_m4
```

Render the stage report:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py report \
  --work-root data/triangular_hmc_production/l6_n1e3_u1_stage \
  --summary-kind stage
```

Open the notebook template against the same summary root for interactive inspection:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC/test
jupyter notebook hmc_report_template.ipynb
```

The notebook now understands `stage`, `tune`, and `benchmark` summaries and rewrites
the image paths from `report.md` so PNG figures display correctly inside the notebook.

If a long stage is still in progress and `collect` cannot yet summarize it, render a live trace report with:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
MPLBACKEND=Agg /home/yyk/conda/envs/notebook/bin/python test/render_hmc_live_progress.py \
  --work-root data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag2048 \
  --thermal-cut 1024 \
  --output-dir data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag2048/live_progress
```

Example: blocker-focused preconditioned probe on `L=6, Nbos=1e4, U2=1`:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_hmc.py stage \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 10000 \
  --u2-values 1 \
  --beta 32 \
  --dtau 0.01 \
  --bins 512 \
  --thermal-cut 256 \
  --warm 512 \
  --repeats 2 \
  --hmc-nfrog 20 \
  --hmc-dt 0.008 \
  --hmc-jitter 2 \
  --hmc-mass 16 \
  --hmc-mass-spatial-uniform 1 \
  --stage-label l6_n1e4_u1_beta32_dtau1em2_m16_mu1_nf20_dt0p008_diag512 \
  --work-root data/triangular_hmc_production/l6_n1e4_u1_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p008_diag512
```

The `stage` report writes:

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

The stage workflow is intended for the gradual production ladder.
The first judgement is whether `squareOcc` and `IPR` visibly stabilize across sample index, seeds, and initial-state choices.
Acceptance and `ESS/sec` are secondary diagnostics.
If you want one command to launch several initial-state families, use
`--ini-type-values` and/or `--ini-ampl-values`; the driver will emit separate case names
for each family so convergence across seeds and initial states can be judged from one report.

Current local workstation examples under `data/triangular_hmc_production/` already include:

- a healthy baseline projector ladder on the same physics point:
  - `l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015`
  - `l6_n1e3_u1_beta64_dtau1em2_stage_m4_nf8_dt0p02`
  - `l6_n1e3_u1_beta96_dtau8em3_stage_m4_nf8_dt0p02`
  - `l6_n1e3_u1_beta128_dtau5em3_stage_m4_nf8_dt0p02`
  - `l6_n1e3_u1_beta160_dtau4em3_stage_m4_nf8_dt0p02`
- a stage-healthy preconditioned replacement for the old `Nbos=1e4, U2=1` blocker:
  - `l6_n1e4_u1_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p008_diag512`
- a stronger-coupling exploratory rung that is already stage-healthy with a smaller step size:
  - `l6_n1e4_u1e1_beta32_dtau1em2_stage_m4_nf16_dt0p006`
- a fully healthy higher-`U2` preconditioned rung:
  - `l6_n1e4_u3e1_beta32_dtau1em2_stage_m16_mu1_nf12_dt0p006_diag1024`
- current higher-`U2` preconditioned tuning / stage roots:
  - `l6_n1e4_u3e1_beta32_dtau1em2_tune_m16_mu1`
  - `l6_n1e4_u1e2_beta32_dtau1em2_tune_m16_mu1_conservative`
  - `l6_n1e4_u1e2_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p001_diag1024`
  - `l6_n1e4_u3e2_beta32_dtau1em2_tune_m16_mu1_conservative`
  - `l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag1024`

On the current representative baseline ladder, `squareOcc` and `IPR` only move by about
`7.6e-4` and `7.7e-4` across `beta=32 -> 192`, `dtau=0.01 -> 0.002`.
That is why the current workstation priority is now the higher-`U2` ramp rather than
deeper projector scans on the easy baseline point.

The current next production rungs are:

- `L=6, Nbos=1e4, U2=100`, `mass=16`, `uniform_mass=1`, `nfrog=20`, `dt=0.001`
- `L=6, Nbos=1e4, U2=300`, `mass=16`, `uniform_mass=1`, `nfrog=28`, `dt=0.0005`
- `L=6, Nbos=1e4, U2=1000` is the current workstation blocker:
  - the old scalar/uniform-only reference geometry `mass=16`, `uniform_mass=1`, `shell1_mass=0`, `nfrog=20`, `dt=0.0004`
    is now a completed `strong_drift` stage
  - `shell1_mass=1` has already been ruled out as too aggressive
  - shell1-only `shell1_mass=4` and `shell1_mass=8` exploratory stages have now both finished as `strong_drift`
  - the first `shell2` family `shell1_mass=8`, `shell2_mass=4` has also now failed on completed stages
  - the lighter `shell2` follow-up `shell1_mass=4`, `shell2_mass=2` also failed to produce
    a selective tune window; all cheap candidates stayed at `acceptance = 1` with huge positive `DeltaH`
  - the active replacement path is now the broader combined low-|k| split:
    - `mass=16`, `uniform_mass=1`, `lowk_mass=4`
    - short-tune grid:
      - `12 x 0.0003`
      - `16 x 0.00025`
      - `20 x 0.0002`
      - `24 x 0.00015`
    - current recommended candidate: `16 x 0.00025`
    - current backup candidate: `20 x 0.0002`
    - active live traces:
      - `l6_n1e4_u1e3_beta32_dtau1em2_stage_m16_mu1_lowk4_nf16_dt2p5em4_diag512_warm0`
      - `l6_n1e4_u1e3_beta32_dtau1em2_stage_m16_mu1_lowk4_nf20_dt2em4_diag512_warm0`

See `HMC-REFERENCE.md` for the running interpretation of these stage results and blockers.
For a single merged entry point across all production directories, open:

- `data/triangular_hmc_production/overview/report.md`

It is generated by:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/render_hmc_production_overview.py \
  --root data/triangular_hmc_production \
  --output-dir data/triangular_hmc_production/overview
```

This merged report is the preferred place to inspect:

- `squareOcc` / `IPR` trends versus `beta` and `dtau`
- recommended tune parameters versus rung
- sample-index traces across different parameter points
- per-case sections that keep different `Nbos` / `U2` scales separate
- per-case relative-change plots that make small `beta/dtau` shifts visible even when the
  absolute observable scales differ strongly across physics points
- `bins`, `thermal_cut`, `warm`, and `samples/post` coverage, so a short trace cannot
  be mistaken for a deep thermalization run
- note that the dashed marker in each trace is only the configured `thermal_cut`
  and the overlaid lines are independent stage runs, not one continued chain
- note that the production stage gate is now judged on the worst retained repeat as well as the mean drift
- note that the merged stage tables now also expose `repeats`, `max drift`, and `repeat span`
  so cross-repeat plateau mismatch is visible from the overview page

For explicit convergence scans at fixed `beta` or fixed `dtau`, use:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_convergence_scan.py \
  --mode beta \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 1000 \
  --u2-values 1 \
  --dtau 0.01 \
  --beta-values 32,64,96,128
```

or

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/production_convergence_scan.py \
  --mode dtau \
  --lattice-type triangular \
  --l-values 6 \
  --nbos-values 1000 \
  --u2-values 1 \
  --beta 128 \
  --dtau-values 0.01,0.008,0.005,0.004
```

The convergence wrapper now defaults to a less conservative sample count:

- `bins=1024`
- `thermal_cut=512`
- `warm=512`

The canonical `test/production_hmc.py stage` / `collect` workflow now uses the same
`1024 / 512 / 512` defaults unless overridden explicitly.
The SLURM helper `test/dqmc_production` mirrors this split:

- `PROD_COMMAND=tune` / `collect-tune` default to `64 / 32 / 32`
- `PROD_COMMAND=stage` / `collect` default to `1024 / 512 / 512`

### Archived Development Tune Scan

The old development-only tune helper is still available as an archived reference:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/archive_small_benchmark/tune_hmc.py --set triangular_weak_u2 --grid '8:0.01,10:0.012,12:0.015' --bins 80 --sweeps 4 --warm 20
```

### Archived Development Benchmark

The old development benchmark helper is also archived:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/archive_small_benchmark/benchmark_hmc.py --sets triangular_weak_u2 --repeats 2
```

The benchmark reports:

- bin-averaged observables after thermal cut
- combined-error `z` scores
- HMC acceptance
- `tau_int(doubleOcc)`
- `ESS/sec`

### Small Triangular Correctness Workflow

The repository has a dedicated small-parameter triangular correctness campaign:

- `L = 6`
- `beta = 32`
- `dtau = 0.01`
- `U1 = 0`
- `Nbos = 10, 100, 1000`
- `U2 = 1e-2, 1e-1, 1e0, 1e1, 1e2`

This benchmark campaign is now treated as closed for correctness validation.
The active entry point is the final rendered report under `data/triangular_hmc_small_benchmark/full_grid_progress_v3/`.
Raw benchmark/tuning directories were moved under `data/triangular_hmc_small_benchmark/archive_20260319/` so future production work is not buried under old development runs.

Use this section as an archived reference, not as the active development target.

Archived benchmark tooling now lives under:

- `test/archive_small_benchmark/`
- `test/README.md` summarizes the current active-vs-archived `test/` layout

The final correctness report is:

- `data/triangular_hmc_small_benchmark/full_grid_progress_v3/report.md`

To re-render the full 3x5 visual report from the archived benchmark directories:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/archive_small_benchmark/render_small_hmc_full_grid.py \
  --python /home/yyk/conda/envs/notebook/bin/python \
  --output-dir data/triangular_hmc_small_benchmark/full_grid_progress_v3
```

This full-grid report intentionally mixes strict benchmark directories and single-repeat visual probes.
The renderer labels them separately:

- `PASS` / `FAIL` only apply to cases with at least two repeats
- `PROBE` means the point is present for visual comparison, but its `z-score` is not treated as a strict statistic
- for `PROBE` cases, the report emphasizes `|diff| / span`, absolute differences, and sample-index traces instead of repeat-based `z-score`

If an archived long strict benchmark already finished its `runs/` directories but the summary
`json/csv` files were not written, rebuild them without rerunning:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/archive_small_benchmark/small_hmc_benchmark.py collect-benchmark \
  --work-root data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n100_u1e1_v3_m4 \
  --nbos-values 100 \
  --u2-values 1e1 \
  --bins 640 \
  --thermal-cut 384 \
  --warm 0 \
  --repeats 4 \
  --hmc-nfrog 32 \
  --hmc-dt 0.04 \
  --hmc-jitter 2 \
  --hmc-mass 4
```

This writes:

- `stage_cases.csv`
- `stage_observables.csv`
- `stage_samples.csv`
- `stage_summary.json`
- `overview.png`
- `observable_*_n*.png`
- `trace_*.png`
- `tune_*.png`
- `report.md`

For notebook-based review of the rendered report, open:

- `test/archive_small_benchmark/small_hmc_stage_report.ipynb`

The report renderer and notebook expect `pandas` and `matplotlib`; in this repository the
recommended interpreter is `/home/yyk/conda/envs/notebook/bin/python`.
The notebook reads the CSV and PNG files from a rendered report directory. Set `REPORT_DIR` inside the notebook before running its cells.
`stage_samples.csv` now stores the key observable traces for every repeat, not only `repeat=0`, so you can do repeated thermal-cut checks directly inside the notebook or in a post-processing script.
The notebook exposes both `NBOS_SELECT` and `TRACE_REPEAT`, which lets you inspect one particle-number slice and one benchmark repeat without regenerating the report.
The archived staged checkpoint with 9 strict-passing health points is `data/triangular_hmc_small_benchmark/stage_report_v6a`.
The fixed full-grid visual report is rendered with `test/archive_small_benchmark/render_small_hmc_full_grid.py`; its default output target is `data/triangular_hmc_small_benchmark/full_grid_progress_v3`.
The current full-grid rendered report is `data/triangular_hmc_small_benchmark/full_grid_progress_v3`.
The archived raw benchmark campaign lives under `data/triangular_hmc_small_benchmark/archive_20260319`.

Practical conclusion from this benchmark campaign:

- the HMC action/force implementation is considered validated on the small-parameter triangular testbed
- the final correctness reference is the rendered full-grid report, not the individual raw run directories
- for large `Nbos` and `U2`, local and HMC can both remain unthermalized on workstation-length runs; those points are not useful correctness references
- future work in this branch should focus on production-grade HMC thermalization, preconditioning, and HPC workflows rather than further local-vs-HMC benchmarking
For the running handoff note covering validated health points, active production rungs,
unresolved slow modes, tuning trends, and HPC-oriented lessons, see `HMC-REFERENCE.md`.
For the cleaned `test/` directory layout and archive locations, see `test/README.md`.
For unresolved slow-mode diagnostics, the same renderer can also be used on non-passing benchmark
directories. A current example is `data/triangular_hmc_small_benchmark/progress_report_n1000_u1e0_v2`,
which shows per-repeat `sample index` traces including `SF_Gamma`.
Another current example is `data/triangular_hmc_small_benchmark/progress_report_n100_u1e1_v3`,
which tracks the strong-coupling `Nbos=100, U2=10` point across multiple benchmark variants.

### Production Tune Workflow

The new production driver is `test/production_hmc.py`.

Tune HMC on a triangular production grid:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/production_hmc.py tune \
  --lattice-type triangular \
  --l-values 12 \
  --nbos-values 100000,1000000,10000000 \
  --u2-values 100,1000,10000 \
  --beta 256 \
  --dtau 0.001 \
  --nwrap 32 \
  --bins 64 \
  --sweeps 1 \
  --thermal-cut 32 \
  --warm 32 \
  --grid '8:0.002,10:0.003,12:0.004,16:0.005' \
  --work-root /tmp/bpqmc_prod_tri_L12
```

This writes:

- `production_tune.json`
- `production_tune.csv`
- `recommended_hmc.json`

### Production Benchmark Workflow

Benchmark local vs HMC on the same parameter grid:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/production_hmc.py benchmark \
  --lattice-type triangular \
  --l-values 12 \
  --nbos-values 100000,1000000,10000000 \
  --u2-values 100,1000,10000 \
  --beta 256 \
  --dtau 0.001 \
  --nwrap 32 \
  --bins 64 \
  --sweeps 1 \
  --thermal-cut 32 \
  --warm 32 \
  --repeats 2 \
  --hmc-json /tmp/bpqmc_prod_tri_L12/recommended_hmc.json \
  --work-root /tmp/bpqmc_prod_tri_L12
```

This writes:

- `production_benchmark.json`
- `production_benchmark_cases.csv`
- `production_benchmark_observables.csv`

### Render Figures and Report

Create PNG figures and a Markdown report:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/render_hmc_report.py /tmp/bpqmc_prod_tri_L12/production_benchmark.json
```

This creates:

- `acceptance.png`
- `ess_per_sec.png`
- `tau_int.png`
- `zscore_heatmap.png`
- `report.md`

For notebook-based post-processing, start from:

- `test/hmc_report_template.ipynb`

## Physics-to-Code Mapping

The code stores one rank-1 bosonic orbital but represents a two-flavor weight through a modulus square. The main mapping is:

- `Nbos` ↔ `N_b`
- `Beta` ↔ `2\theta`
- `Prop%UUR(:,1)` ↔ `P_R(\tau)=B(\tau,0)P`
- `Prop%UUL(1,:)` ↔ `P_L^\dagger(\tau)=P^\dagger B(2\theta,\tau)`
- `Prop%overlap` ↔ `P_L^\dagger(\tau) P_R(\tau)`
- `Conf%phi_list(nf,ii,nt)` stores both auxiliary-field flavors
- `OperatorHubbard%alpha` stores the Hubbard-Stratonovich coefficient for each flavor

The local-update ratio is implemented through:

```fortran
ratio_abs = abs(ratio_exp * ratio_Pfa * dconjg(ratio_Pfa))
```

The HMC effective action is

```math
H[\phi,\pi] = K[\pi] + S[\phi], \qquad
K[\pi] = \frac12 \sum_{\tau,i,n_f}\pi_{n_f,i}^2(\tau),
```

```math
S[\phi] = \frac12 \sum_{\tau,i,n_f}\phi_{n_f,i}^2(\tau)
  - N_b \ln \left| P^\dagger B(2\theta,0) P \right|^2.
```

The force uses the same sign convention as the local Metropolis ratio:

```math
F_{n_f,i}(\tau) =
  -\phi_{n_f,i}(\tau) + 2\,\mathrm{Re}\!\left[\alpha_{n_f}\,\bar{G}_{ii}^{(n_f)}(\tau)\right].
```

This README uses `\mathrm{Re}` for the real-part operator so GitHub math rendering works.

## Rank-1 Formulation

All bosons occupy one orbital:

```math
|\Phi_T\rangle = \frac{1}{\sqrt{N_b!}} (a^\dagger P)^{N_b} |0\rangle.
```

The Gaussian propagator acts directly on the orbital vector:

```math
e^{-a^\dagger T a} (a^\dagger P)^{N_b} |0\rangle
= (a^\dagger e^{-T} P)^{N_b} |0\rangle.
```

This is why the Monte Carlo loop stores only:

- the right vector
- the left vector
- their overlap
- wrap/log-normalization data

Equal-time Green's functions are reconstructed on demand from the normalized vectors:

```math
\bar{G}_{ij} = N_b \frac{[B(\theta,0)P]_i [P^\dagger B(2\theta,\theta)]_j}{P^\dagger B(2\theta,0)P}.
```

## Numerical Stabilization

The code periodically rescales left and right propagated vectors and accumulates the discarded logarithmic norms. HMC uses these log norms to evaluate

```math
\ln \left| P^\dagger B(2\theta,0) P \right|^2
= 2\ln Z_R + 2\ln Z_L + \ln |P_L^\dagger P_R|^2.
```

This is required for stable long trajectories and large `beta`.

## Current Status of This Branch

This branch now includes:

- runtime `kagome/triangular` switching
- HMC wired into the main flow
- triangular `iniHam = 5`
- production tuning and benchmark scripts
- report rendering and notebook template
- reduced HMC memory footprint by removing one full-size force buffer

The production-scale `L=12` and `L=21` runs with `beta=256`, `Dtau=0.001`, `Nbos=1e5..1e7`, `U2=1e2..1e4` are intended to be driven by `test/production_hmc.py`.
