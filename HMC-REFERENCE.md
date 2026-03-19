# HMC Reference Notes

This document is the working handoff note for the triangular-lattice HMC campaign.
It records what has already been validated, what is still unresolved, and which tuning
trends look robust enough to reuse when moving toward larger lattices and HPC runs.
For the current `test/` directory layout, see `test/README.md`; archived small-benchmark
helpers referenced here now live under `test/archive_small_benchmark/`.

## Current Priority

The small-parameter local-vs-HMC correctness campaign is closed.
The active task is now production-grade HMC bring-up:

- move gradually from easy workstation points toward the production physics line
- judge stages primarily by `squareOcc` and `IPR` thermalization
- treat acceptance and `ESS/sec` as secondary diagnostics
- keep every stage resumable and reportable for later HPC handoff

The active production CLI is `test/production_hmc.py`:

- `tune`: short candidate scans
- `stage`: long HMC-only runs with `squareOcc`/`IPR` traces
- `collect`: rebuild stage summaries from completed `runs/`
- `report`: render Markdown + PNG summaries

The default output root for this new workflow should be under:

- `data/triangular_hmc_production/`

## Current Stage / Next Stage / Blocker

Current healthy ladder:

- `L=6`, `Nbos=1e3`, `U2=1`
- `beta=32`, `dtau=0.01`
- `beta=64`, `dtau=0.01`
- `beta=96`, `dtau=0.008`
- `beta=128`, `dtau=0.005`

Current next rung:

- `L=6`, `Nbos=1e3`, `U2=1`, `beta=160`, `dtau=0.004`
- start from scalar-mass full-field HMC with `mass=4`, jittered trajectories, and a tune scan before the long stage

Current main blocker:

- `L=6`, `Nbos=1e4`, `U2=1`, `beta=32`, `dtau=0.01`
- deeper thermal cuts alone do not remove the residual `squareOcc/IPR` slow drift
- longer trajectories with the same scalar mass also did not beat the `16 x 0.01` baseline
- this rung is the first serious candidate for stronger preconditioning rather than more blind scalar-mass scans

## Current Production Checkpoint

As of the current workstation bring-up round, the first production-style stage data under
`data/triangular_hmc_production/` already separates four regimes clearly.

### Healthy projector ladder on the baseline physics point

- [l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015/report/report.md)
- [l6_n1e3_u1_beta64_dtau1em2_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta64_dtau1em2_stage_m4_nf8_dt0p02/report/report.md)
- [l6_n1e3_u1_beta96_dtau8em3_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta96_dtau8em3_stage_m4_nf8_dt0p02/report/report.md)
- [l6_n1e3_u1_beta128_dtau5em3_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta128_dtau5em3_stage_m4_nf8_dt0p02/report/report.md)
- parameters:
  - `L=6`
  - `Nbos=1e3`
  - `U2=1`
  - scalar mass `m=4`
  - jitter `=2`
- observed stage status:
  - `beta=32`, `dtau=0.01`: `stable_window`, `squareOcc/IPR drift/span ~ 0.141`
  - `beta=64`, `dtau=0.01`: `stable_window`, `squareOcc/IPR drift/span ~ 0.030`
  - `beta=96`, `dtau=0.008`: `stable_window`, `squareOcc/IPR drift/span ~ 0.153`
  - `beta=128`, `dtau=0.005`: `stable_window`, `squareOcc/IPR drift/span ~ 0.068`

Interpretation:

- The baseline triangular rung is now healthy through three projector settings.
- The baseline triangular rung is now healthy through four projector settings.
- This is the first real evidence that the production ramp can proceed in projector depth
  rather than only at one easy point.
- The `beta=96` rung is slower in wall-clock time than `beta=64`, but its traces still pass
  the current `squareOcc/IPR` gate cleanly.
- The `beta=128` rung continues the same trend: wall-clock time rises again, but the retained
  `squareOcc/IPR` window is still convincingly stable.

### Healthy reference rung

- [l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015/report/report.md)
- parameters:
  - `L=6`
  - `Nbos=1e3`
  - `U2=1`
  - `beta=32`
  - `dtau=0.01`
  - `mass=4`
  - `nfrog=12`
  - `dt=0.015`
  - `jitter=2`
- observed stage status:
  - `stable_window`
  - `squareOcc drift/span ~ 0.141`
  - `IPR drift/span ~ 0.141`

Interpretation:

- This is the first production-style HMC rung that already looks thermalized by the new
  `squareOcc/IPR` gate.
- It is the current baseline checkpoint for future promotions.

### Slow-drift rung

- [l6_n1e4_u1_beta32_dtau1em2_stage_m4_nf16_dt0p01](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1_beta32_dtau1em2_stage_m4_nf16_dt0p01/report/report.md)
- parameters:
  - `L=6`
  - `Nbos=1e4`
  - `U2=1`
  - `beta=32`
  - `dtau=0.01`
  - `mass=4`
  - `nfrog=16`
  - `dt=0.01`
  - `jitter=2`
- observed stage status:
  - `slow_drift`
  - `squareOcc drift/span ~ 0.331`
  - `IPR drift/span ~ 0.331`

Interpretation:

- The chain is no longer stuck, but this rung has not fully stabilized.
- This is currently the clearest example of a point that likely needs either:
  - longer thermal cut / deeper stage runs
  - better proposal geometry or preconditioning

### Stronger-coupling exploratory rung

- [l6_n1e4_u1e1_beta32_dtau1em2_stage_m4_nf16_dt0p006](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e1_beta32_dtau1em2_stage_m4_nf16_dt0p006/report/report.md)
- parameters:
  - `L=6`
  - `Nbos=1e4`
  - `U2=10`
  - `beta=32`
  - `dtau=0.01`
  - `mass=4`
  - `nfrog=16`
  - `dt=0.006`
  - `jitter=2`
- observed stage status:
  - `stable_window`
  - `squareOcc drift/span ~ 0.218`
  - `IPR drift/span ~ 0.218`
  - repeat acceptance roughly `0.63` to `0.77`

Interpretation:

- Stronger coupling does not necessarily fail first.
- With a sufficiently smaller step size, this point already passes the current
  `squareOcc/IPR` stage gate even though acceptance is below the old ideal band.
- This is direct evidence that, in the current production phase, acceptance should remain
  secondary to thermalization.

### Immediate tuning lesson from the `Nbos=1e4, U2=10` scan

- [l6_n1e4_u1e1_beta32_dtau1em2_tune_m4](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e1_beta32_dtau1em2_tune_m4/report/report.md)

Relevant outcomes:

- `8 x 0.01` was unusable:
  - acceptance collapsed to zero
  - `DeltaH` blew up by many orders of magnitude
- `12 x 0.008` was mobile but still rough:
  - acceptance around `0.53`
  - large `DeltaH`
- `16 x 0.006` was the best short-scan compromise:
  - acceptance around `0.67`
  - much smaller `tau_int`
  - best `ESS/sec`

Interpretation:

- On this stronger rung, the stability wall is controlled by the step size much more than
  by the nominal trajectory length.
- The practical recipe was:
  - keep `mass=4`
  - reduce `dt`
  - let acceptance fall below `0.70` if the traces improve

## Scope

- Lattice: `triangular`
- Development benchmark size: `L=6`
- Development physics baseline: `beta=32`, `dtau=0.01`, `U1=0`
- Main small-parameter scan:
  - `Nbos = 10, 100, 1000`
  - `U2 = 1e-2, 1e-1, 1e0, 1e1, 1e2`

The current staged visual report is:

- [stage_report_v6a](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/archive_20260319/stage_report_v6a/report.md)

An example unresolved progress report is:

- [progress_report_n1000_u1e0_v2](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/archive_20260319/progress_report_n1000_u1e0_v2/report.md)
- [progress_report_n100_u1e1_v3](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/archive_20260319/progress_report_n100_u1e1_v3/report.md)

The notebook entry point is:

- [small_hmc_stage_report.ipynb](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/test/archive_small_benchmark/small_hmc_stage_report.ipynb)

The fixed full-grid renderer is:

- [render_small_hmc_full_grid.py](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/test/archive_small_benchmark/render_small_hmc_full_grid.py)

The current full-grid visual report is:

- [full_grid_progress_v3](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/full_grid_progress_v3/report.md)

The archived raw benchmark and tuning directories are now under:

- [archive_20260319](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/archive_20260319)

## Campaign Status

The small-parameter triangular correctness benchmark is considered complete enough to close.

Final interpretation:

- The HMC action and force implementation passed the ratio and finite-difference checks.
- The weak/moderate-coupling triangular benchmarks established a nontrivial strict-pass set.
- The final 15-point visual report is sufficient as the correctness handoff artifact for this campaign.
- On large `Nbos` and large `U2`, workstation-length local and HMC runs can both remain unthermalized, so those points are not useful as correctness references.
- Future work should stop centering local-vs-HMC comparison and instead focus on production-grade HMC thermalization, preconditioning, and HPC execution.

## Production Ramp

The recommended production ramp is intentionally conservative.
Do not jump directly to the final target.

Current preferred order:

1. `L=6`, easy projector, gradual `(Nbos, U2)` increase
2. `L=6`, same physics point, gradual `(beta, dtau)` refinement
3. `L=8 -> 10 -> 12`, with a reset to easier projector settings at each new `L`
4. `L=12`, push the target physics line toward `beta=256`, `dtau=0.001`
5. `L=21` spot checks only after `L=12` is believable

The current production stage gate is:

- `squareOcc` trace stabilizes across sample index
- `IPR` trace stabilizes across sample index
- different seeds or init-state families approach the same band

Do not promote a rung simply because acceptance looks good.

## What Is Already Confirmed

### 1. The HMC action and force are self-consistent

The strongest evidence is the ratio and finite-difference checks run on actual benchmark
configurations, including strong-coupling cases. Relevant diagnostics live under:

- [diag_n10_u1e2_ratio](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/diag_n10_u1e2_ratio)

Representative results:

- `ratio_scan_mean_reldiff ~ 1e-12`
- `force_fd_mean_reldiff ~ 1e-6`

Interpretation:

- If a benchmark mismatch remains after this, the default suspicion should be mixing,
  thermalization, or proposal design, not a sign error in the force.

### 2. There is now a genuine strict-pass health set

The current 9-point strict-pass set is rendered in `stage_report_v6a`.
Those cases are:

- `Nbos=10`: `U2=1e-2, 1e-1, 1e0, 1e1`
- `Nbos=100`: `U2=1e-2, 1e-1, 1e0`
- `Nbos=1000`: `U2=1e-2, 1e-1`

This is enough to support the statement that the triangular HMC implementation is
already correct over a nontrivial region of parameter space.

### 3. For some health points, HMC is already clearly better than local

Representative examples from `stage_report_v6a`:

- `Nbos=100, U2=1e-2`: `tau_local ~ 5.84`, `tau_hmc ~ 0.77`
- `Nbos=1000, U2=1e-2`: `tau_local ~ 6.23`, `tau_hmc ~ 0.81`

Interpretation:

- HMC is not merely reproducing the same means.
- In parts of the weak-coupling region it is already decorrelating faster than local.

## Tuning Rules That Already Look Reliable

### Acceptance is a diagnostic, not the objective

For this code, acceptance alone is not a safe optimizer.
Some candidates with very high acceptance still have poor `tau_int`.

Practical rule:

- First optimize for correctness and visible thermalization.
- Then compare `tau_int(doubleOcc)` and `ESS/sec`.

### Short warm=0 scans are useful but not sufficient

The short `tune` scans are still worth doing because they quickly separate:

- obviously stuck candidates
- numerically unstable candidates
- candidates with at least some useful motion

But the short scan ranking is not final.
For example, on `Nbos=1000, U2=1`, the short scan favored `mass=16, 32x0.04`,
yet long strict benchmarks still failed on `SF_Gamma`.

Practical rule:

- Always re-rank with a long direct local-vs-HMC benchmark before declaring a point healthy.

### Thermal cut really matters

The current workflow has already hit multiple cases where a candidate failed only because
the cut was too shallow.

Concrete example:

- `Nbos=100, U2=1e-2` initially failed only `SF_Gamma`
- raw-run repeated cut scan showed that `thermal_cut=448` removed the mismatch
- that point then passed the strict benchmark as
  [strict_healthy_n100_u1em2_v2](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/strict_healthy_n100_u1em2_v2)

Practical rule:

- If only one observable fails and the action/force checks are already clean, do not stop.
- First do a repeated cut scan from the raw run directories.

### Full trace storage is worth keeping

`small_benchmark_samples.csv` now stores the key series for all repeats, not only `repeat=0`.
This is the correct design for future production debugging because it enables:

- repeated thermal-cut scans
- repeat-to-repeat trace inspection
- quick notebook-based diagnosis without reopening every run manually

The staged report tooling now also writes per-repeat thermalization traces and includes
`SF_Gamma` directly in the default trace panel, because this observable has become the
main diagnostic for the unresolved `Nbos=1000, U2=1` mode.
If a long benchmark finishes its `runs/` directories but the summary files are missing,
`test/archive_small_benchmark/small_hmc_benchmark.py collect-benchmark` can now rebuild the `json/csv` tables
without rerunning the QMC jobs.

If you want one command for the current 15-point visual report,
`test/archive_small_benchmark/render_small_hmc_full_grid.py` wraps the staged renderer with the repository's
current best directory selection.

## Current Unresolved Modes

### `Nbos=1000, U2=1`

This is the current main unresolved correctness point.

Tested long candidates:

- [strict_healthy_n1000_u1e0_v1](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/strict_healthy_n1000_u1e0_v1)
  - `mass=16, nfrog=20, dt=0.08`
  - failed only `SF_Gamma`
- [strict_healthy_n1000_u1e0_v2](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/strict_healthy_n1000_u1e0_v2)
  - `mass=16, nfrog=32, dt=0.04`
  - failed only `SF_Gamma`
- [strict_healthy_n1000_u1e0_v3_w512](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/strict_healthy_n1000_u1e0_v3_w512)
  - same HMC parameters as `v2`, but with `warm=512`
  - still under evaluation during this stage of the campaign

Observed pattern:

- Most observables agree.
- `SF_Gamma` remains systematically lower in HMC than in local.
- Segment means suggest the HMC chain is not just noisy; it is spending too much time in a
  lower-`SF_Gamma` region.

Interpretation:

- This point is likely limited by a slow collective mode that the current global HMC
  proposal still does not decorrelate well enough.

### `Nbos=100, U2=10`

This is the current strong-coupling tuning point.

Useful tuning directory:

- [tune_n100_u1e1_massscan_v1](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/tune_n100_u1e1_massscan_v1)

Current trend:

- `mass=1` can be made mobile only with much smaller steps, but `tau` remains poor
- `mass=4` improves robustness and acceptance
- the first promising strong-coupling candidate from this scan is
  `mass=64, nfrog=24, dt=0.06`

A strict benchmark has been launched from that candidate:

- [strict_healthy_n100_u1e1_v2](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/strict_healthy_n100_u1e1_v2)

## Preconditioning Trends

### Uniform mass scaling already helps

Two cases where mass preconditioning clearly matters:

- `Nbos=1000, U2=1`
- `Nbos=100, U2=10`

Empirical lesson:

- Increasing the mass can move a point from "completely stuck" to "mobile and benchmarkable"
- but the best long-run choice is not necessarily the one that looks best in a short scan

### Acceptance can become too high

At large masses the code often reaches `acceptance > 0.95`.
This is not automatically good.

Empirical lesson:

- Very high acceptance frequently means the trajectory is too conservative
- that can still leave `tau_int` large, especially for the hard collective modes

## Why Local Can Sometimes Mix Better Than the Current HMC

This has now shown up clearly enough in the small-parameter campaign that it should be
recorded explicitly.

Observed pattern:

- local updates are not always faster in wall-clock time
- but on some unresolved points they do wander away from their current state more easily
- the current HMC implementation can keep high acceptance while still exploring only a narrow
  part of configuration space

The likely reason is structural:

- local Metropolis updates are stochastic and irreversible at the one-site/one-time-slice level
- even if each accepted move is small, many such moves can slowly diffuse along directions that
  the current HMC mass matrix does not match well
- the present HMC uses one uniform virtual mass and one leapfrog scale for all modes
- when the hard modes set the stability bound, the soft collective modes can end up under-driven
- that produces the practical symptom of "high acceptance but poor global mixing"

Interpretation:

- if local and HMC disagree at a point, and the action/force checks are already clean,
  do not assume the force is wrong
- first ask whether HMC is trapped on a poorly preconditioned slow manifold

This is the main motivation for the next preconditioning ideas:

- mode-dependent masses
- Fourier acceleration
- multi-timescale integration
- hybrid HMC + local interleaving

## Next Preconditioning Ideas Worth Trying

These are now justified by actual data, not only by theory.

### 1. Mode-dependent masses

Motivation:

- the unresolved `SF_Gamma` mismatch at `Nbos=1000, U2=1`
- the strong-coupling slowdown at `Nbos=100, U2=10`

Likely design:

- assign larger virtual masses to the stiff short-wavelength modes
- keep softer masses for the slow collective modes

### 2. Fourier acceleration

Motivation:

- the problematic modes already look collective rather than site-local
- a lattice-momentum basis should be more natural for preconditioning than uniform real-space mass

### 3. Multi-timescale integration

Motivation:

- the current leapfrog has to obey the hard-mode stability bound
- that makes the soft-mode motion inefficient

### 4. Hybrid HMC + local updates

Motivation:

- local updates can still help shuffle the modes that HMC moves poorly
- this is especially attractive for the unresolved `SF_Gamma` problem

## Recommended Workflow For Future HPC Campaigns

The future production targets are expected to look like:

- `L >= 21`
- `beta = 256`
- `dtau = 0.001`
- `Ltrot = 256000`
- `Nbos = 1e5, 1e6, 1e7`
- `U2 = 1e2, 1e3, 1e4`

At those scales, the local development lessons should be reused as follows:

### 1. Separate correctness from efficiency

- First prove that a candidate reproduces the local benchmark on the smallest tractable point
- Only then optimize for `ESS/sec`

### 2. Keep staged reports, not only pass/fail tables

Even when a point is not yet strict-pass, keep:

- tune CSV
- strict benchmark CSV
- trace PNG
- raw run directory

The trend itself is often informative enough to guide the next scan.

### 3. Expect much deeper thermal cuts

The small campaign already shows that `thermal_cut` is observable-dependent.
On HPC production points, assume that:

- `warm`
- `thermal_cut`
- and even measurement window length

may all need to be increased together.

### 4. Do not trust one ranking source

Use all three:

- short tune scan
- long strict benchmark
- raw-run trace inspection

### 5. Keep a production candidate table

For each physical point, keep a small table of:

- sampler mode
- `nfrog`
- `dt`
- `mass`
- `jitter`
- acceptance
- `tau_int`
- `ESS/sec`
- pass/fail status
- notes on which observable was hardest to match

## Commands Worth Reusing

Small strong-coupling tune:

```bash
/home/yyk/conda/envs/notebook/bin/python test/archive_small_benchmark/small_hmc_benchmark.py tune \
  --nbos-values 100 \
  --u2-values 1e1 \
  --warm 0 \
  --bins 160 \
  --sweeps 1 \
  --grid '24:0.06,32:0.04,40:0.03,48:0.025,64:0.02' \
  --hmc-jitter 2 \
  --hmc-mass-grid '1,4,16,64' \
  --repeats 2 \
  --confin-root data/triangular_hmc_small_benchmark/local_seed_grid_v2 \
  --work-root data/triangular_hmc_small_benchmark/tune_n100_u1e1_massscan_v1
```

Strict benchmark from a tuned candidate:

```bash
/home/yyk/conda/envs/notebook/bin/python test/archive_small_benchmark/small_hmc_benchmark.py benchmark \
  --nbos-values 100 \
  --u2-values 1e1 \
  --bins 640 \
  --thermal-cut 384 \
  --warm 0 \
  --repeats 4 \
  --hmc-nfrog 24 \
  --hmc-dt 0.06 \
  --hmc-jitter 2 \
  --hmc-mass 64 \
  --confin-root data/triangular_hmc_small_benchmark/local_seed_grid_v2 \
  --work-root data/triangular_hmc_small_benchmark/strict_healthy_n100_u1e1_v2
```

Rebuild a staged report from the known strict-pass set:

```bash
/home/yyk/conda/envs/notebook/bin/python test/archive_small_benchmark/render_small_hmc_stage_report.py \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n10_u1em2_v2 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n10_u1em1_v1 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n10_u1e0_v1 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n10_u1e1_v1 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n100_u1em2_v2 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n100_u1em1_v1 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n100_u1e0_v1 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n1000_u1em2_v1 \
  data/triangular_hmc_small_benchmark/archive_20260319/strict_healthy_n1000_u1em1_v1 \
  --tune-csv data/triangular_hmc_small_benchmark/archive_20260319/full_global_tune_grid_v2/small_tune.csv \
  --output-dir data/triangular_hmc_small_benchmark/archive_20260319/stage_report_v6a
```
