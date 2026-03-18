# HMC Reference Notes

This document is the working handoff note for the triangular-lattice HMC campaign.
It records what has already been validated, what is still unresolved, and which tuning
trends look robust enough to reuse when moving toward larger lattices and HPC runs.

## Scope

- Lattice: `triangular`
- Development benchmark size: `L=6`
- Development physics baseline: `beta=32`, `dtau=0.01`, `U1=0`
- Main small-parameter scan:
  - `Nbos = 10, 100, 1000`
  - `U2 = 1e-2, 1e-1, 1e0, 1e1, 1e2`

The current staged visual report is:

- [stage_report_v6a](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_small_benchmark/stage_report_v6a/report.md)

The notebook entry point is:

- [small_hmc_stage_report.ipynb](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/test/small_hmc_stage_report.ipynb)

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
/home/yyk/conda/envs/notebook/bin/python test/small_hmc_benchmark.py tune \
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
/home/yyk/conda/envs/notebook/bin/python test/small_hmc_benchmark.py benchmark \
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
/home/yyk/conda/envs/notebook/bin/python test/render_small_hmc_stage_report.py \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1em2_v2 \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1em1_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1e0_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1e1_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n100_u1em2_v2 \
  data/triangular_hmc_small_benchmark/strict_healthy_n100_u1em1_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n100_u1e0_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n1000_u1em2_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n1000_u1em1_v1 \
  --tune-csv data/triangular_hmc_small_benchmark/full_global_tune_grid_v2/small_tune.csv \
  --output-dir data/triangular_hmc_small_benchmark/stage_report_v6a
```
