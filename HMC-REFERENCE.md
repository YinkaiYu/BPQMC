# HMC Reference Notes

This document is the working handoff note for the triangular-lattice HMC campaign.
It records what has already been validated, what is still unresolved, and which tuning
trends look robust enough to reuse when moving toward larger lattices and HPC runs.
For the current `test/` directory layout, see `test/README.md`; archived small-benchmark
helpers referenced here now live under `test/archive_small_benchmark/`.
For the persistent action checklist that should survive long tuning sessions, see
`PRODUCTION-RAMP.md`.

## Current Priority

The small-parameter local-vs-HMC correctness campaign is closed.
The active task is now production-grade HMC bring-up:

- move gradually from easy workstation points toward the production physics line
- judge stages primarily by `squareOcc` and `IPR` thermalization
- treat acceptance and `ESS/sec` as secondary diagnostics
- keep every stage resumable and reportable for later HPC handoff

The active production CLI is `test/production_hmc.py`:

- `tune`: short candidate scans
- `collect-tune`: rebuild tune summaries from finished or partially finished tune `runs/`
- `stage`: long HMC-only runs with `squareOcc`/`IPR` traces
- `collect`: rebuild stage summaries from completed or partially completed stage `runs/`
- `report`: render Markdown + PNG summaries
- `--hmc-mass-spatial-uniform`: optional per-time-slice spatial-uniform-mode mass split
  on top of the residual `--hmc-mass`; `0` disables it
- `--hmc-mass-spatial-shell1`: optional triangular lowest-|k|-shell mass split
  on top of both the residual `--hmc-mass` and the spatially uniform split;
  `0` disables it

For seed/init-state convergence checks, the production driver can now expand multiple
initial-state families in one invocation via:

- `--ini-type-values`
- `--ini-ampl-values`

The default output root for this new workflow should be under:

- `data/triangular_hmc_production/`

The unified entry point for the accumulated production results is now:

- [overview/report.md](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/overview/report.md)
- the same `report.md` now starts with `## Current Representative Stage Per Case`,
  so one current retained-window result per fixed `(L, Nbos, U2)` can be reviewed first
- the same `report.md` now also contains a `## Live Progress` section for in-flight large-`U2` stages
- explicit fixed-`beta` or fixed-`dtau` scans can be launched with
  [production_convergence_scan.py](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/test/production_convergence_scan.py)
- the overview report now shows `bins`, `thermal_cut`, `warm`, and `samples/post` per stage,
  so late slow-mode drift is easier to spot without opening each rung directory separately
- the overview report now also renders per-case relative-change plots and a representative
  one-stage-per-`(beta, dtau)` convergence table, so `beta/dtau` trends can be judged without
  mixing together different `Nbos`/`U2` scales or old failed tuning attempts
- in the merged trace panels, the dashed line is only the configured `thermal_cut`, and the
  overlaid lines are independent stage runs rather than one continued chain
- the production stage gate now checks the worst retained repeat as well as the mean drift,
  so one drifting chain can no longer be hidden by averaging over healthier repeats
- the production stage gate also checks cross-repeat retained-window mismatch, reported in the
  overview as `repeat span`, so different seeds settling on different plateaus are no longer hidden

## Current Stage / Next Stage / Blocker

Current healthy ladder:

- `L=6`, `Nbos=1e3`, `U2=1`
- `beta=32`, `dtau=0.01`
- `beta=64`, `dtau=0.01`
- `beta=96`, `dtau=0.008`
- `beta=128`, `dtau=0.005`
- `beta=160`, `dtau=0.004`

Current next rung:

- the workstation priority has shifted from pushing `beta/dtau` deeper on the easy baseline point
  to pushing `U2` upward at fixed `beta=32`, `dtau=0.01`
- current active `U2` ladder:
  - `L=6`, `Nbos=1e4`, `U2=30`
    - long-stage winner:
      - residual mass `m=16`, uniform mass `mu=1`
      - `nfrog=12`, `dt=0.006`, jitter `=2`
      - `overall_status = healthy`
    - completed long stage:
      [l6_n1e4_u3e1_beta32_dtau1em2_stage_m16_mu1_nf12_dt0p006_diag1024](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e1_beta32_dtau1em2_stage_m16_mu1_nf12_dt0p006_diag1024)
  - `L=6`, `Nbos=1e4`, `U2=100`
    - conservative short-tune winner: residual mass `m=16`, uniform mass `mu=1`,
      `nfrog=20`, `dt=0.001`, jitter `=2`
    - completed long stage:
      [l6_n1e4_u1e2_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p001_diag1024](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e2_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p001_diag1024)
    - current result:
      - `overall_status = healthy`
      - `Accept_HMC ≈ 0.980 / 0.981`
      - `squareOcc/IPR drift/span max ≈ 0.0135`
      - retained repeat-span ratio `≈ 0.051`
  - `L=6`, `Nbos=1e4`, `U2=300`
    - conservative short-tune winner: residual mass `m=16`, uniform mass `mu=1`,
      `nfrog=28`, `dt=0.0005`, jitter `=2`
    - current long stage:
      [l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag1024](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag1024)
    - current result:
      - `overall_status = needs_review`
      - case status `slow_drift`
      - `Accept_HMC ≈ 0.986 / 0.979`
      - `squareOcc/IPR drift/span max ≈ 0.401`
      - retained repeat-span ratio `≈ 0.055`
    - deeper formal stage:
      [l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag2048](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag2048)
    - current deeper-stage reading:
      - the first completed repeat at `2048 / 1024 / 1024` already gives a formal
        `stable_window` summary, with `squareOcc/IPR drift/span max ≈ 0.019`
      - the active task is now to finish the missing repeat and check whether the
        cross-repeat retained-window mismatch stays small
  - `L=6`, `Nbos=1e4`, `U2=1000`
    - exploratory tune reports:
      - [m16 probe](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e3_beta32_dtau1em2_tune_m16_mu1_probe/report/report.md)
      - [m16 aggressive](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e3_beta32_dtau1em2_tune_m16_mu1_aggressive/report/report.md)
      - [m32 probe](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e3_beta32_dtau1em2_tune_m32_mu1_probe/report/report.md)
    - completed reference stage:
      [l6_n1e4_u1e3_beta32_dtau1em2_stage_m16_mu1_nf20_dt4em4_diag1024_probe/report/report.md](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e3_beta32_dtau1em2_stage_m16_mu1_nf20_dt4em4_diag1024_probe/report/report.md)
    - current reading:
      - none of the scalar/uniform-only short scans removes the slow mode yet; all remain at `tau_int(doubleOcc) ~ 25`
      - the mildly more aggressive `m=16`, `mu=1`, `20 x 0.0004` geometry has now been ruled only a reference failure:
        its completed long stage is `strong_drift`, with `squareOcc/IPR drift/span ≈ 0.660`
      - the active next step is now the new lowest-shell preconditioner, not another scalar/uniform-only rerun
    - new shell1 in-flight tune:
      [l6_n1e4_u1e3_beta32_dtau1em2_tune_m16_mu1_mk1_probe](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e3_beta32_dtau1em2_tune_m16_mu1_mk1_probe)
      - residual mass `m=16`
      - uniform mass `mu=1`
      - shell1 mass `mk=1`
      - grid: `12x5e-4`, `16x4e-4`, `20x4e-4`, `24x3e-4`, `28x2.5e-4`

Current main blocker:

- the scalar-mass `L=6`, `Nbos=1e4`, `U2=1`, `beta=32`, `dtau=0.01` rung is no longer the
  main blocker; it has now been converted into a usable preconditioned base geometry
- current production blocker is the higher-`U2` ladder itself:
  - how far the spatial-uniform preconditioner can be pushed before the retained
    `squareOcc/IPR` windows split again across seeds
  - whether the same mass split remains useful when `U2` is increased toward `1e2`, `1e3`,
    and later `Nbos` is also raised toward `1e5`
  - the immediate blocker rung is now `L=6`, `Nbos=1e4`, `U2=1000`
    - the old scalar/uniform-only reference geometry `m=16`, `mu=1`, `20 x 0.0004` is formally `strong_drift`
    - the new active idea is to split the triangular lowest nonzero momentum shell away from the residual modes
  - `L=6`, `Nbos=1e4`, `U2=300` is still active, but its deeper `2048 / 1024 / 1024` rerun already looks much healthier;
    the main remaining question there is cross-repeat agreement, not obvious single-trace drift
- trace-based tuning lesson from the original `Nbos=1e4, U2=1` blocker point:
  - a short traced run with `m=16`, `mu=1`, `20 x 0.008` gives an approximate uniform-mode period
    `T_md ~ 0.16` from `phi_mean_f2`
  - this means `Nfrog * dt ~ 0.16` is closer to a full cycle than to the earlier
    quarter-period heuristic target
  - the next short scan was therefore centered near `Nfrog * dt ~ 0.04`
  - that quarter-period-informed short scan favored `6 x 0.008`, jitter `=1`
    on `ESS/sec(doubleOcc)` alone, but the longer stage later showed it was still `slow_drift`
  - the decisive longer-stage winner for this physics point is instead:
    - residual mass `m=16`
    - spatial-uniform mass `mu=1`
    - `nfrog=20`, `dt=0.008`, jitter `=2`
    - work root:
      [l6_n1e4_u1_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p008_diag512](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p008_diag512)
    - completed-stage result:
      - `Accept_HMC ≈ 0.985`
      - `squareOcc/IPR drift/span max ≈ 0.092`
      - retained repeat-span ratio `≈ 0.022`
      - final stage status `stable_window`

## Current Production Checkpoint

As of the current workstation bring-up round, the first production-style stage data under
`data/triangular_hmc_production/` already separates four regimes clearly.

### Healthy projector ladder on the baseline physics point

- [l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta32_dtau1em2_stage_m4_nf12_dt0p015/report/report.md)
- [l6_n1e3_u1_beta64_dtau1em2_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta64_dtau1em2_stage_m4_nf8_dt0p02/report/report.md)
- [l6_n1e3_u1_beta96_dtau8em3_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta96_dtau8em3_stage_m4_nf8_dt0p02/report/report.md)
- [l6_n1e3_u1_beta128_dtau5em3_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta128_dtau5em3_stage_m4_nf8_dt0p02/report/report.md)
- [l6_n1e3_u1_beta160_dtau4em3_stage_m4_nf8_dt0p02](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e3_u1_beta160_dtau4em3_stage_m4_nf8_dt0p02/report/report.md)
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
  - `beta=160`, `dtau=0.004`: `stable_window`, `squareOcc/IPR drift/span ~ 0.094`

Interpretation:

- The baseline triangular rung is now healthy through three projector settings.
- The baseline triangular rung is now healthy through five projector settings.
- This is the first real evidence that the production ramp can proceed in projector depth
  rather than only at one easy point.
- The `beta=96` rung is slower in wall-clock time than `beta=64`, but its traces still pass
  the current `squareOcc/IPR` gate cleanly.
- The `beta=128` rung continues the same trend: wall-clock time rises again, but the retained
  `squareOcc/IPR` window is still convincingly stable.
- The `beta=160` rung continues to pass with the same scalar-mass setting (`m=4`, `8 x 0.02`,
  jitter `=2`). The main cost increase is wall-clock time rather than a sudden thermalization failure.

Current convergence interpretation for the healthy baseline ladder:

- `squareOcc` and `IPR` are already numerically close across the healthy `beta=96 -> 128 -> 160`
  rungs.
- With the current representative overview selection, the baseline `L=6, Nbos=1e3, U2=1`
  ladder now shows:
  - `squareOcc` relative spread across `beta=32 -> 192`, `dtau=0.01 -> 0.002`:
    about `7.56e-4`
  - `IPR` relative spread across the same ladder:
    about `7.74e-4`
- This is still not a strict fixed-`beta` / fixed-`dtau` convergence proof because the ladder
  changes both parameters together, but it is already strong enough to justify moving the main
  workstation effort from deeper projector scans to higher-`U2` scans.
- The correct next step is:
  - if observables stop changing as `beta` increases, do not push to larger `beta` just because it is possible
  - if observables stop changing as `dtau` decreases, do not push to smaller `dtau` just because it is possible
  - otherwise, continue the ladder or run explicit fixed-`beta` / fixed-`dtau` comparison scans
  - in the current branch state, deeper projector scans are temporarily secondary to the
    `U2=30 -> 100 -> ...` production ramp
- For those explicit convergence scans, treat `O(10^3)` post-warm samples as normal rather than exceptional.
  The convenience wrapper now defaults to `bins=1024`, `thermal_cut=512`, `warm=512`.
- The canonical `test/production_hmc.py stage` / `collect` workflow now uses the same
  `1024 / 512 / 512` defaults unless overridden explicitly.
- The SLURM helper `test/dqmc_production` mirrors that split:
  - `PROD_COMMAND=tune` / `collect-tune` default to `64 / 32 / 32`
  - `PROD_COMMAND=stage` / `collect` default to `1024 / 512 / 512`
- The new spatial-uniform mass split is intentionally script-driven rather than baked into the
  parameter file format:
  - `hmc_mass` remains the residual-mode mass recorded in `paramC_sets.txt`
  - `--hmc-mass-spatial-uniform` is forwarded by the Python tooling as
    `BPQMC_HMC_MASS_SPATIAL_UNIFORM`
  - `--hmc-mass-spatial-shell1` is forwarded by the Python tooling as
    `BPQMC_HMC_MASS_SPATIAL_SHELL1`
  - on the triangular lattice it builds a real orthonormal basis of the first nonzero
    reciprocal-shell cosine/sine modes and assigns them their own leapfrog mass
  - `info.txt` now records both masses so archived runs remain self-describing

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
- The older `thermal_cut=192` setting is visibly too shallow on the `rep0` trace.
  In the `384`-sample stage, the clear upward drift continues until roughly `sample index 320~350`
  before the trace approaches the later plateau.
- This is exactly why the newer long rerun was promoted to `bins=1024`, `thermal_cut=512`, `warm=512`:
  the short/deeper stages were still underestimating how late the slow mode settles.
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

### Higher-`U2` tuning updates on the `L=6, Nbos=1e4` ladder

- `U2=30`
  - scalar short-tune winner:
    - `m=4`, `nfrog=20`, `dt=0.003`, jitter `=2`
    - `acceptance ≈ 0.859`
    - `tau_int(doubleOcc) ≈ 1.94`
    - `ESS/sec(doubleOcc) ≈ 0.574`
  - spatial-uniform preconditioned short-tune winner:
    - `m=16`, `mu=1`, `nfrog=12`, `dt=0.006`, jitter `=2`
    - `acceptance ≈ 0.844`
    - `tau_int(doubleOcc) ≈ 0.757`
    - `ESS/sec(doubleOcc) ≈ 2.34`
  - takeaway:
    - on this rung the preconditioner is no longer just a stability aid; it is already the
      better efficiency geometry as well
  - completed long stage:
    - [l6_n1e4_u3e1_beta32_dtau1em2_stage_m16_mu1_nf12_dt0p006_diag1024](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e1_beta32_dtau1em2_stage_m16_mu1_nf12_dt0p006_diag1024/report/report.md)
    - `2/2` repeats complete
    - `overall_status = healthy`
    - `acceptance ≈ 0.790`
    - `squareOcc/IPR drift/span max ≈ 0.038`
    - retained repeat-span ratio `≈ 0.0068`
- `U2=100`
  - first short-tune / long-stage lesson:
    - `16 x 0.003` looked viable on the shallow `64 / 32 / 32` tune
    - the corresponding long stage froze immediately with `Accept_HMC = 0`
  - current conservative short-tune winner:
    - `m=16`, `mu=1`, `nfrog=20`, `dt=0.001`, jitter `=2`
    - `acceptance ≈ 0.992`
    - `tau_int(doubleOcc) ≈ 5.88`
    - `ESS/sec(doubleOcc) ≈ 0.193`
  - current long-stage checkpoint:
    - [l6_n1e4_u1e2_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p001_diag1024](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u1e2_beta32_dtau1em2_stage_m16_mu1_nf20_dt0p001_diag1024/report/report.md)
    - `rep0` complete, `rep1` running
    - partial status: `stable_window`
    - `acceptance ≈ 0.979`
    - retained-window `squareOcc/IPR drift/span ≈ 0.0135`
  - takeaway:
    - by `U2=100`, the allowed step-size window has narrowed enough that shallow short scans
      are no longer trustworthy; the tuning workflow itself must become more conservative
    - the rung is now mobile again, so the blocker has moved from “can it move at all?” to
      “will multiple repeats land on the same retained plateau?”
- `U2=300`
  - first short-tune grid:
    - [l6_n1e4_u3e2_beta32_dtau1em2_tune_m16_mu1](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_tune_m16_mu1/report/report.md)
    - no viable candidate
  - current conservative short-tune winner:
    - [l6_n1e4_u3e2_beta32_dtau1em2_tune_m16_mu1_conservative](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_tune_m16_mu1_conservative/report/report.md)
    - `m=16`, `mu=1`, `nfrog=28`, `dt=0.0005`, jitter `=2`
    - `acceptance ≈ 0.992`
    - `tau_int(doubleOcc) ≈ 1.59`
    - `ESS/sec(doubleOcc) ≈ 0.524`
  - active long stage:
    - [l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag1024](/mnt/c/users/newton/documents/ligroupiop/2408_bosonsignproblem/code_bpqmc/data/triangular_hmc_production/l6_n1e4_u3e2_beta32_dtau1em2_stage_m16_mu1_nf28_dt5em4_diag1024)
  - takeaway:
    - the strong-coupling ladder still opens if `dt` is reduced aggressively enough; the
      window has not closed yet, but it is becoming narrow quickly

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
