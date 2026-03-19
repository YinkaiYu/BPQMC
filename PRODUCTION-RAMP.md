# Production Ramp Plan

This note is the persistent execution plan for the triangular-lattice production HMC ramp.
It is intentionally shorter and more action-oriented than `HMC-REFERENCE.md`.
`HMC-REFERENCE.md` remains the main narrative handoff note; this file is the live checklist.

## Final Target

- lattice: `triangular`
- physics target: `Nbos = 1e5`, `U2 = 1e3`, `U1 = 0`
- size target: `L >= 21`
- projector target: `beta = 256`, `dtau = 0.001`
- key observables: `squareOcc`, `IPR`

## Decision Rules

- Do not promote to a harder rung until `squareOcc` and `IPR` look stable across retained windows.
- Prefer deeper `U2` / `Nbos` ramps over deeper `beta/dtau` ramps when the representative baseline
  ladder already changes only weakly with projector settings.
- Treat short `tune` scans at large `U2` as provisional only.
  Any candidate that looks good on `64 / 32 / 32` must still be checked by a long `stage`.
- If every tune candidate is stuck, do not create a usable recommendation.
  Shrink `dt`, increase `nfrog`, or deepen the tune before launching a long stage.

## Current Baseline Conclusion

- Representative baseline: `L = 6`, `Nbos = 1e3`, `U2 = 1`
- Healthy ladder through:
  - `beta = 32`, `dtau = 0.01`
  - `beta = 64`, `dtau = 0.01`
  - `beta = 96`, `dtau = 0.008`
  - `beta = 128`, `dtau = 0.005`
  - `beta = 160`, `dtau = 0.004`
- The representative overview now shows only about `7.6e-4` relative movement in `squareOcc`
  and `7.7e-4` in `IPR` across `beta = 32 -> 192`, `dtau = 0.01 -> 0.002`.
- Therefore the active workstation priority is the higher-`U2` ladder, not deeper projector scans.

## Current Healthy Rungs

- `L = 6`, `Nbos = 1e4`, `U2 = 1`
  - preconditioned stage winner:
    - `mass = 16`
    - `uniform_mass = 1`
    - `shell1_mass = 0`
    - `nfrog = 20`
    - `dt = 0.008`
- `L = 6`, `Nbos = 1e4`, `U2 = 10`
  - scalar-mass stage winner:
    - `mass = 4`
    - `nfrog = 16`
    - `dt = 0.006`
- `L = 6`, `Nbos = 1e4`, `U2 = 30`
  - preconditioned stage winner:
    - `mass = 16`
    - `uniform_mass = 1`
    - `shell1_mass = 0`
    - `nfrog = 12`
    - `dt = 0.006`
- `L = 6`, `Nbos = 1e4`, `U2 = 100`
  - conservative preconditioned stage winner:
    - `mass = 16`
    - `uniform_mass = 1`
    - `shell1_mass = 0`
    - `nfrog = 20`
    - `dt = 0.001`

## Active Rungs

- `L = 6`, `Nbos = 1e4`, `U2 = 300`
  - first tune grid had no viable candidate
  - current conservative long stage:
    - `mass = 16`
    - `uniform_mass = 1`
    - `shell1_mass = 0`
    - `nfrog = 28`
    - `dt = 0.0005`
  - current expectation:
    - the old `1024 / 512 / 512` stage is no longer the one to trust
    - the deeper `2048 / 1024 / 1024` rerun has already produced a much healthier formal
      `1/2` retained-window `stable_window` summary
    - the missing repeat still needs to be completed before `U2 = 300` can be treated as a fully
      settled representative rung
    - promotion to `U2 = 1e3` should wait for that deeper `2/2` result
- `L = 6`, `Nbos = 1e4`, `U2 = 1e3`
  - the old scalar/uniform-only reference rung
    `m = 16`, `mu = 1`, `mk = 0`, `nfrog = 20`, `dt = 0.0004`
    is now closed as `strong_drift`
  - shell1-only replacement path:
    - `mk = 1` is closed as too aggressive
    - completed `mk = 4` and `mk = 8` `1024 / 512 / 512` stages are both still `strong_drift`
  - current active replacement path:
    - keep `mass = 16`
    - keep `uniform_mass = 1`
    - keep `shell1_mass = 8`
    - add the new `shell2_mass` split
    - current active short-tune probe:
      - `shell2_mass = 4`
      - grid:
        - `12 x 0.0003`
        - `16 x 0.00025`
        - `20 x 0.0002`
        - `24 x 0.00015`
      - current winner:
        - `12 x 0.0003`, jitter `=2`
      - matching warm=`0` long trace is now running

## Immediate Next Steps

1. Finish the deeper `U2 = 300` long stage to a full `2/2` repeat set and re-check
   `squareOcc` / `IPR` retained-window agreement in the unified overview report.
2. Finish reading the active `U2 = 1e3` shell2 warm=`0` long trace:
   - `mk1 = 8`, `mk2 = 4`
   - winner from the first short shell2 probe:
     - `12 x 0.0003`
   If the early trace is materially flatter than the old shell1-only lines,
   promote it immediately to a deeper `1024 / 512 / 512` stage.
3. Compare the shell2 probe against the closed bad references:
   - `mk = 0`, `20 x 0.0004`
   - shell1-only `mk = 4`, `24 x 0.0002`
   - shell1-only `mk = 8`, `24 x 0.00025`
   If shell2 still does not flatten materially faster, widen the low-|k| basis again rather than
   continuing to rescan the same scalar/uniform/shell1-only family.
4. If `U2 = 300` stays healthy at `2048 / 1024 / 1024`, keep that geometry as the
   representative preconditioned rung and use it as the production-side reference.
5. Only after `U2 = 1e3` has at least a partially healthy geometry should the ramp
   move to larger `Nbos` or larger `L`.

## Current Early-Trace Reading

- `U2 = 1e3`, old closed bad reference:
  - `mk = 0`, `20 x 0.0004`
  - completed long-stage `squareOcc` head-to-tail drift/span is about `0.60`
- `U2 = 1e3`, shell1-only completed stages:
  - `mk = 4`, `24 x 0.0002`
  - `mk = 8`, `24 x 0.00025`
  - both still end as `strong_drift`
- `U2 = 1e3`, current shell2 probe:
  - `mk1 = 8`, `mk2 = 4`
  - grid:
    - `12 x 0.0003`
    - `16 x 0.00025`
    - `20 x 0.0002`
    - `24 x 0.00015`
- Working interpretation:
  - shell1-only preconditioning helped the very early trace shape but did not solve the retained-window drift
  - the next meaningful question is whether widening the low-|k| split to shell2 can convert that
    early improvement into a genuinely stable retained window

## Medium-Term Ladder

Use this order unless a rung clearly blocks:

1. `L = 6`, `Nbos = 1e4`, `U2 = 1e2`
2. `L = 6`, `Nbos = 1e4`, `U2 = 3e2`
3. `L = 6`, `Nbos = 1e4`, `U2 = 1e3`
4. `L = 6`, `Nbos = 1e5`, `U2 = 1e3`
5. `L = 8`, easier anchor + main target line
6. `L = 10`, easier anchor + main target line
7. `L = 12`, easier anchor + main target line
8. only then revisit whether deeper `beta/dtau` scans are still needed on the target line

## HPC Handoff Notes

- Every stage should remain resumable with `collect` / `report`.
- Keep using the unified production entry point:
  - `data/triangular_hmc_production/overview/report.md`
- For future HPC work, do not assume the last tuned workstation parameters are final.
  Re-run at least a short `tune` on the new rung before committing a long batch campaign.
