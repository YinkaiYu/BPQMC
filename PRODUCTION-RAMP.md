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
    - `nfrog = 12`
    - `dt = 0.006`
- `L = 6`, `Nbos = 1e4`, `U2 = 100`
  - conservative preconditioned stage winner:
    - `mass = 16`
    - `uniform_mass = 1`
    - `nfrog = 20`
    - `dt = 0.001`

## Active Rungs

- `L = 6`, `Nbos = 1e4`, `U2 = 300`
  - first tune grid had no viable candidate
  - current conservative long stage:
    - `mass = 16`
    - `uniform_mass = 1`
    - `nfrog = 28`
    - `dt = 0.0005`
  - current expectation:
    - the current `2/2` stage is no longer stuck, but it is still `slow_drift`
    - the next step is a deeper or smaller-`dt` retune, not immediate promotion
- `L = 6`, `Nbos = 1e4`, `U2 = 1e3`
  - this is the next physics rung once the `U2 = 300` geometry is made reliable enough
  - first action should be a conservative preconditioned tune, not a blind reuse of `U2 = 300`

## Immediate Next Steps

1. Re-render the merged production overview so the latest `U2 = 100` healthy stage and
   `U2 = 300` slow-drift stage appear in the unified entry point.
2. Launch a conservative retune for `L = 6`, `Nbos = 1e4`, `U2 = 300` around smaller `dt`
   than `28 x 0.0005`, because the current retained window still drifts.
3. If the retuned `U2 = 300` geometry becomes stable across repeats, promote to `U2 = 1e3`.
4. Only after `U2 = 1e3` has at least a partially healthy geometry should the ramp move to
   larger `Nbos` or larger `L`.

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
