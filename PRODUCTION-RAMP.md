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
    - shell-by-shell probes are now effectively closed:
      - first completed shell2 family:
        - `shell1_mass = 8`
        - `shell2_mass = 4`
        - completed `warm=0` and deeper stages are both still `strong_drift`
      - lighter shell2 follow-up:
        - `shell1_mass = 4`
        - `shell2_mass = 2`
        - grid:
          - `12 x 0.0003`
          - `16 x 0.00025`
          - `20 x 0.0002`
          - `24 x 0.00015`
        - all candidates still kept `acceptance = 1` with huge positive `DeltaH`
    - current active probe:
      - broader combined low-|k| split:
        - `lowk_mass = 4`
        - grid:
          - `12 x 0.0003`
          - `16 x 0.00025`
          - `20 x 0.0002`
          - `24 x 0.00015`
        - completed `warm=0` stages are now both `strong_drift`
    - current completed broader-family leads:
      - `lowk_mass = 2`
        - tune winner: `16 x 0.0002`
        - moving-family stages on `51001` / `52001` still remain `strong_drift`
        - frozen-family stage on `50001` remains `stuck_or_invalid`
      - `lowk_mass = 1`
        - tune winner: `20 x 0.00012`
        - completed moving-family stages at `51001` / `52001` still remain `strong_drift`
        - partial three-family `256 / 128 / 32` stage also remains `strong_drift`
    - current active broader-family probe:
      - `lowk_mass = 0.5`
      - fixed-family short-tune winner:
        - `20 x 0.00008`
        - `acceptance ≈ 0.995`
        - `tau_int(doubleOcc) ≈ 4.41`
        - `ESS/sec ≈ 0.591`
      - active retained-window follow-up:
        - `512 / 256 / 32`
        - `20 x 0.00008`
        - explicit seed family block `50001, 51001, 52001`

## Immediate Next Steps

1. Finish the deeper `U2 = 300` long stage to a full `2/2` repeat set and re-check
   `squareOcc` / `IPR` retained-window agreement in the unified overview report.
2. Move past the now-closed three-shell `m_lowk = 0.5` family and start the widened
   grouped low-|k| probe:
   - `m = 16`, `mu = 1`, `m_lowk = 0.5`, `lowk_shells = 4`
   - explicit seed-family block `50001, 51001, 52001`
   - first target grid:
     - `16 x 0.0001`
     - `20 x 0.00008`
     - `24 x 0.00006`
     - `28 x 0.00005`
   - immediate check:
     - whether widening the Fourier-accelerated low-|k| basis reduces the
       worst-repeat retained drift below the closed three-shell
       `m_lowk = 0.5` / `1` / `2` / `4` families
3. Compare the active broader low-|k| probe against the closed bad references:
   - `mk = 0`, `20 x 0.0004`
   - shell1-only `mk = 4`, `24 x 0.0002`
   - shell1-only `mk = 8`, `24 x 0.00025`
   - first shell2 family `mk1 = 8`, `mk2 = 4`
   - lighter shell2 follow-up `mk1 = 4`, `mk2 = 2`
   - combined low-|k| `m_lowk = 4`
   - combined low-|k| `m_lowk = 2`
   - combined low-|k| `m_lowk = 1`
   - closed three-shell combined low-|k| `m_lowk = 0.5`
   - active widened combined low-|k| `m_lowk = 0.5`, `lowk_shells = 4`
   If the widened `lowk_shells = 4` family still does not flatten materially faster,
   widen the Fourier-accelerated low-|k| basis again and keep multi-seed-family diagnosis mandatory.
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
- `U2 = 1e3`, first shell2 family:
  - `mk1 = 8`, `mk2 = 4`
  - short-tune winner: `12 x 0.0003`
  - completed `warm=0` stage: `strong_drift`, `drift/span ≈ 0.762`
  - completed deeper stage: `strong_drift`, `drift/span ≈ 0.744`
  - aggressive follow-up still keeps `acceptance = 1`
- `U2 = 1e3`, lighter shell2 follow-up:
  - `mk1 = 4`, `mk2 = 2`
  - grid:
    - `12 x 0.0003`
    - `16 x 0.00025`
    - `20 x 0.0002`
    - `24 x 0.00015`
  - all candidates still kept `acceptance = 1` with huge positive `DeltaH`
- `U2 = 1e3`, current combined low-|k| probe:
  - first family `m_lowk = 4`
  - `12 x 0.0003` already dies with `acceptance = 0`
  - `16 x 0.00025` and `20 x 0.0002` both later end as `strong_drift`
  - second family `m_lowk = 2`
  - `12 x 0.00025` dies with `acceptance = 0`
  - `16 x 0.0002` is the closed moving-family reference
  - completed moving-family stages still remain `strong_drift`
  - third family `m_lowk = 1`
  - tune winner: `20 x 0.00012`
  - completed moving-family stages at `51001` / `52001` still remain `strong_drift`
  - partial three-family `256 / 128 / 32` stage also remains `strong_drift`
  - current active probe: `m_lowk = 0.5`
    - `16 x 0.0001`
    - `20 x 0.00008`
    - `24 x 0.00006`
    - `28 x 0.00005`
  - current fixed-family partial tune result:
    - recommended candidate: `20 x 0.00008`
    - `acceptance ≈ 0.995`
    - `tau_int(doubleOcc) ≈ 4.41`
    - `ESS/sec ≈ 0.591`
  - current retained-window follow-up:
    - `512 / 256 / 32`
    - `20 x 0.00008`
    - `50001,51001,52001` seed family block
- Working interpretation:
  - shell1-only preconditioning helped the very early trace shape but did not solve the retained-window drift
  - the first shell2 family also failed on completed stages
  - the lighter shell2 follow-up also failed to produce a genuinely selective tune window
  - the broader combined low-|k| family is the first path that produces both
    unstable points and nearby sane-`DeltaH` points
  - however, `m_lowk = 4`, `2`, and the first completed `m_lowk = 1` stages are all still too drifty
    on retained windows
  - the main question has therefore shifted from mere movement to whether a lighter broad low-|k|
    mass can reduce worst-repeat retained drift without falling back into frozen or zero-accept basins

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
