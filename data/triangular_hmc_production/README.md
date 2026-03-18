# Triangular HMC Production Data

This directory is reserved for active triangular-lattice HMC production bring-up.

Recommended layout:

- one subdirectory per stage rung
- each stage rung contains:
  - `runs/`
  - `production_stage.json`
  - `production_stage_cases.csv`
  - `production_stage_observables.csv`
  - `production_stage_runs.csv`
  - `production_stage_samples.csv`
  - `report/`

Typical examples:

- `l6_n1e3_u1_stage/`
- `l6_n1e5_u1e3_beta64_stage/`
- `l12_n1e5_u1e3_beta256_dtau1em3_stage/`

The active judgement for these stages is thermalization of `squareOcc` and `IPR`,
not local-vs-HMC agreement.
