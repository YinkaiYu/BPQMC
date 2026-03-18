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

### Development Tune Scan

Scan HMC parameters on the built-in development sets:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/tune_hmc.py --set triangular_weak_u2 --grid '8:0.01,10:0.012,12:0.015' --bins 80 --sweeps 4 --warm 20
```

### Development Benchmark

Compare local and HMC on a development set:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
python3 test/benchmark_hmc.py --sets triangular_weak_u2 --repeats 2
```

The benchmark reports:

- bin-averaged observables after thermal cut
- combined-error `z` scores
- HMC acceptance
- `tau_int(doubleOcc)`
- `ESS/sec`

### Small Triangular Correctness Workflow

The repository now has a dedicated small-parameter triangular workflow for correctness-first HMC validation:

- `L = 6`
- `beta = 32`
- `dtau = 0.01`
- `U1 = 0`
- `Nbos = 10, 100, 1000`
- `U2 = 1e-2, 1e-1, 1e0, 1e1, 1e2`

All data from this workflow is meant to live under `data/triangular_hmc_small_benchmark/`.

1. Build the executable:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
make -C src
```

2. Generate validated local warm-start seeds for all 15 points:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/small_hmc_benchmark.py local-seed \
  --seed-warm 256 \
  --seed-bins 160 \
  --work-root data/triangular_hmc_small_benchmark/local_seed_grid_v2
```

3. Run a full-field HMC tune scan on the small grid:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/small_hmc_benchmark.py tune \
  --bins 64 \
  --warm 0 \
  --sweeps 1 \
  --grid '8:0.15,12:0.12,16:0.10,20:0.08,24:0.06' \
  --hmc-mass-grid 1 \
  --hmc-block-grid 0 \
  --hmc-site-block-grid 0 \
  --repeats 1 \
  --hmc-jitter 2 \
  --confin-root data/triangular_hmc_small_benchmark/local_seed_grid_v2 \
  --work-root data/triangular_hmc_small_benchmark/full_global_tune_grid_v2
```

4. Run a strict local-vs-HMC benchmark for one health point. Example: `Nbos=10, U2=1e-2`:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/small_hmc_benchmark.py benchmark \
  --nbos-values 10 \
  --u2-values 1e-2 \
  --bins 640 \
  --thermal-cut 384 \
  --warm 0 \
  --repeats 4 \
  --hmc-nfrog 16 \
  --hmc-dt 0.1 \
  --hmc-jitter 2 \
  --hmc-mass 1.0 \
  --hmc-block-tau 0 \
  --hmc-block-sites 0 \
  --confin-root data/triangular_hmc_small_benchmark/local_seed_grid_v2 \
  --work-root data/triangular_hmc_small_benchmark/strict_healthy_n10_u1em2_v2
```

5. Render a stage report from the strict benchmark directories that already pass:

```bash
cd /mnt/c/Users/Newton/Documents/LigroupIOP/2408_bosonSignProblem/code_BPQMC
/home/yyk/conda/envs/notebook/bin/python test/render_small_hmc_stage_report.py \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1em2_v2 \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1em1_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1e0_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n10_u1e1_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n100_u1em1_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n100_u1e0_v1 \
  data/triangular_hmc_small_benchmark/strict_healthy_n1000_u1em1_v1 \
  --tune-csv data/triangular_hmc_small_benchmark/full_global_tune_grid_v2/small_tune.csv \
  --output-dir data/triangular_hmc_small_benchmark/stage_report_v6a
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

- `test/small_hmc_stage_report.ipynb`

The notebook reads the CSV and PNG files from a rendered report directory. Set `REPORT_DIR` inside the notebook before running its cells.
`stage_samples.csv` now stores the key observable traces for every repeat, not only `repeat=0`, so you can do repeated thermal-cut checks directly inside the notebook or in a post-processing script.
The notebook exposes both `NBOS_SELECT` and `TRACE_REPEAT`, which lets you inspect one particle-number slice and one benchmark repeat without regenerating the report.
The current staged checkpoint with 9 strict-passing health points is `data/triangular_hmc_small_benchmark/stage_report_v6a`.
For the running handoff note covering validated health points, unresolved modes, tuning trends,
and production-oriented lessons, see `HMC-REFERENCE.md`.

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
