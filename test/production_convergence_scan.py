#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path


def parse_list(text: str) -> list[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def run_command(cmd: list[str], cwd: Path) -> None:
    print("+", " ".join(shlex.quote(part) for part in cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch fixed-beta or fixed-dtau production HMC convergence scans.")
    parser.add_argument("--mode", choices=("beta", "dtau"), required=True)
    parser.add_argument("--lattice-type", default="triangular")
    parser.add_argument("--l-values", default="6")
    parser.add_argument("--nbos-values", default="1000")
    parser.add_argument("--u2-values", default="1")
    parser.add_argument("--beta", type=float, default=128.0, help="Used in dtau-scan mode.")
    parser.add_argument("--dtau", type=float, default=0.005, help="Used in beta-scan mode.")
    parser.add_argument("--beta-values", default="32,64,96,128,160,192")
    parser.add_argument("--dtau-values", default="0.01,0.008,0.005,0.004,0.002")
    parser.add_argument("--bins", type=int, default=1024)
    parser.add_argument("--thermal-cut", type=int, default=512)
    parser.add_argument("--warm", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--hmc-nfrog", type=int, default=8)
    parser.add_argument("--hmc-dt", type=float, default=0.02)
    parser.add_argument("--hmc-jitter", type=int, default=2)
    parser.add_argument("--hmc-mass", type=float, default=4.0)
    parser.add_argument("--hmc-mass-spatial-uniform", type=float, default=0.0)
    parser.add_argument("--hmc-mass-spatial-lowk", type=float, default=0.0)
    parser.add_argument("--hmc-mass-spatial-shell1", type=float, default=0.0)
    parser.add_argument("--hmc-mass-spatial-shell2", type=float, default=0.0)
    parser.add_argument("--root", default="data/triangular_hmc_production/convergence")
    parser.add_argument("--binary", default="src/BPQMC.out")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    root = (repo_root / args.root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    production_hmc = repo_root / "test" / "production_hmc.py"
    overview = repo_root / "test" / "render_hmc_production_overview.py"

    sweep_values = parse_list(args.beta_values if args.mode == "beta" else args.dtau_values)
    work_dirs: list[Path] = []
    for value in sweep_values:
        if args.mode == "beta":
            beta = value
            dtau = f"{args.dtau:.6g}"
            label = f"beta_{value}".replace(".", "p").replace("-", "m")
        else:
            beta = f"{args.beta:.6g}"
            dtau = value
            label = f"dtau_{value}".replace(".", "p").replace("-", "m")
        work_root = root / label
        work_dirs.append(work_root)
        cmd = [
            "/home/yyk/conda/envs/notebook/bin/python",
            str(production_hmc),
            "stage",
            "--lattice-type", args.lattice_type,
            "--l-values", args.l_values,
            "--nbos-values", args.nbos_values,
            "--u2-values", args.u2_values,
            "--beta", str(beta),
            "--dtau", str(dtau),
            "--bins", str(args.bins),
            "--thermal-cut", str(args.thermal_cut),
            "--warm", str(args.warm),
            "--repeats", str(args.repeats),
            "--hmc-nfrog", str(args.hmc_nfrog),
            "--hmc-dt", str(args.hmc_dt),
            "--hmc-jitter", str(args.hmc_jitter),
            "--hmc-mass", str(args.hmc_mass),
            "--hmc-mass-spatial-uniform", str(args.hmc_mass_spatial_uniform),
            "--hmc-mass-spatial-lowk", str(args.hmc_mass_spatial_lowk),
            "--hmc-mass-spatial-shell1", str(args.hmc_mass_spatial_shell1),
            "--hmc-mass-spatial-shell2", str(args.hmc_mass_spatial_shell2),
            "--work-root", str(work_root),
            "--binary", str((repo_root / args.binary).resolve()),
            "--stage-label", label,
        ]
        run_command(cmd, repo_root)

    overview_dir = root / "overview"
    cmd = [
        "/home/yyk/conda/envs/notebook/bin/python",
        str(overview),
        "--root", str(root),
        "--output-dir", str(overview_dir),
    ]
    run_command(cmd, repo_root)
    print(f"overview: {overview_dir / 'report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
