#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from hmc_tools import (
    DEFAULT_BINARY,
    default_parameter_sets,
    integrated_autocorr_time,
    lag1_autocorr,
    parse_info_metrics,
    prepare_run_dir,
    read_scalar_series,
    run_case,
)


DEFAULT_GRID = [
    (4, 0.004),
    (4, 0.006),
    (8, 0.004),
    (8, 0.006),
    (8, 0.008),
    (10, 0.006),
]


def parse_grid(text: str) -> list[tuple[int, float]]:
    grid: list[tuple[int, float]] = []
    for item in text.split(","):
        nfrog_text, dt_text = item.split(":")
        grid.append((int(nfrog_text), float(dt_text)))
    return grid


def main() -> int:
    parser = argparse.ArgumentParser(description="Offline HMC tuning scan.")
    parser.add_argument("--set", dest="set_name", default="weak_u2", choices=sorted(default_parameter_sets().keys()))
    parser.add_argument("--grid", default="", help="Comma-separated Nfrog:dt pairs, for example 8:0.006,8:0.008")
    parser.add_argument("--bins", type=int, default=32)
    parser.add_argument("--sweeps", type=int, default=6)
    parser.add_argument("--warm", type=int, default=8)
    parser.add_argument("--seed", type=int, default=314159)
    parser.add_argument("--work-root", default=str(Path("/tmp") / "bpqmc_tune"))
    parser.add_argument("--binary", default=str(DEFAULT_BINARY))
    args = parser.parse_args()

    cfg = default_parameter_sets()[args.set_name].with_sampling(nbin=args.bins, nsweep=args.sweeps, is_warm=True, nwarm=args.warm)
    grid = DEFAULT_GRID if not args.grid else parse_grid(args.grid)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    rows = []
    for idx, (nfrog, dt) in enumerate(grid):
        run_dir = work_root / f"{cfg.name}_nf{nfrog}_dt{dt:.4f}".replace(".", "p")
        prepare_run_dir(run_dir, cfg, is_global=True, nfrog=nfrog, hmc_dt=dt, seed=args.seed + 1000 * idx, binary=binary)
        run_case(run_dir)
        info = parse_info_metrics(run_dir)
        acc = info["Accept_HMC"]
        double_occ = read_scalar_series(run_dir, "doubleOcc")
        tau_int = integrated_autocorr_time(double_occ)
        rho1 = lag1_autocorr(double_occ)
        rows.append((nfrog, dt, acc, tau_int, rho1))

    print(f"# tuning set: {cfg.name}")
    print("Nfrog dt acceptance tau_int(doubleOcc) lag1(doubleOcc)")
    for nfrog, dt, acc, tau_int, rho1 in rows:
        print(f"{nfrog:5d} {dt:7.4f} {acc:10.3f} {tau_int:18.3f} {rho1:16.3f}")

    candidates = [row for row in rows if 0.70 <= row[2] <= 0.85]
    if candidates:
        candidates.sort(key=lambda row: (row[3], row[0] * row[1]))
        best = candidates[0]
        print("\n# recommended candidate")
        print(f"Nfrog={best[0]} dt={best[1]:.4f} acceptance={best[2]:.3f} tau_int={best[3]:.3f} lag1={best[4]:.3f}")
    else:
        print("\n# no candidate hit the target acceptance window [0.70, 0.85]")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
