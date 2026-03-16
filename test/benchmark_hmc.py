#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from hmc_tools import (
    DEFAULT_BINARY,
    RunConfig,
    default_parameter_sets,
    parse_info_metrics,
    prepare_run_dir,
    read_complex_series_real,
    read_scalar_series,
    run_case,
    sample_stderr,
    series_mean,
)


OBSERVABLES = {
    "density_up": read_scalar_series,
    "density_do": read_scalar_series,
    "kinetic": read_scalar_series,
    "doubleOcc": read_scalar_series,
    "squareOcc": read_scalar_series,
    "nearestOcc": read_scalar_series,
    "IPR": read_scalar_series,
    "num_up": read_scalar_series,
    "num_do": read_scalar_series,
    "numsquare_up": read_scalar_series,
    "numsquare_do": read_scalar_series,
    "C3breaking": read_scalar_series,
    "C3breaking_up": read_scalar_series,
    "SF_Gamma": read_complex_series_real,
    "SF_K": read_complex_series_real,
    "PF_Gamma": read_complex_series_real,
    "C3_Gamma": read_complex_series_real,
    "dentot_Gamma": read_complex_series_real,
    "denden_Gamma": read_complex_series_real,
}


TUNED_HMC = {
    "weak_u2": (8, 0.010),
    "mixed_u1_u2": (12, 0.400),
    "strong_u2": (10, 0.004),
}


def compare_modes(summary: dict[str, dict[str, list[float]]]) -> tuple[bool, list[str]]:
    ok = True
    lines: list[str] = []
    for obs_name in OBSERVABLES:
        local_values = summary["local"][obs_name]
        hmc_values = summary["hmc"][obs_name]
        local_mean = series_mean(local_values)
        hmc_mean = series_mean(hmc_values)
        local_err = sample_stderr(local_values)
        hmc_err = sample_stderr(hmc_values)
        combined = (local_err ** 2 + hmc_err ** 2) ** 0.5
        diff = abs(local_mean - hmc_mean)
        passed = diff <= combined if combined > 0.0 else diff == 0.0
        ok = ok and passed
        status = "PASS" if passed else "FAIL"
        lines.append(
            f"{status:4s} {obs_name:14s} "
            f"local={local_mean: .6e}±{local_err:.2e} "
            f"hmc={hmc_mean: .6e}±{hmc_err:.2e} "
            f"|diff|={diff:.2e} combined_err={combined:.2e}"
        )
    return ok, lines


def collect_run_means(run_dir: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    for obs_name, reader in OBSERVABLES.items():
        result[obs_name] = series_mean(reader(run_dir, obs_name))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark local updates against HMC.")
    parser.add_argument("--sets", default="weak_u2,mixed_u1_u2,strong_u2", help="Comma-separated parameter-set names")
    parser.add_argument("--repeats", type=int, default=4)
    parser.add_argument("--bins", type=int, default=24)
    parser.add_argument("--sweeps", type=int, default=6)
    parser.add_argument("--thermal-cut", type=int, default=-1)
    parser.add_argument("--seed-base", type=int, default=12001)
    parser.add_argument("--work-root", default=str(Path("/tmp") / "bpqmc_benchmark"))
    parser.add_argument("--binary", default=str(DEFAULT_BINARY))
    args = parser.parse_args()

    parameter_sets = default_parameter_sets()
    chosen_sets = [parameter_sets[name] for name in args.sets.split(",")]
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    overall_ok = True
    for cfg in chosen_sets:
        thermal_cut = args.bins // 2 if args.thermal_cut < 0 else args.thermal_cut
        cfg = cfg.with_sampling(nbin=args.bins, nsweep=args.sweeps, is_warm=False, nwarm=0)
        cfg = replace(cfg, nthermal=thermal_cut)
        summary: dict[str, dict[str, list[float]]] = {
            "local": defaultdict(list),
            "hmc": defaultdict(list),
        }
        hmc_accept = []
        nfrog, dt = TUNED_HMC[cfg.name]
        print(f"# set: {cfg.name}")
        print(f"# HMC parameters: Nfrog={nfrog} dt={dt:.4f}")
        print(f"# thermal cut: {thermal_cut} bins")
        for repeat in range(args.repeats):
            seed = args.seed_base + 100 * repeat
            for mode_name, is_global in (("local", False), ("hmc", True)):
                run_dir = work_root / f"{cfg.name}_{mode_name}_rep{repeat}"
                prepare_run_dir(run_dir, cfg, is_global=is_global, nfrog=nfrog, hmc_dt=dt, seed=seed, binary=binary)
                run_case(run_dir)
                means = collect_run_means(run_dir)
                for obs_name, value in means.items():
                    summary[mode_name][obs_name].append(value)
                if is_global:
                    info = parse_info_metrics(run_dir)
                    hmc_accept.append(info["Accept_HMC"])

        ok, lines = compare_modes(summary)
        overall_ok = overall_ok and ok
        for line in lines:
            print(line)
        accept_mean = series_mean(hmc_accept)
        accept_err = sample_stderr(hmc_accept)
        print(f"ACPT HMC acceptance mean={accept_mean:.3f} stderr={accept_err:.3f}")
        print(f"RESULT {'PASS' if ok else 'FAIL'}\n")

    return 0 if overall_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
