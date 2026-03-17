#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import asdict, replace
from pathlib import Path

from hmc_tools import (
    DEFAULT_BINARY,
    ess_per_second,
    integrated_autocorr_time,
    lag1_autocorr,
    parse_info_metrics,
    prepare_run_dir,
    production_parameter_sets,
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


def parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(item.strip()) for item in text.split(",") if item.strip())


def parse_float_list(text: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in text.split(",") if item.strip())


def parse_grid(text: str) -> list[tuple[int, float]]:
    grid: list[tuple[int, float]] = []
    for item in text.split(","):
        nfrog_text, dt_text = item.split(":")
        grid.append((int(nfrog_text), float(dt_text)))
    return grid


def format_dt_tag(dt: float) -> str:
    return f"{dt:.9e}".replace("+", "").replace("-", "m").replace(".", "p")


def collect_run_means(run_dir: Path, thermal_cut: int) -> dict[str, float]:
    means: dict[str, float] = {}
    for obs_name, reader in OBSERVABLES.items():
        values = reader(run_dir, obs_name)
        means[obs_name] = series_mean(values[thermal_cut:])
    return means


def summarize_perf(run_dir: Path, thermal_cut: int) -> dict[str, float]:
    values = read_scalar_series(run_dir, "doubleOcc")[thermal_cut:]
    info = parse_info_metrics(run_dir)
    summary = {
        "tau_int_doubleOcc": integrated_autocorr_time(values),
        "lag1_doubleOcc": lag1_autocorr(values),
        "ess_per_sec_doubleOcc": ess_per_second(values, info["Tot_CPU_time"]),
        "cpu_time": info["Tot_CPU_time"],
        "span_doubleOcc": max(values) - min(values) if values else 0.0,
    }
    if "HMC_DeltaH_mean" in info:
        summary["hmc_deltaH_mean"] = info["HMC_DeltaH_mean"]
    if "HMC_DeltaH_abs_max" in info:
        summary["hmc_deltaH_abs_max"] = info["HMC_DeltaH_abs_max"]
    return summary


def aggregate_tune_repeat_rows(repeat_rows: list[dict[str, float | int | str]]) -> dict[str, float | int | str]:
    summary: dict[str, float | int | str] = {
        "name": repeat_rows[0]["name"],
        "nfrog": repeat_rows[0]["nfrog"],
        "hmc_dt": repeat_rows[0]["hmc_dt"],
        "hmc_jitter": repeat_rows[0]["hmc_jitter"],
        "repeats": len(repeat_rows),
        "stuck_repeats": sum(int(row["stuck"]) for row in repeat_rows),
    }
    metric_keys = [
        "acceptance",
        "tau_int_doubleOcc",
        "lag1_doubleOcc",
        "ess_per_sec_doubleOcc",
        "cpu_time",
    ]
    for key in metric_keys:
        values = [float(row[key]) for row in repeat_rows]
        summary[f"{key}_mean"] = series_mean(values)
        summary[f"{key}_stderr"] = sample_stderr(values)
    if "hmc_deltaH_mean" in repeat_rows[0]:
        values = [float(row["hmc_deltaH_mean"]) for row in repeat_rows]
        summary["hmc_deltaH_mean_mean"] = series_mean(values)
        summary["hmc_deltaH_mean_stderr"] = sample_stderr(values)
    if "hmc_deltaH_abs_max" in repeat_rows[0]:
        values = [float(row["hmc_deltaH_abs_max"]) for row in repeat_rows]
        summary["hmc_deltaH_abs_max"] = max(values)
    return summary


def compare_mode_stats(summary: dict[str, dict[str, list[float]]]) -> tuple[bool, list[dict[str, float | str | bool]]]:
    rows: list[dict[str, float | str | bool]] = []
    overall_ok = True
    abs_tol_floor = 1.0e-12
    rel_tol = 1.0e-12
    sigma_cut = 2.0
    for obs_name in OBSERVABLES:
        local_values = summary["local"][obs_name]
        hmc_values = summary["hmc"][obs_name]
        local_mean = series_mean(local_values)
        hmc_mean = series_mean(hmc_values)
        local_err = sample_stderr(local_values)
        hmc_err = sample_stderr(hmc_values)
        abs_tol = max(abs_tol_floor, rel_tol * max(abs(local_mean), abs(hmc_mean)))
        combined_err = max((local_err ** 2 + hmc_err ** 2) ** 0.5, abs_tol)
        diff = abs(local_mean - hmc_mean)
        z_score = diff / combined_err if combined_err > 0.0 else 0.0
        passed = z_score <= sigma_cut
        overall_ok = overall_ok and passed
        rows.append(
            {
                "observable": obs_name,
                "local_mean": local_mean,
                "local_err": local_err,
                "hmc_mean": hmc_mean,
                "hmc_err": hmc_err,
                "abs_diff": diff,
                "combined_err": combined_err,
                "z_score": z_score,
                "passed": passed,
            }
        )
    return overall_ok, rows


def get_hmc_params(case_name: str, args: argparse.Namespace, tuned_map: dict[str, dict[str, float]]) -> tuple[int, float, int]:
    if case_name in tuned_map:
        params = tuned_map[case_name]
        return int(params["nfrog"]), float(params["hmc_dt"]), int(params.get("hmc_jitter", 0))
    return args.hmc_nfrog, args.hmc_dt, args.hmc_jitter


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_configs(args: argparse.Namespace) -> dict[str, object]:
    return production_parameter_sets(
        lattice_type=args.lattice_type,
        l_values=parse_int_list(args.l_values),
        nbos_values=parse_int_list(args.nbos_values),
        u2_values=parse_float_list(args.u2_values),
        rt=args.rt,
        ru1=args.ru1,
        beta=args.beta,
        dtau=args.dtau,
        nwrap=args.nwrap,
        nbin=args.bins,
        nsweep=args.sweeps,
        nthermal=args.thermal_cut,
        ini_type=args.ini_type,
        ini_ampl=args.ini_ampl,
        ini_ham=args.ini_ham,
        ini_twist=args.ini_twist,
        imbalance=args.imbalance,
    )


def run_tune(args: argparse.Namespace) -> int:
    configs = build_configs(args)
    grid = parse_grid(args.grid)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    tuned_cases = []
    csv_rows = []
    repeat_csv_rows = []
    for cfg in configs.values():
        cfg = cfg.with_sampling(is_warm=args.warm > 0, nwarm=max(args.warm, 0))
        case_rows = []
        for idx, (nfrog, dt) in enumerate(grid):
            repeat_rows = []
            for repeat in range(args.repeats):
                seed = args.seed_base + args.seed_step * idx + args.repeat_seed_step * repeat
                run_dir = work_root / "runs" / cfg.name / f"nf{nfrog}_dt{format_dt_tag(dt)}" / f"rep{repeat}"
                prepare_run_dir(
                    run_dir,
                    cfg,
                    is_global=True,
                    nfrog=nfrog,
                    hmc_dt=dt,
                    hmc_jitter=args.hmc_jitter,
                    seed=seed,
                    binary=binary,
                )
                run_case(run_dir, np_ranks=args.np)
                info = parse_info_metrics(run_dir)
                values = read_scalar_series(run_dir, "doubleOcc")
                row = {
                    "name": cfg.name,
                    "nfrog": nfrog,
                    "hmc_dt": dt,
                    "hmc_jitter": args.hmc_jitter,
                    "repeat": repeat,
                    "seed": seed,
                    "acceptance": info["Accept_HMC"],
                    "tau_int_doubleOcc": integrated_autocorr_time(values),
                    "lag1_doubleOcc": lag1_autocorr(values),
                    "ess_per_sec_doubleOcc": ess_per_second(values, info["Tot_CPU_time"]),
                    "cpu_time": info["Tot_CPU_time"],
                }
                row["stuck"] = int(info["Accept_HMC"] <= args.min_run_accept or row["ess_per_sec_doubleOcc"] <= 0.0)
                if "HMC_DeltaH_mean" in info:
                    row["hmc_deltaH_mean"] = info["HMC_DeltaH_mean"]
                if "HMC_DeltaH_abs_max" in info:
                    row["hmc_deltaH_abs_max"] = info["HMC_DeltaH_abs_max"]
                repeat_rows.append(row)
                repeat_csv_rows.append(row)

            summary_row = aggregate_tune_repeat_rows(repeat_rows)
            case_rows.append(summary_row)
            csv_rows.append(summary_row)

        healthy_rows = [row for row in case_rows if row["stuck_repeats"] == 0]
        candidates = [row for row in healthy_rows if args.accept_min <= row["acceptance_mean"] <= args.accept_max]
        if candidates:
            candidates.sort(key=lambda row: (-row["ess_per_sec_doubleOcc_mean"], row["tau_int_doubleOcc_mean"], row["nfrog"] * row["hmc_dt"]))
            best = candidates[0]
        elif healthy_rows:
            healthy_rows.sort(key=lambda row: (abs(row["acceptance_mean"] - 0.775), -row["ess_per_sec_doubleOcc_mean"], row["tau_int_doubleOcc_mean"]))
            best = healthy_rows[0]
        else:
            case_rows.sort(key=lambda row: (row["stuck_repeats"], abs(row["acceptance_mean"] - 0.775), -row["ess_per_sec_doubleOcc_mean"], row["tau_int_doubleOcc_mean"]))
            best = case_rows[0]
        tuned_cases.append(
            {
                "name": cfg.name,
                "config": asdict(cfg),
                "rows": case_rows,
                "recommended": best,
            }
        )

    summary = {
        "mode": "tune",
        "grid": [{"nfrog": nfrog, "hmc_dt": dt} for nfrog, dt in grid],
        "cases": tuned_cases,
    }
    (work_root / "production_tune.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        work_root / "production_tune.csv",
        csv_rows,
        [
            "name",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "repeats",
            "stuck_repeats",
            "acceptance_mean",
            "acceptance_stderr",
            "tau_int_doubleOcc_mean",
            "tau_int_doubleOcc_stderr",
            "lag1_doubleOcc_mean",
            "lag1_doubleOcc_stderr",
            "ess_per_sec_doubleOcc_mean",
            "ess_per_sec_doubleOcc_stderr",
            "cpu_time_mean",
            "cpu_time_stderr",
            "hmc_deltaH_mean_mean",
            "hmc_deltaH_mean_stderr",
            "hmc_deltaH_abs_max",
        ],
    )
    write_csv(
        work_root / "production_tune_repeats.csv",
        repeat_csv_rows,
        [
            "name",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "repeat",
            "seed",
            "stuck",
            "acceptance",
            "tau_int_doubleOcc",
            "lag1_doubleOcc",
            "ess_per_sec_doubleOcc",
            "cpu_time",
            "hmc_deltaH_mean",
            "hmc_deltaH_abs_max",
        ],
    )
    tuned_map = {
        case["name"]: {
            "nfrog": case["recommended"]["nfrog"],
            "hmc_dt": case["recommended"]["hmc_dt"],
            "hmc_jitter": case["recommended"]["hmc_jitter"],
        }
        for case in tuned_cases
    }
    (work_root / "recommended_hmc.json").write_text(json.dumps(tuned_map, indent=2), encoding="utf-8")
    return 0


def run_benchmark(args: argparse.Namespace) -> int:
    configs = build_configs(args)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    tuned_map: dict[str, dict[str, float]] = {}
    if args.hmc_json:
        tuned_map = json.loads(Path(args.hmc_json).read_text(encoding="utf-8"))

    case_rows = []
    observable_rows = []
    cases = []
    overall_ok = True

    for cfg in configs.values():
        if args.warm >= 0:
            cfg = cfg.with_sampling(is_warm=args.warm > 0, nwarm=max(args.warm, 0))
        nfrog, hmc_dt, jitter = get_hmc_params(cfg.name, args, tuned_map)
        summary: dict[str, dict[str, list[float]]] = {"local": defaultdict(list), "hmc": defaultdict(list)}
        perf: dict[str, dict[str, list[float]]] = {"local": defaultdict(list), "hmc": defaultdict(list)}
        hmc_accept = []
        hmc_stuck_repeats = 0
        local_stuck_repeats = 0
        for repeat in range(args.repeats):
            seed = args.seed_base + 100 * repeat
            for mode_name, is_global in (("local", False), ("hmc", True)):
                run_dir = work_root / "runs" / cfg.name / f"{mode_name}_rep{repeat}"
                prepare_run_dir(
                    run_dir,
                    cfg,
                    is_global=is_global,
                    nfrog=nfrog,
                    hmc_dt=hmc_dt,
                    hmc_jitter=jitter if is_global else 0,
                    seed=seed,
                    binary=binary,
                )
                run_case(run_dir, np_ranks=args.np)
                means = collect_run_means(run_dir, cfg.nthermal)
                for obs_name, value in means.items():
                    summary[mode_name][obs_name].append(value)
                perf_stats = summarize_perf(run_dir, cfg.nthermal)
                for key, value in perf_stats.items():
                    perf[mode_name][key].append(value)
                if is_global:
                    info = parse_info_metrics(run_dir)
                    hmc_accept.append(info["Accept_HMC"])
                    if info["Accept_HMC"] <= args.min_run_accept or perf_stats["ess_per_sec_doubleOcc"] <= 0.0:
                        hmc_stuck_repeats += 1
                else:
                    if perf_stats["ess_per_sec_doubleOcc"] <= 0.0:
                        local_stuck_repeats += 1

        ok, obs_rows = compare_mode_stats(summary)
        if hmc_stuck_repeats > 0 or local_stuck_repeats > 0:
            ok = False
        overall_ok = overall_ok and ok
        accept_mean = series_mean(hmc_accept)
        accept_err = sample_stderr(hmc_accept)
        local_speed = series_mean(perf["local"]["ess_per_sec_doubleOcc"])
        hmc_speed = series_mean(perf["hmc"]["ess_per_sec_doubleOcc"])
        speed_ratio = hmc_speed / max(local_speed, 1.0e-12)
        case_summary = {
            "name": cfg.name,
            "config": asdict(cfg),
            "hmc": {"nfrog": nfrog, "hmc_dt": hmc_dt, "hmc_jitter": jitter},
            "acceptance_mean": accept_mean,
            "acceptance_stderr": accept_err,
            "perf": {
                "local": {key: series_mean(values) for key, values in perf["local"].items()},
                "hmc": {key: series_mean(values) for key, values in perf["hmc"].items()},
                "speed_ratio_hmc_over_local": speed_ratio,
            },
            "local_stuck_repeats": local_stuck_repeats,
            "hmc_stuck_repeats": hmc_stuck_repeats,
            "observables": obs_rows,
            "passed": ok,
        }
        cases.append(case_summary)
        case_rows.append(
            {
                "name": cfg.name,
                "lattice_type": cfg.lattice_type,
                "Lx": cfg.nlx,
                "Ly": cfg.nly,
                "Nbos": cfg.nbos,
                "U2": cfg.ru2,
                "beta": cfg.beta,
                "dtau": cfg.beta / cfg.ltrot,
                "nfrog": nfrog,
                "hmc_dt": hmc_dt,
                "hmc_jitter": jitter,
                "acceptance_mean": accept_mean,
                "acceptance_stderr": accept_err,
                "local_stuck_repeats": local_stuck_repeats,
                "hmc_stuck_repeats": hmc_stuck_repeats,
                "local_tau_int_doubleOcc": case_summary["perf"]["local"]["tau_int_doubleOcc"],
                "hmc_tau_int_doubleOcc": case_summary["perf"]["hmc"]["tau_int_doubleOcc"],
                "hmc_deltaH_mean": case_summary["perf"]["hmc"].get("hmc_deltaH_mean", 0.0),
                "hmc_deltaH_abs_max": case_summary["perf"]["hmc"].get("hmc_deltaH_abs_max", 0.0),
                "local_ess_per_sec_doubleOcc": local_speed,
                "hmc_ess_per_sec_doubleOcc": hmc_speed,
                "speed_ratio_hmc_over_local": speed_ratio,
                "passed": ok,
            }
        )
        for row in obs_rows:
            observable_rows.append({"name": cfg.name, **row})

    summary = {"mode": "benchmark", "cases": cases, "passed": overall_ok}
    (work_root / "production_benchmark.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        work_root / "production_benchmark_cases.csv",
        case_rows,
        [
            "name",
            "lattice_type",
            "Lx",
            "Ly",
            "Nbos",
            "U2",
            "beta",
            "dtau",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "acceptance_mean",
            "acceptance_stderr",
            "local_stuck_repeats",
            "hmc_stuck_repeats",
            "local_tau_int_doubleOcc",
            "hmc_tau_int_doubleOcc",
            "hmc_deltaH_mean",
            "hmc_deltaH_abs_max",
            "local_ess_per_sec_doubleOcc",
            "hmc_ess_per_sec_doubleOcc",
            "speed_ratio_hmc_over_local",
            "passed",
        ],
    )
    write_csv(
        work_root / "production_benchmark_observables.csv",
        observable_rows,
        ["name", "observable", "local_mean", "local_err", "hmc_mean", "hmc_err", "abs_diff", "combined_err", "z_score", "passed"],
    )
    return 0 if overall_ok else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Production tuning and benchmark driver for HMC/local PQMC runs.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--lattice-type", default="triangular")
        subparser.add_argument("--l-values", default="12")
        subparser.add_argument("--nbos-values", default="100000,1000000,10000000")
        subparser.add_argument("--u2-values", default="100,1000,10000")
        subparser.add_argument("--rt", type=float, default=1.0)
        subparser.add_argument("--ru1", type=float, default=0.0)
        subparser.add_argument("--beta", type=float, default=256.0)
        subparser.add_argument("--dtau", type=float, default=1.0e-3)
        subparser.add_argument("--nwrap", type=int, default=32)
        subparser.add_argument("--bins", type=int, default=64)
        subparser.add_argument("--sweeps", type=int, default=1)
        subparser.add_argument("--thermal-cut", type=int, default=32)
        subparser.add_argument("--warm", type=int, default=32)
        subparser.add_argument("--ini-type", type=int, default=2)
        subparser.add_argument("--ini-ampl", type=float, default=0.1)
        subparser.add_argument("--ini-ham", type=int, default=5)
        subparser.add_argument("--ini-twist", type=float, default=1.0e-4)
        subparser.add_argument("--imbalance", type=float, default=0.0)
        subparser.add_argument("--np", type=int, default=1)
        subparser.add_argument("--seed-base", type=int, default=50001)
        subparser.add_argument("--binary", default=str(DEFAULT_BINARY))
        subparser.add_argument("--work-root", default=str(Path("/tmp") / "bpqmc_production"))

    tune = subparsers.add_parser("tune", help="Scan HMC parameters on a production parameter grid.")
    add_common(tune)
    tune.add_argument("--grid", required=True, help="Comma-separated Nfrog:dt pairs.")
    tune.add_argument("--hmc-jitter", type=int, default=0)
    tune.add_argument("--repeats", type=int, default=1)
    tune.add_argument("--seed-step", type=int, default=0, help="Increment added to the base seed for each grid point.")
    tune.add_argument("--repeat-seed-step", type=int, default=100, help="Increment added to the seed for each repeat of the same grid point.")
    tune.add_argument("--accept-min", type=float, default=0.70)
    tune.add_argument("--accept-max", type=float, default=0.85)
    tune.add_argument("--min-run-accept", type=float, default=0.05)
    tune.set_defaults(func=run_tune)

    bench = subparsers.add_parser("benchmark", help="Run direct local-vs-HMC benchmarks on a production parameter grid.")
    add_common(bench)
    bench.add_argument("--repeats", type=int, default=2)
    bench.add_argument("--hmc-json", default="", help="JSON file from the tune step with per-case HMC parameters.")
    bench.add_argument("--hmc-nfrog", type=int, default=10)
    bench.add_argument("--hmc-dt", type=float, default=0.012)
    bench.add_argument("--hmc-jitter", type=int, default=0)
    bench.add_argument("--min-run-accept", type=float, default=0.05)
    bench.set_defaults(func=run_benchmark)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
