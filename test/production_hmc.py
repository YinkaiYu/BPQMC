#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import asdict
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

TRACE_OBSERVABLES = ("squareOcc", "IPR", "doubleOcc", "nearestOcc")
GATE_OBSERVABLES = ("squareOcc", "IPR")
PLOT_OBSERVABLES = ("kinetic", "doubleOcc", "SF_Gamma", "SF_K", "PF_Gamma", "denden_Gamma")


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


def resolve_seed_bases(
    args: argparse.Namespace,
    *,
    requested_repeats: int,
    base_offset: int = 0,
) -> list[int]:
    if args.seed_base_values:
        return list(parse_int_list(args.seed_base_values))
    return [args.seed_base + base_offset + repeat * args.repeat_seed_step for repeat in range(requested_repeats)]


def format_dt_tag(dt: float) -> str:
    return f"{dt:.9e}".replace("+", "").replace("-", "m").replace(".", "p")


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_configs(args: argparse.Namespace) -> dict[str, object]:
    ini_type_values = parse_int_list(args.ini_type_values) if args.ini_type_values else ()
    ini_ampl_values = parse_float_list(args.ini_ampl_values) if args.ini_ampl_values else ()
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
        ini_type_values=ini_type_values,
        ini_ampl_values=ini_ampl_values,
        ini_ham=args.ini_ham,
        ini_twist=args.ini_twist,
        imbalance=args.imbalance,
    )


def load_hmc_map(path: str) -> dict[str, dict[str, float]]:
    if not path:
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


def get_hmc_params(
    case_name: str,
    args: argparse.Namespace,
    tuned_map: dict[str, dict[str, float]],
) -> tuple[int, float, int, float, float, float, int, float, float]:
    if case_name in tuned_map:
        params = tuned_map[case_name]
        return (
            int(params["nfrog"]),
            float(params["hmc_dt"]),
            int(params.get("hmc_jitter", 0)),
            float(params.get("hmc_mass", 1.0)),
            float(params.get("hmc_mass_spatial_uniform", 0.0)),
            float(params.get("hmc_mass_spatial_lowk", 0.0)),
            int(params.get("hmc_mass_spatial_lowk_shells", 3)),
            float(params.get("hmc_mass_spatial_shell1", 0.0)),
            float(params.get("hmc_mass_spatial_shell2", 0.0)),
        )
    return (
        args.hmc_nfrog,
        args.hmc_dt,
        args.hmc_jitter,
        args.hmc_mass,
        args.hmc_mass_spatial_uniform,
        args.hmc_mass_spatial_lowk,
        args.hmc_mass_spatial_lowk_shells,
        args.hmc_mass_spatial_shell1,
        args.hmc_mass_spatial_shell2,
    )


def hmc_env_overrides(
    hmc_mass_spatial_uniform: float,
    hmc_mass_spatial_lowk: float,
    hmc_mass_spatial_lowk_shells: int,
    hmc_mass_spatial_shell1: float,
    hmc_mass_spatial_shell2: float,
) -> dict[str, str]:
    return {
        "BPQMC_HMC_MASS_SPATIAL_UNIFORM": f"{hmc_mass_spatial_uniform:.16g}",
        "BPQMC_HMC_MASS_SPATIAL_LOWK": f"{hmc_mass_spatial_lowk:.16g}",
        "BPQMC_HMC_MASS_SPATIAL_LOWK_SHELLS": str(hmc_mass_spatial_lowk_shells),
        "BPQMC_HMC_MASS_SPATIAL_SHELL1": f"{hmc_mass_spatial_shell1:.16g}",
        "BPQMC_HMC_MASS_SPATIAL_SHELL2": f"{hmc_mass_spatial_shell2:.16g}",
    }


def collect_run_means(run_dir: Path, thermal_cut: int) -> dict[str, float]:
    means: dict[str, float] = {}
    for obs_name, reader in OBSERVABLES.items():
        values = reader(run_dir, obs_name)
        trimmed = values[thermal_cut:]
        means[obs_name] = series_mean(trimmed) if trimmed else values[-1]
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


def trace_window_stats(values: list[float], thermal_cut: int, window_frac: float, min_window: int) -> dict[str, float]:
    trimmed = values[thermal_cut:] if thermal_cut < len(values) else []
    window_pool = trimmed if trimmed else values
    if not window_pool:
        return {
            "n_samples": 0,
            "mean": 0.0,
            "stderr": 0.0,
            "start_mean": 0.0,
            "end_mean": 0.0,
            "drift": 0.0,
            "span": 0.0,
            "drift_over_span": 0.0,
        }
    window = min(len(window_pool), max(min_window, int(round(window_frac * len(window_pool)))))
    start = window_pool[:window]
    end = window_pool[-window:]
    overall_mean = series_mean(window_pool)
    start_mean = series_mean(start)
    end_mean = series_mean(end)
    span = max(window_pool) - min(window_pool)
    drift = end_mean - start_mean
    scale = max(span, 1.0e-12 * max(1.0, abs(overall_mean)))
    return {
        "n_samples": len(window_pool),
        "mean": overall_mean,
        "stderr": sample_stderr(window_pool),
        "start_mean": start_mean,
        "end_mean": end_mean,
        "drift": drift,
        "span": span,
        "drift_over_span": abs(drift) / scale,
    }


def aggregate_tune_repeat_rows(repeat_rows: list[dict[str, float | int | str]]) -> dict[str, float | int | str]:
    summary: dict[str, float | int | str] = {
        "name": repeat_rows[0]["name"],
        "nfrog": repeat_rows[0]["nfrog"],
        "hmc_dt": repeat_rows[0]["hmc_dt"],
        "hmc_jitter": repeat_rows[0]["hmc_jitter"],
        "hmc_mass": repeat_rows[0]["hmc_mass"],
        "hmc_mass_spatial_uniform": repeat_rows[0]["hmc_mass_spatial_uniform"],
        "hmc_mass_spatial_lowk": repeat_rows[0]["hmc_mass_spatial_lowk"],
        "hmc_mass_spatial_lowk_shells": repeat_rows[0]["hmc_mass_spatial_lowk_shells"],
        "hmc_mass_spatial_shell1": repeat_rows[0]["hmc_mass_spatial_shell1"],
        "hmc_mass_spatial_shell2": repeat_rows[0]["hmc_mass_spatial_shell2"],
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


def summarize_stage_observables(
    case_name: str,
    repeat_observables: dict[str, list[dict[str, float | int | str]]],
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for observable, obs_rows in repeat_observables.items():
        mean_values = [float(row["mean"]) for row in obs_rows]
        span_values = [float(row["span"]) for row in obs_rows]
        span_scale = max(series_mean(span_values), 1.0e-12 * max(1.0, abs(series_mean(mean_values))))
        rows.append(
            {
                "name": case_name,
                "observable": observable,
                "repeats": len(obs_rows),
                "mean": series_mean(mean_values),
                "stderr": sample_stderr(mean_values),
                "start_mean": series_mean([float(row["start_mean"]) for row in obs_rows]),
                "start_mean_stderr": sample_stderr([float(row["start_mean"]) for row in obs_rows]),
                "end_mean": series_mean([float(row["end_mean"]) for row in obs_rows]),
                "end_mean_stderr": sample_stderr([float(row["end_mean"]) for row in obs_rows]),
                "drift_mean": series_mean([float(row["drift"]) for row in obs_rows]),
                "drift_stderr": sample_stderr([float(row["drift"]) for row in obs_rows]),
                "span_mean": series_mean([float(row["span"]) for row in obs_rows]),
                "span_stderr": sample_stderr([float(row["span"]) for row in obs_rows]),
                "drift_over_span_mean": series_mean([float(row["drift_over_span"]) for row in obs_rows]),
                "drift_over_span_max": max(float(row["drift_over_span"]) for row in obs_rows),
                "repeat_mean_span": max(mean_values) - min(mean_values),
                "repeat_mean_span_over_span": (max(mean_values) - min(mean_values)) / span_scale,
                "n_samples_mean": series_mean([float(row["n_samples"]) for row in obs_rows]),
            }
        )
    return rows


def classify_stage_case(
    observable_rows: list[dict[str, float | int | str]],
    acceptance_mean: float,
    min_run_accept: float,
    stable_ratio: float,
    warning_ratio: float,
    stuck_repeats: int,
) -> tuple[str, str]:
    if stuck_repeats > 0 or acceptance_mean <= min_run_accept:
        return "stuck_or_invalid", "At least one repeat is stuck or acceptance is too low."
    gate_rows = [row for row in observable_rows if row["observable"] in GATE_OBSERVABLES]
    if not gate_rows:
        return "insufficient_gate_data", "Missing squareOcc/IPR traces."
    worst_mean_ratio = max(float(row["drift_over_span_mean"]) for row in gate_rows)
    worst_repeat_ratio = max(float(row["drift_over_span_max"]) for row in gate_rows)
    cross_repeat_ratio = max(float(row["repeat_mean_span_over_span"]) for row in gate_rows)
    if (
        worst_mean_ratio <= stable_ratio
        and worst_repeat_ratio <= stable_ratio
        and cross_repeat_ratio <= stable_ratio
    ):
        return "stable_window", "squareOcc/IPR drift is small and the retained repeat means agree."
    if (
        worst_mean_ratio <= warning_ratio
        and worst_repeat_ratio <= warning_ratio
        and cross_repeat_ratio <= warning_ratio
    ):
        return "slow_drift", "squareOcc/IPR is improving but at least one repeat still shows residual drift or plateau mismatch."
    return "strong_drift", "squareOcc/IPR still drifts strongly or different repeats settle onto clearly different retained windows."


def collect_stage_case(
    cfg,
    args: argparse.Namespace,
    binary: Path,
    tuned_map: dict[str, dict[str, float]],
    *,
    execute: bool,
    case_index: int,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    if args.warm >= 0:
        cfg = cfg.with_sampling(is_warm=args.warm > 0, nwarm=max(args.warm, 0))
    nfrog, hmc_dt, jitter, hmc_mass, hmc_mass_spatial_uniform, hmc_mass_spatial_lowk, hmc_mass_spatial_lowk_shells, hmc_mass_spatial_shell1, hmc_mass_spatial_shell2 = get_hmc_params(cfg.name, args, tuned_map)
    repeat_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    repeat_observables: dict[str, list[dict[str, float | int | str]]] = defaultdict(list)

    seed_bases = resolve_seed_bases(args, requested_repeats=args.repeats, base_offset=case_index * args.seed_step)
    for repeat, seed in enumerate(seed_bases):
        run_dir = Path(args.work_root).resolve() / "runs" / cfg.name / f"hmc_rep{repeat}"
        print(
            "  [{mode}] repeat={repeat} seed={seed} nfrog={nfrog} dt={dt:.6g} mass={mass:g} uniform_mass={uniform_mass:g} lowk_mass={lowk_mass:g} lowk_shells={lowk_shells} shell1_mass={shell1_mass:g} shell2_mass={shell2_mass:g}".format(
                mode="run" if execute else "reuse",
                repeat=repeat,
                seed=seed,
                nfrog=nfrog,
                dt=hmc_dt,
                mass=hmc_mass,
                uniform_mass=hmc_mass_spatial_uniform,
                lowk_mass=hmc_mass_spatial_lowk,
                lowk_shells=hmc_mass_spatial_lowk_shells,
                shell1_mass=hmc_mass_spatial_shell1,
                shell2_mass=hmc_mass_spatial_shell2,
            ),
            flush=True,
        )
        if execute:
            prepare_run_dir(
                run_dir,
                cfg,
                is_global=True,
                nfrog=nfrog,
                hmc_dt=hmc_dt,
                hmc_jitter=jitter,
                hmc_mass=hmc_mass,
                seed=seed,
                binary=binary,
            )
            run_case(
                run_dir,
                np_ranks=args.np,
                env_overrides=hmc_env_overrides(
                    hmc_mass_spatial_uniform,
                    hmc_mass_spatial_lowk,
                    hmc_mass_spatial_lowk_shells,
                    hmc_mass_spatial_shell1,
                    hmc_mass_spatial_shell2,
                ),
            )
        elif not run_dir.exists():
            print("    [skip] missing run directory", flush=True)
            continue

        try:
            info = parse_info_metrics(run_dir)
            perf_stats = summarize_perf(run_dir, cfg.nthermal)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            print(f"    [skip] incomplete run: {exc}", flush=True)
            continue
        repeat_row = {
            "name": cfg.name,
            "repeat": repeat,
            "seed": seed,
            "nfrog": nfrog,
            "hmc_dt": hmc_dt,
            "hmc_jitter": jitter,
            "hmc_mass": hmc_mass,
            "hmc_mass_spatial_uniform": hmc_mass_spatial_uniform,
            "hmc_mass_spatial_lowk": hmc_mass_spatial_lowk,
            "hmc_mass_spatial_lowk_shells": hmc_mass_spatial_lowk_shells,
            "hmc_mass_spatial_shell1": hmc_mass_spatial_shell1,
            "hmc_mass_spatial_shell2": hmc_mass_spatial_shell2,
            "acceptance": info["Accept_HMC"],
            "tau_int_doubleOcc": perf_stats["tau_int_doubleOcc"],
            "lag1_doubleOcc": perf_stats["lag1_doubleOcc"],
            "ess_per_sec_doubleOcc": perf_stats["ess_per_sec_doubleOcc"],
            "cpu_time": perf_stats["cpu_time"],
            "span_doubleOcc": perf_stats["span_doubleOcc"],
            "hmc_deltaH_mean": perf_stats.get("hmc_deltaH_mean", 0.0),
            "hmc_deltaH_abs_max": perf_stats.get("hmc_deltaH_abs_max", 0.0),
            "stuck": int(info["Accept_HMC"] <= args.min_run_accept or perf_stats["ess_per_sec_doubleOcc"] <= 0.0),
        }
        repeat_rows.append(repeat_row)

        for observable in TRACE_OBSERVABLES:
            values = OBSERVABLES[observable](run_dir, observable)
            obs_summary = trace_window_stats(values, cfg.nthermal, args.trace_window_frac, args.trace_min_window)
            repeat_observables[observable].append(
                {
                    "name": cfg.name,
                    "repeat": repeat,
                    "observable": observable,
                    **obs_summary,
                }
            )
            for sample_index, value in enumerate(values):
                sample_rows.append(
                    {
                        "name": cfg.name,
                        "repeat": repeat,
                        "observable": observable,
                        "sample_index": sample_index,
                        "post_thermal": int(sample_index >= cfg.nthermal),
                        "value": value,
                    }
                )

    if not repeat_rows:
        raise FileNotFoundError(f"No completed stage repeats found for case: {cfg.name}")

    observable_rows = summarize_stage_observables(cfg.name, repeat_observables)
    acceptance_values = [float(row["acceptance"]) for row in repeat_rows]
    tau_values = [float(row["tau_int_doubleOcc"]) for row in repeat_rows]
    ess_values = [float(row["ess_per_sec_doubleOcc"]) for row in repeat_rows]
    delta_h_values = [float(row["hmc_deltaH_mean"]) for row in repeat_rows]
    stuck_repeats = sum(int(row["stuck"]) for row in repeat_rows)
    status, status_detail = classify_stage_case(
        observable_rows,
        acceptance_mean=series_mean(acceptance_values),
        min_run_accept=args.min_run_accept,
        stable_ratio=args.stable_drift_ratio,
        warning_ratio=args.warning_drift_ratio,
        stuck_repeats=stuck_repeats,
    )

    gate_summary = {row["observable"]: row for row in observable_rows if row["observable"] in GATE_OBSERVABLES}
    case_row = {
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
        "hmc_mass": hmc_mass,
        "hmc_mass_spatial_uniform": hmc_mass_spatial_uniform,
        "hmc_mass_spatial_lowk": hmc_mass_spatial_lowk,
        "hmc_mass_spatial_lowk_shells": hmc_mass_spatial_lowk_shells,
        "hmc_mass_spatial_shell1": hmc_mass_spatial_shell1,
        "hmc_mass_spatial_shell2": hmc_mass_spatial_shell2,
        "repeats": len(repeat_rows),
        "requested_repeats": len(seed_bases),
        "missing_repeats": max(len(seed_bases) - len(repeat_rows), 0),
        "acceptance_mean": series_mean(acceptance_values),
        "acceptance_stderr": sample_stderr(acceptance_values),
        "tau_int_doubleOcc_mean": series_mean(tau_values),
        "tau_int_doubleOcc_stderr": sample_stderr(tau_values),
        "ess_per_sec_doubleOcc_mean": series_mean(ess_values),
        "ess_per_sec_doubleOcc_stderr": sample_stderr(ess_values),
        "hmc_deltaH_mean": series_mean(delta_h_values),
        "hmc_deltaH_abs_max": max(float(row["hmc_deltaH_abs_max"]) for row in repeat_rows),
        "stuck_repeats": stuck_repeats,
        "squareOcc_drift_ratio": float(gate_summary.get("squareOcc", {}).get("drift_over_span_mean", 0.0)),
        "squareOcc_drift_ratio_max": float(gate_summary.get("squareOcc", {}).get("drift_over_span_max", 0.0)),
        "squareOcc_repeat_span_ratio": float(gate_summary.get("squareOcc", {}).get("repeat_mean_span_over_span", 0.0)),
        "IPR_drift_ratio": float(gate_summary.get("IPR", {}).get("drift_over_span_mean", 0.0)),
        "IPR_drift_ratio_max": float(gate_summary.get("IPR", {}).get("drift_over_span_max", 0.0)),
        "IPR_repeat_span_ratio": float(gate_summary.get("IPR", {}).get("repeat_mean_span_over_span", 0.0)),
        "status": status,
        "status_detail": status_detail,
    }
    case_summary = {
        "name": cfg.name,
        "config": asdict(cfg),
        "hmc": {
            "nfrog": nfrog,
            "hmc_dt": hmc_dt,
            "hmc_jitter": jitter,
            "hmc_mass": hmc_mass,
            "hmc_mass_spatial_uniform": hmc_mass_spatial_uniform,
            "hmc_mass_spatial_lowk": hmc_mass_spatial_lowk,
            "hmc_mass_spatial_lowk_shells": hmc_mass_spatial_lowk_shells,
            "hmc_mass_spatial_shell1": hmc_mass_spatial_shell1,
            "hmc_mass_spatial_shell2": hmc_mass_spatial_shell2,
        },
        "requested_repeats": len(seed_bases),
        "missing_repeats": max(len(seed_bases) - len(repeat_rows), 0),
        "repeat_runs": repeat_rows,
        "observables": observable_rows,
        "status": status,
        "status_detail": status_detail,
    }
    return case_summary, [case_row], observable_rows, repeat_rows + sample_rows


def run_tune_or_collect(args: argparse.Namespace, *, execute: bool) -> int:
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
        print(
            f"[{'tune' if execute else 'collect-tune'}] case={cfg.name} mass={args.hmc_mass:g} "
            f"uniform_mass={args.hmc_mass_spatial_uniform:g} lowk_mass={args.hmc_mass_spatial_lowk:g} lowk_shells={args.hmc_mass_spatial_lowk_shells} shell1_mass={args.hmc_mass_spatial_shell1:g} "
            f"shell2_mass={args.hmc_mass_spatial_shell2:g} "
            f"jitter={args.hmc_jitter}",
            flush=True,
        )
        for idx, (nfrog, dt) in enumerate(grid):
            repeat_rows = []
            seed_bases = resolve_seed_bases(args, requested_repeats=args.repeats, base_offset=args.seed_step * idx)
            for repeat, seed in enumerate(seed_bases):
                run_dir = work_root / "runs" / cfg.name / f"nf{nfrog}_dt{format_dt_tag(dt)}" / f"rep{repeat}"
                print(
                    f"  [{'run' if execute else 'reuse'}] nfrog={nfrog} dt={dt:g} repeat={repeat} seed={seed}",
                    flush=True,
                )
                if execute:
                    prepare_run_dir(
                        run_dir,
                        cfg,
                        is_global=True,
                        nfrog=nfrog,
                        hmc_dt=dt,
                        hmc_jitter=args.hmc_jitter,
                        hmc_mass=args.hmc_mass,
                        seed=seed,
                        binary=binary,
                    )
                    run_case(
                        run_dir,
                        np_ranks=args.np,
                        env_overrides=hmc_env_overrides(
                            args.hmc_mass_spatial_uniform,
                            args.hmc_mass_spatial_lowk,
                            args.hmc_mass_spatial_lowk_shells,
                            args.hmc_mass_spatial_shell1,
                            args.hmc_mass_spatial_shell2,
                        ),
                    )
                elif not run_dir.exists():
                    print("    [skip] missing run directory", flush=True)
                    continue
                try:
                    info = parse_info_metrics(run_dir)
                    values = read_scalar_series(run_dir, "doubleOcc")
                except (FileNotFoundError, KeyError, ValueError) as exc:
                    print(f"    [skip] incomplete run: {exc}", flush=True)
                    continue
                if "Accept_HMC" not in info or "Tot_CPU_time" not in info:
                    missing_keys = [key for key in ("Accept_HMC", "Tot_CPU_time") if key not in info]
                    print(f"    [skip] incomplete run: missing {','.join(missing_keys)}", flush=True)
                    continue
                row = {
                    "name": cfg.name,
                    "nfrog": nfrog,
                    "hmc_dt": dt,
                    "hmc_jitter": args.hmc_jitter,
                    "hmc_mass": args.hmc_mass,
                    "hmc_mass_spatial_uniform": args.hmc_mass_spatial_uniform,
                    "hmc_mass_spatial_lowk": args.hmc_mass_spatial_lowk,
                    "hmc_mass_spatial_lowk_shells": args.hmc_mass_spatial_lowk_shells,
                    "hmc_mass_spatial_shell1": args.hmc_mass_spatial_shell1,
                    "hmc_mass_spatial_shell2": args.hmc_mass_spatial_shell2,
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
            if not repeat_rows:
                continue

            summary_row = aggregate_tune_repeat_rows(repeat_rows)
            summary_row["requested_repeats"] = len(seed_bases)
            summary_row["missing_repeats"] = max(len(seed_bases) - len(repeat_rows), 0)
            case_rows.append(summary_row)
            csv_rows.append(summary_row)

        if not case_rows:
            raise FileNotFoundError(f"No completed tune rows found for case: {cfg.name}")
        healthy_rows = [row for row in case_rows if row["stuck_repeats"] == 0]
        candidates = healthy_rows if healthy_rows else list(case_rows)
        candidates.sort(
            key=lambda row: (
                int(row.get("missing_repeats", 0)),
                row["stuck_repeats"],
                -row["ess_per_sec_doubleOcc_mean"],
                row["tau_int_doubleOcc_mean"],
                0 if args.accept_min <= row["acceptance_mean"] <= args.accept_max else 1,
                abs(row["acceptance_mean"] - 0.80),
                row["nfrog"] * row["hmc_dt"],
            )
        )
        best = candidates[0]
        recommended_viable = bool(healthy_rows)
        if recommended_viable:
            print(
                "  [recommended] nfrog={nfrog} dt={dt:.6g} acc={acc:.3f} tau={tau:.3f} ess={ess:.3f}".format(
                    nfrog=int(best["nfrog"]),
                    dt=float(best["hmc_dt"]),
                    acc=float(best["acceptance_mean"]),
                    tau=float(best["tau_int_doubleOcc_mean"]),
                    ess=float(best["ess_per_sec_doubleOcc_mean"]),
                ),
                flush=True,
            )
        else:
            print(
                "  [recommended] none viable; best fallback nfrog={nfrog} dt={dt:.6g} acc={acc:.3f} tau={tau:.3f} ess={ess:.3f}".format(
                    nfrog=int(best["nfrog"]),
                    dt=float(best["hmc_dt"]),
                    acc=float(best["acceptance_mean"]),
                    tau=float(best["tau_int_doubleOcc_mean"]),
                    ess=float(best["ess_per_sec_doubleOcc_mean"]),
                ),
                flush=True,
            )
        tuned_cases.append(
            {
                "name": cfg.name,
                "config": asdict(cfg),
                "rows": case_rows,
                "recommended": best,
                "recommended_viable": recommended_viable,
            }
        )

    summary = {
        "mode": "tune",
        "summary_source": "run" if execute else "collect",
        "grid": [
            {
                "nfrog": nfrog,
                "hmc_dt": dt,
                "hmc_mass": args.hmc_mass,
                "hmc_mass_spatial_uniform": args.hmc_mass_spatial_uniform,
                "hmc_mass_spatial_lowk": args.hmc_mass_spatial_lowk,
                "hmc_mass_spatial_lowk_shells": args.hmc_mass_spatial_lowk_shells,
                "hmc_mass_spatial_shell1": args.hmc_mass_spatial_shell1,
                "hmc_mass_spatial_shell2": args.hmc_mass_spatial_shell2,
            }
            for nfrog, dt in grid
        ],
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
            "hmc_mass",
            "hmc_mass_spatial_uniform",
            "hmc_mass_spatial_lowk",
            "hmc_mass_spatial_lowk_shells",
            "hmc_mass_spatial_shell1",
            "hmc_mass_spatial_shell2",
            "repeats",
            "requested_repeats",
            "missing_repeats",
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
            "hmc_mass",
            "hmc_mass_spatial_uniform",
            "hmc_mass_spatial_lowk",
            "hmc_mass_spatial_lowk_shells",
            "hmc_mass_spatial_shell1",
            "hmc_mass_spatial_shell2",
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
            "hmc_mass": case["recommended"]["hmc_mass"],
            "hmc_mass_spatial_uniform": case["recommended"]["hmc_mass_spatial_uniform"],
            "hmc_mass_spatial_lowk": case["recommended"].get("hmc_mass_spatial_lowk", 0.0),
            "hmc_mass_spatial_lowk_shells": case["recommended"].get("hmc_mass_spatial_lowk_shells", 3),
            "hmc_mass_spatial_shell1": case["recommended"].get("hmc_mass_spatial_shell1", 0.0),
            "hmc_mass_spatial_shell2": case["recommended"].get("hmc_mass_spatial_shell2", 0.0),
        }
        for case in tuned_cases
        if case.get("recommended_viable", True)
    }
    (work_root / "recommended_hmc.json").write_text(json.dumps(tuned_map, indent=2), encoding="utf-8")
    return 0


def run_tune(args: argparse.Namespace) -> int:
    return run_tune_or_collect(args, execute=True)


def collect_tune(args: argparse.Namespace) -> int:
    return run_tune_or_collect(args, execute=False)


def run_benchmark(args: argparse.Namespace) -> int:
    configs = build_configs(args)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    tuned_map = load_hmc_map(args.hmc_json)

    case_rows = []
    observable_rows = []
    cases = []
    overall_ok = True

    for case_index, cfg in enumerate(configs.values()):
        if args.warm >= 0:
            cfg = cfg.with_sampling(is_warm=args.warm > 0, nwarm=max(args.warm, 0))
        nfrog, hmc_dt, jitter, hmc_mass, hmc_mass_spatial_uniform, hmc_mass_spatial_lowk, hmc_mass_spatial_lowk_shells, hmc_mass_spatial_shell1, hmc_mass_spatial_shell2 = get_hmc_params(cfg.name, args, tuned_map)
        summary: dict[str, dict[str, list[float]]] = {"local": defaultdict(list), "hmc": defaultdict(list)}
        perf: dict[str, dict[str, list[float]]] = {"local": defaultdict(list), "hmc": defaultdict(list)}
        hmc_accept = []
        hmc_stuck_repeats = 0
        local_stuck_repeats = 0
        for repeat in range(args.repeats):
            seed = args.seed_base + case_index * args.seed_step + repeat * args.repeat_seed_step
            for mode_name, is_global in (("local", False), ("hmc", True)):
                run_dir = work_root / "runs" / cfg.name / f"{mode_name}_rep{repeat}"
                prepare_run_dir(
                    run_dir,
                    cfg,
                    is_global=is_global,
                    nfrog=nfrog,
                    hmc_dt=hmc_dt,
                    hmc_jitter=jitter if is_global else 0,
                    hmc_mass=hmc_mass if is_global else 1.0,
                    seed=seed,
                    binary=binary,
                )
                run_case(
                    run_dir,
                    np_ranks=args.np,
                    env_overrides=hmc_env_overrides(
                        hmc_mass_spatial_uniform if is_global else 0.0,
                        hmc_mass_spatial_lowk if is_global else 0.0,
                        hmc_mass_spatial_lowk_shells if is_global else 3,
                        hmc_mass_spatial_shell1 if is_global else 0.0,
                        hmc_mass_spatial_shell2 if is_global else 0.0,
                    ),
                )
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
            "hmc": {
                "nfrog": nfrog,
                "hmc_dt": hmc_dt,
                "hmc_jitter": jitter,
                "hmc_mass": hmc_mass,
                "hmc_mass_spatial_uniform": hmc_mass_spatial_uniform,
                "hmc_mass_spatial_lowk": hmc_mass_spatial_lowk,
                "hmc_mass_spatial_lowk_shells": hmc_mass_spatial_lowk_shells,
                "hmc_mass_spatial_shell1": hmc_mass_spatial_shell1,
                "hmc_mass_spatial_shell2": hmc_mass_spatial_shell2,
            },
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
                "hmc_mass": hmc_mass,
                "hmc_mass_spatial_uniform": hmc_mass_spatial_uniform,
                "hmc_mass_spatial_lowk": hmc_mass_spatial_lowk,
                "hmc_mass_spatial_lowk_shells": hmc_mass_spatial_lowk_shells,
                "hmc_mass_spatial_shell1": hmc_mass_spatial_shell1,
                "hmc_mass_spatial_shell2": hmc_mass_spatial_shell2,
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
            "hmc_mass",
            "hmc_mass_spatial_uniform",
            "hmc_mass_spatial_lowk",
            "hmc_mass_spatial_lowk_shells",
            "hmc_mass_spatial_shell1",
            "hmc_mass_spatial_shell2",
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


def run_stage_or_collect(args: argparse.Namespace, *, execute: bool) -> int:
    configs = build_configs(args)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    tuned_map = load_hmc_map(args.hmc_json)

    cases = []
    case_rows: list[dict[str, object]] = []
    observable_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []

    for case_index, cfg in enumerate(configs.values()):
        print(
            f"[{'stage' if execute else 'collect'}] case={cfg.name}",
            flush=True,
        )
        case_summary, case_row_list, case_obs_rows, row_bundle = collect_stage_case(
            cfg,
            args,
            binary,
            tuned_map,
            execute=execute,
            case_index=case_index,
        )
        cases.append(case_summary)
        case_rows.extend(case_row_list)
        observable_rows.extend(case_obs_rows)
        split_index = len(case_summary["repeat_runs"])
        run_rows.extend(row_bundle[:split_index])
        sample_rows.extend(row_bundle[split_index:])

    overall_status = "healthy" if all(case["status"] == "stable_window" and int(case.get("missing_repeats", 0)) == 0 for case in cases) else "needs_review"
    summary = {
        "mode": "stage",
        "stage_label": args.stage_label,
        "overall_status": overall_status,
        "gate_observables": list(GATE_OBSERVABLES),
        "trace_observables": list(TRACE_OBSERVABLES),
        "config_summary": {
            "lattice_type": args.lattice_type,
            "l_values": args.l_values,
            "nbos_values": args.nbos_values,
            "u2_values": args.u2_values,
            "beta": args.beta,
            "dtau": args.dtau,
            "bins": args.bins,
            "thermal_cut": args.thermal_cut,
            "warm": args.warm,
            "hmc_mass_spatial_uniform": args.hmc_mass_spatial_uniform,
            "hmc_mass_spatial_lowk": args.hmc_mass_spatial_lowk,
            "hmc_mass_spatial_lowk_shells": args.hmc_mass_spatial_lowk_shells,
            "hmc_mass_spatial_shell1": args.hmc_mass_spatial_shell1,
            "hmc_mass_spatial_shell2": args.hmc_mass_spatial_shell2,
        },
        "cases": cases,
    }
    (work_root / "production_stage.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        work_root / "production_stage_cases.csv",
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
            "hmc_mass",
            "hmc_mass_spatial_uniform",
            "hmc_mass_spatial_lowk",
            "hmc_mass_spatial_lowk_shells",
            "hmc_mass_spatial_shell1",
            "hmc_mass_spatial_shell2",
            "repeats",
            "requested_repeats",
            "missing_repeats",
            "acceptance_mean",
            "acceptance_stderr",
            "tau_int_doubleOcc_mean",
            "tau_int_doubleOcc_stderr",
            "ess_per_sec_doubleOcc_mean",
            "ess_per_sec_doubleOcc_stderr",
            "hmc_deltaH_mean",
            "hmc_deltaH_abs_max",
            "stuck_repeats",
            "squareOcc_drift_ratio",
            "squareOcc_drift_ratio_max",
            "squareOcc_repeat_span_ratio",
            "IPR_drift_ratio",
            "IPR_drift_ratio_max",
            "IPR_repeat_span_ratio",
            "status",
            "status_detail",
        ],
    )
    write_csv(
        work_root / "production_stage_observables.csv",
        observable_rows,
        [
            "name",
            "observable",
            "repeats",
            "mean",
            "stderr",
            "start_mean",
            "start_mean_stderr",
            "end_mean",
            "end_mean_stderr",
            "drift_mean",
            "drift_stderr",
            "span_mean",
            "span_stderr",
            "drift_over_span_mean",
            "drift_over_span_max",
            "repeat_mean_span",
            "repeat_mean_span_over_span",
            "n_samples_mean",
        ],
    )
    write_csv(
        work_root / "production_stage_runs.csv",
        run_rows,
        [
            "name",
            "repeat",
            "seed",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "hmc_mass",
            "hmc_mass_spatial_uniform",
            "hmc_mass_spatial_lowk",
            "hmc_mass_spatial_lowk_shells",
            "hmc_mass_spatial_shell1",
            "hmc_mass_spatial_shell2",
            "acceptance",
            "tau_int_doubleOcc",
            "lag1_doubleOcc",
            "ess_per_sec_doubleOcc",
            "cpu_time",
            "span_doubleOcc",
            "hmc_deltaH_mean",
            "hmc_deltaH_abs_max",
            "stuck",
        ],
    )
    write_csv(
        work_root / "production_stage_samples.csv",
        sample_rows,
        ["name", "repeat", "observable", "sample_index", "post_thermal", "value"],
    )
    return 0


def run_stage(args: argparse.Namespace) -> int:
    return run_stage_or_collect(args, execute=True)


def collect_stage(args: argparse.Namespace) -> int:
    return run_stage_or_collect(args, execute=False)


def render_summary(args: argparse.Namespace) -> int:
    from render_hmc_report import render_summary_file

    work_root = Path(args.work_root).resolve()
    if args.summary_json:
        summary_json = Path(args.summary_json).resolve()
    else:
        candidate = work_root / f"production_{args.summary_kind}.json"
        if not candidate.exists():
            raise FileNotFoundError(f"Could not find summary JSON: {candidate}")
        summary_json = candidate
    output_dir = Path(args.output_dir).resolve() if args.output_dir else work_root / "report"
    render_summary_file(summary_json, output_dir)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Production tuning, stage runs, and reporting for HMC/local PQMC.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(
        subparser: argparse.ArgumentParser,
        *,
        bins_default: int = 64,
        thermal_cut_default: int = 32,
        warm_default: int = 32,
    ) -> None:
        subparser.add_argument("--lattice-type", default="triangular")
        subparser.add_argument("--l-values", default="12")
        subparser.add_argument("--nbos-values", default="100000,1000000,10000000")
        subparser.add_argument("--u2-values", default="100,1000,10000")
        subparser.add_argument("--rt", type=float, default=1.0)
        subparser.add_argument("--ru1", type=float, default=0.0)
        subparser.add_argument("--beta", type=float, default=256.0)
        subparser.add_argument("--dtau", type=float, default=1.0e-3)
        subparser.add_argument("--nwrap", type=int, default=32)
        subparser.add_argument("--bins", type=int, default=bins_default)
        subparser.add_argument("--sweeps", type=int, default=1)
        subparser.add_argument("--thermal-cut", type=int, default=thermal_cut_default)
        subparser.add_argument("--warm", type=int, default=warm_default)
        subparser.add_argument("--ini-type", type=int, default=2)
        subparser.add_argument("--ini-ampl", type=float, default=0.1)
        subparser.add_argument("--ini-type-values", default="", help="Optional comma-separated initial-state type list.")
        subparser.add_argument("--ini-ampl-values", default="", help="Optional comma-separated initial-state amplitude list.")
        subparser.add_argument("--ini-ham", type=int, default=5)
        subparser.add_argument("--ini-twist", type=float, default=1.0e-4)
        subparser.add_argument("--imbalance", type=float, default=0.0)
        subparser.add_argument("--np", type=int, default=1)
        subparser.add_argument("--seed-base", type=int, default=50001)
        subparser.add_argument("--seed-base-values", default="", help="Optional comma-separated explicit seed-base list used verbatim for repeats; overrides --seed-base/--seed-step/--repeat-seed-step family generation.")
        subparser.add_argument("--seed-step", type=int, default=1000)
        subparser.add_argument("--repeat-seed-step", type=int, default=100)
        subparser.add_argument("--binary", default=str(DEFAULT_BINARY))
        subparser.add_argument("--work-root", default=str(Path("/tmp") / "bpqmc_production"))

    tune = subparsers.add_parser("tune", help="Scan HMC parameters on a production parameter grid.")
    add_common(tune)
    tune.add_argument("--grid", required=True, help="Comma-separated Nfrog:dt pairs.")
    tune.add_argument("--hmc-jitter", type=int, default=0)
    tune.add_argument("--hmc-mass", type=float, default=1.0)
    tune.add_argument("--hmc-mass-spatial-uniform", type=float, default=0.0)
    tune.add_argument("--hmc-mass-spatial-lowk", type=float, default=0.0)
    tune.add_argument("--hmc-mass-spatial-lowk-shells", type=int, default=3)
    tune.add_argument("--hmc-mass-spatial-shell1", type=float, default=0.0)
    tune.add_argument("--hmc-mass-spatial-shell2", type=float, default=0.0)
    tune.add_argument("--repeats", type=int, default=1)
    tune.add_argument("--accept-min", type=float, default=0.70)
    tune.add_argument("--accept-max", type=float, default=0.85)
    tune.add_argument("--min-run-accept", type=float, default=0.05)
    tune.set_defaults(func=run_tune)

    collect_tune_parser = subparsers.add_parser("collect-tune", help="Rebuild tune summaries from completed tune run directories.")
    add_common(collect_tune_parser)
    collect_tune_parser.add_argument("--grid", required=True, help="Comma-separated Nfrog:dt pairs.")
    collect_tune_parser.add_argument("--hmc-jitter", type=int, default=0)
    collect_tune_parser.add_argument("--hmc-mass", type=float, default=1.0)
    collect_tune_parser.add_argument("--hmc-mass-spatial-uniform", type=float, default=0.0)
    collect_tune_parser.add_argument("--hmc-mass-spatial-lowk", type=float, default=0.0)
    collect_tune_parser.add_argument("--hmc-mass-spatial-lowk-shells", type=int, default=3)
    collect_tune_parser.add_argument("--hmc-mass-spatial-shell1", type=float, default=0.0)
    collect_tune_parser.add_argument("--hmc-mass-spatial-shell2", type=float, default=0.0)
    collect_tune_parser.add_argument("--repeats", type=int, default=1)
    collect_tune_parser.add_argument("--accept-min", type=float, default=0.70)
    collect_tune_parser.add_argument("--accept-max", type=float, default=0.85)
    collect_tune_parser.add_argument("--min-run-accept", type=float, default=0.05)
    collect_tune_parser.set_defaults(func=collect_tune)

    bench = subparsers.add_parser("benchmark", help="Run direct local-vs-HMC benchmarks on a production parameter grid.")
    add_common(bench)
    bench.add_argument("--repeats", type=int, default=2)
    bench.add_argument("--hmc-json", default="", help="JSON file from the tune step with per-case HMC parameters.")
    bench.add_argument("--hmc-nfrog", type=int, default=10)
    bench.add_argument("--hmc-dt", type=float, default=0.012)
    bench.add_argument("--hmc-jitter", type=int, default=0)
    bench.add_argument("--hmc-mass", type=float, default=1.0)
    bench.add_argument("--hmc-mass-spatial-uniform", type=float, default=0.0)
    bench.add_argument("--hmc-mass-spatial-lowk", type=float, default=0.0)
    bench.add_argument("--hmc-mass-spatial-lowk-shells", type=int, default=3)
    bench.add_argument("--hmc-mass-spatial-shell1", type=float, default=0.0)
    bench.add_argument("--hmc-mass-spatial-shell2", type=float, default=0.0)
    bench.add_argument("--min-run-accept", type=float, default=0.05)
    bench.set_defaults(func=run_benchmark)

    def add_stage_args(subparser: argparse.ArgumentParser) -> None:
        add_common(subparser, bins_default=1024, thermal_cut_default=512, warm_default=512)
        subparser.add_argument("--repeats", type=int, default=2)
        subparser.add_argument("--stage-label", default="")
        subparser.add_argument("--hmc-json", default="", help="JSON file with per-case HMC parameters, usually recommended_hmc.json.")
        subparser.add_argument("--hmc-nfrog", type=int, default=10)
        subparser.add_argument("--hmc-dt", type=float, default=0.012)
        subparser.add_argument("--hmc-jitter", type=int, default=0)
        subparser.add_argument("--hmc-mass", type=float, default=1.0)
        subparser.add_argument("--hmc-mass-spatial-uniform", type=float, default=0.0)
        subparser.add_argument("--hmc-mass-spatial-lowk", type=float, default=0.0)
        subparser.add_argument("--hmc-mass-spatial-lowk-shells", type=int, default=3)
        subparser.add_argument("--hmc-mass-spatial-shell1", type=float, default=0.0)
        subparser.add_argument("--hmc-mass-spatial-shell2", type=float, default=0.0)
        subparser.add_argument("--min-run-accept", type=float, default=0.05)
        subparser.add_argument("--trace-window-frac", type=float, default=0.25)
        subparser.add_argument("--trace-min-window", type=int, default=16)
        subparser.add_argument("--stable-drift-ratio", type=float, default=0.25)
        subparser.add_argument("--warning-drift-ratio", type=float, default=0.50)

    stage = subparsers.add_parser("stage", help="Run long HMC stage jobs and store production trace diagnostics.")
    add_stage_args(stage)
    stage.set_defaults(func=run_stage)

    collect = subparsers.add_parser("collect", help="Rebuild production stage summaries from completed run directories.")
    add_stage_args(collect)
    collect.set_defaults(func=collect_stage)

    report = subparsers.add_parser("report", help="Render figures and Markdown from a production summary JSON.")
    report.add_argument("--work-root", default=str(Path("/tmp") / "bpqmc_production"))
    report.add_argument("--summary-json", default="")
    report.add_argument("--summary-kind", choices=("stage", "tune", "benchmark"), default="stage")
    report.add_argument("--output-dir", default="")
    report.set_defaults(func=render_summary)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
