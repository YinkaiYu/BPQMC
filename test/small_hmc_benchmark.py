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
    RunConfig,
    ess_per_second,
    expected_confin_lines,
    integrated_autocorr_time,
    lag1_autocorr,
    parse_info_metrics,
    prepare_run_dir,
    read_complex_series_real,
    read_scalar_series,
    run_case,
    sample_stderr,
    series_mean,
    validate_confin_file,
)
from production_hmc import OBSERVABLES, compare_mode_stats


DEFAULT_NBOS = (10, 100, 1000)
DEFAULT_U2 = (1.0e-2, 1.0e-1, 1.0, 10.0, 100.0)
SERIES_OBSERVABLES = (
    "kinetic",
    "doubleOcc",
    "squareOcc",
    "nearestOcc",
    "IPR",
    "SF_Gamma",
    "SF_K",
    "PF_Gamma",
    "C3_Gamma",
    "dentot_Gamma",
    "denden_Gamma",
)


def parse_int_list(text: str) -> tuple[int, ...]:
    return tuple(int(item.strip()) for item in text.split(",") if item.strip())


def parse_float_list(text: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in text.split(",") if item.strip())


def parse_grid(text: str) -> list[tuple[int, float]]:
    rows: list[tuple[int, float]] = []
    for item in text.split(","):
        nfrog_text, dt_text = item.split(":")
        rows.append((int(nfrog_text), float(dt_text)))
    return rows


def parse_mass_grid(text: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in text.split(",") if item.strip())


def parse_block_grid(text: str) -> tuple[int, ...]:
    return tuple(int(item.strip()) for item in text.split(",") if item.strip())


def parse_site_block_grid(text: str) -> tuple[int, ...]:
    return tuple(int(item.strip()) for item in text.split(",") if item.strip())


def format_u2_tag(value: float) -> str:
    return f"{value:.0e}".replace("+", "")


def format_dt_tag(dt: float) -> str:
    return f"{dt:.9e}".replace("+", "").replace("-", "m").replace(".", "p")


def make_small_configs(args: argparse.Namespace) -> dict[str, RunConfig]:
    ltrot = int(round(args.beta / args.dtau))
    configs: dict[str, RunConfig] = {}
    for nbos in parse_int_list(args.nbos_values):
        for u2 in parse_float_list(args.u2_values):
            name = f"triangular_L{args.l_value}_N{nbos}_U2_{format_u2_tag(u2)}"
            configs[name] = RunConfig(
                name=name,
                lattice_type="triangular",
                rt=args.rt,
                ru1=args.ru1,
                ru2=u2,
                nbos=nbos,
                nlx=args.l_value,
                nly=args.l_value,
                ltrot=ltrot,
                beta=args.beta,
                nwrap=args.nwrap,
                nbin=args.bins,
                nsweep=args.sweeps,
                nthermal=args.thermal_cut,
                is_tau=False,
                is_warm=args.warm > 0,
                nwarm=max(args.warm, 0),
                ini_type=args.ini_type,
                ini_ampl=args.ini_ampl,
                ini_ham=args.ini_ham,
                ini_twist=args.ini_twist,
                imbalance=args.imbalance,
            )
    return configs


def load_series(run_dir: Path, obs_name: str) -> list[float]:
    reader = OBSERVABLES.get(obs_name, read_scalar_series)
    return reader(run_dir, obs_name)


def collect_means(run_dir: Path, thermal_cut: int) -> dict[str, float]:
    means: dict[str, float] = {}
    for obs_name, reader in OBSERVABLES.items():
        values = reader(run_dir, obs_name)
        means[obs_name] = series_mean(values[thermal_cut:])
    return means


def collect_perf(run_dir: Path, thermal_cut: int) -> dict[str, float]:
    values = read_scalar_series(run_dir, "doubleOcc")[thermal_cut:]
    info = parse_info_metrics(run_dir)
    metrics = {
        "tau_int_doubleOcc": integrated_autocorr_time(values),
        "lag1_doubleOcc": lag1_autocorr(values),
        "ess_per_sec_doubleOcc": ess_per_second(values, info["Tot_CPU_time"]),
        "cpu_time": info["Tot_CPU_time"],
    }
    if "HMC_DeltaH_mean" in info:
        metrics["hmc_deltaH_mean"] = info["HMC_DeltaH_mean"]
    if "HMC_DeltaH_abs_max" in info:
        metrics["hmc_deltaH_abs_max"] = info["HMC_DeltaH_abs_max"]
    return metrics


def run_or_collect_benchmark_repeat(
    args: argparse.Namespace,
    cfg: RunConfig,
    run_dir: Path,
    is_global: bool,
    repeat: int,
    seed: int,
    hmc_params: dict[str, object],
    binary: Path,
    run_missing: bool,
) -> tuple[dict[str, float], dict[str, float], dict[str, float], list[dict[str, object]]]:
    if run_missing:
        prepare_run_dir(
            run_dir,
            cfg,
            is_global=is_global,
            nfrog=int(hmc_params["nfrog"]),
            hmc_dt=float(hmc_params["hmc_dt"]),
            hmc_jitter=int(hmc_params.get("hmc_jitter", 0)) if is_global else 0,
            hmc_mass=float(hmc_params.get("hmc_mass", 1.0)) if is_global else 1.0,
            hmc_block_tau=int(hmc_params.get("hmc_block_tau", 0)) if is_global else 0,
            hmc_block_sites=int(hmc_params.get("hmc_block_sites", 0)) if is_global else 0,
            seed=seed,
            binary=binary,
            confin_from=case_confin_path(args.confin_root, cfg),
        )
        run_case(run_dir, np_ranks=args.np)
    elif not run_dir.exists():
        raise FileNotFoundError(f"missing benchmark run directory: {run_dir}")

    run_means = collect_means(run_dir, cfg.nthermal)
    perf = collect_perf(run_dir, cfg.nthermal)
    info = parse_info_metrics(run_dir)
    sample_rows: list[dict[str, object]] = []
    mode_name = "hmc" if is_global else "local"
    for obs_name in SERIES_OBSERVABLES:
        values = load_series(run_dir, obs_name)
        for sample_idx, value in enumerate(values):
            sample_rows.append(
                {
                    "name": cfg.name,
                    "Nbos": cfg.nbos,
                    "U2": cfg.ru2,
                    "thermal_cut": cfg.nthermal,
                    "mode": mode_name,
                    "repeat": repeat,
                    "observable": obs_name,
                    "sample_index": sample_idx,
                    "value": value,
                }
            )
    return run_means, perf, info, sample_rows


def aggregate_perf(rows: list[dict[str, float]]) -> dict[str, float]:
    keys = rows[0].keys()
    return {key: series_mean([row[key] for row in rows]) for key in keys}


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def case_confin_path(confin_root: str, cfg: RunConfig) -> Path | None:
    if not confin_root:
        return None
    path = Path(confin_root).resolve() / cfg.name / "confout.txt"
    ok, detail = validate_confin_file(path, cfg)
    if ok:
        return path
    print(
        f"  [warm-start] skip invalid confout for {cfg.name}: {detail}; "
        f"need {expected_confin_lines(cfg)} lines",
        flush=True,
    )
    return None


def run_local_seed(args: argparse.Namespace) -> int:
    configs = make_small_configs(args)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    rows = []
    for cfg in configs.values():
        print(f"[local-seed] case={cfg.name}", flush=True)
        cfg = cfg.with_sampling(nbin=args.seed_bins, nsweep=args.sweeps, is_warm=args.seed_warm > 0, nwarm=max(args.seed_warm, 0))
        run_dir = work_root / cfg.name
        prepare_run_dir(
            run_dir,
            cfg,
            is_global=False,
            nfrog=4,
            hmc_dt=1.0e-3,
            seed=args.seed_base,
            binary=binary,
        )
        run_case(run_dir, np_ranks=args.np)
        info = parse_info_metrics(run_dir)
        rows.append(
            {
                "name": cfg.name,
                "seed_bins": args.seed_bins,
                "seed_warm": args.seed_warm,
                "cpu_time": info["Tot_CPU_time"],
                "confout": str(run_dir / "confout.txt"),
            }
        )
        print(f"  [seeded] cpu={info['Tot_CPU_time']:.3f} confout={run_dir / 'confout.txt'}", flush=True)

    write_csv(work_root / "local_seed_runs.csv", rows, ["name", "seed_bins", "seed_warm", "cpu_time", "confout"])
    return 0


def run_tune(args: argparse.Namespace) -> int:
    configs = make_small_configs(args)
    grid = parse_grid(args.grid)
    masses = parse_mass_grid(args.hmc_mass_grid)
    blocks = parse_block_grid(args.hmc_block_grid)
    site_blocks = parse_site_block_grid(args.hmc_site_block_grid)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    case_rows = []
    repeat_rows = []
    cases = []
    for cfg in configs.values():
        print(f"[tune] case={cfg.name}", flush=True)
        scan_rows = []
        for block_tau in blocks:
            for block_sites in site_blocks:
                for mass in masses:
                    for nfrog, dt in grid:
                        print(
                            f"  [candidate] nfrog={nfrog} dt={dt:.6g} jitter={args.hmc_jitter} "
                            f"mass={mass:g} block_tau={block_tau} block_sites={block_sites}",
                            flush=True,
                        )
                        rows_for_choice = []
                        for repeat in range(args.repeats):
                            seed = (
                                args.seed_base
                                + 1000 * repeat
                                + 17 * nfrog
                                + int(round(1.0e8 * dt))
                                + int(round(10 * mass))
                                + 7919 * block_tau
                                + 1543 * block_sites
                            )
                            run_dir = (
                                work_root
                                / "runs"
                                / cfg.name
                                / f"nf{nfrog}_dt{format_dt_tag(dt)}_m{mass:g}_b{block_tau}_s{block_sites}_j{args.hmc_jitter}"
                                / f"rep{repeat}"
                            )
                            prepare_run_dir(
                                run_dir,
                                cfg,
                                is_global=True,
                                nfrog=nfrog,
                                hmc_dt=dt,
                                hmc_jitter=args.hmc_jitter,
                                hmc_mass=mass,
                                hmc_block_tau=block_tau,
                                hmc_block_sites=block_sites,
                                seed=seed,
                                binary=binary,
                                confin_from=case_confin_path(args.confin_root, cfg),
                            )
                            run_case(run_dir, np_ranks=args.np)
                            info = parse_info_metrics(run_dir)
                            values = read_scalar_series(run_dir, "doubleOcc")
                            row = {
                                "name": cfg.name,
                                "nfrog": nfrog,
                                "hmc_dt": dt,
                                "hmc_jitter": args.hmc_jitter,
                                "hmc_mass": mass,
                                "hmc_block_tau": block_tau,
                                "hmc_block_sites": block_sites,
                                "repeat": repeat,
                                "seed": seed,
                                "acceptance": info["Accept_HMC"],
                                "tau_int_doubleOcc": integrated_autocorr_time(values),
                                "lag1_doubleOcc": lag1_autocorr(values),
                                "ess_per_sec_doubleOcc": ess_per_second(values, info["Tot_CPU_time"]),
                                "cpu_time": info["Tot_CPU_time"],
                            }
                            row["stuck"] = int(row["acceptance"] <= args.min_run_accept or row["ess_per_sec_doubleOcc"] <= 0.0)
                            if "HMC_DeltaH_mean" in info:
                                row["hmc_deltaH_mean"] = info["HMC_DeltaH_mean"]
                            if "HMC_DeltaH_abs_max" in info:
                                row["hmc_deltaH_abs_max"] = info["HMC_DeltaH_abs_max"]
                            repeat_rows.append(row)
                            rows_for_choice.append(row)
                            print(
                                "    [repeat] repeat={repeat} accept={accept:.3f} tau={tau:.3f} ess/sec={ess:.4g} cpu={cpu:.3f}".format(
                                    repeat=repeat,
                                    accept=row["acceptance"],
                                    tau=row["tau_int_doubleOcc"],
                                    ess=row["ess_per_sec_doubleOcc"],
                                    cpu=row["cpu_time"],
                                ),
                                flush=True,
                            )
                        summary = {
                            "name": cfg.name,
                            "nfrog": nfrog,
                            "hmc_dt": dt,
                            "hmc_jitter": args.hmc_jitter,
                            "hmc_mass": mass,
                            "hmc_block_tau": block_tau,
                            "hmc_block_sites": block_sites,
                            "repeats": len(rows_for_choice),
                            "stuck_repeats": sum(int(row["stuck"]) for row in rows_for_choice),
                        }
                        for key in ("acceptance", "tau_int_doubleOcc", "lag1_doubleOcc", "ess_per_sec_doubleOcc", "cpu_time"):
                            values = [float(row[key]) for row in rows_for_choice]
                            summary[f"{key}_mean"] = series_mean(values)
                            summary[f"{key}_stderr"] = sample_stderr(values)
                        scan_rows.append(summary)
                        case_rows.append(summary)
                        print(
                            "  [summary] accept={accept:.3f} tau={tau:.3f} ess/sec={ess:.4g} stuck={stuck}".format(
                                accept=summary["acceptance_mean"],
                                tau=summary["tau_int_doubleOcc_mean"],
                                ess=summary["ess_per_sec_doubleOcc_mean"],
                                stuck=summary["stuck_repeats"],
                            ),
                            flush=True,
                        )
        scan_rows.sort(
            key=lambda row: (
                row["stuck_repeats"],
                0 if row["acceptance_mean"] > args.min_run_accept else 1,
                0 if row["ess_per_sec_doubleOcc_mean"] > 0.0 else 1,
                abs(row["lag1_doubleOcc_mean"]),
                row["tau_int_doubleOcc_mean"],
                -row["acceptance_mean"],
                -row["ess_per_sec_doubleOcc_mean"],
                row["hmc_block_tau"],
                row["hmc_block_sites"],
                row["nfrog"] * row["hmc_dt"],
            )
        )
        cases.append({"name": cfg.name, "config": asdict(cfg), "rows": scan_rows, "recommended": scan_rows[0]})
        print(
            "  [recommended] nfrog={nfrog} dt={dt:.6g} mass={mass:g} block_tau={block_tau} block_sites={block_sites} "
            "accept={accept:.3f} tau={tau:.3f} ess/sec={ess:.4g}".format(
                nfrog=scan_rows[0]["nfrog"],
                dt=scan_rows[0]["hmc_dt"],
                mass=scan_rows[0]["hmc_mass"],
                block_tau=scan_rows[0]["hmc_block_tau"],
                block_sites=scan_rows[0]["hmc_block_sites"],
                accept=scan_rows[0]["acceptance_mean"],
                tau=scan_rows[0]["tau_int_doubleOcc_mean"],
                ess=scan_rows[0]["ess_per_sec_doubleOcc_mean"],
            ),
            flush=True,
        )

    summary = {
        "mode": "small_tune",
        "cases": cases,
        "grid": [{"nfrog": nfrog, "hmc_dt": dt} for nfrog, dt in grid],
        "masses": masses,
        "blocks": blocks,
        "site_blocks": site_blocks,
    }
    (work_root / "small_tune.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        work_root / "small_tune.csv",
        case_rows,
        [
            "name",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "hmc_mass",
            "hmc_block_tau",
            "hmc_block_sites",
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
        ],
    )
    write_csv(
        work_root / "small_tune_repeats.csv",
        repeat_rows,
        [
            "name",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "hmc_mass",
            "hmc_block_tau",
            "hmc_block_sites",
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
            "hmc_block_tau": case["recommended"]["hmc_block_tau"],
            "hmc_block_sites": case["recommended"]["hmc_block_sites"],
        }
        for case in cases
    }
    (work_root / "recommended_hmc.json").write_text(json.dumps(tuned_map, indent=2), encoding="utf-8")
    return 0


def benchmark_core(args: argparse.Namespace, run_missing: bool) -> int:
    configs = make_small_configs(args)
    binary = Path(args.binary).resolve()
    work_root = Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    tuned_map = {}
    if args.hmc_json:
        tuned_map = json.loads(Path(args.hmc_json).read_text(encoding="utf-8"))

    cases = []
    case_rows = []
    observable_rows = []
    sample_rows = []
    overall_ok = True

    for cfg in configs.values():
        print(f"[benchmark] case={cfg.name}", flush=True)
        hmc_params = tuned_map.get(
            cfg.name,
            {
                "nfrog": args.hmc_nfrog,
                "hmc_dt": args.hmc_dt,
                "hmc_jitter": args.hmc_jitter,
                "hmc_mass": args.hmc_mass,
                "hmc_block_tau": args.hmc_block_tau,
                "hmc_block_sites": args.hmc_block_sites,
            },
        )
        summary: dict[str, dict[str, list[float]]] = {"local": defaultdict(list), "hmc": defaultdict(list)}
        perf_rows: dict[str, list[dict[str, float]]] = {"local": [], "hmc": []}
        hmc_accept = []
        mode_runs: dict[str, list[dict[str, object]]] = {"local": [], "hmc": []}
        local_stuck_repeats = 0
        hmc_stuck_repeats = 0

        for repeat in range(args.repeats):
            seed = args.seed_base + 100 * repeat
            for mode_name, is_global in (("local", False), ("hmc", True)):
                run_dir = work_root / "runs" / cfg.name / f"{mode_name}_rep{repeat}"
                run_means, perf, info, run_samples = run_or_collect_benchmark_repeat(
                    args,
                    cfg,
                    run_dir,
                    is_global,
                    repeat,
                    seed,
                    hmc_params,
                    binary,
                    run_missing,
                )
                for obs_name, value in run_means.items():
                    summary[mode_name][obs_name].append(value)
                perf_rows[mode_name].append(perf)
                if is_global:
                    hmc_accept.append(info["Accept_HMC"])
                    if info["Accept_HMC"] <= args.min_run_accept or perf["ess_per_sec_doubleOcc"] <= 0.0:
                        hmc_stuck_repeats += 1
                elif perf["ess_per_sec_doubleOcc"] <= 0.0:
                    local_stuck_repeats += 1
                mode_runs[mode_name].append(
                    {
                        "repeat": repeat,
                        "seed": seed,
                        "run_dir": str(run_dir),
                        "perf": perf,
                        "info": info,
                    }
                )
                sample_rows.extend(run_samples)
                print(
                    "  [repeat] mode={mode} repeat={repeat} accept={accept} tau={tau:.3f} ess/sec={ess:.4g} cpu={cpu:.3f}".format(
                        mode=mode_name,
                        repeat=repeat,
                        accept=f"{info['Accept_HMC']:.3f}" if is_global else "n/a",
                        tau=perf["tau_int_doubleOcc"],
                        ess=perf["ess_per_sec_doubleOcc"],
                        cpu=perf["cpu_time"],
                    ),
                    flush=True,
                )

        ok, obs_rows = compare_mode_stats(summary)
        if local_stuck_repeats > 0 or hmc_stuck_repeats > 0:
            ok = False
        overall_ok = overall_ok and ok
        perf_summary = {mode: aggregate_perf(rows) for mode, rows in perf_rows.items()}
        accept_mean = series_mean(hmc_accept)
        accept_err = sample_stderr(hmc_accept)
        speed_ratio = perf_summary["hmc"]["ess_per_sec_doubleOcc"] / max(perf_summary["local"]["ess_per_sec_doubleOcc"], 1.0e-12)
        case = {
            "name": cfg.name,
            "config": asdict(cfg),
            "hmc": hmc_params,
            "acceptance_mean": accept_mean,
            "acceptance_stderr": accept_err,
            "perf": {
                "local": perf_summary["local"],
                "hmc": perf_summary["hmc"],
                "speed_ratio_hmc_over_local": speed_ratio,
            },
            "local_stuck_repeats": local_stuck_repeats,
            "hmc_stuck_repeats": hmc_stuck_repeats,
            "runs": mode_runs,
            "observables": obs_rows,
            "passed": ok,
        }
        cases.append(case)
        case_rows.append(
            {
                "name": cfg.name,
                "L": cfg.nlx,
                "Nbos": cfg.nbos,
                "U2": cfg.ru2,
                "beta": cfg.beta,
                "dtau": cfg.beta / cfg.ltrot,
                "thermal_cut": cfg.nthermal,
                "warm": cfg.nwarm,
                "nfrog": hmc_params["nfrog"],
                "hmc_dt": hmc_params["hmc_dt"],
                "hmc_jitter": hmc_params.get("hmc_jitter", 0),
                "hmc_mass": hmc_params.get("hmc_mass", 1.0),
                "hmc_block_tau": hmc_params.get("hmc_block_tau", 0),
                "hmc_block_sites": hmc_params.get("hmc_block_sites", 0),
                "acceptance_mean": accept_mean,
                "acceptance_stderr": accept_err,
                "local_stuck_repeats": local_stuck_repeats,
                "hmc_stuck_repeats": hmc_stuck_repeats,
                "local_tau_int_doubleOcc": perf_summary["local"]["tau_int_doubleOcc"],
                "hmc_tau_int_doubleOcc": perf_summary["hmc"]["tau_int_doubleOcc"],
                "local_ess_per_sec_doubleOcc": perf_summary["local"]["ess_per_sec_doubleOcc"],
                "hmc_ess_per_sec_doubleOcc": perf_summary["hmc"]["ess_per_sec_doubleOcc"],
                "speed_ratio_hmc_over_local": speed_ratio,
                "passed": ok,
            }
        )
        for row in obs_rows:
            observable_rows.append({"name": cfg.name, "Nbos": cfg.nbos, "U2": cfg.ru2, **row})
        failing = [row["observable"] for row in obs_rows if not row["passed"]]
        print(
            "  [result] pass={passed} accept={accept:.3f} tau_local={tau_l:.3f} tau_hmc={tau_h:.3f} fail_obs={fail_obs}".format(
                passed=ok,
                accept=accept_mean,
                tau_l=perf_summary["local"]["tau_int_doubleOcc"],
                tau_h=perf_summary["hmc"]["tau_int_doubleOcc"],
                fail_obs=",".join(failing) if failing else "none",
            ),
            flush=True,
        )

    summary = {
        "mode": "small_benchmark",
        "description": "Triangular L=6 small-parameter local-vs-HMC correctness study",
        "cases": cases,
        "passed": overall_ok,
    }
    (work_root / "small_benchmark.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        work_root / "small_benchmark_cases.csv",
        case_rows,
        [
            "name",
            "L",
            "Nbos",
            "U2",
            "beta",
            "dtau",
            "thermal_cut",
            "warm",
            "nfrog",
            "hmc_dt",
            "hmc_jitter",
            "hmc_mass",
            "hmc_block_tau",
            "hmc_block_sites",
            "acceptance_mean",
            "acceptance_stderr",
            "local_stuck_repeats",
            "hmc_stuck_repeats",
            "local_tau_int_doubleOcc",
            "hmc_tau_int_doubleOcc",
            "local_ess_per_sec_doubleOcc",
            "hmc_ess_per_sec_doubleOcc",
            "speed_ratio_hmc_over_local",
            "passed",
        ],
    )
    write_csv(
        work_root / "small_benchmark_observables.csv",
        observable_rows,
        ["name", "Nbos", "U2", "observable", "local_mean", "local_err", "hmc_mean", "hmc_err", "abs_diff", "combined_err", "z_score", "passed"],
    )
    write_csv(
        work_root / "small_benchmark_samples.csv",
        sample_rows,
        ["name", "Nbos", "U2", "thermal_cut", "mode", "repeat", "observable", "sample_index", "value"],
    )
    return 0 if overall_ok else 1


def run_benchmark(args: argparse.Namespace) -> int:
    return benchmark_core(args, run_missing=True)


def collect_benchmark(args: argparse.Namespace) -> int:
    return benchmark_core(args, run_missing=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Small triangular HMC/local benchmark workflow.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--l-value", type=int, default=6)
        subparser.add_argument("--nbos-values", default="10,100,1000")
        subparser.add_argument("--u2-values", default="1e-2,1e-1,1e0,1e1,1e2")
        subparser.add_argument("--rt", type=float, default=1.0)
        subparser.add_argument("--ru1", type=float, default=0.0)
        subparser.add_argument("--beta", type=float, default=32.0)
        subparser.add_argument("--dtau", type=float, default=0.01)
        subparser.add_argument("--nwrap", type=int, default=16)
        subparser.add_argument("--bins", type=int, default=320)
        subparser.add_argument("--sweeps", type=int, default=1)
        subparser.add_argument("--thermal-cut", type=int, default=192)
        subparser.add_argument("--warm", type=int, default=256)
        subparser.add_argument("--ini-type", type=int, default=2)
        subparser.add_argument("--ini-ampl", type=float, default=0.1)
        subparser.add_argument("--ini-ham", type=int, default=5)
        subparser.add_argument("--ini-twist", type=float, default=1.0e-4)
        subparser.add_argument("--imbalance", type=float, default=0.0)
        subparser.add_argument("--np", type=int, default=1)
        subparser.add_argument("--seed-base", type=int, default=70001)
        subparser.add_argument("--binary", default=str(DEFAULT_BINARY))
        subparser.add_argument("--work-root", default="data/triangular_hmc_small_benchmark")
        subparser.add_argument("--confin-root", default="", help="Optional root directory with per-case confout.txt warm-starts.")

    seed_local = subparsers.add_parser("local-seed", help="Generate per-case local thermalized confout seeds.")
    add_common(seed_local)
    seed_local.add_argument("--seed-warm", type=int, default=256)
    seed_local.add_argument("--seed-bins", type=int, default=160)
    seed_local.set_defaults(func=run_local_seed)

    tune = subparsers.add_parser("tune", help="Scan HMC parameters for the L=6 triangular study.")
    add_common(tune)
    tune.add_argument("--grid", required=True, help="Comma-separated Nfrog:dt pairs.")
    tune.add_argument("--hmc-jitter", type=int, default=0)
    tune.add_argument("--hmc-mass-grid", default="1.0")
    tune.add_argument("--hmc-block-grid", default="0")
    tune.add_argument("--hmc-site-block-grid", default="0")
    tune.add_argument("--repeats", type=int, default=2)
    tune.add_argument("--min-run-accept", type=float, default=0.05)
    tune.set_defaults(func=run_tune)

    bench = subparsers.add_parser("benchmark", help="Run strict local-vs-HMC benchmarks for the L=6 triangular study.")
    add_common(bench)
    bench.add_argument("--repeats", type=int, default=6)
    bench.add_argument("--hmc-json", default="")
    bench.add_argument("--hmc-nfrog", type=int, default=8)
    bench.add_argument("--hmc-dt", type=float, default=6.0e-5)
    bench.add_argument("--hmc-jitter", type=int, default=0)
    bench.add_argument("--hmc-mass", type=float, default=1.0)
    bench.add_argument("--hmc-block-tau", type=int, default=0)
    bench.add_argument("--hmc-block-sites", type=int, default=0)
    bench.add_argument("--min-run-accept", type=float, default=0.0)
    bench.set_defaults(func=run_benchmark)

    collect = subparsers.add_parser(
        "collect-benchmark",
        help="Rebuild small benchmark summaries from an existing work-root without rerunning.",
    )
    add_common(collect)
    collect.add_argument("--repeats", type=int, default=6)
    collect.add_argument("--hmc-json", default="")
    collect.add_argument("--hmc-nfrog", type=int, default=8)
    collect.add_argument("--hmc-dt", type=float, default=6.0e-5)
    collect.add_argument("--hmc-jitter", type=int, default=0)
    collect.add_argument("--hmc-mass", type=float, default=1.0)
    collect.add_argument("--hmc-block-tau", type=int, default=0)
    collect.add_argument("--hmc-block-sites", type=int, default=0)
    collect.add_argument("--min-run-accept", type=float, default=0.0)
    collect.set_defaults(func=collect_benchmark)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
