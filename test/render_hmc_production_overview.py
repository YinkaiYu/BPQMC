#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TRACE_OBSERVABLES = ("squareOcc", "IPR", "doubleOcc", "nearestOcc")
TREND_OBSERVABLES = ("squareOcc", "IPR")
GATE_OBSERVABLES = ("squareOcc", "IPR")
STATUS_RANK = {"stable_window": 0, "slow_drift": 1, "strong_drift": 2}


def scan_summaries(root: Path) -> tuple[list[Path], list[Path]]:
    stage_jsons = sorted(root.glob("*/production_stage.json"))
    tune_jsons = sorted(root.glob("*/production_tune.json"))
    return stage_jsons, tune_jsons


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_tag(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in value)


def case_sort_key(row: dict[str, object]) -> tuple[int, float, float, float, str]:
    return (
        int(row.get("Lx", 0)),
        float(row.get("Nbos", 0)),
        float(row.get("U2", 0)),
        float(row.get("beta", 0)),
        str(row.get("label", "")),
    )


def rel_link(target: Path, base_dir: Path) -> str:
    return os.path.relpath(target.resolve(), base_dir.resolve())


def linked_label(work_root: str, label: str, base_dir: Path) -> str:
    root = Path(work_root)
    for rel in ("report/report.md", "live_progress/report.md"):
        candidate = root / rel
        if candidate.exists():
            return f"[{label}]({rel_link(candidate, base_dir)})"
    return label


def load_stage_rows(stage_jsons: list[Path]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    case_rows: list[dict[str, object]] = []
    obs_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for path in stage_jsons:
        data = json.loads(path.read_text(encoding="utf-8"))
        work_root = path.parent
        samples_path = work_root / "production_stage_samples.csv"
        config_summary = data.get("config_summary", {})
        case = data["cases"][0]
        cfg = case["config"]
        case_key = case["name"]
        beta = float(cfg["beta"])
        dtau = float(cfg["beta"]) / float(cfg["ltrot"])
        gate_rows = [obs for obs in case["observables"] if obs["observable"] in GATE_OBSERVABLES]
        trace_counts: dict[int, int] = defaultdict(int)
        post_counts: dict[int, int] = defaultdict(int)
        if samples_path.exists():
            with samples_path.open("r", encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                for row in reader:
                    if row["observable"] not in TRACE_OBSERVABLES:
                        continue
                    repeat = int(row["repeat"])
                    post_thermal = int(row["post_thermal"])
                    sample_rows.append(
                        {
                            "label": work_root.name,
                            "stage_label": data["stage_label"],
                            "case": row["name"],
                            "repeat": repeat,
                            "observable": row["observable"],
                            "sample_index": int(row["sample_index"]),
                            "post_thermal": post_thermal,
                            "thermal_cut": int(config_summary.get("thermal_cut", 0)),
                            "warm": int(config_summary.get("warm", 0)),
                            "value": float(row["value"]),
                        }
                    )
                    if row["observable"] == "squareOcc":
                        trace_counts[repeat] += 1
                        post_counts[repeat] += post_thermal
        trace_count_values = list(trace_counts.values())
        post_count_values = list(post_counts.values())
        case_rows.append(
            {
                "label": work_root.name,
                "stage_label": data["stage_label"],
                "work_root": str(work_root),
                "case": case_key,
                "Lx": cfg["nlx"],
                "Ly": cfg["nly"],
                "Nbos": cfg["nbos"],
                "U2": cfg["ru2"],
                "beta": beta,
                "dtau": dtau,
                "bins": int(config_summary.get("bins", 0)),
                "thermal_cut": int(config_summary.get("thermal_cut", 0)),
                "warm": int(config_summary.get("warm", 0)),
                "nfrog": case["hmc"]["nfrog"],
                "hmc_dt": case["hmc"]["hmc_dt"],
                "hmc_mass": case["hmc"]["hmc_mass"],
                "hmc_mass_spatial_uniform": case["hmc"].get("hmc_mass_spatial_uniform", 0.0),
                "hmc_mass_spatial_lowk": case["hmc"].get("hmc_mass_spatial_lowk", 0.0),
                "hmc_mass_spatial_lowk_shells": case["hmc"].get("hmc_mass_spatial_lowk_shells", 3),
                "hmc_mass_spatial_midk": case["hmc"].get("hmc_mass_spatial_midk", 0.0),
                "hmc_mass_spatial_midk_shells": case["hmc"].get("hmc_mass_spatial_midk_shells", 0),
                "hmc_mass_spatial_shell1": case["hmc"].get("hmc_mass_spatial_shell1", 0.0),
                "hmc_mass_spatial_shell2": case["hmc"].get("hmc_mass_spatial_shell2", 0.0),
                "completed_repeats": len(case["repeat_runs"]),
                "requested_repeats": int(case.get("requested_repeats", len(case["repeat_runs"]))),
                "missing_repeats": int(case.get("missing_repeats", 0)),
                "trace_samples_mean": (sum(trace_count_values) / len(trace_count_values)) if trace_count_values else 0.0,
                "post_thermal_samples_mean": (sum(post_count_values) / len(post_count_values)) if post_count_values else 0.0,
                "trace_samples_max": max(trace_count_values) if trace_count_values else 0,
                "gate_drift_ratio_max": max(float(obs.get("drift_over_span_max", 0.0)) for obs in gate_rows) if gate_rows else 0.0,
                "gate_repeat_span_ratio": max(float(obs.get("repeat_mean_span_over_span", 0.0)) for obs in gate_rows) if gate_rows else 0.0,
                "acceptance_mean": sum(float(item["acceptance"]) for item in case["repeat_runs"]) / len(case["repeat_runs"]),
                "tau_int_doubleOcc_mean": sum(float(item["tau_int_doubleOcc"]) for item in case["repeat_runs"]) / len(case["repeat_runs"]),
                "ess_per_sec_doubleOcc_mean": sum(float(item["ess_per_sec_doubleOcc"]) for item in case["repeat_runs"]) / len(case["repeat_runs"]),
                "status": case["status"],
            }
        )
        for obs in case["observables"]:
            obs_rows.append(
                {
                    "label": work_root.name,
                    "stage_label": data["stage_label"],
                    "case": case_key,
                    "Lx": cfg["nlx"],
                    "Ly": cfg["nly"],
                    "Nbos": cfg["nbos"],
                    "U2": cfg["ru2"],
                    "beta": beta,
                    "dtau": dtau,
                    "observable": obs["observable"],
                    "mean": obs["mean"],
                    "stderr": obs["stderr"],
                    "drift_over_span": obs["drift_over_span_mean"],
                    "status": case["status"],
                }
            )
    return case_rows, obs_rows, sample_rows


def load_tune_rows(tune_jsons: list[Path]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    candidate_rows: list[dict[str, object]] = []
    recommended_rows: list[dict[str, object]] = []
    for path in tune_jsons:
        data = json.loads(path.read_text(encoding="utf-8"))
        work_root = path.parent
        case = data["cases"][0]
        cfg = case["config"]
        beta = float(cfg["beta"])
        dtau = float(cfg["beta"]) / float(cfg["ltrot"])
        for row in case["rows"]:
            candidate_rows.append(
                {
                    "label": work_root.name,
                    "work_root": str(work_root),
                    "case": case["name"],
                    "beta": beta,
                    "dtau": dtau,
                    "nfrog": row["nfrog"],
                    "hmc_dt": row["hmc_dt"],
                    "hmc_mass": row["hmc_mass"],
                    "hmc_mass_spatial_uniform": row.get("hmc_mass_spatial_uniform", 0.0),
                    "hmc_mass_spatial_lowk": row.get("hmc_mass_spatial_lowk", 0.0),
                    "hmc_mass_spatial_lowk_shells": row.get("hmc_mass_spatial_lowk_shells", 3),
                    "hmc_mass_spatial_midk": row.get("hmc_mass_spatial_midk", 0.0),
                    "hmc_mass_spatial_midk_shells": row.get("hmc_mass_spatial_midk_shells", 0),
                    "hmc_mass_spatial_shell1": row.get("hmc_mass_spatial_shell1", 0.0),
                    "hmc_mass_spatial_shell2": row.get("hmc_mass_spatial_shell2", 0.0),
                    "acceptance_mean": row["acceptance_mean"],
                    "tau_int_doubleOcc_mean": row["tau_int_doubleOcc_mean"],
                    "ess_per_sec_doubleOcc_mean": row["ess_per_sec_doubleOcc_mean"],
                    "hmc_deltaH_abs_max": row.get("hmc_deltaH_abs_max", 0.0),
                }
            )
        rec = case["recommended"]
        recommended_rows.append(
            {
                "label": work_root.name,
                "work_root": str(work_root),
                "case": case["name"],
                "beta": beta,
                "dtau": dtau,
                "recommended_viable": int(case.get("recommended_viable", True)),
                "nfrog": rec["nfrog"],
                "hmc_dt": rec["hmc_dt"],
                "hmc_mass": rec["hmc_mass"],
                "hmc_mass_spatial_uniform": rec.get("hmc_mass_spatial_uniform", 0.0),
                "hmc_mass_spatial_lowk": rec.get("hmc_mass_spatial_lowk", 0.0),
                "hmc_mass_spatial_lowk_shells": rec.get("hmc_mass_spatial_lowk_shells", 3),
                "hmc_mass_spatial_midk": rec.get("hmc_mass_spatial_midk", 0.0),
                "hmc_mass_spatial_midk_shells": rec.get("hmc_mass_spatial_midk_shells", 0),
                "hmc_mass_spatial_shell1": rec.get("hmc_mass_spatial_shell1", 0.0),
                "hmc_mass_spatial_shell2": rec.get("hmc_mass_spatial_shell2", 0.0),
                "acceptance_mean": rec["acceptance_mean"],
                "tau_int_doubleOcc_mean": rec["tau_int_doubleOcc_mean"],
                "ess_per_sec_doubleOcc_mean": rec["ess_per_sec_doubleOcc_mean"],
                "hmc_deltaH_abs_max": rec.get("hmc_deltaH_abs_max", 0.0),
            }
        )
    return candidate_rows, recommended_rows


def render_trend_plot(rows: list[dict[str, object]], observable: str, xkey: str, output_path: Path) -> None:
    obs_rows = [row for row in rows if row["observable"] == observable]
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in obs_rows:
        grouped[str(row["case"])].append(row)
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for case, items in sorted(grouped.items()):
        items.sort(key=lambda row: float(row[xkey]))
        ax.errorbar(
            [float(row[xkey]) for row in items],
            [float(row["mean"]) for row in items],
            yerr=[float(row["stderr"]) for row in items],
            marker="o",
            linewidth=1.5,
            capsize=3,
            label=case,
        )
    ax.set_xlabel(xkey)
    ax.set_ylabel(observable)
    ax.set_title(f"{observable} vs {xkey}")
    ax.grid(alpha=0.25)
    if len(grouped) > 1:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def render_case_trend_plot(
    rows: list[dict[str, object]],
    case: str,
    observable: str,
    xkey: str,
    output_path: Path,
) -> None:
    obs_rows = [row for row in rows if row["observable"] == observable and str(row["case"]) == case]
    if not obs_rows:
        return
    obs_rows.sort(key=lambda row: float(row[xkey]))
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.errorbar(
        [float(row[xkey]) for row in obs_rows],
        [float(row["mean"]) for row in obs_rows],
        yerr=[float(row["stderr"]) for row in obs_rows],
        marker="o",
        linewidth=1.6,
        capsize=3,
    )
    ax.set_xlabel(xkey)
    ax.set_ylabel(observable)
    ax.set_title(f"{case}: {observable} vs {xkey}")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def render_case_relative_trend_plot(
    rows: list[dict[str, object]],
    case: str,
    observable: str,
    xkey: str,
    output_path: Path,
) -> None:
    obs_rows = [row for row in rows if row["observable"] == observable and str(row["case"]) == case]
    if not obs_rows:
        return
    obs_rows.sort(key=lambda row: float(row[xkey]))
    baseline = float(obs_rows[0]["mean"])
    if baseline == 0.0:
        return
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    xvals = [float(row[xkey]) for row in obs_rows]
    yvals = [(float(row["mean"]) - baseline) / baseline for row in obs_rows]
    yerrs = [float(row["stderr"]) / abs(baseline) for row in obs_rows]
    ax.errorbar(
        xvals,
        yvals,
        yerr=yerrs,
        marker="o",
        linewidth=1.6,
        capsize=3,
    )
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
    ax.set_xlabel(xkey)
    ax.set_ylabel(f"relative change in {observable}")
    ax.set_title(f"{case}: relative {observable} vs {xkey}")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def render_trace_plot(sample_rows: list[dict[str, object]], case: str, observable: str, output_path: Path) -> None:
    rows = [row for row in sample_rows if row["case"] == case and row["observable"] == observable and int(row["repeat"]) == 0]
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["stage_label"])].append(row)
    if not grouped:
        return
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for stage_label, items in sorted(grouped.items()):
        items.sort(key=lambda row: int(row["sample_index"]))
        first = items[0]
        legend_label = f"{stage_label} [warm={int(first['warm'])}, cut={int(first['thermal_cut'])}]"
        line = ax.plot(
            [int(row["sample_index"]) for row in items],
            [float(row["value"]) for row in items],
            linewidth=1.2,
            label=legend_label,
        )[0]
        post_rows = [row for row in items if int(row["post_thermal"]) == 1]
        if post_rows:
            cut_index = min(int(row["sample_index"]) for row in post_rows)
            ax.axvline(cut_index, color=line.get_color(), linestyle="--", linewidth=0.9, alpha=0.45)
    ax.set_xlabel("sample index")
    ax.set_ylabel(observable)
    ax.set_title(f"{case}: {observable} trace (repeat 0, dashed = configured post-thermal start)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def render_recommended_plot(rows: list[dict[str, object]], metric: str, xkey: str, output_path: Path) -> None:
    rows = [row for row in rows if int(row.get("recommended_viable", 1)) == 1]
    if not rows:
        return
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["case"])].append(row)
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for case, items in sorted(grouped.items()):
        items.sort(key=lambda row: float(row[xkey]))
        ax.plot(
            [float(row[xkey]) for row in items],
            [float(row[metric]) for row in items],
            marker="o",
            linewidth=1.5,
            label=case,
        )
    ax.set_xlabel(xkey)
    ax.set_ylabel(metric)
    ax.set_title(f"recommended {metric} vs {xkey}")
    ax.grid(alpha=0.25)
    if len(grouped) > 1:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def select_representative_stage_rows(
    stage_cases: list[dict[str, object]],
    stage_obs: list[dict[str, object]],
) -> list[dict[str, object]]:
    best_by_key: dict[tuple[str, float, float], dict[str, object]] = {}
    for row in stage_cases:
        key = (str(row["case"]), float(row["beta"]), float(row["dtau"]))
        current = best_by_key.get(key)
        candidate_score = (
            STATUS_RANK.get(str(row["status"]), 99),
            int(row["missing_repeats"]),
            -int(row["completed_repeats"]),
            float(row["gate_repeat_span_ratio"]),
            float(row["gate_drift_ratio_max"]),
            str(row["label"]),
        )
        if current is None:
            best_by_key[key] = row
            continue
        current_score = (
            STATUS_RANK.get(str(current["status"]), 99),
            int(current["missing_repeats"]),
            -int(current["completed_repeats"]),
            float(current["gate_repeat_span_ratio"]),
            float(current["gate_drift_ratio_max"]),
            str(current["label"]),
        )
        if candidate_score < current_score:
            best_by_key[key] = row
    keep_labels = {str(row["label"]) for row in best_by_key.values()}
    return [row for row in stage_obs if str(row["label"]) in keep_labels]


def select_representative_stage_cases(stage_cases: list[dict[str, object]]) -> list[dict[str, object]]:
    best_by_case: dict[str, dict[str, object]] = {}
    for row in stage_cases:
        key = str(row["case"])
        current = best_by_case.get(key)
        candidate_score = (
            STATUS_RANK.get(str(row["status"]), 99),
            int(row["missing_repeats"]),
            float(row["gate_repeat_span_ratio"]),
            float(row["gate_drift_ratio_max"]),
            -float(row["post_thermal_samples_mean"]),
            -float(row["beta"]),
            float(row["dtau"]),
            str(row["label"]),
        )
        if current is None:
            best_by_case[key] = row
            continue
        current_score = (
            STATUS_RANK.get(str(current["status"]), 99),
            int(current["missing_repeats"]),
            float(current["gate_repeat_span_ratio"]),
            float(current["gate_drift_ratio_max"]),
            -float(current["post_thermal_samples_mean"]),
            -float(current["beta"]),
            float(current["dtau"]),
            str(current["label"]),
        )
        if candidate_score < current_score:
            best_by_case[key] = row
    return sorted(best_by_case.values(), key=case_sort_key)


def build_case_convergence_rows(
    stage_cases: list[dict[str, object]],
    stage_obs: list[dict[str, object]],
) -> list[dict[str, object]]:
    representative_obs = select_representative_stage_rows(stage_cases, stage_obs)
    case_meta: dict[str, dict[str, object]] = {}
    for row in stage_cases:
        case = str(row["case"])
        current = case_meta.get(case)
        if current is None:
            case_meta[case] = {
                "Lx": row["Lx"],
                "Ly": row["Ly"],
                "Nbos": row["Nbos"],
                "U2": row["U2"],
                "beta_min": float(row["beta"]),
                "beta_max": float(row["beta"]),
                "dtau_min": float(row["dtau"]),
                "dtau_max": float(row["dtau"]),
            }
            continue
        current["beta_min"] = min(float(current["beta_min"]), float(row["beta"]))
        current["beta_max"] = max(float(current["beta_max"]), float(row["beta"]))
        current["dtau_min"] = min(float(current["dtau_min"]), float(row["dtau"]))
        current["dtau_max"] = max(float(current["dtau_max"]), float(row["dtau"]))
    out: list[dict[str, object]] = []
    for case in sorted({str(row["case"]) for row in representative_obs}):
        case_rows = [row for row in representative_obs if str(row["case"]) == case and row["observable"] in TREND_OBSERVABLES]
        for observable in TREND_OBSERVABLES:
            obs_rows = [row for row in case_rows if row["observable"] == observable]
            if not obs_rows:
                continue
            obs_rows.sort(key=lambda row: (float(row["beta"]), float(row["dtau"]), str(row["label"])))
            first = obs_rows[0]
            last = obs_rows[-1]
            means = [float(row["mean"]) for row in obs_rows]
            rel_spread = 0.0
            if means and means[0] != 0.0:
                rel_spread = (max(means) - min(means)) / abs(means[0])
            rel_last = 0.0
            if float(first["mean"]) != 0.0:
                rel_last = (float(last["mean"]) - float(first["mean"])) / abs(float(first["mean"]))
            meta = case_meta.get(case, {})
            out.append(
                {
                    "case": case,
                    "Lx": meta.get("Lx", ""),
                    "Ly": meta.get("Ly", ""),
                    "Nbos": meta.get("Nbos", ""),
                    "U2": meta.get("U2", ""),
                    "observable": observable,
                    "n_points": len(obs_rows),
                    "beta_min": meta.get("beta_min", float(first["beta"])),
                    "beta_max": meta.get("beta_max", float(last["beta"])),
                    "dtau_min": meta.get("dtau_min", min(float(item["dtau"]) for item in obs_rows)),
                    "dtau_max": meta.get("dtau_max", max(float(item["dtau"]) for item in obs_rows)),
                    "first_label": first["label"],
                    "last_label": last["label"],
                    "first_mean": float(first["mean"]),
                    "last_mean": float(last["mean"]),
                    "max_mean": max(means),
                    "min_mean": min(means),
                    "rel_spread": rel_spread,
                    "rel_last_minus_first": rel_last,
                }
            )
    return out


def scan_live_reports(root: Path) -> list[dict[str, str]]:
    def read_trace(path: Path) -> list[float]:
        if not path.exists():
            return []
        values: list[float] = []
        with path.open("r", encoding="utf-8", errors="ignore") as stream:
            for line in stream:
                stripped = line.strip()
                if stripped:
                    values.append(float(stripped.split()[0]))
        return values

    def partial_trace_stats(values: list[float], thermal_cut: int) -> tuple[int, str, float]:
        if not values:
            return 0, "empty", 0.0
        if thermal_cut > 0 and len(values) > thermal_cut:
            tail = values[thermal_cut:]
            mode = "post-cut"
        else:
            tail = values
            mode = "pre-cut"
        window = max(16, len(tail) // 4)
        start = tail[:window]
        end = tail[-window:]
        span = max(tail) - min(tail) if len(tail) > 1 else 0.0
        drift = (sum(end) / len(end)) - (sum(start) / len(start))
        ratio = abs(drift) / span if span > 0.0 else 0.0
        return len(values), mode, ratio

    def summarize_run_dirs(run_dirs: list[Path], thermal_cut: int) -> tuple[int, str, float, float]:
        longest_trace = 0
        window_modes: set[str] = set()
        square_ratio = 0.0
        ipr_ratio = 0.0
        for run_dir in run_dirs:
            run_square = read_trace(run_dir / "squareOcc")
            run_ipr = read_trace(run_dir / "IPR")
            run_len, run_mode, run_square_ratio = partial_trace_stats(run_square, thermal_cut)
            _, run_ipr_mode, run_ipr_ratio = partial_trace_stats(run_ipr, thermal_cut)
            longest_trace = max(longest_trace, run_len)
            window_modes.add(run_mode)
            window_modes.add(run_ipr_mode)
            square_ratio = max(square_ratio, run_square_ratio)
            ipr_ratio = max(ipr_ratio, run_ipr_ratio)
        if not window_modes:
            window_mode = "empty"
        elif len(window_modes) == 1:
            window_mode = next(iter(window_modes))
        else:
            window_mode = "mixed"
        return longest_trace, window_mode, square_ratio, ipr_ratio

    def parse_info_value(info_path: Path, key: str) -> int:
        if not info_path.exists():
            return 0
        for line in info_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if ":" not in line:
                continue
            lhs, rhs = line.split(":", 1)
            if lhs.strip() != key:
                continue
            try:
                return int(float(rhs.strip()))
            except ValueError:
                return 0
        return 0

    rows: list[dict[str, str]] = []
    for report_path in sorted(root.glob("*/live_progress/report.md")):
        stage_root = report_path.parent.parent
        case_name = stage_root.name
        pngs = sorted(report_path.parent.glob("live_trace_*.png"))
        stage_json = stage_root / "production_stage.json"
        stage_config_json = stage_root / "stage_config.json"
        live_summary_json = report_path.parent / "summary.json"
        thermal_cut = 0
        if stage_json.exists():
            data = json.loads(stage_json.read_text(encoding="utf-8"))
            thermal_cut = int(data.get("config_summary", {}).get("thermal_cut", 0))
        elif stage_config_json.exists():
            data = json.loads(stage_config_json.read_text(encoding="utf-8"))
            thermal_cut = int(data.get("thermal_cut", 0))
        elif live_summary_json.exists():
            live_data = json.loads(live_summary_json.read_text(encoding="utf-8"))
            thermal_cut = int(live_data.get("thermal_cut", 0))
        run_dirs = sorted((stage_root / "runs").glob("*/*"))
        longest_trace, window_mode, square_ratio, ipr_ratio = summarize_run_dirs(run_dirs, thermal_cut)
        rows.append(
            {
                "stage_root": str(stage_root.resolve()),
                "label": case_name,
                "trace_path": str(pngs[0].resolve()) if pngs else "",
                "longest_trace": str(longest_trace),
                "window_mode": window_mode,
                "square_ratio": f"{square_ratio:.3f}",
                "ipr_ratio": f"{ipr_ratio:.3f}",
                "updated_epoch": report_path.stat().st_mtime,
            }
        )
    return rows


def write_report(
    output_dir: Path,
    stage_cases: list[dict[str, object]],
    stage_obs: list[dict[str, object]],
    tune_recommended: list[dict[str, object]],
    convergence_rows: list[dict[str, object]],
    live_reports: list[dict[str, str]],
) -> None:
    def fill_hmc_defaults(row: dict[str, object]) -> dict[str, object]:
        row.setdefault("hmc_mass_spatial_uniform", 0.0)
        row.setdefault("hmc_mass_spatial_lowk", 0.0)
        row.setdefault("hmc_mass_spatial_lowk_shells", 3)
        row.setdefault("hmc_mass_spatial_midk", 0.0)
        row.setdefault("hmc_mass_spatial_midk_shells", 0)
        row.setdefault("hmc_mass_spatial_shell1", 0.0)
        row.setdefault("hmc_mass_spatial_shell2", 0.0)
        return row

    recent_live_reports = sorted(live_reports, key=lambda row: float(row.get("updated_epoch", 0.0)), reverse=True)
    if recent_live_reports:
        latest_epoch = float(recent_live_reports[0].get("updated_epoch", 0.0))
        recent_live_reports = [
            row for row in recent_live_reports
            if latest_epoch - float(row.get("updated_epoch", 0.0)) <= 48.0 * 3600.0
        ] or recent_live_reports[:6]
        recent_live_reports = recent_live_reports[:6]
    lines = [
        "# HMC Production Overview",
        "",
        "This is the unified entry point for the current `data/triangular_hmc_production/` campaign.",
        "It merges stage and tune summaries so the production ladder can be reviewed from one report.",
        "Short traces can look stable before a late slow-mode drift becomes visible, so always read the `samples/post` coverage before trusting a stage.",
        "The dashed marker in each trace is only the configured `thermal_cut`; it is not an inferred optimal cut.",
        "Each overlaid trace belongs to an independent stage run at the same physics point, not to a single continued Markov chain.",
        "",
        "## Current Representative Stage Per Case",
        "",
        "This table picks one current best/most representative stage for each fixed `(L, Nbos, U2)` case.",
        "The selection prefers healthier retained-window behavior first, then lower cross-repeat mismatch, then deeper retained windows.",
        "",
        "| case | representative stage | beta | dtau | repeats | bins | cut | warm | samples/post | max drift | repeat span | nfrog | dt | mass | m_uniform | m_lowk | lowk_shells | m_midk | midk_shells | m_shell1 | m_shell2 | acceptance | tau_int(doubleOcc) | ESS/sec | status |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in select_representative_stage_cases(stage_cases):
        row = fill_hmc_defaults(dict(row))
        lines.append(
            "| {case} | {label} | {beta:.6g} | {dtau:.6g} | {completed_repeats}/{requested_repeats} | {bins} | {thermal_cut} | {warm} | {trace_samples_mean:.0f}/{post_thermal_samples_mean:.0f} | {gate_drift_ratio_max:.3f} | {gate_repeat_span_ratio:.3f} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {hmc_mass_spatial_lowk:.6g} | {hmc_mass_spatial_lowk_shells:.0f} | {hmc_mass_spatial_midk:.6g} | {hmc_mass_spatial_midk_shells:.0f} | {hmc_mass_spatial_shell1:.6g} | {hmc_mass_spatial_shell2:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
        "## Stage Summary",
        "",
        "| label | case | beta | dtau | repeats | bins | cut | warm | samples/post | max drift | repeat span | nfrog | dt | mass | m_uniform | m_lowk | lowk_shells | m_midk | midk_shells | m_shell1 | m_shell2 | acceptance | tau_int(doubleOcc) | ESS/sec | status |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in sorted(stage_cases, key=case_sort_key):
        row = fill_hmc_defaults(dict(row))
        lines.append(
            "| {label} | {case} | {beta:.6g} | {dtau:.6g} | {completed_repeats}/{requested_repeats} | {bins} | {thermal_cut} | {warm} | {trace_samples_mean:.0f}/{post_thermal_samples_mean:.0f} | {gate_drift_ratio_max:.3f} | {gate_repeat_span_ratio:.3f} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {hmc_mass_spatial_lowk:.6g} | {hmc_mass_spatial_lowk_shells:.0f} | {hmc_mass_spatial_midk:.6g} | {hmc_mass_spatial_midk_shells:.0f} | {hmc_mass_spatial_shell1:.6g} | {hmc_mass_spatial_shell2:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Recommended Tune Summary",
            "",
            "`viable=0` means every candidate in that tune sweep was flagged as stuck or below the minimum run acceptance.",
            "",
            "| label | case | beta | dtau | viable | nfrog | dt | mass | m_uniform | m_lowk | lowk_shells | m_midk | midk_shells | m_shell1 | m_shell2 | acceptance | tau_int(doubleOcc) | ESS/sec | DeltaH_abs_max |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in sorted(tune_recommended, key=lambda item: (str(item["case"]), float(item["Nbos"]) if "Nbos" in item else 0.0, float(item["U2"]) if "U2" in item else 0.0, float(item["beta"]), str(item["label"]))):
        row = fill_hmc_defaults(dict(row))
        lines.append(
            "| {label} | {case} | {beta:.6g} | {dtau:.6g} | {recommended_viable} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {hmc_mass_spatial_lowk:.6g} | {hmc_mass_spatial_lowk_shells:.0f} | {hmc_mass_spatial_midk:.6g} | {hmc_mass_spatial_midk_shells:.0f} | {hmc_mass_spatial_shell1:.6g} | {hmc_mass_spatial_shell2:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {hmc_deltaH_abs_max:.3f} |".format(
                **row
            )
        )
    if recent_live_reports:
        lines.extend(
            [
                "",
                "## Live Progress",
                "",
                "These are the currently active in-flight stage traces. They are shown inline here so the production campaign can be reviewed from this one report.",
                "If `window` is `pre-cut`, the drift ratios are computed on the currently available trace because the configured `thermal_cut` has not been crossed yet.",
                "",
                "| label | longest trace | window | squareOcc drift/span | IPR drift/span |",
                "| --- | ---: | --- | ---: | ---: |",
            ]
        )
        for row in recent_live_reports:
            lines.append(
                f"| {row['label']} | {row['longest_trace']} | {row['window_mode']} | {row['square_ratio']} | {row['ipr_ratio']} |"
            )
        for row in recent_live_reports:
            if not row.get("trace_path"):
                continue
            trace_rel = rel_link(Path(str(row["trace_path"])), output_dir)
            lines.extend(
                [
                    "",
                    f"### Live Trace: {row['label']}",
                    "",
                    f"![{row['label']} live trace]({trace_rel})",
                ]
            )
    lines.extend(
        [
            "",
            "## Observable Trends",
            "",
            "These plots are the main place to judge whether `squareOcc` and `IPR` are already converged in `beta` or `dtau`.",
            "The combined plots are only a quick map. The per-case sections below are the ones to trust when different parameter points sit on very different absolute scales.",
            "",
            "### squareOcc",
            "",
            "![squareOcc vs beta](squareOcc_vs_beta.png)",
            "",
            "![squareOcc vs dtau](squareOcc_vs_dtau.png)",
            "",
            "### IPR",
            "",
            "![IPR vs beta](IPR_vs_beta.png)",
            "",
            "![IPR vs dtau](IPR_vs_dtau.png)",
            "",
            "## Recommended Tune Trends",
            "",
            "![recommended tau vs beta](recommended_tau_int_doubleOcc_mean_vs_beta.png)",
            "",
            "![recommended ess vs beta](recommended_ess_per_sec_doubleOcc_mean_vs_beta.png)",
            "",
            "## Per-Case Convergence Summary",
            "",
            "This table compresses each fixed `(L, Nbos, U2)` ladder into relative changes of `squareOcc` and `IPR` across the available `beta/dtau` points.",
            "`rel spread` is `(max(mean) - min(mean)) / |first mean|`; `last-first` is `(last mean - first mean) / |first mean|`.",
            "",
            "| case | observable | n points | beta range | dtau range | rel spread | last-first |",
            "| --- | --- | ---: | --- | --- | ---: | ---: |",
        ]
    )
    for row in convergence_rows:
        lines.append(
            "| {case} | {observable} | {n_points} | {beta_min:.6g} -> {beta_max:.6g} | {dtau_max:.6g} -> {dtau_min:.6g} | {rel_spread:.6e} | {rel_last_minus_first:.6e} |".format(
                **row
            )
        )
    for case in sorted({str(row["case"]) for row in stage_cases}):
        case_tag = safe_tag(case)
        lines.extend(
            [
                "",
                f"## Case Summary: {case}",
                "",
                "| label | beta | dtau | repeats | bins | cut | warm | samples/post | max drift | repeat span | nfrog | dt | mass | m_uniform | m_lowk | lowk_shells | m_midk | midk_shells | m_shell1 | m_shell2 | acceptance | tau_int(doubleOcc) | ESS/sec | status |",
                "| --- | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
            ]
        )
        case_rows = [row for row in stage_cases if str(row["case"]) == case]
        for row in sorted(case_rows, key=lambda item: (float(item["beta"]), float(item["dtau"]), str(item["label"]))):
            row = fill_hmc_defaults(dict(row))
            lines.append(
                "| {label} | {beta:.6g} | {dtau:.6g} | {completed_repeats}/{requested_repeats} | {bins} | {thermal_cut} | {warm} | {trace_samples_mean:.0f}/{post_thermal_samples_mean:.0f} | {gate_drift_ratio_max:.3f} | {gate_repeat_span_ratio:.3f} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {hmc_mass_spatial_lowk:.6g} | {hmc_mass_spatial_lowk_shells:.0f} | {hmc_mass_spatial_midk:.6g} | {hmc_mass_spatial_midk_shells:.0f} | {hmc_mass_spatial_shell1:.6g} | {hmc_mass_spatial_shell2:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {status} |".format(
                    **row
                )
            )
        lines.extend(
            [
                "",
                f"## Trends: {case}",
                "",
                f"![{case} squareOcc vs beta](trend_{case_tag}_squareOcc_vs_beta.png)",
                "",
                f"![{case} squareOcc vs dtau](trend_{case_tag}_squareOcc_vs_dtau.png)",
                "",
                f"![{case} relative squareOcc vs beta](trend_{case_tag}_squareOcc_relative_vs_beta.png)",
                "",
                f"![{case} relative squareOcc vs dtau](trend_{case_tag}_squareOcc_relative_vs_dtau.png)",
                "",
                f"![{case} IPR vs beta](trend_{case_tag}_IPR_vs_beta.png)",
                "",
                f"![{case} IPR vs dtau](trend_{case_tag}_IPR_vs_dtau.png)",
                "",
                f"![{case} relative IPR vs beta](trend_{case_tag}_IPR_relative_vs_beta.png)",
                "",
                f"![{case} relative IPR vs dtau](trend_{case_tag}_IPR_relative_vs_dtau.png)",
                "",
                f"## Traces: {case}",
                "",
                f"![{case} squareOcc trace](trace_{case_tag}_squareOcc.png)",
                "",
                f"![{case} IPR trace](trace_{case_tag}_IPR.png)",
            ]
        )
    (output_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render one unified report for data/triangular_hmc_production.")
    parser.add_argument("--root", default="data/triangular_hmc_production")
    parser.add_argument("--output-dir", default="data/triangular_hmc_production/overview")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    stage_jsons, tune_jsons = scan_summaries(root)
    live_reports = scan_live_reports(root)
    stage_cases, stage_obs, stage_samples = load_stage_rows(stage_jsons)
    representative_stage_obs = select_representative_stage_rows(stage_cases, stage_obs)
    tune_candidates, tune_recommended = load_tune_rows(tune_jsons)
    convergence_rows = build_case_convergence_rows(stage_cases, stage_obs)

    write_csv(
        output_dir / "stage_cases.csv",
        stage_cases,
        [
            "label",
            "stage_label",
            "work_root",
            "case",
            "Lx",
            "Ly",
            "Nbos",
            "U2",
            "beta",
            "dtau",
            "bins",
            "thermal_cut",
            "warm",
            "nfrog",
            "hmc_dt",
            "hmc_mass",
            "hmc_mass_spatial_uniform",
            "hmc_mass_spatial_lowk",
            "hmc_mass_spatial_lowk_shells",
            "hmc_mass_spatial_midk",
            "hmc_mass_spatial_midk_shells",
            "hmc_mass_spatial_shell1",
            "hmc_mass_spatial_shell2",
            "completed_repeats",
            "requested_repeats",
            "missing_repeats",
            "trace_samples_mean",
            "post_thermal_samples_mean",
            "trace_samples_max",
            "gate_drift_ratio_max",
            "gate_repeat_span_ratio",
            "acceptance_mean",
            "tau_int_doubleOcc_mean",
            "ess_per_sec_doubleOcc_mean",
            "status",
        ],
    )
    write_csv(
        output_dir / "stage_observables.csv",
        representative_stage_obs,
        ["label", "stage_label", "case", "Lx", "Ly", "Nbos", "U2", "beta", "dtau", "observable", "mean", "stderr", "drift_over_span", "status"],
    )
    write_csv(
        output_dir / "case_convergence.csv",
        convergence_rows,
        ["case", "Lx", "Ly", "Nbos", "U2", "observable", "n_points", "beta_min", "beta_max", "dtau_min", "dtau_max", "first_label", "last_label", "first_mean", "last_mean", "min_mean", "max_mean", "rel_spread", "rel_last_minus_first"],
    )
    write_csv(
        output_dir / "stage_samples.csv",
        stage_samples,
        ["label", "stage_label", "case", "repeat", "observable", "sample_index", "post_thermal", "thermal_cut", "warm", "value"],
    )
    write_csv(
        output_dir / "tune_candidates.csv",
        tune_candidates,
        ["label", "work_root", "case", "beta", "dtau", "nfrog", "hmc_dt", "hmc_mass", "hmc_mass_spatial_uniform", "hmc_mass_spatial_lowk", "hmc_mass_spatial_lowk_shells", "hmc_mass_spatial_midk", "hmc_mass_spatial_midk_shells", "hmc_mass_spatial_shell1", "hmc_mass_spatial_shell2", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "hmc_deltaH_abs_max"],
    )
    write_csv(
        output_dir / "tune_recommended.csv",
        tune_recommended,
        ["label", "work_root", "case", "beta", "dtau", "recommended_viable", "nfrog", "hmc_dt", "hmc_mass", "hmc_mass_spatial_uniform", "hmc_mass_spatial_lowk", "hmc_mass_spatial_lowk_shells", "hmc_mass_spatial_midk", "hmc_mass_spatial_midk_shells", "hmc_mass_spatial_shell1", "hmc_mass_spatial_shell2", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "hmc_deltaH_abs_max"],
    )

    for observable in TREND_OBSERVABLES:
        render_trend_plot(representative_stage_obs, observable, "beta", output_dir / f"{observable}_vs_beta.png")
        render_trend_plot(representative_stage_obs, observable, "dtau", output_dir / f"{observable}_vs_dtau.png")
    render_recommended_plot(tune_recommended, "tau_int_doubleOcc_mean", "beta", output_dir / "recommended_tau_int_doubleOcc_mean_vs_beta.png")
    render_recommended_plot(tune_recommended, "ess_per_sec_doubleOcc_mean", "beta", output_dir / "recommended_ess_per_sec_doubleOcc_mean_vs_beta.png")
    for case in sorted({str(row["case"]) for row in stage_cases}):
        case_tag = safe_tag(case)
        for observable in TREND_OBSERVABLES:
            render_case_trend_plot(representative_stage_obs, case, observable, "beta", output_dir / f"trend_{case_tag}_{observable}_vs_beta.png")
            render_case_trend_plot(representative_stage_obs, case, observable, "dtau", output_dir / f"trend_{case_tag}_{observable}_vs_dtau.png")
            render_case_relative_trend_plot(representative_stage_obs, case, observable, "beta", output_dir / f"trend_{case_tag}_{observable}_relative_vs_beta.png")
            render_case_relative_trend_plot(representative_stage_obs, case, observable, "dtau", output_dir / f"trend_{case_tag}_{observable}_relative_vs_dtau.png")
        render_trace_plot(stage_samples, case, "squareOcc", output_dir / f"trace_{case_tag}_squareOcc.png")
        render_trace_plot(stage_samples, case, "IPR", output_dir / f"trace_{case_tag}_IPR.png")
    write_report(output_dir, stage_cases, stage_obs, tune_recommended, convergence_rows, live_reports)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
