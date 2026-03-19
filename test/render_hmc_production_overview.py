#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TRACE_OBSERVABLES = ("squareOcc", "IPR", "doubleOcc", "nearestOcc")
TREND_OBSERVABLES = ("squareOcc", "IPR")
GATE_OBSERVABLES = ("squareOcc", "IPR")


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
                "nfrog": rec["nfrog"],
                "hmc_dt": rec["hmc_dt"],
                "hmc_mass": rec["hmc_mass"],
                "hmc_mass_spatial_uniform": rec.get("hmc_mass_spatial_uniform", 0.0),
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


def write_report(
    output_dir: Path,
    stage_cases: list[dict[str, object]],
    stage_obs: list[dict[str, object]],
    tune_recommended: list[dict[str, object]],
) -> None:
    lines = [
        "# HMC Production Overview",
        "",
        "This is the unified entry point for the current `data/triangular_hmc_production/` campaign.",
        "It merges stage and tune summaries so the production ladder can be reviewed from one report.",
        "Short traces can look stable before a late slow-mode drift becomes visible, so always read the `samples/post` coverage before trusting a stage.",
        "The dashed marker in each trace is only the configured `thermal_cut`; it is not an inferred optimal cut.",
        "Each overlaid trace belongs to an independent stage run at the same physics point, not to a single continued Markov chain.",
        "",
        "## Stage Summary",
        "",
        "| label | case | beta | dtau | repeats | bins | cut | warm | samples/post | max drift | repeat span | nfrog | dt | mass | m_uniform | acceptance | tau_int(doubleOcc) | ESS/sec | status |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in sorted(stage_cases, key=lambda item: (str(item["case"]), float(item["beta"]), str(item["label"]))):
        lines.append(
            "| {label} | {case} | {beta:.6g} | {dtau:.6g} | {completed_repeats}/{requested_repeats} | {bins} | {thermal_cut} | {warm} | {trace_samples_mean:.0f}/{post_thermal_samples_mean:.0f} | {gate_drift_ratio_max:.3f} | {gate_repeat_span_ratio:.3f} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Recommended Tune Summary",
            "",
            "| label | case | beta | dtau | nfrog | dt | mass | m_uniform | acceptance | tau_int(doubleOcc) | ESS/sec | DeltaH_abs_max |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in sorted(tune_recommended, key=lambda item: (str(item["case"]), float(item["beta"]), str(item["label"]))):
        lines.append(
            "| {label} | {case} | {beta:.6g} | {dtau:.6g} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {hmc_deltaH_abs_max:.3f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Observable Trends",
            "",
            "These plots are the main place to judge whether `squareOcc` and `IPR` are already converged in `beta` or `dtau`.",
            "The combined plots are only for a quick overview; the per-case plots below are the ones to trust when different parameter points have very different scales.",
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
        ]
    )
    for case in sorted({str(row["case"]) for row in stage_cases}):
        case_tag = safe_tag(case)
        lines.extend(
            [
                "",
                f"## Case Summary: {case}",
                "",
                "| label | beta | dtau | repeats | bins | cut | warm | samples/post | max drift | repeat span | nfrog | dt | mass | m_uniform | acceptance | tau_int(doubleOcc) | ESS/sec | status |",
                "| --- | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
            ]
        )
        case_rows = [row for row in stage_cases if str(row["case"]) == case]
        for row in sorted(case_rows, key=lambda item: (float(item["beta"]), str(item["label"]))):
            lines.append(
                "| {label} | {beta:.6g} | {dtau:.6g} | {completed_repeats}/{requested_repeats} | {bins} | {thermal_cut} | {warm} | {trace_samples_mean:.0f}/{post_thermal_samples_mean:.0f} | {gate_drift_ratio_max:.3f} | {gate_repeat_span_ratio:.3f} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {hmc_mass_spatial_uniform:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {status} |".format(
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
                f"![{case} IPR vs beta](trend_{case_tag}_IPR_vs_beta.png)",
                "",
                f"![{case} IPR vs dtau](trend_{case_tag}_IPR_vs_dtau.png)",
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
    stage_cases, stage_obs, stage_samples = load_stage_rows(stage_jsons)
    tune_candidates, tune_recommended = load_tune_rows(tune_jsons)

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
        stage_obs,
        ["label", "stage_label", "case", "beta", "dtau", "observable", "mean", "stderr", "drift_over_span", "status"],
    )
    write_csv(
        output_dir / "stage_samples.csv",
        stage_samples,
        ["label", "stage_label", "case", "repeat", "observable", "sample_index", "post_thermal", "thermal_cut", "warm", "value"],
    )
    write_csv(
        output_dir / "tune_candidates.csv",
        tune_candidates,
        ["label", "work_root", "case", "beta", "dtau", "nfrog", "hmc_dt", "hmc_mass", "hmc_mass_spatial_uniform", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "hmc_deltaH_abs_max"],
    )
    write_csv(
        output_dir / "tune_recommended.csv",
        tune_recommended,
        ["label", "work_root", "case", "beta", "dtau", "nfrog", "hmc_dt", "hmc_mass", "hmc_mass_spatial_uniform", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "hmc_deltaH_abs_max"],
    )

    for observable in TREND_OBSERVABLES:
        render_trend_plot(stage_obs, observable, "beta", output_dir / f"{observable}_vs_beta.png")
        render_trend_plot(stage_obs, observable, "dtau", output_dir / f"{observable}_vs_dtau.png")
    render_recommended_plot(tune_recommended, "tau_int_doubleOcc_mean", "beta", output_dir / "recommended_tau_int_doubleOcc_mean_vs_beta.png")
    render_recommended_plot(tune_recommended, "ess_per_sec_doubleOcc_mean", "beta", output_dir / "recommended_ess_per_sec_doubleOcc_mean_vs_beta.png")
    for case in sorted({str(row["case"]) for row in stage_cases}):
        case_tag = safe_tag(case)
        for observable in TREND_OBSERVABLES:
            render_case_trend_plot(stage_obs, case, observable, "beta", output_dir / f"trend_{case_tag}_{observable}_vs_beta.png")
            render_case_trend_plot(stage_obs, case, observable, "dtau", output_dir / f"trend_{case_tag}_{observable}_vs_dtau.png")
        render_trace_plot(stage_samples, case, "squareOcc", output_dir / f"trace_{case_tag}_squareOcc.png")
        render_trace_plot(stage_samples, case, "IPR", output_dir / f"trace_{case_tag}_IPR.png")
    write_report(output_dir, stage_cases, stage_obs, tune_recommended)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
