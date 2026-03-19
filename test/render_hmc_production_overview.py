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
        case = data["cases"][0]
        cfg = case["config"]
        case_key = case["name"]
        beta = float(cfg["beta"])
        dtau = float(cfg["beta"]) / float(cfg["ltrot"])
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
                "nfrog": case["hmc"]["nfrog"],
                "hmc_dt": case["hmc"]["hmc_dt"],
                "hmc_mass": case["hmc"]["hmc_mass"],
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
        if samples_path.exists():
            with samples_path.open("r", encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                for row in reader:
                    if row["observable"] not in TRACE_OBSERVABLES:
                        continue
                    sample_rows.append(
                        {
                            "stage_label": data["stage_label"],
                            "case": row["name"],
                            "repeat": int(row["repeat"]),
                            "observable": row["observable"],
                            "sample_index": int(row["sample_index"]),
                            "post_thermal": int(row["post_thermal"]),
                            "value": float(row["value"]),
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
        ax.plot(
            [int(row["sample_index"]) for row in items],
            [float(row["value"]) for row in items],
            linewidth=1.2,
            label=stage_label,
        )
    ax.set_xlabel("sample index")
    ax.set_ylabel(observable)
    ax.set_title(f"{case}: {observable} trace (repeat 0)")
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
        "",
        "## Stage Summary",
        "",
        "| label | case | beta | dtau | nfrog | dt | mass | acceptance | tau_int(doubleOcc) | ESS/sec | status |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in sorted(stage_cases, key=lambda item: (str(item["case"]), float(item["beta"]), str(item["label"]))):
        lines.append(
            "| {label} | {case} | {beta:.6g} | {dtau:.6g} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Recommended Tune Summary",
            "",
            "| label | case | beta | dtau | nfrog | dt | mass | acceptance | tau_int(doubleOcc) | ESS/sec | DeltaH_abs_max |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in sorted(tune_recommended, key=lambda item: (str(item["case"]), float(item["beta"]), str(item["label"]))):
        lines.append(
            "| {label} | {case} | {beta:.6g} | {dtau:.6g} | {nfrog} | {hmc_dt:.6g} | {hmc_mass:.6g} | {acceptance_mean:.3f} | {tau_int_doubleOcc_mean:.3f} | {ess_per_sec_doubleOcc_mean:.3f} | {hmc_deltaH_abs_max:.3f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Observable Trends",
            "",
            "These plots are the main place to judge whether `squareOcc` and `IPR` are already converged in `beta` or `dtau`.",
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
        ["label", "stage_label", "work_root", "case", "Lx", "Ly", "Nbos", "U2", "beta", "dtau", "nfrog", "hmc_dt", "hmc_mass", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "status"],
    )
    write_csv(
        output_dir / "stage_observables.csv",
        stage_obs,
        ["label", "stage_label", "case", "beta", "dtau", "observable", "mean", "stderr", "drift_over_span", "status"],
    )
    write_csv(
        output_dir / "stage_samples.csv",
        stage_samples,
        ["stage_label", "case", "repeat", "observable", "sample_index", "post_thermal", "value"],
    )
    write_csv(
        output_dir / "tune_candidates.csv",
        tune_candidates,
        ["label", "work_root", "case", "beta", "dtau", "nfrog", "hmc_dt", "hmc_mass", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "hmc_deltaH_abs_max"],
    )
    write_csv(
        output_dir / "tune_recommended.csv",
        tune_recommended,
        ["label", "work_root", "case", "beta", "dtau", "nfrog", "hmc_dt", "hmc_mass", "acceptance_mean", "tau_int_doubleOcc_mean", "ess_per_sec_doubleOcc_mean", "hmc_deltaH_abs_max"],
    )

    for observable in TREND_OBSERVABLES:
        render_trend_plot(stage_obs, observable, "beta", output_dir / f"{observable}_vs_beta.png")
        render_trend_plot(stage_obs, observable, "dtau", output_dir / f"{observable}_vs_dtau.png")
    render_recommended_plot(tune_recommended, "tau_int_doubleOcc_mean", "beta", output_dir / "recommended_tau_int_doubleOcc_mean_vs_beta.png")
    render_recommended_plot(tune_recommended, "ess_per_sec_doubleOcc_mean", "beta", output_dir / "recommended_ess_per_sec_doubleOcc_mean_vs_beta.png")
    for case in sorted({str(row["case"]) for row in stage_cases}):
        case_tag = safe_tag(case)
        render_trace_plot(stage_samples, case, "squareOcc", output_dir / f"trace_{case_tag}_squareOcc.png")
        render_trace_plot(stage_samples, case, "IPR", output_dir / f"trace_{case_tag}_IPR.png")
    write_report(output_dir, stage_cases, stage_obs, tune_recommended)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
