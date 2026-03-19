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
import numpy as np


PLOT_OBSERVABLES = ["kinetic", "doubleOcc", "SF_Gamma", "SF_K", "PF_Gamma", "denden_Gamma"]
TRACE_OBSERVABLES = ["squareOcc", "IPR", "doubleOcc", "nearestOcc"]


def load_summary(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _bar_error_plot(labels: list[str], values: list[float], errors: list[float], *, ylabel: str, title: str, output: Path, target_band: tuple[float, float] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.5))
    x = np.arange(len(labels))
    ax.bar(x, values, color="#3b6fb6")
    ax.errorbar(x, values, yerr=errors, fmt="none", ecolor="black", capsize=4)
    if target_band is not None:
        ax.axhspan(target_band[0], target_band[1], color="#d7ebff", alpha=0.8)
        ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x, labels, rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def save_benchmark_acceptance_plot(summary: dict, out_dir: Path) -> None:
    cases = summary["cases"]
    labels = [case["name"] for case in cases]
    values = [case["acceptance_mean"] for case in cases]
    errors = [case["acceptance_stderr"] for case in cases]
    _bar_error_plot(
        labels,
        values,
        errors,
        ylabel="Acceptance",
        title="HMC Acceptance by Case",
        output=out_dir / "acceptance.png",
        target_band=(0.70, 0.85),
    )


def save_benchmark_perf_plot(summary: dict, out_dir: Path) -> None:
    cases = summary["cases"]
    labels = [case["name"] for case in cases]
    local_values = [case["perf"]["local"]["ess_per_sec_doubleOcc"] for case in cases]
    hmc_values = [case["perf"]["hmc"]["ess_per_sec_doubleOcc"] for case in cases]
    x = np.arange(len(labels))
    width = 0.38
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.5))
    ax.bar(x - width / 2, local_values, width, label="Local", color="#999999")
    ax.bar(x + width / 2, hmc_values, width, label="HMC", color="#d55e00")
    ax.set_ylabel("ESS / sec")
    ax.set_title("Sampling Efficiency")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "ess_per_sec.png", dpi=180)
    plt.close(fig)


def save_benchmark_tau_plot(summary: dict, out_dir: Path) -> None:
    cases = summary["cases"]
    labels = [case["name"] for case in cases]
    local_values = [case["perf"]["local"]["tau_int_doubleOcc"] for case in cases]
    hmc_values = [case["perf"]["hmc"]["tau_int_doubleOcc"] for case in cases]
    x = np.arange(len(labels))
    width = 0.38
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.5))
    ax.bar(x - width / 2, local_values, width, label="Local", color="#999999")
    ax.bar(x + width / 2, hmc_values, width, label="HMC", color="#009e73")
    ax.set_ylabel("Integrated Autocorrelation Time")
    ax.set_title("tau_int(doubleOcc)")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "tau_int.png", dpi=180)
    plt.close(fig)


def save_benchmark_zscore_heatmap(summary: dict, out_dir: Path) -> None:
    cases = summary["cases"]
    labels = [case["name"] for case in cases]
    matrix = np.zeros((len(PLOT_OBSERVABLES), len(cases)))
    for jdx, case in enumerate(cases):
        by_name = {row["observable"]: row for row in case["observables"]}
        for idx, obs in enumerate(PLOT_OBSERVABLES):
            matrix[idx, jdx] = by_name[obs]["z_score"]
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.5))
    im = ax.imshow(matrix, aspect="auto", cmap="magma")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(PLOT_OBSERVABLES)), PLOT_OBSERVABLES)
    ax.set_title("Local vs HMC z-scores")
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(out_dir / "zscore_heatmap.png", dpi=180)
    plt.close(fig)


def write_benchmark_markdown(summary: dict, out_dir: Path) -> None:
    lines = [
        "# HMC Benchmark Report",
        "",
        f"Overall result: `{'PASS' if summary['passed'] else 'FAIL'}`",
        "",
        "| Case | Acceptance | tau_int local | tau_int HMC | ESS/sec local | ESS/sec HMC | Speed ratio | Result |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for case in summary["cases"]:
        lines.append(
            "| {name} | {accept:.3f} | {tau_l:.3f} | {tau_h:.3f} | {ess_l:.3f} | {ess_h:.3f} | {ratio:.3f} | {result} |".format(
                name=case["name"],
                accept=case["acceptance_mean"],
                tau_l=case["perf"]["local"]["tau_int_doubleOcc"],
                tau_h=case["perf"]["hmc"]["tau_int_doubleOcc"],
                ess_l=case["perf"]["local"]["ess_per_sec_doubleOcc"],
                ess_h=case["perf"]["hmc"]["ess_per_sec_doubleOcc"],
                ratio=case["perf"]["speed_ratio_hmc_over_local"],
                result="PASS" if case["passed"] else "FAIL",
            )
        )
    lines.extend(
        [
            "",
            "## Figures",
            "",
            "![Acceptance](acceptance.png)",
            "",
            "![ESS per sec](ess_per_sec.png)",
            "",
            "![tau_int](tau_int.png)",
            "",
            "![z-score heatmap](zscore_heatmap.png)",
            "",
        ]
    )
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def render_benchmark_summary(summary: dict, out_dir: Path) -> None:
    save_benchmark_acceptance_plot(summary, out_dir)
    save_benchmark_perf_plot(summary, out_dir)
    save_benchmark_tau_plot(summary, out_dir)
    save_benchmark_zscore_heatmap(summary, out_dir)
    write_benchmark_markdown(summary, out_dir)


def save_tune_summary(summary: dict, out_dir: Path) -> None:
    lines = ["# HMC Tune Report", ""]
    for case in summary["cases"]:
        rows = case["rows"]
        rows = sorted(
            rows,
            key=lambda row: (
                row["hmc_mass"],
                row.get("hmc_mass_spatial_uniform", 0.0),
                row.get("hmc_mass_spatial_shell1", 0.0),
                row.get("hmc_mass_spatial_shell2", 0.0),
                row["nfrog"] * row["hmc_dt"],
            ),
        )
        traj = [row["nfrog"] * row["hmc_dt"] for row in rows]
        acc = [row["acceptance_mean"] for row in rows]
        tau = [row["tau_int_doubleOcc_mean"] for row in rows]
        ess = [row["ess_per_sec_doubleOcc_mean"] for row in rows]
        labels = [
            (
                f"{int(row['nfrog'])}x{row['hmc_dt']:g}\n"
                f"m={row['hmc_mass']:g}\n"
                f"mu={row.get('hmc_mass_spatial_uniform', 0.0):g}\n"
                f"mk1={row.get('hmc_mass_spatial_shell1', 0.0):g}\n"
                f"mk2={row.get('hmc_mass_spatial_shell2', 0.0):g}"
            )
            for row in rows
        ]

        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
        axes[0].scatter(traj, acc, s=70, color="#3b6fb6")
        axes[0].axhspan(0.70, 0.85, color="#d7ebff", alpha=0.8)
        axes[0].set_ylabel("Acceptance")
        axes[0].set_xlabel("Nfrog * dt")
        axes[0].set_title("Acceptance")
        axes[1].scatter(traj, tau, s=70, color="#009e73")
        axes[1].set_ylabel("tau_int(doubleOcc)")
        axes[1].set_xlabel("Nfrog * dt")
        axes[1].set_title("Autocorrelation")
        axes[2].scatter(traj, ess, s=70, color="#d55e00")
        axes[2].set_ylabel("ESS / sec")
        axes[2].set_xlabel("Nfrog * dt")
        axes[2].set_title("Efficiency")
        series_list = [acc, tau, ess]
        for axis_index, axis in enumerate(axes):
            for xpos, ypos, label in zip(traj, series_list[axis_index], labels):
                axis.annotate(label, (xpos, ypos), textcoords="offset points", xytext=(4, 4), fontsize=8)
        fig.tight_layout()
        filename = f"tune_{case['name']}.png"
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)

        rec = case["recommended"]
        viable = case.get("recommended_viable", True)
        lines.extend(
            [
                f"## {case['name']}",
                "",
                "Recommended candidate:" if viable else "No viable recommended candidate in this scan.",
                "",
                "- `nfrog={}`".format(int(rec["nfrog"])),
                "- `dt={}`".format(rec["hmc_dt"]),
                "- `jitter={}`".format(int(rec["hmc_jitter"])),
                "- `mass={}`".format(rec["hmc_mass"]),
                "- `uniform_mass={}`".format(rec.get("hmc_mass_spatial_uniform", 0.0)),
                "- `shell1_mass={}`".format(rec.get("hmc_mass_spatial_shell1", 0.0)),
                "- `shell2_mass={}`".format(rec.get("hmc_mass_spatial_shell2", 0.0)),
                "- `acceptance={:.3f} ± {:.3f}`".format(rec["acceptance_mean"], rec["acceptance_stderr"]),
                "- `tau_int(doubleOcc)={:.3f} ± {:.3f}`".format(rec["tau_int_doubleOcc_mean"], rec["tau_int_doubleOcc_stderr"]),
                "- `ESS/sec={:.3f} ± {:.3f}`".format(rec["ess_per_sec_doubleOcc_mean"], rec["ess_per_sec_doubleOcc_stderr"]),
                "",
                f"![{case['name']} tune]({filename})",
                "",
            ]
        )
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def save_stage_acceptance_plot(case_rows: list[dict[str, str]], out_dir: Path) -> None:
    labels = [row["name"] for row in case_rows]
    values = [float(row["acceptance_mean"]) for row in case_rows]
    errors = [float(row["acceptance_stderr"]) for row in case_rows]
    _bar_error_plot(
        labels,
        values,
        errors,
        ylabel="Acceptance",
        title="HMC Acceptance by Stage Case",
        output=out_dir / "acceptance.png",
        target_band=None,
    )


def save_stage_perf_plot(case_rows: list[dict[str, str]], out_dir: Path) -> None:
    labels = [row["name"] for row in case_rows]
    values = [float(row["ess_per_sec_doubleOcc_mean"]) for row in case_rows]
    errors = [float(row["ess_per_sec_doubleOcc_stderr"]) for row in case_rows]
    _bar_error_plot(
        labels,
        values,
        errors,
        ylabel="ESS / sec",
        title="HMC ESS/sec by Stage Case",
        output=out_dir / "ess_per_sec.png",
        target_band=None,
    )


def save_stage_tau_plot(case_rows: list[dict[str, str]], out_dir: Path) -> None:
    labels = [row["name"] for row in case_rows]
    values = [float(row["tau_int_doubleOcc_mean"]) for row in case_rows]
    errors = [float(row["tau_int_doubleOcc_stderr"]) for row in case_rows]
    _bar_error_plot(
        labels,
        values,
        errors,
        ylabel="Integrated Autocorrelation Time",
        title="tau_int(doubleOcc) by Stage Case",
        output=out_dir / "tau_int.png",
        target_band=None,
    )


def save_stage_drift_plot(case_rows: list[dict[str, str]], out_dir: Path) -> None:
    labels = [row["name"] for row in case_rows]
    square = [float(row["squareOcc_drift_ratio"]) for row in case_rows]
    ipr = [float(row["IPR_drift_ratio"]) for row in case_rows]
    x = np.arange(len(labels))
    width = 0.38
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.5))
    ax.bar(x - width / 2, square, width, label="squareOcc", color="#0072b2")
    ax.bar(x + width / 2, ipr, width, label="IPR", color="#cc79a7")
    ax.axhline(0.25, color="black", linestyle="--", linewidth=1.0, label="stable guide")
    ax.axhline(0.50, color="black", linestyle=":", linewidth=1.0, label="warning guide")
    ax.set_ylabel("|last window - first window| / span")
    ax.set_title("Gate Observable Drift Ratios")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "drift_ratios.png", dpi=180)
    plt.close(fig)


def save_stage_trace_plots(case_rows: list[dict[str, str]], sample_rows: list[dict[str, str]], out_dir: Path) -> list[str]:
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in sample_rows:
        grouped[row["name"]][row["observable"]].append(row)
    filenames: list[str] = []
    for case_row in case_rows:
        case_name = case_row["name"]
        fig, axes = plt.subplots(len(TRACE_OBSERVABLES), 1, figsize=(10, 2.8 * len(TRACE_OBSERVABLES)), sharex=True)
        if len(TRACE_OBSERVABLES) == 1:
            axes = [axes]
        for ax, observable in zip(axes, TRACE_OBSERVABLES):
            obs_rows = grouped[case_name].get(observable, [])
            by_repeat: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in obs_rows:
                by_repeat[row["repeat"]].append(row)
            for repeat, repeat_rows in sorted(by_repeat.items(), key=lambda item: int(item[0])):
                repeat_rows.sort(key=lambda row: int(row["sample_index"]))
                x = [int(row["sample_index"]) for row in repeat_rows]
                y = [float(row["value"]) for row in repeat_rows]
                ax.plot(x, y, linewidth=1.2, label=f"rep {repeat}")
            ax.set_ylabel(observable)
            ax.grid(alpha=0.2)
            if observable == "squareOcc":
                ax.set_title(case_name)
            if len(by_repeat) <= 6:
                ax.legend(loc="best", fontsize=8)
        axes[-1].set_xlabel("Sample index")
        fig.tight_layout()
        filename = f"trace_{case_name}.png"
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)
        filenames.append(filename)
    return filenames


def write_stage_markdown(summary: dict, case_rows: list[dict[str, str]], trace_files: list[str], out_dir: Path) -> None:
    lines = [
        "# HMC Production Stage Report",
        "",
        f"Stage label: `{summary.get('stage_label', '')}`",
        "",
        f"Overall status: `{summary.get('overall_status', 'needs_review')}`",
        "",
        "Primary gate observables:",
        "",
        "- `squareOcc`",
        "- `IPR`",
        "",
        "Heuristic status meanings:",
        "",
        "- `stable_window`: gate-observable drift is small and different repeats land on consistent retained windows",
        "- `slow_drift`: visible drift or repeat-to-repeat plateau mismatch remains, but the run is no longer obviously stuck",
        "- `strong_drift`: gate observables still drift strongly or different repeats settle onto clearly different retained windows",
        "- `stuck_or_invalid`: acceptance or ESS indicates that at least one repeat is effectively unusable",
        "",
        "| Case | mass | m_uniform | m_shell1 | m_shell2 | Acceptance | tau_int(doubleOcc) | ESS/sec | squareOcc drift/span | IPR drift/span | Status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in case_rows:
        lines.append(
            "| {name} | {mass:.6g} | {uniform:.6g} | {shell1:.6g} | {shell2:.6g} | {accept:.3f} | {tau:.3f} | {ess:.3f} | {sq:.3f} | {ipr:.3f} | {status} |".format(
                name=row["name"],
                mass=float(row["hmc_mass"]),
                uniform=float(row.get("hmc_mass_spatial_uniform", 0.0)),
                shell1=float(row.get("hmc_mass_spatial_shell1", 0.0)),
                shell2=float(row.get("hmc_mass_spatial_shell2", 0.0)),
                accept=float(row["acceptance_mean"]),
                tau=float(row["tau_int_doubleOcc_mean"]),
                ess=float(row["ess_per_sec_doubleOcc_mean"]),
                sq=float(row["squareOcc_drift_ratio"]),
                ipr=float(row["IPR_drift_ratio"]),
                status=row["status"],
            )
        )
    lines.extend(
        [
            "",
            "## Overview",
            "",
            "This report is meant for production bring-up rather than local-vs-HMC correctness benchmarking. "
            "The main question is whether `squareOcc` and `IPR` look thermalized enough to promote the point to a harder rung.",
            "",
            "![Acceptance](acceptance.png)",
            "",
            "![ESS per sec](ess_per_sec.png)",
            "",
            "![tau_int](tau_int.png)",
            "",
            "![Drift ratios](drift_ratios.png)",
            "",
            "## Sample-Index Traces",
            "",
        ]
    )
    for filename in trace_files:
        lines.extend([f"![{filename}]({filename})", ""])
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def render_stage_summary(summary: dict, summary_path: Path, out_dir: Path) -> None:
    case_rows = read_csv_rows(summary_path.parent / "production_stage_cases.csv")
    sample_rows = read_csv_rows(summary_path.parent / "production_stage_samples.csv")
    save_stage_acceptance_plot(case_rows, out_dir)
    save_stage_perf_plot(case_rows, out_dir)
    save_stage_tau_plot(case_rows, out_dir)
    save_stage_drift_plot(case_rows, out_dir)
    trace_files = save_stage_trace_plots(case_rows, sample_rows, out_dir)
    write_stage_markdown(summary, case_rows, trace_files, out_dir)


def render_summary_file(summary_json: Path, out_dir: Path) -> None:
    summary_path = summary_json.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = load_summary(summary_path)
    mode = summary.get("mode", "")
    if mode == "benchmark":
        render_benchmark_summary(summary, out_dir)
        return
    if mode == "tune":
        save_tune_summary(summary, out_dir)
        return
    if mode == "stage":
        render_stage_summary(summary, summary_path, out_dir)
        return
    raise ValueError(f"Unsupported summary mode for rendering: {mode}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render figures and a Markdown report from a production summary JSON.")
    parser.add_argument("summary_json")
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()

    summary_path = Path(args.summary_json).resolve()
    out_dir = Path(args.output_dir).resolve() if args.output_dir else summary_path.parent / "report"
    render_summary_file(summary_path, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
