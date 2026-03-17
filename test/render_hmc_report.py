#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PLOT_OBSERVABLES = ["kinetic", "doubleOcc", "SF_Gamma", "SF_K", "PF_Gamma", "denden_Gamma"]


def load_summary(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_acceptance_plot(summary: dict, out_dir: Path) -> None:
    cases = summary["cases"]
    labels = [case["name"] for case in cases]
    values = [case["acceptance_mean"] for case in cases]
    errors = [case["acceptance_stderr"] for case in cases]
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(labels)), 4.5))
    x = np.arange(len(labels))
    ax.bar(x, values, color="#3b6fb6")
    ax.errorbar(x, values, yerr=errors, fmt="none", ecolor="black", capsize=4)
    ax.axhspan(0.70, 0.85, color="#d7ebff", alpha=0.8)
    ax.set_ylabel("Acceptance")
    ax.set_title("HMC Acceptance by Case")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.set_ylim(0.0, 1.0)
    fig.tight_layout()
    fig.savefig(out_dir / "acceptance.png", dpi=180)
    plt.close(fig)


def save_perf_plot(summary: dict, out_dir: Path) -> None:
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


def save_tau_plot(summary: dict, out_dir: Path) -> None:
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


def save_zscore_heatmap(summary: dict, out_dir: Path) -> None:
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


def write_markdown(summary: dict, out_dir: Path) -> None:
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Render figures and a Markdown report from production_benchmark.json.")
    parser.add_argument("summary_json")
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()

    summary_path = Path(args.summary_json).resolve()
    out_dir = Path(args.output_dir).resolve() if args.output_dir else summary_path.parent / "report"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = load_summary(summary_path)
    save_acceptance_plot(summary, out_dir)
    save_perf_plot(summary, out_dir)
    save_tau_plot(summary, out_dir)
    save_zscore_heatmap(summary, out_dir)
    write_markdown(summary, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
