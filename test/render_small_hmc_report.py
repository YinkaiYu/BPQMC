#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CURVE_OBSERVABLES = [
    "IPR",
    "kinetic",
    "doubleOcc",
    "squareOcc",
    "nearestOcc",
    "SF_Gamma",
    "SF_K",
    "PF_Gamma",
    "C3_Gamma",
    "dentot_Gamma",
    "denden_Gamma",
]
SERIES_OBSERVABLES = ["IPR", "kinetic", "doubleOcc", "nearestOcc", "SF_Gamma", "SF_K"]


def load_summary(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_pass_matrix(cases_df: pd.DataFrame, out_dir: Path) -> None:
    pivot = cases_df.pivot(index="Nbos", columns="U2", values="passed").sort_index().sort_index(axis=1)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    matrix = pivot.to_numpy(dtype=float)
    im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0.0, vmax=1.0)
    ax.set_xticks(np.arange(pivot.shape[1]), [f"{value:g}" for value in pivot.columns])
    ax.set_yticks(np.arange(pivot.shape[0]), [f"{value:g}" for value in pivot.index])
    ax.set_xlabel("U2")
    ax.set_ylabel("Nbos")
    ax.set_title("Strict benchmark pass matrix")
    for (i, j), value in np.ndenumerate(matrix):
        ax.text(j, i, "PASS" if value > 0.5 else "FAIL", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="pass=1 fail=0")
    fig.tight_layout()
    fig.savefig(out_dir / "pass_matrix.png", dpi=180)
    plt.close(fig)


def save_acceptance_tau(cases_df: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    for nbos, group in cases_df.groupby("Nbos", sort=True):
        group = group.sort_values("U2")
        x = group["U2"].to_numpy()
        axes[0].errorbar(x, group["acceptance_mean"], yerr=group["acceptance_stderr"], marker="o", capsize=3, label=f"Nbos={nbos}")
        axes[1].plot(x, group["local_tau_int_doubleOcc"], marker="o", linestyle="-", label=f"Local N={nbos}")
        axes[1].plot(x, group["hmc_tau_int_doubleOcc"], marker="s", linestyle="--", label=f"HMC N={nbos}")
        axes[2].plot(x, group["speed_ratio_hmc_over_local"], marker="o", linestyle="-", label=f"Nbos={nbos}")
    axes[0].set_xscale("log")
    axes[1].set_xscale("log")
    axes[2].set_xscale("log")
    axes[0].set_ylabel("Acceptance")
    axes[1].set_ylabel("tau_int(doubleOcc)")
    axes[2].set_ylabel("ESS/sec ratio (HMC / local)")
    axes[0].set_xlabel("U2")
    axes[1].set_xlabel("U2")
    axes[2].set_xlabel("U2")
    axes[0].set_title("Acceptance")
    axes[1].set_title("Autocorrelation")
    axes[2].set_title("Relative efficiency")
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=7, ncol=2)
    axes[2].legend(fontsize=8)
    fig.savefig(out_dir / "acceptance_tau_efficiency.png", dpi=180)
    plt.close(fig)


def save_observable_curves(obs_df: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for obs_name in CURVE_OBSERVABLES:
        subset = obs_df[obs_df["observable"] == obs_name].copy()
        if subset.empty:
            continue
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        for nbos, group in subset.groupby("Nbos", sort=True):
            group = group.sort_values("U2")
            x = group["U2"].to_numpy()
            ax.errorbar(x, group["local_mean"], yerr=group["local_err"], marker="o", linestyle="-", capsize=3, label=f"Local N={nbos}")
            ax.errorbar(x, group["hmc_mean"], yerr=group["hmc_err"], marker="s", linestyle="--", capsize=3, label=f"HMC N={nbos}")
        ax.set_xscale("log")
        ax.set_xlabel("U2")
        ax.set_ylabel(obs_name)
        ax.set_title(f"{obs_name} vs U2")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(out_dir / f"{obs_name}_vs_u2.png", dpi=180)
        plt.close(fig)


def save_tune_panels(tune_df: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, case_df in tune_df.groupby("name", sort=True):
        case_df = case_df.copy()
        case_df["traj_len"] = case_df["nfrog"] * case_df["hmc_dt"]
        case_df.sort_values(["hmc_mass", "traj_len"], inplace=True)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
        for mass, mass_df in case_df.groupby("hmc_mass", sort=True):
            label = f"mass={mass:g}"
            axes[0].plot(mass_df["traj_len"], mass_df["acceptance_mean"], marker="o", label=label)
            axes[1].plot(mass_df["traj_len"], mass_df["tau_int_doubleOcc_mean"], marker="o", label=label)
            axes[2].plot(mass_df["traj_len"], mass_df["ess_per_sec_doubleOcc_mean"], marker="o", label=label)
        axes[0].set_ylabel("Acceptance")
        axes[1].set_ylabel("tau_int(doubleOcc)")
        axes[2].set_ylabel("ESS / sec")
        for ax in axes:
            ax.set_xlabel("trajectory length = Nfrog * dt")
            ax.set_xscale("log")
            ax.legend(fontsize=8)
        fig.suptitle(name)
        fig.savefig(out_dir / f"{name}_tune.png", dpi=180)
        plt.close(fig)


def save_sample_traces(samples_df: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for (nbos, u2), case_df in samples_df.groupby(["Nbos", "U2"], sort=True):
        fig, axes = plt.subplots(len(SERIES_OBSERVABLES), 1, figsize=(10, 2.7 * len(SERIES_OBSERVABLES)), sharex=True, constrained_layout=True)
        for idx, obs_name in enumerate(SERIES_OBSERVABLES):
            ax = axes[idx]
            subset = case_df[case_df["observable"] == obs_name].sort_values(["mode", "sample_index"])
            for mode, mode_df in subset.groupby("mode"):
                ax.plot(mode_df["sample_index"], mode_df["value"], label=mode.upper(), linewidth=1.1)
            ax.set_ylabel(obs_name)
            ax.legend(fontsize=8, loc="best")
        axes[-1].set_xlabel("sample index after thermal cut")
        fig.suptitle(f"Nbos={nbos}, U2={u2:g}")
        fig.savefig(out_dir / f"trace_N{nbos}_U2_{u2:g}.png", dpi=180)
        plt.close(fig)


def write_markdown(summary: dict, cases_df: pd.DataFrame, out_dir: Path, tune_dir: Path | None) -> None:
    lines = [
        "# Triangular L=6 HMC Small-Parameter Benchmark",
        "",
        f"Overall result: `{'PASS' if summary['passed'] else 'FAIL'}`",
        "",
        "## Summary Table",
        "",
        "| Case | Acceptance | tau_int local | tau_int HMC | ESS/sec local | ESS/sec HMC | Result |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for _, row in cases_df.sort_values(["Nbos", "U2"]).iterrows():
        lines.append(
            "| {name} | {accept:.3f} | {tau_l:.3f} | {tau_h:.3f} | {ess_l:.4f} | {ess_h:.4f} | {result} |".format(
                name=row["name"],
                accept=row["acceptance_mean"],
                tau_l=row["local_tau_int_doubleOcc"],
                tau_h=row["hmc_tau_int_doubleOcc"],
                ess_l=row["local_ess_per_sec_doubleOcc"],
                ess_h=row["hmc_ess_per_sec_doubleOcc"],
                result="PASS" if row["passed"] else "FAIL",
            )
        )
    lines.extend(
        [
            "",
            "## Figures",
            "",
            "![Pass matrix](pass_matrix.png)",
            "",
            "![Acceptance/tau/efficiency](acceptance_tau_efficiency.png)",
            "",
            "### Observable curves",
            "",
        ]
    )
    for obs_name in CURVE_OBSERVABLES:
        png = f"observables/{obs_name}_vs_u2.png"
        if (out_dir / png).exists():
            lines.extend([f"![{obs_name}]({png})", ""])
    if tune_dir is not None and any(tune_dir.glob("*_tune.png")):
        lines.extend(["### Tuning panels", ""])
        for png in sorted(tune_dir.glob("*_tune.png")):
            lines.extend([f"![{png.stem}](tune/{png.name})", ""])
    lines.extend(["### Sample-index traces", ""])
    for png in sorted((out_dir / "traces").glob("trace_*.png")):
        lines.extend([f"![{png.stem}](traces/{png.name})", ""])
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render the triangular L=6 small-parameter HMC benchmark report.")
    parser.add_argument("summary_json")
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()

    summary_path = Path(args.summary_json).resolve()
    root = summary_path.parent
    out_dir = Path(args.output_dir).resolve() if args.output_dir else root / "report"
    out_dir.mkdir(parents=True, exist_ok=True)

    cases_csv = root / "small_benchmark_cases.csv"
    observables_csv = root / "small_benchmark_observables.csv"
    samples_csv = root / "small_benchmark_samples.csv"
    tune_csv = root / "small_tune.csv"
    summary = load_summary(summary_path)
    cases_df = pd.read_csv(cases_csv)
    obs_df = pd.read_csv(observables_csv)
    samples_df = pd.read_csv(samples_csv)
    tune_dir = None
    if tune_csv.exists():
        tune_df = pd.read_csv(tune_csv)
        tune_dir = out_dir / "tune"
        save_tune_panels(tune_df, tune_dir)

    save_pass_matrix(cases_df, out_dir)
    save_acceptance_tau(cases_df, out_dir)
    save_observable_curves(obs_df, out_dir / "observables")
    save_sample_traces(samples_df, out_dir / "traces")
    write_markdown(summary, cases_df, out_dir, tune_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
