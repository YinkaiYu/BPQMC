#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


OBS_PLOTS = ("doubleOcc", "squareOcc", "nearestOcc", "IPR", "PF_Gamma", "kinetic")
TRACE_PLOTS = ("doubleOcc", "squareOcc", "nearestOcc", "IPR")


def load_tables(bench_roots: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cases = []
    observables = []
    samples = []
    for root in bench_roots:
        case_df = pd.read_csv(root / "small_benchmark_cases.csv")
        obs_df = pd.read_csv(root / "small_benchmark_observables.csv")
        sample_df = pd.read_csv(root / "small_benchmark_samples.csv")
        case_df["benchmark_root"] = str(root)
        obs_df["benchmark_root"] = str(root)
        sample_df["benchmark_root"] = str(root)
        cases.append(case_df)
        observables.append(obs_df)
        samples.append(sample_df)
    return pd.concat(cases, ignore_index=True), pd.concat(observables, ignore_index=True), pd.concat(samples, ignore_index=True)


def save_overview(cases: pd.DataFrame, out_dir: Path) -> None:
    cases = cases.sort_values(["Nbos", "U2"]).reset_index(drop=True)
    labels = [f"N={int(row.Nbos)} U2={row.U2:g}" for row in cases.itertuples()]
    x = range(len(cases))

    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

    axes[0].bar(x, cases["acceptance_mean"], color="#3b6fb6")
    axes[0].axhspan(0.65, 0.85, color="#d7ebff", alpha=0.8)
    axes[0].set_ylabel("Acceptance")
    axes[0].set_title("Strict Benchmark Overview")

    width = 0.38
    axes[1].bar([idx - width / 2 for idx in x], cases["local_tau_int_doubleOcc"], width, label="Local", color="#999999")
    axes[1].bar([idx + width / 2 for idx in x], cases["hmc_tau_int_doubleOcc"], width, label="HMC", color="#009e73")
    axes[1].set_ylabel("tau_int(doubleOcc)")
    axes[1].legend()

    axes[2].bar([idx - width / 2 for idx in x], cases["local_ess_per_sec_doubleOcc"], width, label="Local", color="#999999")
    axes[2].bar([idx + width / 2 for idx in x], cases["hmc_ess_per_sec_doubleOcc"], width, label="HMC", color="#d55e00")
    axes[2].set_ylabel("ESS / sec")
    axes[2].legend()
    axes[2].set_xticks(list(x), labels, rotation=20, ha="right")

    fig.tight_layout()
    fig.savefig(out_dir / "overview.png", dpi=180)
    plt.close(fig)


def save_observable_curves(cases: pd.DataFrame, observables: pd.DataFrame, out_dir: Path) -> list[str]:
    saved = []
    merged = observables.merge(cases[["name", "Nbos", "U2"]], on=["name", "Nbos", "U2"], how="left")
    for obs_name in OBS_PLOTS:
        obs_df = merged[merged["observable"] == obs_name].sort_values(["Nbos", "U2"])
        if obs_df.empty:
            continue
        fig, ax = plt.subplots(figsize=(8.5, 5.0))
        for nbos, group in obs_df.groupby("Nbos"):
            ax.errorbar(
                group["U2"],
                group["local_mean"],
                yerr=group["local_err"],
                marker="o",
                linestyle="-",
                label=f"Local N={int(nbos)}",
                color="#666666",
            )
            ax.errorbar(
                group["U2"],
                group["hmc_mean"],
                yerr=group["hmc_err"],
                marker="s",
                linestyle="--",
                label=f"HMC N={int(nbos)}",
                color="#c44e52",
            )
        ax.set_xscale("log")
        ax.set_xlabel("U2")
        ax.set_ylabel(obs_name)
        ax.set_title(f"{obs_name} vs U2")
        ax.legend(ncol=2, fontsize=9)
        fig.tight_layout()
        filename = f"observable_{obs_name}.png"
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)
        saved.append(filename)
    return saved


def save_trace_plots(cases: pd.DataFrame, samples: pd.DataFrame, out_dir: Path) -> list[str]:
    saved = []
    samples = samples[samples["repeat"] == 0]
    for row in cases.itertuples():
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
        axes = axes.flatten()
        case_samples = samples[samples["name"] == row.name]
        if case_samples.empty:
            plt.close(fig)
            continue
        for ax, obs_name in zip(axes, TRACE_PLOTS):
            obs_df = case_samples[case_samples["observable"] == obs_name]
            for mode_name, color in (("local", "#666666"), ("hmc", "#c44e52")):
                mode_df = obs_df[obs_df["mode"] == mode_name]
                ax.plot(mode_df["sample_index"], mode_df["value"], lw=1.0, label=mode_name, color=color)
            ax.axvline(row.thermal_cut, color="#1f77b4", linestyle=":", lw=1.0)
            ax.set_title(obs_name)
        axes[0].legend()
        fig.suptitle(f"Thermalization traces: N={int(row.Nbos)}, U2={row.U2:g}")
        fig.tight_layout()
        filename = f"trace_{row.name}.png"
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)
        saved.append(filename)
    return saved


def save_tune_plots(cases: pd.DataFrame, tune_df: pd.DataFrame | None, out_dir: Path) -> list[str]:
    saved = []
    if tune_df is None:
        return saved
    for row in cases.itertuples():
        case_tune = tune_df[tune_df["name"] == row.name].copy()
        if case_tune.empty:
            continue
        case_tune["traj_len"] = case_tune["nfrog"] * case_tune["hmc_dt"]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

        sc0 = axes[0].scatter(case_tune["traj_len"], case_tune["acceptance_mean"], c=case_tune["ess_per_sec_doubleOcc_mean"], cmap="viridis", s=70)
        for tune_row in case_tune.itertuples():
            axes[0].annotate(f"{int(tune_row.nfrog)}x{tune_row.hmc_dt:g}", (tune_row.traj_len, tune_row.acceptance_mean), fontsize=8)
        axes[0].set_xlabel("Trajectory length = Nfrog * dt")
        axes[0].set_ylabel("Acceptance")
        axes[0].set_title("Tune scan")
        fig.colorbar(sc0, ax=axes[0], label="ESS / sec")

        sc1 = axes[1].scatter(case_tune["traj_len"], case_tune["tau_int_doubleOcc_mean"], c=case_tune["acceptance_mean"], cmap="plasma", s=70)
        for tune_row in case_tune.itertuples():
            axes[1].annotate(f"{int(tune_row.nfrog)}x{tune_row.hmc_dt:g}", (tune_row.traj_len, tune_row.tau_int_doubleOcc_mean), fontsize=8)
        axes[1].set_xlabel("Trajectory length = Nfrog * dt")
        axes[1].set_ylabel("tau_int(doubleOcc)")
        axes[1].set_title("Tune scan")
        fig.colorbar(sc1, ax=axes[1], label="Acceptance")

        fig.suptitle(f"Tuning evidence: {row.name}")
        fig.tight_layout()
        filename = f"tune_{row.name}.png"
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)
        saved.append(filename)
    return saved


def write_markdown(
    cases: pd.DataFrame,
    observable_figs: list[str],
    trace_figs: list[str],
    tune_figs: list[str],
    out_dir: Path,
) -> None:
    passed = int(cases["passed"].sum())
    lines = [
        "# Small Triangular HMC Stage Report",
        "",
        f"Health points passed: `{passed}/{len(cases)}`",
        "",
        "| Case | Acceptance | tau_local | tau_hmc | ESS/sec local | ESS/sec HMC | Speed ratio | Result |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in cases.sort_values(["Nbos", "U2"]).itertuples():
        lines.append(
            "| {name} | {acc:.3f} | {tau_l:.3f} | {tau_h:.3f} | {ess_l:.3f} | {ess_h:.3f} | {ratio:.3f} | {result} |".format(
                name=row.name,
                acc=row.acceptance_mean,
                tau_l=row.local_tau_int_doubleOcc,
                tau_h=row.hmc_tau_int_doubleOcc,
                ess_l=row.local_ess_per_sec_doubleOcc,
                ess_h=row.hmc_ess_per_sec_doubleOcc,
                ratio=row.speed_ratio_hmc_over_local,
                result="PASS" if row.passed else "FAIL",
            )
        )
    lines.extend(
        [
            "",
            "## Overview",
            "",
            "![overview](overview.png)",
            "",
            "## Observable Curves",
            "",
        ]
    )
    for filename in observable_figs:
        lines.extend([f"![{filename}]({filename})", ""])
    lines.extend(["## Thermalization Traces", ""])
    for filename in trace_figs:
        lines.extend([f"![{filename}]({filename})", ""])
    if tune_figs:
        lines.extend(["## Tuning Evidence", ""])
        for filename in tune_figs:
            lines.extend([f"![{filename}]({filename})", ""])
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render a stage report from one or more small benchmark directories.")
    parser.add_argument("bench_roots", nargs="+", help="Benchmark directories containing small_benchmark_*.csv files.")
    parser.add_argument("--tune-csv", default="", help="Optional small_tune.csv for tuning-evidence plots.")
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()

    bench_roots = [Path(item).resolve() for item in args.bench_roots]
    out_dir = Path(args.output_dir).resolve() if args.output_dir else bench_roots[0] / "stage_report"
    out_dir.mkdir(parents=True, exist_ok=True)

    cases, observables, samples = load_tables(bench_roots)
    tune_df = pd.read_csv(args.tune_csv) if args.tune_csv else None

    cases.to_csv(out_dir / "stage_cases.csv", index=False)
    observables.to_csv(out_dir / "stage_observables.csv", index=False)
    samples.to_csv(out_dir / "stage_samples.csv", index=False)
    summary = {
        "cases": json.loads(cases.to_json(orient="records")),
        "bench_roots": [str(path) for path in bench_roots],
        "tune_csv": args.tune_csv,
    }
    (out_dir / "stage_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    save_overview(cases, out_dir)
    observable_figs = save_observable_curves(cases, observables, out_dir)
    trace_figs = save_trace_plots(cases, samples, out_dir)
    tune_figs = save_tune_plots(cases, tune_df, out_dir)
    write_markdown(cases, observable_figs, trace_figs, tune_figs, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
