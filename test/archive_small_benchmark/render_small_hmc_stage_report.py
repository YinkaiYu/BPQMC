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
TRACE_PLOTS = ("doubleOcc", "nearestOcc", "IPR", "SF_Gamma")


def load_tables(bench_roots: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cases = []
    observables = []
    samples = []
    for root in bench_roots:
        case_df = pd.read_csv(root / "small_benchmark_cases.csv")
        obs_df = pd.read_csv(root / "small_benchmark_observables.csv")
        sample_df = pd.read_csv(root / "small_benchmark_samples.csv")
        case_df["case_id"] = case_df["name"].astype(str) + "__" + root.name
        obs_df["case_id"] = obs_df["name"].astype(str) + "__" + root.name
        sample_df["case_id"] = sample_df["name"].astype(str) + "__" + root.name
        case_df["benchmark_root"] = str(root)
        obs_df["benchmark_root"] = str(root)
        sample_df["benchmark_root"] = str(root)
        repeat_rows = []
        for case_name in case_df["name"].astype(str).unique():
            sample_repeats = int(sample_df[sample_df["name"] == case_name]["repeat"].nunique())
            run_root = root / "runs" / case_name
            run_repeats = len(list(run_root.glob("local_rep*"))) if run_root.exists() else 0
            repeat_rows.append({"name": case_name, "num_repeats": max(sample_repeats, run_repeats)})
        repeat_df = pd.DataFrame(repeat_rows)
        case_df = case_df.merge(repeat_df, on="name", how="left")
        cases.append(case_df)
        observables.append(obs_df)
        samples.append(sample_df)
    case_df = pd.concat(cases, ignore_index=True)
    obs_df = pd.concat(observables, ignore_index=True)
    sample_df = pd.concat(samples, ignore_index=True)
    return case_df, obs_df, sample_df


def enrich_tables(cases: pd.DataFrame, observables: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    observables = observables.copy()
    span_rows = []
    for (nbos, obs_name), group in observables.groupby(["Nbos", "observable"]):
        mid = 0.5 * (group["local_mean"] + group["hmc_mean"])
        span = float(mid.max() - mid.min())
        span_rows.append(
            {
                "Nbos": nbos,
                "observable": obs_name,
                "physics_span_same_nbos": span,
            }
        )
    span_df = pd.DataFrame(span_rows)
    observables = observables.merge(span_df, on=["Nbos", "observable"], how="left")
    observables["diff_over_span_same_nbos"] = observables.apply(
        lambda row: 0.0
        if abs(float(row["physics_span_same_nbos"])) <= 1.0e-14
        else float(row["abs_diff"]) / float(row["physics_span_same_nbos"]),
        axis=1,
    )

    case_rows = []
    for case_id, group in observables.groupby("case_id"):
        worst_z = group.sort_values("z_score", ascending=False).iloc[0]
        worst_span = group.sort_values("diff_over_span_same_nbos", ascending=False).iloc[0]
        fail_count = int((~group["passed"]).sum())
        case_rows.append(
            {
                "case_id": case_id,
                "worst_z_observable": worst_z["observable"],
                "worst_z_score": worst_z["z_score"],
                "worst_span_observable": worst_span["observable"],
                "worst_diff_over_span": worst_span["diff_over_span_same_nbos"],
                "failed_observable_count": fail_count,
            }
        )
    case_df = pd.DataFrame(case_rows)
    cases = cases.merge(case_df, on="case_id", how="left")
    return cases, observables


def save_overview(cases: pd.DataFrame, out_dir: Path) -> None:
    cases = cases.sort_values(["Nbos", "U2"]).reset_index(drop=True)
    labels = [f"N={int(row.Nbos)} U2={row.U2:g}" for row in cases.itertuples()]
    x = range(len(cases))

    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

    axes[0].bar(x, cases["acceptance_mean"], color="#3b6fb6")
    axes[0].axhspan(0.65, 0.85, color="#d7ebff", alpha=0.8)
    axes[0].set_ylabel("Acceptance")
    if int(cases["passed"].sum()) == len(cases):
        axes[0].set_title("Strict Benchmark Overview")
    else:
        axes[0].set_title("Benchmark Progress Overview")

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


def save_agreement_overview(cases: pd.DataFrame, out_dir: Path) -> None:
    cases = cases.sort_values(["Nbos", "U2"]).reset_index(drop=True)
    labels = [f"N={int(row.Nbos)} U2={row.U2:g}" for row in cases.itertuples()]
    x = range(len(cases))
    width = 0.38

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axes[0].bar(x, cases["worst_z_score"], color="#7a5195")
    axes[0].axhline(2.0, color="#333333", linestyle="--", linewidth=1.0)
    axes[0].set_ylabel("Max z-score")
    axes[0].set_title("Worst local-vs-HMC mismatch per case")

    axes[1].bar([idx - width / 2 for idx in x], cases["worst_diff_over_span"], width, color="#ef5675", label="Worst |diff| / span")
    axes[1].bar([idx + width / 2 for idx in x], cases["failed_observable_count"], width, color="#ffa600", label="Failed observable count")
    axes[1].set_ylabel("Relative mismatch / count")
    axes[1].legend()
    axes[1].set_xticks(list(x), labels, rotation=20, ha="right")

    fig.tight_layout()
    fig.savefig(out_dir / "agreement_overview.png", dpi=180)
    plt.close(fig)


def save_observable_curves(cases: pd.DataFrame, observables: pd.DataFrame, out_dir: Path) -> list[str]:
    saved = []
    merged = observables.merge(cases[["case_id", "name", "Nbos", "U2"]], on=["case_id", "name", "Nbos", "U2"], how="left")
    nbos_values = sorted(merged["Nbos"].unique())
    color_map = {nbos: plt.get_cmap("tab10")(idx % 10) for idx, nbos in enumerate(nbos_values)}
    for obs_name in OBS_PLOTS:
        obs_df = merged[merged["observable"] == obs_name].sort_values(["Nbos", "U2"])
        if obs_df.empty:
            continue
        for nbos, group in obs_df.groupby("Nbos"):
            fig, ax = plt.subplots(figsize=(7.6, 4.8))
            color = color_map[nbos]
            ax.errorbar(
                group["U2"],
                group["local_mean"],
                yerr=group["local_err"],
                marker="o",
                linestyle="-",
                label="Local",
                color=color,
            )
            ax.errorbar(
                group["U2"],
                group["hmc_mean"],
                yerr=group["hmc_err"],
                marker="s",
                linestyle="--",
                label="HMC",
                color=color,
                alpha=0.95,
            )
            ax.set_xscale("log")
            ax.set_xlabel("U2")
            ax.set_ylabel(obs_name)
            ax.set_title(f"{obs_name} vs U2, Nbos={int(nbos)}")
            ax.legend()
            fig.tight_layout()
            filename = f"observable_{obs_name}_n{int(nbos)}.png"
            fig.savefig(out_dir / filename, dpi=180)
            plt.close(fig)
            saved.append(filename)
    return saved


def save_trace_plots(cases: pd.DataFrame, samples: pd.DataFrame, out_dir: Path) -> list[str]:
    saved = []
    for row in cases.itertuples():
        case_samples = samples[samples["case_id"] == row.case_id]
        if case_samples.empty:
            continue
        for repeat in sorted(case_samples["repeat"].unique()):
            rep_samples = case_samples[case_samples["repeat"] == repeat]
            fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
            axes = axes.flatten()
            for ax, obs_name in zip(axes, TRACE_PLOTS):
                obs_df = rep_samples[rep_samples["observable"] == obs_name]
                for mode_name, color in (("local", "#666666"), ("hmc", "#c44e52")):
                    mode_df = obs_df[obs_df["mode"] == mode_name]
                    ax.plot(mode_df["sample_index"], mode_df["value"], lw=1.0, label=mode_name, color=color)
                ax.axvline(row.thermal_cut, color="#1f77b4", linestyle=":", lw=1.0)
                ax.set_title(obs_name)
            axes[0].legend()
            fig.suptitle(
                f"Thermalization traces: N={int(row.Nbos)}, U2={row.U2:g}, repeat={int(repeat)}, "
                f"run={Path(row.benchmark_root).name}"
            )
            fig.tight_layout()
            filename = f"trace_{row.case_id}_rep{int(repeat)}.png"
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
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

        sc0 = axes[0].scatter(case_tune["traj_len"], case_tune["acceptance_mean"], c=case_tune["hmc_mass"], cmap="viridis", s=80)
        for tune_row in case_tune.itertuples():
            axes[0].annotate(
                f"{int(tune_row.nfrog)}x{tune_row.hmc_dt:g}\nm={tune_row.hmc_mass:g}",
                (tune_row.traj_len, tune_row.acceptance_mean),
                fontsize=7,
            )
        axes[0].set_xlabel("Trajectory length = Nfrog * dt")
        axes[0].set_ylabel("Acceptance")
        axes[0].set_title("Acceptance")
        fig.colorbar(sc0, ax=axes[0], label="Mass")

        sc1 = axes[1].scatter(case_tune["traj_len"], case_tune["tau_int_doubleOcc_mean"], c=case_tune["acceptance_mean"], cmap="plasma", s=80)
        for tune_row in case_tune.itertuples():
            axes[1].annotate(
                f"{int(tune_row.nfrog)}x{tune_row.hmc_dt:g}\nm={tune_row.hmc_mass:g}",
                (tune_row.traj_len, tune_row.tau_int_doubleOcc_mean),
                fontsize=7,
            )
        axes[1].set_xlabel("Trajectory length = Nfrog * dt")
        axes[1].set_ylabel("tau_int(doubleOcc)")
        axes[1].set_title("Autocorrelation")
        fig.colorbar(sc1, ax=axes[1], label="Acceptance")

        sc2 = axes[2].scatter(case_tune["traj_len"], case_tune["ess_per_sec_doubleOcc_mean"], c=case_tune["acceptance_mean"], cmap="cividis", s=80)
        for tune_row in case_tune.itertuples():
            axes[2].annotate(
                f"{int(tune_row.nfrog)}x{tune_row.hmc_dt:g}\nm={tune_row.hmc_mass:g}",
                (tune_row.traj_len, tune_row.ess_per_sec_doubleOcc_mean),
                fontsize=7,
            )
        axes[2].set_xlabel("Trajectory length = Nfrog * dt")
        axes[2].set_ylabel("ESS / sec")
        axes[2].set_title("Efficiency")
        fig.colorbar(sc2, ax=axes[2], label="Acceptance")

        fig.suptitle(f"Tuning evidence: {row.name} ({Path(row.benchmark_root).name})")
        fig.tight_layout()
        filename = f"tune_{row.case_id}.png"
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)
        saved.append(filename)
    return saved


def write_markdown(
    cases: pd.DataFrame,
    observables: pd.DataFrame,
    observable_figs: list[str],
    trace_figs: list[str],
    tune_figs: list[str],
    tune_df: pd.DataFrame | None,
    out_dir: Path,
) -> None:
    strict_mask = cases["num_repeats"].fillna(0).astype(int) >= 2
    strict_total = int(strict_mask.sum())
    strict_passed = int((strict_mask & cases["passed"]).sum())
    probe_total = int((~strict_mask).sum())
    lattice_values = sorted(cases["name"].str.split("_").str[0].unique())
    l_values = ", ".join(str(int(value)) for value in sorted(cases["L"].unique()))
    beta_values = ", ".join(f"{value:g}" for value in sorted(cases["beta"].unique()))
    dtau_values = ", ".join(f"{value:g}" for value in sorted(cases["dtau"].unique()))
    nbos_values = ", ".join(str(int(value)) for value in sorted(cases["Nbos"].unique()))
    u2_values = ", ".join(f"{value:g}" for value in sorted(cases["U2"].unique()))
    overview_cases = cases.sort_values(["Nbos", "U2"])
    lines = [
        "# Small Triangular HMC Stage Report",
        "",
        f"Strict benchmark cases passed: `{strict_passed}/{strict_total}`",
        "",
        f"Single-repeat visual probes: `{probe_total}`",
        "",
        "## Campaign Parameters",
        "",
        f"- Lattice: `{', '.join(lattice_values)}`",
        f"- `L`: `{l_values}`",
        f"- `beta`: `{beta_values}`",
        f"- `dtau`: `{dtau_values}`",
        f"- `U1`: `0`",
        f"- `Nbos` covered in this stage report: `{nbos_values}`",
        f"- `U2` covered in this stage report: `{u2_values}`",
        "",
        "## Case Summary",
        "",
    ]
    if strict_passed < strict_total or probe_total > 0:
        lines.extend(
            [
                "This report includes unresolved or still-failing benchmark variants, plus visual single-repeat probes.",
                "Use the case table and the sample-index traces to see whether the mismatch looks like a slow drift, a thermal-cut issue, or a persistent sampler bias.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "This report only aggregates strict local-vs-HMC benchmark directories that already passed the current health criterion.",
                "Acceptance is shown as a diagnostic, but the actual health criterion is agreement of post-cut observable means within the combined statistical error bars.",
                "",
            ]
        )
    lines.extend(
        [
            "| Case | Repeats | Acceptance | tau_local | tau_hmc | ESS/sec local | ESS/sec HMC | Speed ratio | Result |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in overview_cases.itertuples():
        if int(row.num_repeats) < 2:
            result = "PROBE"
        else:
            result = "PASS" if row.passed else "FAIL"
        lines.append(
            "| {name} ({root}) | {repeats} | {acc:.3f} | {tau_l:.3f} | {tau_h:.3f} | {ess_l:.3f} | {ess_h:.3f} | {ratio:.3f} | {result} |".format(
                name=row.name,
                root=Path(row.benchmark_root).name,
                repeats=int(row.num_repeats),
                acc=row.acceptance_mean,
                tau_l=row.local_tau_int_doubleOcc,
                tau_h=row.hmc_tau_int_doubleOcc,
                ess_l=row.local_ess_per_sec_doubleOcc,
                ess_h=row.hmc_ess_per_sec_doubleOcc,
                ratio=row.speed_ratio_hmc_over_local,
                result=result,
            )
        )
    lines.extend(
        [
            "",
            "## Metric Definitions",
            "",
            "- `tau_int(doubleOcc)`: integrated autocorrelation time estimated from the post-cut `doubleOcc` series. Larger values mean slower mixing.",
            "- `ESS / sec`: effective sample size per wall-clock second, estimated from the same `doubleOcc` series. Larger values mean more statistically independent samples per unit time.",
            "- `Speed ratio`: `ESS/sec(HMC) / ESS/sec(Local)`.",
            "- `z-score`: `|mean_local - mean_hmc| / sqrt(err_local^2 + err_hmc^2)`. Here `err_local` and `err_hmc` are the standard errors of the post-cut repeat means, not naive per-sample errors.",
            "- `|diff| / span`: the observable mismatch divided by the observable span across the available `U2` values at the same `Nbos`, using the midpoint `(local + hmc)/2` at each point. This is a visual/physics-scale diagnostic, not a strict pass criterion.",
            "",
            "## Overview",
            "",
            "The top panel shows the mean HMC acceptance for each strict health point.",
            "The middle panel compares the integrated autocorrelation time of `doubleOcc`; smaller is better.",
            "The bottom panel compares `ESS/sec`; larger is better.",
            "At this stage, the main purpose of these three panels is to separate correctness from efficiency: a point can be correct even when HMC is slower than local.",
            "",
            "![overview](overview.png)",
            "",
            "## Agreement Overview",
            "",
            "The first panel below shows the worst observable z-score for each case.",
            "The dashed line marks the current strict `2 sigma` benchmark cut.",
            "The second panel adds two softer diagnostics: the worst mismatch as a fraction of the same-`Nbos` observable span across `U2`, and the count of observables that fail the strict cut.",
            "This is useful when a case looks visually healthy but still has one or two borderline `z-score` failures.",
            "",
            "![agreement_overview](agreement_overview.png)",
            "",
            "## Case Diagnostics",
            "",
            "For unresolved cases, the tables below list the observables that still sit outside the combined error-bar criterion.",
            "The `z-score` is the absolute local-vs-HMC mean mismatch divided by the combined standard error.",
            "",
        ]
    )
    for row in overview_cases.itertuples():
        case_obs = observables[observables["case_id"] == row.case_id].copy()
        fail_df = case_obs[~case_obs["passed"]].sort_values("z_score", ascending=False)
        lines.extend([f"### {row.name} ({Path(row.benchmark_root).name})", ""])
        is_probe = int(row.num_repeats) < 2
        if is_probe:
            lines.extend(
                [
                    "This is a single-repeat probe. Its trace plots are still useful,",
                    "but the reported standard errors and z-scores are not reliable strict-benchmark statistics.",
                    "",
                ]
            )
        lines.extend(
            [
                "Case summary:",
                "- Benchmark parameters: "
                f"`L={int(row.L)}`, `beta={row.beta:g}`, `dtau={row.dtau:g}`, `Nbos={int(row.Nbos)}`, `U2={row.U2:g}`, "
                f"`nfrog={int(row.nfrog)}`, `dt={row.hmc_dt:g}`, `jitter={int(row.hmc_jitter)}`, `mass={row.hmc_mass:g}`.",
                f"- Worst |diff| / span: `{row.worst_span_observable}` = `{row.worst_diff_over_span:.3f}`.",
                f"- Failed observable count: `{int(row.failed_observable_count)}`.",
                "",
            ]
        )
        if is_probe:
            lines.extend(
                [
                    "- Worst z-score: suppressed for single-repeat probes because the repeat-to-repeat standard error is not defined.",
                    "",
                ]
            )
        else:
            lines.extend([f"- Worst z-score: `{row.worst_z_observable}` = `{row.worst_z_score:.3f}`.", ""])
        if fail_df.empty:
            lines.extend(["All tracked observables pass the current strict benchmark criterion.", ""])
            continue
        if is_probe:
            probe_df = case_obs.sort_values("diff_over_span_same_nbos", ascending=False).head(8)
            lines.extend(
                [
                    "| Observable | |diff| / span | Abs diff | Local mean | HMC mean |",
                    "| --- | ---: | ---: | ---: | ---: |",
                ]
            )
            for fail in probe_df.itertuples():
                lines.append(
                    "| {obs} | {span:.3f} | {diff:.6g} | {local:.6g} | {hmc:.6g} |".format(
                        obs=fail.observable,
                        span=fail.diff_over_span_same_nbos,
                        diff=fail.abs_diff,
                        local=fail.local_mean,
                        hmc=fail.hmc_mean,
                    )
                )
        else:
            lines.extend(
                [
                    "| Observable | z-score | |diff| / span | Local mean | HMC mean | Local err | HMC err |",
                    "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
                ]
            )
            for fail in fail_df.itertuples():
                lines.append(
                    "| {obs} | {z:.3f} | {span:.3f} | {local:.6g} | {hmc:.6g} | {local_err:.3g} | {hmc_err:.3g} |".format(
                        obs=fail.observable,
                        z=fail.z_score,
                        span=fail.diff_over_span_same_nbos,
                        local=fail.local_mean,
                        hmc=fail.hmc_mean,
                        local_err=fail.local_err,
                        hmc_err=fail.hmc_err,
                    )
                )
        lines.extend([""])
    lines.extend(
        [
            "## Observable Curves",
            "",
        ]
    )
    for filename in observable_figs:
        lines.extend([f"![{filename}]({filename})", ""])
    lines.extend(
        [
            "Each observable figure now fixes one `Nbos` and compares only the two samplers on the same axes.",
            "This avoids the scale-compression problem that happens when different `Nbos` values are mixed in one panel.",
            "",
            "## Thermalization Traces",
            "",
            "Each trace figure overlays the local and HMC sample histories for one strict benchmark case.",
            "The vertical dotted line marks the thermal cut used when computing the benchmark means and error bars.",
            "",
        ]
    )
    for filename in trace_figs:
        lines.extend([f"![{filename}]({filename})", ""])
    if tune_figs:
        lines.extend(
            [
                "## Tuning Evidence",
                "",
                "Each tuning panel corresponds to one physical point.",
                "Every marker is one HMC candidate from the short warm=0 scan.",
                "The x-axis is the MD trajectory length `Nfrog * dt`; changing this moves the proposal across different fractions of the effective auxiliary-field oscillation period.",
                "The left panel is acceptance, the middle panel is `tau_int(doubleOcc)`, and the right panel is `ESS/sec(doubleOcc)`.",
                "Mass values are shown both by color and by the text annotations next to the markers.",
                "These plots are meant to show why a candidate was chosen before the long strict benchmark was launched.",
                "",
            ]
        )
        for filename in tune_figs:
            lines.extend([f"![{filename}]({filename})", ""])
        if tune_df is not None:
            lines.extend([""])
            for row in overview_cases.itertuples():
                case_tune = tune_df[tune_df["name"] == row.name].copy()
                if case_tune.empty:
                    continue
                case_tune["traj_len"] = case_tune["nfrog"] * case_tune["hmc_dt"]
                lines.extend(
                    [
                        f"### {row.name} ({Path(row.benchmark_root).name})",
                        "",
                        "| Nfrog | dt | mass | traj_len | acceptance | tau_int | ESS/sec |",
                        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
                    ]
                )
                for tune_row in case_tune.sort_values(["hmc_mass", "traj_len"]).itertuples():
                    lines.append(
                        "| {nfrog} | {dt:.6g} | {mass:.3g} | {traj:.6g} | {acc:.3f} | {tau:.3f} | {ess:.3f} |".format(
                            nfrog=int(tune_row.nfrog),
                            dt=tune_row.hmc_dt,
                            mass=tune_row.hmc_mass,
                            traj=tune_row.traj_len,
                            acc=tune_row.acceptance_mean,
                            tau=tune_row.tau_int_doubleOcc_mean,
                            ess=tune_row.ess_per_sec_doubleOcc_mean,
                        )
                    )
                lines.extend([""])
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

    cases, observables = enrich_tables(cases, observables)

    cases.to_csv(out_dir / "stage_cases.csv", index=False)
    observables.to_csv(out_dir / "stage_observables.csv", index=False)
    samples.to_csv(out_dir / "stage_samples.csv", index=False)
    summary = {
        "cases": json.loads(cases.to_json(orient="records")),
        "observables": json.loads(observables.to_json(orient="records")),
        "bench_roots": [str(path) for path in bench_roots],
        "tune_csv": args.tune_csv,
    }
    (out_dir / "stage_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    save_overview(cases, out_dir)
    save_agreement_overview(cases, out_dir)
    observable_figs = save_observable_curves(cases, observables, out_dir)
    trace_figs = save_trace_plots(cases, samples, out_dir)
    tune_figs = save_tune_plots(cases, tune_df, out_dir)
    write_markdown(cases, observables, observable_figs, trace_figs, tune_figs, tune_df, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
