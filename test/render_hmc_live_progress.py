from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


TRACE_OBSERVABLES = ("squareOcc", "IPR", "doubleOcc", "nearestOcc")
INFO_VALUE_RE = re.compile(r"^\s*(.*?)\s*:\s*([0-9Ee+\-.]+)\s*$")


def parse_info_values(info_path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    if not info_path.exists():
        return values
    for line in info_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = INFO_VALUE_RE.match(line)
        if not match:
            continue
        key = match.group(1).strip()
        try:
            values[key] = float(match.group(2))
        except ValueError:
            continue
    return values


def read_trace(path: Path) -> list[float]:
    if not path.exists():
        return []
    values: list[float] = []
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            values.append(float(line))
    return values


def discover_run_dirs(work_root: Path) -> list[Path]:
    runs_root = work_root / "runs"
    if not runs_root.exists():
        return []
    return sorted(path for path in runs_root.glob("*/*") if path.is_dir())


def partial_trace_stats(values: list[float], thermal_cut: int) -> dict[str, float | int | str]:
    if not values:
        return {
            "longest_trace": 0,
            "window_mode": "empty",
            "span": 0.0,
            "drift": 0.0,
            "drift_ratio": 0.0,
        }
    if thermal_cut > 0 and len(values) > thermal_cut:
        tail = values[thermal_cut:]
        window_mode = "post-cut"
    else:
        tail = values
        window_mode = "pre-cut"
    window = max(16, len(tail) // 4)
    start = tail[:window]
    end = tail[-window:]
    start_mean = sum(start) / len(start)
    end_mean = sum(end) / len(end)
    span = max(tail) - min(tail) if len(tail) > 1 else 0.0
    drift = end_mean - start_mean
    drift_ratio = abs(drift) / span if span > 0.0 else 0.0
    return {
        "longest_trace": len(values),
        "window_mode": window_mode,
        "span": span,
        "drift": drift,
        "drift_ratio": drift_ratio,
    }


def render_case_plot(case_name: str, run_dirs: list[Path], thermal_cut: int, out_dir: Path) -> str:
    fig, axes = plt.subplots(len(TRACE_OBSERVABLES), 1, figsize=(10, 2.8 * len(TRACE_OBSERVABLES)), sharex=True)
    if len(TRACE_OBSERVABLES) == 1:
        axes = [axes]

    max_len = 0
    warm_values: list[int] = []
    for run_dir in run_dirs:
        info_values = parse_info_values(run_dir / "info.txt")
        if "# Warm" in info_values:
            warm_values.append(int(info_values["# Warm"]))
        label = run_dir.name
        for ax, observable in zip(axes, TRACE_OBSERVABLES):
            values = read_trace(run_dir / observable)
            if not values:
                continue
            max_len = max(max_len, len(values))
            ax.plot(range(len(values)), values, linewidth=1.2, label=label)
            ax.set_ylabel(observable)
            ax.grid(alpha=0.2)

    warm_marker = max(warm_values) if warm_values else None
    for i, ax in enumerate(axes):
        if warm_marker is not None:
            ax.axvline(warm_marker, color="#666666", linestyle=":", linewidth=1.0, label="warm" if i == 0 else None)
        if thermal_cut > 0:
            ax.axvline(thermal_cut, color="#cc3311", linestyle="--", linewidth=1.0, label="thermal_cut" if i == 0 else None)
        if i == 0:
            ax.set_title(f"{case_name} live progress")
            ax.legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("Sample index")
    fig.tight_layout()
    filename = f"live_trace_{case_name}.png"
    fig.savefig(out_dir / filename, dpi=180)
    plt.close(fig)
    return filename


def write_markdown(
    work_root: Path,
    cases: dict[str, list[Path]],
    thermal_cut: int,
    out_dir: Path,
    trace_files: list[tuple[str, str]],
    case_summaries: dict[str, dict[str, object]],
) -> None:
    lines = [
        "# HMC Live Progress",
        "",
        f"Work root: `{work_root}`",
        "",
        "This report is intended for in-flight stage runs that do not yet have a completed repeat.",
        "The dashed red line marks the configured `thermal_cut`; the dotted gray line marks the configured warm-up bins when it can be read from `info.txt`.",
        "If the longest trace has not crossed the configured cut yet, the drift ratios below are computed on the currently available pre-cut samples.",
        "",
        "| Case | Runs seen | Longest trace | Window | squareOcc drift/span | IPR drift/span |",
        "| --- | ---: | ---: | --- | ---: | ---: |",
    ]
    for case_name, run_dirs in cases.items():
        summary = case_summaries[case_name]
        lines.append(
            "| {case_name} | {runs_seen} | {longest_trace} | {window_mode} | {square_ratio:.3f} | {ipr_ratio:.3f} |".format(
                case_name=case_name,
                runs_seen=len(run_dirs),
                longest_trace=int(summary["longest_trace"]),
                window_mode=str(summary["window_mode"]),
                square_ratio=float(summary["square_ratio"]),
                ipr_ratio=float(summary["ipr_ratio"]),
            )
        )
    lines.extend(["", "## Traces", ""])
    for case_name, filename in trace_files:
        lines.extend([f"### {case_name}", "", f"![{case_name}]({filename})", ""])
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render live trace plots from in-flight HMC stage runs.")
    parser.add_argument("--work-root", required=True, help="Stage work root containing runs/.")
    parser.add_argument("--output-dir", default="", help="Output directory for live-progress figures and Markdown.")
    parser.add_argument("--thermal-cut", type=int, default=0, help="Configured thermal-cut marker to draw.")
    args = parser.parse_args()

    work_root = Path(args.work_root).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else work_root / "live_progress"
    output_dir.mkdir(parents=True, exist_ok=True)

    run_dirs = discover_run_dirs(work_root)
    if not run_dirs:
        raise FileNotFoundError(f"No run directories found under {work_root / 'runs'}")

    cases: dict[str, list[Path]] = defaultdict(list)
    for run_dir in run_dirs:
        cases[run_dir.parent.name].append(run_dir)

    trace_files: list[tuple[str, str]] = []
    case_summaries: dict[str, dict[str, object]] = {}
    for case_name, case_run_dirs in sorted(cases.items()):
        filename = render_case_plot(case_name, case_run_dirs, args.thermal_cut, output_dir)
        trace_files.append((case_name, filename))
        best_square: list[float] = []
        best_ipr: list[float] = []
        longest = -1
        for run_dir in case_run_dirs:
            square_vals = read_trace(run_dir / "squareOcc")
            if len(square_vals) > longest:
                longest = len(square_vals)
                best_square = square_vals
                best_ipr = read_trace(run_dir / "IPR")
        square_stats = partial_trace_stats(best_square, args.thermal_cut)
        ipr_stats = partial_trace_stats(best_ipr, args.thermal_cut)
        case_summaries[case_name] = {
            "longest_trace": square_stats["longest_trace"],
            "window_mode": square_stats["window_mode"],
            "square_ratio": square_stats["drift_ratio"],
            "ipr_ratio": ipr_stats["drift_ratio"],
        }
    write_markdown(work_root, cases, args.thermal_cut, output_dir, trace_files, case_summaries)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
