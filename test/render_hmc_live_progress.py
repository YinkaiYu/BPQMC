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


def write_markdown(work_root: Path, cases: dict[str, list[Path]], thermal_cut: int, out_dir: Path, trace_files: list[tuple[str, str]]) -> None:
    lines = [
        "# HMC Live Progress",
        "",
        f"Work root: `{work_root}`",
        "",
        "This report is intended for in-flight stage runs that do not yet have a completed repeat.",
        "The dashed red line marks the configured `thermal_cut`; the dotted gray line marks the configured warm-up bins when it can be read from `info.txt`.",
        "",
        "| Case | Runs seen | Longest trace |",
        "| --- | ---: | ---: |",
    ]
    for case_name, run_dirs in cases.items():
        longest = 0
        for run_dir in run_dirs:
            longest = max(longest, len(read_trace(run_dir / "squareOcc")))
        lines.append(f"| {case_name} | {len(run_dirs)} | {longest} |")
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
    for case_name, case_run_dirs in sorted(cases.items()):
        filename = render_case_plot(case_name, case_run_dirs, args.thermal_cut, output_dir)
        trace_files.append((case_name, filename))
    write_markdown(work_root, cases, args.thermal_cut, output_dir, trace_files)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
