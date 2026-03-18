#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / "data" / "triangular_hmc_small_benchmark"
ARCHIVE_ROOT = DATA_ROOT / "archive_20260319"
RENDER_SCRIPT = Path(__file__).resolve().with_name("render_small_hmc_stage_report.py")

DEFAULT_CASE_ROOTS = (
    "strict_healthy_n10_u1em2_v2",
    "strict_healthy_n10_u1em1_v1",
    "strict_healthy_n10_u1e0_v1",
    "strict_healthy_n10_u1e1_v1",
    "strict_healthy_n10_u1e2_v3",
    "strict_healthy_n100_u1em2_v2",
    "strict_healthy_n100_u1em1_v1",
    "strict_healthy_n100_u1e0_v1",
    "strict_seeded_n100_u1e1_v5_localseed",
    "probe_n100_u1e2_v3",
    "strict_healthy_n1000_u1em2_v1",
    "strict_healthy_n1000_u1em1_v1",
    "strict_seeded_n1000_u1e0_v4_localseed",
    "probe_n1000_u1e1_v2",
    "probe_n1000_u1e2_v2",
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render the fixed 3x5 triangular small-benchmark full-grid report.")
    parser.add_argument(
        "--output-dir",
        default=str(DATA_ROOT / "full_grid_progress_v3"),
        help="Output directory for the rendered report.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to run render_small_hmc_stage_report.py.",
    )
    parser.add_argument(
        "--case-root",
        action="append",
        default=[],
        help="Optional extra or replacement benchmark roots under data/triangular_hmc_small_benchmark.",
    )
    args = parser.parse_args()

    case_roots = []
    for item in (args.case_root or DEFAULT_CASE_ROOTS):
        path = Path(item)
        if not path.is_absolute():
            direct = DATA_ROOT / item
            archived = ARCHIVE_ROOT / item
            if (direct / "small_benchmark_cases.csv").exists():
                path = direct
            else:
                path = archived
        case_roots.append(path)
    missing = [path for path in case_roots if not (path / "small_benchmark_cases.csv").exists()]
    if missing:
        print("Missing benchmark directories:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 2

    cmd = [args.python, str(RENDER_SCRIPT), *(str(path) for path in case_roots), "--output-dir", str(Path(args.output_dir).resolve())]
    return subprocess.run(cmd, cwd=ROOT, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
