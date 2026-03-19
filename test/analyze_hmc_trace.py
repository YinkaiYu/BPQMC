#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path


TRACE_KEYS_OLD = (
    "proposal",
    "stage",
    "lf_step",
    "nsteps",
    "md_time",
    "action",
    "phi_f1",
    "phi_f2",
    "mom_f1",
    "mom_f2",
    "force_f1",
    "force_f2",
    "phi_rms_f1",
    "phi_rms_f2",
)

TRACE_KEYS_NEW = (
    "proposal",
    "stage",
    "lf_step",
    "nsteps",
    "md_time",
    "action",
    "phi_f1",
    "phi_f2",
    "phi_mean_f1",
    "phi_mean_f2",
    "mom_f1",
    "mom_f2",
    "force_f1",
    "force_f2",
    "phi_rms_f1",
    "phi_rms_f2",
)


def load_trace(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == len(TRACE_KEYS_OLD):
            trace_keys = TRACE_KEYS_OLD
        elif len(parts) == len(TRACE_KEYS_NEW):
            trace_keys = TRACE_KEYS_NEW
        else:
            raise SystemExit(f"Unexpected column count in {path}: {line}")
        row: dict[str, float | int] = {}
        for key, value in zip(trace_keys, parts):
            if key in {"proposal", "stage", "lf_step", "nsteps"}:
                row[key] = int(value)
            else:
                row[key] = float(value)
        if "phi_mean_f1" not in row:
            row["phi_mean_f1"] = 0.0
            row["phi_mean_f2"] = 0.0
        rows.append(row)
    return rows


def select_rows(
    rows: list[dict[str, float | int]],
    *,
    stage: int | None,
    proposal: int | None,
) -> list[dict[str, float | int]]:
    filtered = rows
    if stage is not None:
        filtered = [row for row in filtered if int(row["stage"]) == stage]
    if proposal is None:
        if not filtered:
            return []
        proposal = int(filtered[0]["proposal"])
    return [row for row in filtered if int(row["proposal"]) == proposal]


def demean(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    return [value - mean for value in values]


def dominant_fft_periods(values: list[float], max_modes: int) -> list[tuple[float, int, float]]:
    npts = len(values)
    peaks: list[tuple[float, int, float]] = []
    for kfreq in range(1, min(npts // 2, 80)):
        re_part = 0.0
        im_part = 0.0
        for idx, value in enumerate(values):
            angle = 2.0 * math.pi * kfreq * idx / npts
            re_part += value * math.cos(angle)
            im_part -= value * math.sin(angle)
        amplitude = math.sqrt(re_part * re_part + im_part * im_part)
        peaks.append((amplitude, kfreq, npts / kfreq))
    peaks.sort(reverse=True)
    return peaks[:max_modes]


def zero_crossings(values: list[float]) -> list[float]:
    zeros: list[float] = []
    for idx in range(1, len(values)):
        prev = values[idx - 1]
        cur = values[idx]
        if prev == 0.0:
            zeros.append(float(idx - 1))
        elif prev * cur < 0.0:
            frac = abs(prev) / (abs(prev) + abs(cur))
            zeros.append((idx - 1) + frac)
    return zeros


def estimate_period_from_zero_cross(values: list[float], dt: float) -> tuple[float, float] | None:
    crossings = zero_crossings(values)
    if len(crossings) < 3:
        return None
    periods_steps = [crossings[idx] - crossings[idx - 2] for idx in range(2, len(crossings))]
    avg_steps = sum(periods_steps) / len(periods_steps)
    return avg_steps, avg_steps * dt


def linear_fit(xs: list[float], ys: list[float]) -> tuple[float, float] | None:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denom = sum((x - mean_x) ** 2 for x in xs)
    if denom == 0.0:
        return None
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denom
    intercept = mean_y - slope * mean_x
    return intercept, slope


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze a trajectory-level HMC trace.")
    parser.add_argument("trace", type=Path, help="Path to hmc_trace.dat")
    parser.add_argument(
        "--series",
        choices=("phi_f1", "phi_f2", "phi_mean_f1", "phi_mean_f2", "mom_f1", "mom_f2", "force_f1", "force_f2"),
        default="phi_f2",
    )
    parser.add_argument("--stage", type=int, default=None, help="Optional stage filter: 0=warm, 1=measurement.")
    parser.add_argument("--proposal", type=int, default=None, help="Proposal index to analyze. Defaults to the first matching proposal.")
    parser.add_argument("--top-k", type=int, default=5, help="Number of FFT peaks to print.")
    args = parser.parse_args()

    rows = load_trace(args.trace)
    proposal_rows = select_rows(rows, stage=args.stage, proposal=args.proposal)
    if len(proposal_rows) < 8:
        raise SystemExit("Not enough trace rows for analysis.")

    proposal = int(proposal_rows[0]["proposal"])
    stage = int(proposal_rows[0]["stage"])
    nsteps = int(proposal_rows[0]["nsteps"])
    md_times = [float(row["md_time"]) for row in proposal_rows]
    values = [float(row[args.series]) for row in proposal_rows]
    centered = demean(values)
    dt = 0.0 if len(md_times) < 2 else md_times[1] - md_times[0]

    peaks = dominant_fft_periods(centered, args.top_k)
    zero_cross_period = estimate_period_from_zero_cross(centered, dt)
    fit_result = None
    if args.series in {"phi_f1", "phi_f2"}:
        force_name = "force_f1" if args.series == "phi_f1" else "force_f2"
        forces = [float(row[force_name]) for row in proposal_rows]
        fit_result = linear_fit(values, forces)

    print(f"proposal={proposal} stage={stage} rows={len(proposal_rows)} nsteps={nsteps}")
    print(f"series={args.series} md_dt={dt:.16e}")
    print(f"start={values[0]:.16e} end={values[-1]:.16e}")
    print(f"min={min(values):.16e} max={max(values):.16e}")
    if zero_cross_period is None:
        print("zero_cross_period: unavailable")
    else:
        steps, md_period = zero_cross_period
        print(f"zero_cross_period_steps={steps:.6f} zero_cross_period_md={md_period:.16e}")
        print(f"quarter_period_md={0.25 * md_period:.16e}")
    if fit_result is not None:
        intercept, slope = fit_result
        print(f"force_fit: force = {intercept:.16e} + ({slope:.16e}) * phi")
        if slope < 0.0:
            k_eff = -slope
            omega_eff = math.sqrt(k_eff)
            period_eff = 2.0 * math.pi / omega_eff
            print(f"k_eff={k_eff:.16e} omega_eff={omega_eff:.16e} period_eff_md={period_eff:.16e} quarter_eff_md={0.25 * period_eff:.16e}")
        else:
            print("k_eff: unavailable (non-restoring linear fit)")
    print("top_fft_peaks:")
    for amplitude, kfreq, period_steps in peaks:
        period_md = period_steps * dt
        print(
            f"  k={kfreq:3d} amp={amplitude:.16e} "
            f"period_steps={period_steps:.6f} period_md={period_md:.16e} "
            f"quarter_md={0.25 * period_md:.16e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
