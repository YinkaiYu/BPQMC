#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path


def load_monitor(path: Path, stage: int) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line or line.startswith("#"):
            continue
        step, row_stage, accepted, nsteps, delta_h, action, phi_f1, phi_f2, phi_rms_f1, phi_rms_f2 = line.split()
        if int(row_stage) != stage:
            continue
        rows.append(
            {
                "step": int(step),
                "accepted": int(accepted),
                "nsteps": int(nsteps),
                "delta_h": float(delta_h),
                "action": float(action),
                "phi_f1": float(phi_f1),
                "phi_f2": float(phi_f2),
                "phi_rms_f1": float(phi_rms_f1),
                "phi_rms_f2": float(phi_rms_f2),
            }
        )
    return rows


def linear_detrend(values: list[float]) -> tuple[float, float, list[float]]:
    npts = len(values)
    xs = list(range(npts))
    mean_x = sum(xs) / npts
    mean_y = sum(values) / npts
    denom = sum((x - mean_x) * (x - mean_x) for x in xs)
    slope = 0.0 if denom == 0.0 else sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, values)) / denom
    intercept = mean_y - slope * mean_x
    detrended = [y - (intercept + slope * x) for x, y in zip(xs, values)]
    return slope, intercept, detrended


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


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze HMC monitor output.")
    parser.add_argument("monitor", type=Path, help="Path to hmc_monitor.dat")
    parser.add_argument("--dt", type=float, required=True, help="Leapfrog step size used in the run.")
    parser.add_argument("--stage", type=int, default=1, help="Stage to analyze: 0=warm, 1=measurement.")
    parser.add_argument("--series", choices=("phi_f1", "phi_f2", "phi_rms_f1", "phi_rms_f2", "action"), default="phi_f2")
    parser.add_argument("--top-k", type=int, default=5, help="Number of FFT peaks to print.")
    args = parser.parse_args()

    rows = load_monitor(args.monitor, args.stage)
    if len(rows) < 8:
        raise SystemExit("Not enough monitor rows for analysis.")

    values = [float(row[args.series]) for row in rows]
    slope, intercept, detrended = linear_detrend(values)
    rms_detrended = math.sqrt(sum(value * value for value in detrended) / len(detrended))
    avg_nsteps = sum(int(row["nsteps"]) for row in rows) / len(rows)
    proposal_md_time = avg_nsteps * args.dt
    peaks = dominant_fft_periods(detrended, args.top_k)

    print(f"rows={len(rows)} stage={args.stage} series={args.series}")
    print(f"start={values[0]:.16e} end={values[-1]:.16e}")
    print(f"slope_per_proposal={slope:.16e}")
    print(f"rms_detrended={rms_detrended:.16e}")
    print(f"avg_nsteps={avg_nsteps:.6f} proposal_md_time={proposal_md_time:.16e}")
    print("top_fft_peaks:")
    for amplitude, kfreq, period_steps in peaks:
        print(
            f"  k={kfreq:3d} amp={amplitude:.16e} "
            f"period_proposals={period_steps:.6f} period_md={period_steps * proposal_md_time:.16e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
