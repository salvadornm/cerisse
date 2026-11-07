#!/usr/bin/env python3
"""Compare Cerisse probe.csv files for a 1-D acoustic-pulse NSCBC test.

Example
-------
python compare_nsbc_probe.py \
  --case LODI=run_lodi/probe.csv \
  --case NSCBC=run_nscbc/probe.csv \
  --case Extrapolation=run_extrap/probe.csv \
  --normalise-time --output nsbc_comparison

python compare_nsbc_probe.py --case SKEW2=runskew2/probe.csv --case SKEW2_LODI=runskew2_lodi/probe.csv  

Example output
    
    method  R_peak_characteristic  R_L2_characteristic  R_energy  R_peak_pressure  residual_pressure_over_dp0
    SKEW2               0.002519             0.002891  0.000008         0.257426                    1.016531
    SKEW2_LODI          0.008320             0.007361  0.000054         0.257426                   20.231927

"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class Result:
    label: str
    path: Path
    t: np.ndarray
    tx: np.ndarray
    pp: np.ndarray
    up: np.ndarray
    wp: np.ndarray
    wm: np.ndarray
    ea: np.ndarray
    ip: np.ndarray
    im: np.ndarray
    penv: np.ndarray | None
    inc_win: tuple[float, float]
    ref_win: tuple[float, float]
    metrics: dict[str, float]


def key(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    for kwargs in (
        dict(comment="#"),
        dict(sep=r"\s+", comment="#", engine="python"),
    ):
        try:
            df = pd.read_csv(path, **kwargs)
            if len(df.columns) > 1:
                return df
        except Exception:
            pass
    raise ValueError(f"Could not parse {path} as CSV or whitespace table")


def find_col(df: pd.DataFrame, names: list[str], required: bool = True) -> str | None:
    cmap = {key(c): c for c in df.columns}
    for name in names:
        if key(name) in cmap:
            return cmap[key(name)]
    for name in names:
        k = key(name)
        matches = [c for nk, c in cmap.items() if nk.startswith(k) or nk.endswith(k)]
        if len(matches) == 1:
            return matches[0]
    if required:
        raise KeyError(f"Could not find {names}; columns are {list(df.columns)}")
    return None


def parse_case(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        label, filename = spec.split("=", 1)
        return label.strip(), Path(filename.strip())
    p = Path(spec)
    return p.stem, p


def window_mask(t: np.ndarray, w: tuple[float, float]) -> np.ndarray:
    return (t >= w[0]) & (t <= w[1])


# def integrate(y: np.ndarray, t: np.ndarray, m: np.ndarray) -> float:
#     if np.count_nonzero(m) < 2:
#         return float("nan")
#     return float(np.trapezoid(y[m], t[m]))

def integrate(y, t, m):
    if np.count_nonzero(m) < 2:
        return np.nan

    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y[m], t[m]))   # NumPy >=2.0
    else:
        return float(np.trapz(y[m], t[m]))       # NumPy <2.0


def auto_windows(args, t: np.ndarray, c0: float):
    cp = args.u0 + c0
    cm = c0 - args.u0
    if cp <= 0 or cm <= 0:
        raise ValueError("Automatic windows require subsonic rightward mean flow")
    tinc = (args.probe_x - args.pulse_x0) / cp
    tout = (args.outlet_x - args.pulse_x0) / cp
    tref = tout + (args.outlet_x - args.probe_x) / cm
    sx = args.L0 / (math.sqrt(2.0) * args.pulse_B)
    wi = args.window_width * sx / cp
    wr = args.window_width * sx / cm
    inc = (max(t[0], tinc - wi), min(t[-1], tinc + wi))
    ref = (max(t[0], tref - wr), min(t[-1], tref + wr))
    return inc, ref

#time, pressure_max((7)(121)), pressure_min((7)(121)), pressure_average((102)(102)), velocity_average((102)(102))


def analyse(label: str, path: Path, args) -> Result:
    df = read_table(path)    
    ct = find_col(df, ["time"])
    cp = find_col(df, ["pressure_average"])
    #cp = find_col(df, ["Pprobe", "Pfluc", "pressureprobe", "pressure"])
    cu = find_col(df, ["velocity_average"])
    #cu = find_col(df, ["Uprobe", "Ufluc", "xvelocity", "velocity0", "velocity"])        
    cmax = find_col(df, ["pressure_max"], required=False)
    cmin = find_col(df, ["pressure_min"], required=False)
    #cmax = find_col(df, ["PressureMAX", "Pmax"], False)
    #cmin = find_col(df, ["PressureMIN", "Pmin"], False)
    

    cols = [ct, cp, cu] + ([cmax] if cmax else []) + ([cmin] if cmin else [])
    d = df[cols].copy()
    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=[ct, cp, cu]).sort_values(ct).drop_duplicates(ct)
    t = d[ct].to_numpy(float)
    if len(t) < 5:
        raise ValueError(f"{path}: too few valid samples")

    p = d[cp].to_numpy(float)
    u = d[cu].to_numpy(float)
    pp = p - args.p0
    up = u - args.u0

    c0 = math.sqrt(args.gamma * args.p0 / args.rho0)
    z0 = args.rho0 * c0
    wp = up + pp / z0
    wm = up - pp / z0

    ea = pp**2 / (2.0 * args.rho0 * c0**2) + 0.5 * args.rho0 * up**2
    ip = 0.25 * args.rho0 * c0 * wp**2
    im = 0.25 * args.rho0 * c0 * wm**2

    penv = None
    if cmax and cmin:
        pmax = d[cmax].to_numpy(float)
        pmin = d[cmin].to_numpy(float)
        penv = np.maximum(np.abs(pmax - args.p0), np.abs(pmin - args.p0))

    inc_auto, ref_auto = auto_windows(args, t, c0)
    inc = tuple(args.incident_window) if args.incident_window else inc_auto
    ref = tuple(args.reflected_window) if args.reflected_window else ref_auto
    mi = window_mask(t, inc)
    mr = window_mask(t, ref)
    if np.count_nonzero(mi) < 2 or np.count_nonzero(mr) < 2:
        raise ValueError(f"{label}: incident/reflected window has too few samples")

    wi_peak = float(np.max(np.abs(wp[mi])))
    wr_peak = float(np.max(np.abs(wm[mr])))
    wi2 = integrate(wp**2, t, mi)
    wr2 = integrate(wm**2, t, mr)
    ein = integrate(ip, t, mi)
    eref = integrate(im, t, mr)
    pin = float(np.max(np.abs(pp[mi])))
    pref = float(np.max(np.abs(pp[mr])))

    R_peak = wr_peak / wi_peak
    R_L2   = np.sqrt(wr2 / wi2)
    R_energy = eref / ein

    # -------------------------------------------------------
    # Total reflected energy after the incident pulse
    # -------------------------------------------------------

    incident_window = (0.00060, 0.00120)

    # Everything after the incident window
    mall = t >= incident_window[1]

    wr2_total = integrate(wm**2, t, mall)

    R_total_energy = wr2_total / wi2
    R_total_amplitude = np.sqrt(R_total_energy)

    post_exit = (args.outlet_x - args.pulse_x0) / (args.u0 + c0)
    post = t >= post_exit
    residual = float(np.max(penv[post])) if penv is not None and np.any(post) else float("nan")

    metrics = {
        "R_peak_characteristic": wr_peak / wi_peak if wi_peak > 0 else float("nan"),
        "R_L2_characteristic": R_L2,
        "R_energy": R_energy,
        "R_total_amplitude": R_total_amplitude,
        "R_total_energy": R_total_energy,
        "R_peak_pressure": pref / pin if pin > 0 else float("nan"),
        "incident_pressure_peak_Pa": pin,
        "reflected_pressure_peak_Pa": pref,
        "incident_energy_per_area_J_m2": ein,
        "reflected_energy_per_area_J_m2": eref,
        "residual_pressure_peak_Pa": residual,
        "residual_pressure_over_dp0": residual / args.dp_ref if np.isfinite(residual) else float("nan"),
        "incident_window_start_s": inc[0],
        "incident_window_end_s": inc[1],
        "reflected_window_start_s": ref[0],
        "reflected_window_end_s": ref[1],
    }
    t0 = args.L0 / (args.u0 + c0)
    tx = t / t0 if args.normalise_time else t
    return Result(label, path, t, tx, pp, up, wp, wm, ea, ip, im, penv, inc, ref, metrics)


def save(fig, out: Path, name: str, dpi: int):
    fig.tight_layout()
    fig.savefig(out / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def shade(ax, r: Result, args):
    c0 = math.sqrt(args.gamma * args.p0 / args.rho0)
    t0 = args.L0 / (args.u0 + c0)
    f = (lambda x: x / t0) if args.normalise_time else (lambda x: x)
    ax.axvspan(f(r.inc_win[0]), f(r.inc_win[1]), alpha=0.08)
    ax.axvspan(f(r.ref_win[0]), f(r.ref_win[1]), alpha=0.08)


def make_plots(results: list[Result], args, out: Path):
    xlabel = r"$t/t_0$" if args.normalise_time else "Time [s]"
    c0 = math.sqrt(args.gamma * args.p0 / args.rho0)
    wref = args.dp_ref / (args.rho0 * c0)
    eref = args.dp_ref**2 / (args.rho0 * c0**2)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    for r in results:
        ax.plot(r.tx, r.pp / args.dp_ref, label=r.label)
    shade(ax, results[0], args)
    ax.set(xlabel=xlabel, ylabel=r"$p'/\Delta p_0$", title="Pressure perturbation at probe")
    ax.grid(True, alpha=0.3); ax.legend()
    save(fig, out, "01_probe_pressure", args.dpi)

    for i, r in enumerate(results, 1):
        fig, ax = plt.subplots(figsize=(8, 4.8))
        ax.plot(r.tx, r.wp / wref, label=r"$w^+$ right-running")
        ax.plot(r.tx, r.wm / wref, label=r"$w^-$ left-running")
        shade(ax, r, args)
        ax.set(xlabel=xlabel, ylabel=r"$w^\pm/[\Delta p_0/(\rho_0c_0)]$", title=f"Characteristics: {r.label}")
        ax.grid(True, alpha=0.3); ax.legend()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", r.label)
        save(fig, out, f"02_characteristics_{i:02d}_{safe}", args.dpi)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    for r in results:
        ax.plot(r.tx, r.ea / eref, label=r.label)
    shade(ax, results[0], args)
    ax.set(xlabel=xlabel, ylabel=r"$e_a/[\Delta p_0^2/(\rho_0c_0^2)]$", title="Local acoustic energy density")
    ax.set_yscale("log"); ax.grid(True, alpha=0.3, which="both"); ax.legend()
    save(fig, out, "03_acoustic_energy_density", args.dpi)

    avail = [r for r in results if r.penv is not None]
    if avail:
        fig, ax = plt.subplots(figsize=(8, 4.8))
        for r in avail:
            ax.plot(r.tx, r.penv / args.dp_ref, label=r.label)
        shade(ax, avail[0], args)
        ax.set(xlabel=xlabel, ylabel=r"$\|p'\|_\infty/\Delta p_0$", title="Interior residual-pressure envelope")
        ax.set_yscale("log"); ax.grid(True, alpha=0.3, which="both"); ax.legend()
        save(fig, out, "04_residual_pressure_envelope", args.dpi)

    x = np.arange(len(results)); width = 0.25
    fig, ax = plt.subplots(figsize=(max(7, 1.25 * len(results) + 3), 4.8))
    ax.bar(x-width, [r.metrics["R_peak_characteristic"] for r in results], width, label=r"$R_{peak}$")
    ax.bar(x, [r.metrics["R_L2_characteristic"] for r in results], width, label=r"$R_2$")
    ax.bar(x+width, [r.metrics["R_energy"] for r in results], width, label=r"$R_E$")
    ax.set_xticks(x, [r.label for r in results], rotation=20, ha="right")
    ax.set(ylabel="Reflection coefficient", title="Boundary-condition reflection")
    ax.grid(True, alpha=0.3, axis="y"); ax.legend()
    save(fig, out, "05_reflection_coefficients", args.dpi)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    for r in results:
        ax.plot(r.tx, np.abs(r.wm) / wref, label=r.label)
    shade(ax, results[0], args)
    ax.set(xlabel=xlabel, ylabel=r"$|w^-|/[\Delta p_0/(\rho_0c_0)]$", title="Reflected characteristic")
    ax.set_yscale("log"); ax.grid(True, alpha=0.3, which="both"); ax.legend()
    save(fig, out, "06_reflected_characteristic", args.dpi)


def parser():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--case", action="append", required=True, metavar="LABEL=FILE")
    p.add_argument("--output", default="nsbc_comparison")
    p.add_argument("--p0", type=float, default=101325.0)
    p.add_argument("--rho0", type=float, default=1.18)
    p.add_argument("--u0", type=float, default=1.0)
    p.add_argument("--gamma", type=float, default=1.4)
    p.add_argument("--L0", type=float, default=1.0)
    p.add_argument("--dp-ref", type=float, default=None, help="Default is 1e-3*p0")
    p.add_argument("--probe-x", type=float, default=0.8)
    p.add_argument("--outlet-x", type=float, default=1.0)
    p.add_argument("--pulse-x0", type=float, default=0.5)
    p.add_argument("--pulse-B", type=float, default=10.0)
    p.add_argument("--window-width", type=float, default=4.0)
    p.add_argument("--incident-window", nargs=2, type=float, metavar=("START", "END"))
    p.add_argument("--reflected-window", nargs=2, type=float, metavar=("START", "END"))
    p.add_argument("--normalise-time", action="store_true")
    p.add_argument("--dpi", type=int, default=220)
    return p


def main() -> int:
    args = parser().parse_args()
    if args.dp_ref is None:
        args.dp_ref = 1e-3 * args.p0
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    try:
        results = [analyse(*parse_case(spec), args) for spec in args.case]
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    rows = []
    for r in results:
        row = {"method": r.label, "file": str(r.path)}; row.update(r.metrics); rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "nsbc_metrics_summary.csv", index=False)
    display_cols = [
    "method",
    "R_peak_characteristic",
    "R_L2_characteristic",
    "R_total_amplitude",
    "R_energy",
    "R_total_energy",
    "R_peak_pressure",
    "residual_pressure_over_dp0",
    ]
    #print(summary[["method", "R_peak_characteristic", "R_L2_characteristic", "R_energy", "R_peak_pressure", "residual_pressure_over_dp0"]].to_string(index=False))
    print(summary[display_cols].to_string(index=False))
    make_plots(results, args, out)
    print(f"\nWrote plots and summary to {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
