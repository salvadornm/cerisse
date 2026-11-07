#!/usr/bin/env python3
"""
compare_vortex_probes.py

Compare 2-D convecting-vortex outlet tests from Cerisse probe CSV files.

Expected probe names from the recommended input section:
    PGlobalMAX, PGlobalMIN
    POutletMAX, POutletMIN
    PCentre, VCentre
    PUpper,  VUpper
    PLower,  VLower

The code tolerates Cerisse headers such as
    pressure_max((13,0)(115,127))
    velocity_average((101,63)(102,64))

and duplicate velocity columns produced for vector components.

Example
-------
python compare_vortex_probes.py \
    --case Outflow=run_outflow/probe_vortex.csv \
    --case LODI=run_lodi/probe_vortex.csv \
    --case Transverse=run_transverse/probe_vortex.csv \
    --output vortex_comparison

Outputs
-------
* pressure_envelope.png/pdf
* outlet_pressure_envelope.png/pdf
* centre_probe_pressure.png/pdf
* centre_probe_velocity.png/pdf
* upper_lower_transverse_velocity.png/pdf
* incoming_acoustic_proxy.png/pdf
* vortex_metrics_summary.csv
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
class Case:
    label: str
    path: Path
    t: np.ndarray
    cols: dict[str, list[np.ndarray]]
    metrics: dict[str, float]


def parse_case(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        label, filename = spec.split("=", 1)
        return label.strip(), Path(filename.strip())
    path = Path(spec)
    return path.stem, path


def canonical(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def split_cerisse_header(line: str) -> list[str]:
    """
    Split a Cerisse probe header while ignoring commas inside index tuples.

    Example:
      pressure_max((13,0)(114,127))
    must remain one column name.
    """
    names: list[str] = []
    token: list[str] = []
    depth = 0

    for char in line.strip():
        if char == "(":
            depth += 1
            token.append(char)
        elif char == ")":
            depth = max(depth - 1, 0)
            token.append(char)
        elif char == "," and depth == 0:
            names.append("".join(token).strip())
            token = []
        else:
            token.append(char)

    if token:
        names.append("".join(token).strip())

    return names


def make_unique(names: list[str]) -> list[str]:
    """Make duplicate probe names unique while preserving their order."""
    counts: dict[str, int] = {}
    unique: list[str] = []

    for name in names:
        n = counts.get(name, 0)
        unique.append(name if n == 0 else f"{name}.{n}")
        counts[name] = n + 1

    return unique


def read_probe(path: Path) -> pd.DataFrame:
    """
    Read a Cerisse probe CSV.

    Cerisse headers contain unquoted commas inside AMReX IntVect boxes, e.g.
      pressure_max((13,0)(114,127))
    so pandas.read_csv cannot parse the header directly. Read and split the
    header manually, then read numerical rows with header=None.
    """
    if not path.exists():
        raise FileNotFoundError(path)

    with path.open("r", encoding="utf-8") as handle:
        header_line = handle.readline()

    if not header_line:
        raise ValueError(f"{path} is empty")

    names = make_unique(split_cerisse_header(header_line))

    df = pd.read_csv(
        path,
        skiprows=1,
        header=None,
        names=names,
        comment="#",
        skipinitialspace=True,
    )

    if len(df.columns) != len(names):
        raise ValueError(
            f"{path}: parsed {len(df.columns)} data columns but "
            f"{len(names)} header names"
        )

    return df


def find_time_column(df: pd.DataFrame) -> str:
    for col in df.columns:
        if canonical(col) in {"time", "t", "simtime", "simulationtime"}:
            return col
    raise KeyError(f"No time column found. Columns: {list(df.columns)}")


def columns_matching(df: pd.DataFrame, *tokens: str) -> list[str]:
    """
    Return columns containing all canonical tokens.

    This works both with user probe names and generated reduction names.
    """
    keys = [canonical(t) for t in tokens]
    result = []
    for col in df.columns:
        c = canonical(col)
        if all(k in c for k in keys):
            result.append(col)
    return result


def first_existing_group(df: pd.DataFrame, alternatives: list[tuple[str, ...]]) -> list[str]:
    for tokens in alternatives:
        found = columns_matching(df, *tokens)
        if found:
            return found
    return []


def extract_groups(df: pd.DataFrame) -> dict[str, list[np.ndarray]]:
    """
    Decode the exact Cerisse header ordering used by the current vortex input.

    Expected order:
      pressure_max     : global, outlet
      pressure_min     : global, outlet
      pressure_average : centre, upper, lower
      velocity_average : centre-u, centre-v,
                         upper-u, upper-v,
                         lower-u, lower-v

    Pandas renames duplicate CSV headers by appending .1, .2, ...; the original
    column order is nevertheless preserved.
    """
    groups: dict[str, list[np.ndarray]] = {}

    def values(col: str) -> np.ndarray:
        return pd.to_numeric(df[col], errors="coerce").to_numpy(float)

    pmax_cols = columns_matching(df, "pressure", "max")
    pmin_cols = columns_matching(df, "pressure", "min")
    pavg_cols = columns_matching(df, "pressure", "average")
    vavg_cols = columns_matching(df, "velocity", "average")

    if len(pmax_cols) < 2:
        raise KeyError(
            f"Expected two pressure_max columns (global, outlet), found {pmax_cols}"
        )
    if len(pmin_cols) < 2:
        raise KeyError(
            f"Expected two pressure_min columns (global, outlet), found {pmin_cols}"
        )
    if len(pavg_cols) < 3:
        raise KeyError(
            f"Expected three pressure_average columns (centre, upper, lower), "
            f"found {pavg_cols}"
        )
    if len(vavg_cols) < 6:
        raise KeyError(
            f"Expected six velocity_average columns "
            f"(centre-u/v, upper-u/v, lower-u/v), found {vavg_cols}"
        )

    groups["PGlobalMAX"] = [values(pmax_cols[0])]
    groups["PGlobalMIN"] = [values(pmin_cols[0])]
    groups["POutletMAX"] = [values(pmax_cols[1])]
    groups["POutletMIN"] = [values(pmin_cols[1])]

    groups["PCentre"] = [values(pavg_cols[0])]
    groups["PUpper"]  = [values(pavg_cols[1])]
    groups["PLower"]  = [values(pavg_cols[2])]

    groups["VCentre"] = [values(vavg_cols[0]), values(vavg_cols[1])]
    groups["VUpper"]  = [values(vavg_cols[2]), values(vavg_cols[3])]
    groups["VLower"]  = [values(vavg_cols[4]), values(vavg_cols[5])]

    return groups


def finite_mask(t: np.ndarray, groups: dict[str, list[np.ndarray]]) -> np.ndarray:
    mask = np.isfinite(t)
    for arrays in groups.values():
        for a in arrays:
            mask &= np.isfinite(a)
    return mask


def safe_max_abs(a: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(a)
    if np.count_nonzero(valid) == 0:
        return float("nan")
    return float(np.max(np.abs(a[valid])))


def trapz(y: np.ndarray, x: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(valid) < 2:
        return float("nan")
    # np.trapz supports NumPy 1.x and 2.x.
    return float(np.trapz(y[valid], x[valid]))


def analyse(label: str, path: Path, args: argparse.Namespace) -> Case:
    df = read_probe(path)
    tcol = find_time_column(df)
    t = pd.to_numeric(df[tcol], errors="coerce").to_numpy(float)
    groups = extract_groups(df)

    if "PGlobalMAX" not in groups or "PGlobalMIN" not in groups:
        raise KeyError(
            f"{label}: could not identify global pressure max/min columns.\n"
            f"Available columns:\n{list(df.columns)}"
        )

    # Keep rows for which time and the two mandatory global-pressure
    # diagnostics are finite. Do not discard a complete row merely because an
    # optional velocity or local probe contains NaN.
    mask = (
        np.isfinite(t)
        & np.isfinite(groups["PGlobalMAX"][0])
        & np.isfinite(groups["PGlobalMIN"][0])
    )

    if not np.any(mask):
        nan_counts = {
            key: [int(np.count_nonzero(~np.isfinite(a))) for a in arrays]
            for key, arrays in groups.items()
        }
        raise ValueError(
            f"{label}: no usable rows remain after reading {path}. "
            f"NaN counts by probe: {nan_counts}. "
            "Check that the CSV contains numerical data rows and the same "
            "number of values as header columns."
        )

    t = t[mask]
    for key in list(groups):
        groups[key] = [a[mask] for a in groups[key]]

    order = np.argsort(t)
    t = t[order]
    for key in groups:
        groups[key] = [a[order] for a in groups[key]]

    c0 = math.sqrt(args.gamma * args.Rair * args.T0)
    u0 = args.Ma * c0
    Rv = args.Rv_fraction * args.Lx
    t_exit = (args.outlet_x - args.xc) / u0

    pressure_centre = args.p0 * math.exp(
        -0.5 * args.gamma * (args.Cv / (c0 * Rv)) ** 2
    )
    dp_vortex = args.p0 - pressure_centre

    pmax = groups["PGlobalMAX"][0]
    pmin = groups["PGlobalMIN"][0]
    p_envelope = np.maximum(np.abs(pmax - args.p0), np.abs(pmin - args.p0))

    late_start = args.late_start if args.late_start is not None else 1.4 * t_exit
    late_end = args.late_end if args.late_end is not None else float(t[-1])
    mlate = (t >= late_start) & (t <= late_end)

    metrics: dict[str, float] = {
        "t_exit_s": t_exit,
        "dp_vortex_Pa": dp_vortex,
        "late_pressure_peak_over_dpv": safe_max_abs(p_envelope, mlate) / dp_vortex,
        "late_pressure_L2_over_dpv": math.sqrt(
            trapz(p_envelope**2, t, mlate)
            / max(late_end - late_start, np.finfo(float).eps)
        )
        / dp_vortex,
    }

    if "POutletMAX" in groups and "POutletMIN" in groups:
        po_max = groups["POutletMAX"][0]
        po_min = groups["POutletMIN"][0]
        po_env = np.maximum(np.abs(po_max - args.p0), np.abs(po_min - args.p0))
        metrics["late_outlet_pressure_peak_over_dpv"] = (
            safe_max_abs(po_env, mlate) / dp_vortex
        )

    if "PCentre" in groups and "VCentre" in groups:
        p = groups["PCentre"][0]
        vel = groups["VCentre"]
        up = vel[0] - u0
        pp = p - args.p0

        # Linear incoming acoustic proxy at a probe upstream of the outlet.
        wminus = up - pp / (args.rho0 * c0)

        t_probe = (args.probe_x - args.xc) / u0
        inc_halfwidth = args.incident_halfwidth * Rv / u0
        minc = (t >= t_probe - inc_halfwidth) & (t <= t_probe + inc_halfwidth)
        mafter = t >= t_probe + inc_halfwidth

        if np.count_nonzero(minc) >= 2:
            wplus = up + pp / (args.rho0 * c0)
            inc_energy = trapz(wplus**2, t, minc)
            incoming_after = trapz(wminus**2, t, mafter)
            if inc_energy > 0.0:
                metrics["incoming_acoustic_proxy_amplitude"] = math.sqrt(
                    incoming_after / inc_energy
                )
                metrics["incoming_acoustic_proxy_energy"] = incoming_after / inc_energy

        metrics["late_centre_pressure_peak_over_dpv"] = (
            safe_max_abs(pp, mlate) / dp_vortex
        )
        metrics["late_centre_u_peak_over_u0"] = safe_max_abs(up, mlate) / u0

        if len(vel) >= 2:
            vp = vel[1]
            metrics["late_centre_v_peak_over_u0"] = safe_max_abs(vp, mlate) / u0

    return Case(label=label, path=path, t=t, cols=groups, metrics=metrics)


def save(fig: plt.Figure, out: Path, name: str, dpi: int) -> None:
    fig.tight_layout()
    fig.savefig(out / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_cases(cases: list[Case], args: argparse.Namespace, out: Path) -> None:
    c0 = math.sqrt(args.gamma * args.Rair * args.T0)
    u0 = args.Ma * c0
    Rv = args.Rv_fraction * args.Lx
    pcentre = args.p0 * math.exp(
        -0.5 * args.gamma * (args.Cv / (c0 * Rv)) ** 2
    )
    dpv = args.p0 - pcentre
    t_exit = (args.outlet_x - args.xc) / u0

    def xcoord(case: Case) -> np.ndarray:
        return case.t / t_exit if args.normalise_time else case.t

    xlabel = r"$t/t_{\rm exit}$" if args.normalise_time else "Time [s]"

    fig, ax = plt.subplots(figsize=(8.3, 4.8))
    for case in cases:
        pmax = case.cols["PGlobalMAX"][0]
        pmin = case.cols["PGlobalMIN"][0]
        env = np.maximum(np.abs(pmax - args.p0), np.abs(pmin - args.p0))
        ax.plot(xcoord(case), env / dpv, label=case.label)
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\|p-p_0\|_\infty/\Delta p_v$")
    ax.set_title("Interior pressure-disturbance envelope")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    save(fig, out, "01_pressure_envelope", args.dpi)

    available = [
        c for c in cases if "POutletMAX" in c.cols and "POutletMIN" in c.cols
    ]
    if available:
        fig, ax = plt.subplots(figsize=(8.3, 4.8))
        for case in available:
            pmax = case.cols["POutletMAX"][0]
            pmin = case.cols["POutletMIN"][0]
            env = np.maximum(np.abs(pmax - args.p0), np.abs(pmin - args.p0))
            ax.plot(xcoord(case), env / dpv, label=case.label)
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\|p-p_0\|_{\infty,\rm outlet}/\Delta p_v$")
        ax.set_title("Near-outlet pressure envelope")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        save(fig, out, "02_outlet_pressure_envelope", args.dpi)

    available = [c for c in cases if "PCentre" in c.cols]
    if available:
        fig, ax = plt.subplots(figsize=(8.3, 4.8))
        for case in available:
            ax.plot(
                xcoord(case),
                (case.cols["PCentre"][0] - args.p0) / dpv,
                label=case.label,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$(p-p_0)/\Delta p_v$")
        ax.set_title("Centreline pressure probe")
        ax.grid(True, alpha=0.3)
        ax.legend()
        save(fig, out, "03_centre_probe_pressure", args.dpi)

    available = [c for c in cases if "VCentre" in c.cols]
    if available:
        fig, ax = plt.subplots(figsize=(8.3, 4.8))
        for case in available:
            vel = case.cols["VCentre"]
            ax.plot(xcoord(case), (vel[0] - u0) / u0, label=f"{case.label}: u")
            if len(vel) >= 2:
                ax.plot(xcoord(case), vel[1] / u0, linestyle="--", label=f"{case.label}: v")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Velocity perturbation / $u_0$")
        ax.set_title("Centreline velocity probe")
        ax.grid(True, alpha=0.3)
        ax.legend()
        save(fig, out, "04_centre_probe_velocity", args.dpi)

    available = [
        c for c in cases if "VUpper" in c.cols and "VLower" in c.cols
    ]
    if available:
        fig, ax = plt.subplots(figsize=(8.3, 4.8))
        for case in available:
            vu = case.cols["VUpper"]
            vl = case.cols["VLower"]
            # Component 1 is transverse velocity when vector output contains u,v.
            idx = 1 if len(vu) >= 2 and len(vl) >= 2 else 0
            ax.plot(xcoord(case), vu[idx] / u0, label=f"{case.label}: upper")
            ax.plot(
                xcoord(case), vl[idx] / u0, linestyle="--", label=f"{case.label}: lower"
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$v/u_0$" if idx == 1 else "Recorded velocity / $u_0$")
        ax.set_title("Upper/lower vortex probe velocities")
        ax.grid(True, alpha=0.3)
        ax.legend()
        save(fig, out, "05_upper_lower_velocity", args.dpi)

    available = [
        c for c in cases if "PCentre" in c.cols and "VCentre" in c.cols
    ]
    if available:
        fig, ax = plt.subplots(figsize=(8.3, 4.8))
        for case in available:
            p = case.cols["PCentre"][0] - args.p0
            up = case.cols["VCentre"][0] - u0
            wm = up - p / (args.rho0 * c0)
            scale = max(args.Cv / Rv, np.finfo(float).eps)
            ax.plot(xcoord(case), np.abs(wm) / scale, label=case.label)
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$|w^-|/(C_v/R_v)$")
        ax.set_title("Incoming acoustic characteristic proxy")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        save(fig, out, "06_incoming_acoustic_proxy", args.dpi)


def write_summary(cases: list[Case], out: Path) -> None:
    rows = []
    for case in cases:
        row = {"method": case.label, "file": str(case.path)}
        row.update(case.metrics)
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "vortex_metrics_summary.csv", index=False)

    preferred = [
        "method",
        "late_pressure_peak_over_dpv",
        "late_outlet_pressure_peak_over_dpv",
        "incoming_acoustic_proxy_amplitude",
        "incoming_acoustic_proxy_energy",
        "late_centre_v_peak_over_u0",
    ]
    display = [c for c in preferred if c in summary.columns]
    shown = summary[display].copy()
    for col in shown.columns:
        if col != "method":
            shown[col] = shown[col].map(
                lambda x: f"{x:.6e}" if pd.notna(x) else "nan"
            )

    print("\n2-D vortex outlet comparison")
    print(shown.to_string(index=False))


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compare Cerisse 2-D vortex outlet probe histories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--case",
        action="append",
        required=True,
        metavar="LABEL=FILE",
        help="Repeat for each case.",
    )
    p.add_argument("--output", default="vortex_comparison")
    p.add_argument("--p0", type=float, default=101325.0)
    p.add_argument("--T0", type=float, default=300.0)
    p.add_argument("--Rair", type=float, default=287.0)
    p.add_argument("--gamma", type=float, default=1.4)
    p.add_argument("--rho0", type=float, default=101325.0 / (287.0 * 300.0))
    p.add_argument("--Ma", type=float, default=0.575)
    p.add_argument("--Lx", type=float, default=0.013)
    p.add_argument("--outlet-x", type=float, default=0.013)
    p.add_argument("--xc", type=float, default=0.0065)
    p.add_argument("--probe-x", type=float, default=0.0104)
    p.add_argument("--Rv-fraction", type=float, default=0.1)
    p.add_argument("--Cv", type=float, default=0.005)
    p.add_argument(
        "--incident-halfwidth",
        type=float,
        default=4.0,
        help="Incident window half-width in Rv/u0.",
    )
    p.add_argument(
        "--late-start",
        type=float,
        default=None,
        help="Start of late residual window [s]. Default 1.4*t_exit.",
    )
    p.add_argument("--late-end", type=float, default=None)
    p.add_argument("--normalise-time", action="store_true")
    p.add_argument("--dpi", type=int, default=220)
    return p


def main() -> int:
    args = parser().parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    try:
        cases = [analyse(*parse_case(spec), args) for spec in args.case]
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    write_summary(cases, out)
    plot_cases(cases, args, out)
    print(f"\nWrote outputs to {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
