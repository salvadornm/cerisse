#!/usr/bin/env python3
"""
Quantify whether domain boundaries are affecting a series of AMReX/Cerisse plotfiles.

The script reads a sequence of plt directories with yt, extracts thin strips near one or
more boundaries, compares them to an interior reference strip, and reports:

- backflow fraction
- net / incoming / outgoing mass-flux density
- relative pressure / density / normal-velocity jump between boundary and interior strips
- characteristic-wave amplitudes based on a local 1D normal decomposition

For an outflow boundary, a small incoming acoustic amplitude relative to the outgoing one
is a practical indicator that the boundary is not reflecting strongly back into the domain.

Examples
--------
  python tools/boundary_impact.py --plotdir wrk/axis2d/plot --boundaries xhi yhi
  python tools/boundary_impact.py --plotdir wrk/axis2d/plot --boundaries yhi --stride 500
  python tools/boundary_impact.py --plotdir wrk/axis2d/plot --boundaries xhi yhi \
      --strip-cells 8 --ref-gap-cells 24 --output wrk/axis2d/boundary_impact.csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from collections import defaultdict
from dataclasses import dataclass

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt


BOUNDARIES = ("xlo", "xhi", "ylo", "yhi")


@dataclass(frozen=True)
class BoundarySlices:
    boundary: tuple[slice, slice]
    reference: tuple[slice, slice]
    first_line: tuple[slice, slice]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnose boundary influence from a series of AMReX plotfiles."
    )
    parser.add_argument("--plotdir", default="plot", help="Directory containing plt* folders")
    parser.add_argument("--step-min", type=int, default=None, help="Minimum step to include")
    parser.add_argument("--step-max", type=int, default=None, help="Maximum step to include")
    parser.add_argument(
        "--stride",
        type=int,
        default=0,
        help="Step stride for subsampling; 0 keeps all selected plotfiles",
    )
    parser.add_argument(
        "--boundaries",
        nargs="+",
        default=["xhi", "yhi"],
        choices=BOUNDARIES,
        help="Domain faces to diagnose",
    )
    parser.add_argument(
        "--strip-cells",
        type=int,
        default=8,
        help="Boundary strip thickness in cells",
    )
    parser.add_argument(
        "--ref-gap-cells",
        type=int,
        default=16,
        help="Gap between boundary strip and interior reference strip",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=1.4,
        help="Specific heat ratio used for acoustic characteristics",
    )
    parser.add_argument(
        "--output",
        default="boundary_impact.csv",
        help="CSV file written with one row per plotfile and boundary",
    )
    parser.add_argument(
        "--penetration-depths",
        default="8,16,32,64",
        help="Comma-separated inward offsets (cells) for multi-depth penetration analysis",
    )
    parser.add_argument(
        "--penetration-output",
        default="",
        help="CSV for penetration rows (default: <output stem>_penetration.csv)",
    )
    parser.add_argument(
        "--focus-boundary",
        choices=BOUNDARIES,
        default=None,
        help="If set, write focused time-series figure for this boundary",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable auto-generated plots",
    )
    return parser.parse_args()


def parse_depths(text: str) -> list[int]:
    if not text.strip():
        return []
    depths: list[int] = []
    for token in text.split(","):
        value = int(token.strip())
        if value < 0:
            raise ValueError("--penetration-depths values must be >= 0")
        depths.append(value)
    return sorted(set(depths))


def collect_plotfiles(plotdir: str, step_min: int | None, step_max: int | None, stride: int) -> list[tuple[int, str]]:
    candidates = sorted(
        [path for path in glob.glob(os.path.join(plotdir, "plt*")) if os.path.isdir(path)],
        key=lambda path: int(os.path.basename(path).replace("plt", "")),
    )
    selected: list[tuple[int, str]] = []
    for path in candidates:
        step = int(os.path.basename(path).replace("plt", ""))
        if step_min is not None and step < step_min:
            continue
        if step_max is not None and step > step_max:
            continue
        if stride and step % stride != 0:
            continue
        selected.append((step, path))
    return selected


def finest_dims(ds) -> np.ndarray:
    dims = np.array(ds.domain_dimensions, dtype=int)
    max_level = int(getattr(ds.index, "max_level", 0))
    refine_by = int(getattr(ds, "refine_by", 2))
    if max_level > 0:
        dims = dims * (refine_by ** max_level)
    return dims


def load_fields(plotfile: str) -> dict[str, np.ndarray | float]:
    ds = yt.load(plotfile)
    dims = finest_dims(ds)
    level = int(getattr(ds.index, "max_level", 0))
    grid = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)

    fields = {
        "time": float(ds.current_time),
        "rho": np.asarray(grid[("boxlib", "Density")], dtype=np.float64)[..., 0],
        "p": np.asarray(grid[("boxlib", "pressure")], dtype=np.float64)[..., 0],
        "ux": np.asarray(grid[("boxlib", "x_velocity")], dtype=np.float64)[..., 0],
        "uy": np.asarray(grid[("boxlib", "y_velocity")], dtype=np.float64)[..., 0],
    }
    return fields


def build_slices(shape: tuple[int, int], boundary: str, strip_cells: int, ref_gap_cells: int) -> BoundarySlices:
    nx, ny = shape
    needed = 2 * strip_cells + ref_gap_cells
    if boundary in ("xlo", "xhi") and nx <= needed:
        raise ValueError(f"Grid too small in x for {boundary}: nx={nx}, need > {needed}")
    if boundary in ("ylo", "yhi") and ny <= needed:
        raise ValueError(f"Grid too small in y for {boundary}: ny={ny}, need > {needed}")

    if boundary == "xlo":
        return BoundarySlices(
            boundary=(slice(0, strip_cells), slice(None)),
            reference=(slice(strip_cells + ref_gap_cells, 2 * strip_cells + ref_gap_cells), slice(None)),
            first_line=(slice(0, 1), slice(None)),
        )
    if boundary == "xhi":
        return BoundarySlices(
            boundary=(slice(nx - strip_cells, nx), slice(None)),
            reference=(slice(nx - 2 * strip_cells - ref_gap_cells, nx - strip_cells - ref_gap_cells), slice(None)),
            first_line=(slice(nx - 1, nx), slice(None)),
        )
    if boundary == "ylo":
        return BoundarySlices(
            boundary=(slice(None), slice(0, strip_cells)),
            reference=(slice(None), slice(strip_cells + ref_gap_cells, 2 * strip_cells + ref_gap_cells)),
            first_line=(slice(None), slice(0, 1)),
        )
    return BoundarySlices(
        boundary=(slice(None), slice(ny - strip_cells, ny)),
        reference=(slice(None), slice(ny - 2 * strip_cells - ref_gap_cells, ny - strip_cells - ref_gap_cells)),
        first_line=(slice(None), slice(ny - 1, ny)),
    )


def boundary_variables(fields: dict[str, np.ndarray | float], boundary: str, strip_cells: int, ref_gap_cells: int):
    rho = fields["rho"]
    p = fields["p"]
    ux = fields["ux"]
    uy = fields["uy"]
    slices = build_slices(rho.shape, boundary, strip_cells, ref_gap_cells)

    if boundary == "xlo":
        un = -ux
        ut = uy
    elif boundary == "xhi":
        un = ux
        ut = uy
    elif boundary == "ylo":
        un = -uy
        ut = ux
    else:
        un = uy
        ut = ux

    return {
        "rho_b": rho[slices.boundary],
        "rho_r": rho[slices.reference],
        "p_b": p[slices.boundary],
        "p_r": p[slices.reference],
        "un_b": un[slices.boundary],
        "un_r": un[slices.reference],
        "ut_b": ut[slices.boundary],
        "ut_r": ut[slices.reference],
        "rho_face": rho[slices.first_line],
        "un_face": un[slices.first_line],
    }


def boundary_variables_offset(
    fields: dict[str, np.ndarray | float],
    boundary: str,
    strip_cells: int,
    ref_gap_cells: int,
    offset_cells: int,
):
    rho = fields["rho"]
    p = fields["p"]
    ux = fields["ux"]
    uy = fields["uy"]
    nx, ny = rho.shape

    if boundary == "xlo":
        slices = BoundarySlices(
            boundary=(slice(offset_cells, offset_cells + strip_cells), slice(None)),
            reference=(
                slice(offset_cells + strip_cells + ref_gap_cells, offset_cells + 2 * strip_cells + ref_gap_cells),
                slice(None),
            ),
            first_line=(slice(offset_cells, offset_cells + 1), slice(None)),
        )
        un, ut = -ux, uy
    elif boundary == "xhi":
        edge = nx - offset_cells
        slices = BoundarySlices(
            boundary=(slice(edge - strip_cells, edge), slice(None)),
            reference=(slice(edge - 2 * strip_cells - ref_gap_cells, edge - strip_cells - ref_gap_cells), slice(None)),
            first_line=(slice(edge - 1, edge), slice(None)),
        )
        un, ut = ux, uy
    elif boundary == "ylo":
        slices = BoundarySlices(
            boundary=(slice(None), slice(offset_cells, offset_cells + strip_cells)),
            reference=(
                slice(None),
                slice(offset_cells + strip_cells + ref_gap_cells, offset_cells + 2 * strip_cells + ref_gap_cells),
            ),
            first_line=(slice(None), slice(offset_cells, offset_cells + 1)),
        )
        un, ut = -uy, ux
    else:
        edge = ny - offset_cells
        slices = BoundarySlices(
            boundary=(slice(None), slice(edge - strip_cells, edge)),
            reference=(slice(None), slice(edge - 2 * strip_cells - ref_gap_cells, edge - strip_cells - ref_gap_cells)),
            first_line=(slice(None), slice(edge - 1, edge)),
        )
        un, ut = uy, ux

    vars_out = {
        "rho_b": rho[slices.boundary],
        "rho_r": rho[slices.reference],
        "p_b": p[slices.boundary],
        "p_r": p[slices.reference],
        "un_b": un[slices.boundary],
        "un_r": un[slices.reference],
        "ut_b": ut[slices.boundary],
        "ut_r": ut[slices.reference],
        "rho_face": rho[slices.first_line],
        "un_face": un[slices.first_line],
    }
    if min(arr.size for arr in vars_out.values()) == 0:
        raise ValueError(
            f"Offset {offset_cells} is too deep for {boundary} with strip={strip_cells}, gap={ref_gap_cells}"
        )
    return vars_out


def rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(arr * arr)))


def rel_jump(boundary_vals: np.ndarray, reference_vals: np.ndarray) -> float:
    ref_scale = rms(reference_vals)
    return rms(boundary_vals - reference_vals) / max(ref_scale, 1.0e-30)


def summarize_strip(arr: np.ndarray, boundary: str) -> np.ndarray:
    axis = 0 if boundary in ("xlo", "xhi") else 1
    return np.mean(arr, axis=axis)


def decompose_characteristics(
    rho_b: np.ndarray,
    rho_r: np.ndarray,
    p_b: np.ndarray,
    p_r: np.ndarray,
    un_b: np.ndarray,
    un_r: np.ndarray,
    ut_b: np.ndarray,
    ut_r: np.ndarray,
    gamma: float,
    boundary: str,
) -> dict[str, float]:
    rho_bar_b = summarize_strip(rho_b, boundary)
    rho_bar_r = summarize_strip(rho_r, boundary)
    p_bar_b = summarize_strip(p_b, boundary)
    p_bar_r = summarize_strip(p_r, boundary)
    un_bar_b = summarize_strip(un_b, boundary)
    un_bar_r = summarize_strip(un_r, boundary)
    ut_bar_b = summarize_strip(ut_b, boundary)
    ut_bar_r = summarize_strip(ut_r, boundary)

    rho0 = np.maximum(0.5 * (rho_bar_b + rho_bar_r), 1.0e-30)
    p0 = np.maximum(0.5 * (p_bar_b + p_bar_r), 1.0e-30)
    c0 = np.sqrt(gamma * p0 / rho0)
    un0 = 0.5 * (un_bar_b + un_bar_r)

    delta_p = p_bar_b - p_bar_r
    delta_un = un_bar_b - un_bar_r
    delta_ut = ut_bar_b - ut_bar_r
    delta_entropy = (rho_bar_b - rho_bar_r) - delta_p / np.maximum(c0 * c0, 1.0e-30)

    wave_plus = delta_un + delta_p / np.maximum(rho0 * c0, 1.0e-30)
    wave_minus = delta_un - delta_p / np.maximum(rho0 * c0, 1.0e-30)
    lam_plus = un0 + c0
    lam_minus = un0 - c0

    incoming_energy = []
    outgoing_energy = []
    incoming_energy.append(wave_plus[lam_plus < 0.0])
    incoming_energy.append(wave_minus[lam_minus < 0.0])
    outgoing_energy.append(wave_plus[lam_plus >= 0.0])
    outgoing_energy.append(wave_minus[lam_minus >= 0.0])
    incoming = np.concatenate([arr for arr in incoming_energy if arr.size]) if any(arr.size for arr in incoming_energy) else np.empty(0)
    outgoing = np.concatenate([arr for arr in outgoing_energy if arr.size]) if any(arr.size for arr in outgoing_energy) else np.empty(0)

    char_in_rms = rms(incoming)
    char_out_rms = rms(outgoing)
    return {
        "char_in_rms": char_in_rms,
        "char_out_rms": char_out_rms,
        "reflection_coeff": char_in_rms / max(char_out_rms, 1.0e-30),
        "char_entropy_rms": rms(delta_entropy),
        "char_tangential_rms": rms(delta_ut),
    }


def diagnose_boundary(fields: dict[str, np.ndarray | float], boundary: str, strip_cells: int, ref_gap_cells: int, gamma: float) -> dict[str, float | str]:
    vars_by_face = boundary_variables(fields, boundary, strip_cells, ref_gap_cells)
    rho_face = vars_by_face["rho_face"]
    un_face = vars_by_face["un_face"]
    flux_face = rho_face * un_face

    metrics: dict[str, float | str] = {
        "boundary": boundary,
        "backflow_frac": float(np.mean(un_face < 0.0)),
        "net_flux_density": float(np.mean(flux_face)),
        "outgoing_flux_density": float(np.mean(np.where(flux_face > 0.0, flux_face, 0.0))),
        "incoming_flux_density": float(np.mean(np.where(flux_face < 0.0, -flux_face, 0.0))),
    }

    rho_bar_b = summarize_strip(vars_by_face["rho_b"], boundary)
    rho_bar_r = summarize_strip(vars_by_face["rho_r"], boundary)
    p_bar_b = summarize_strip(vars_by_face["p_b"], boundary)
    p_bar_r = summarize_strip(vars_by_face["p_r"], boundary)
    un_bar_b = summarize_strip(vars_by_face["un_b"], boundary)
    un_bar_r = summarize_strip(vars_by_face["un_r"], boundary)

    metrics.update(
        {
            "rho_jump_rel": rel_jump(rho_bar_b, rho_bar_r),
            "p_jump_rel": rel_jump(p_bar_b, p_bar_r),
            "un_jump_rel": rel_jump(un_bar_b, un_bar_r),
            "p_boundary_rms_rel": rms(p_bar_b - np.mean(p_bar_b)) / max(abs(float(np.mean(p_bar_r))), 1.0e-30),
        }
    )
    metrics.update(
        decompose_characteristics(
            vars_by_face["rho_b"],
            vars_by_face["rho_r"],
            vars_by_face["p_b"],
            vars_by_face["p_r"],
            vars_by_face["un_b"],
            vars_by_face["un_r"],
            vars_by_face["ut_b"],
            vars_by_face["ut_r"],
            gamma,
            boundary,
        )
    )
    return metrics


def diagnose_penetration(
    fields: dict[str, np.ndarray | float],
    boundary: str,
    strip_cells: int,
    ref_gap_cells: int,
    gamma: float,
    offset_cells: int,
) -> dict[str, float | str]:
    vars_by_face = boundary_variables_offset(fields, boundary, strip_cells, ref_gap_cells, offset_cells)
    metrics: dict[str, float | str] = {
        "boundary": boundary,
        "depth_cells": float(offset_cells),
    }

    rho_bar_b = summarize_strip(vars_by_face["rho_b"], boundary)
    rho_bar_r = summarize_strip(vars_by_face["rho_r"], boundary)
    p_bar_b = summarize_strip(vars_by_face["p_b"], boundary)
    p_bar_r = summarize_strip(vars_by_face["p_r"], boundary)
    un_bar_b = summarize_strip(vars_by_face["un_b"], boundary)
    un_bar_r = summarize_strip(vars_by_face["un_r"], boundary)

    metrics.update(
        {
            "rho_jump_rel": rel_jump(rho_bar_b, rho_bar_r),
            "p_jump_rel": rel_jump(p_bar_b, p_bar_r),
            "un_jump_rel": rel_jump(un_bar_b, un_bar_r),
        }
    )
    metrics.update(
        decompose_characteristics(
            vars_by_face["rho_b"],
            vars_by_face["rho_r"],
            vars_by_face["p_b"],
            vars_by_face["p_r"],
            vars_by_face["un_b"],
            vars_by_face["un_r"],
            vars_by_face["ut_b"],
            vars_by_face["ut_r"],
            gamma,
            boundary,
        )
    )
    return metrics


def write_csv(rows: list[dict[str, float | str]], output: str) -> None:
    if not rows:
        return
    fieldnames = [
        "step",
        "time",
        "boundary",
        "backflow_frac",
        "net_flux_density",
        "outgoing_flux_density",
        "incoming_flux_density",
        "rho_jump_rel",
        "p_jump_rel",
        "un_jump_rel",
        "p_boundary_rms_rel",
        "char_in_rms",
        "char_out_rms",
        "reflection_coeff",
        "char_entropy_rms",
        "char_tangential_rms",
    ]
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict[str, float | str]]) -> None:
    by_boundary: dict[str, list[dict[str, float | str]]] = {}
    for row in rows:
        by_boundary.setdefault(str(row["boundary"]), []).append(row)

    print("Summary")
    for boundary, entries in by_boundary.items():
        peak_refl = max(entries, key=lambda row: float(row["reflection_coeff"]))
        peak_backflow = max(entries, key=lambda row: float(row["backflow_frac"]))
        peak_jump = max(entries, key=lambda row: float(row["p_jump_rel"]))
        print(
            f"  {boundary}: max reflection={float(peak_refl['reflection_coeff']):.3e} at step {int(peak_refl['step'])}, "
            f"max backflow={float(peak_backflow['backflow_frac']):.3f} at step {int(peak_backflow['step'])}, "
            f"max p_jump_rel={float(peak_jump['p_jump_rel']):.3e} at step {int(peak_jump['step'])}"
        )


def output_stem(path: str) -> str:
    root, _ = os.path.splitext(path)
    return root if root else path


def write_penetration_csv(rows: list[dict[str, float | str]], output: str) -> None:
    if not rows:
        return
    fieldnames = [
        "step",
        "time",
        "boundary",
        "depth_cells",
        "rho_jump_rel",
        "p_jump_rel",
        "un_jump_rel",
        "char_in_rms",
        "char_out_rms",
        "reflection_coeff",
        "char_entropy_rms",
        "char_tangential_rms",
    ]
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _to_array(entries: list[dict[str, float | str]], key: str) -> np.ndarray:
    return np.array([float(entry[key]) for entry in entries], dtype=np.float64)


def plot_timeseries(rows: list[dict[str, float | str]], png_path: str, focus_boundary: str | None) -> None:
    by_boundary: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        by_boundary[str(row["boundary"])].append(row)

    boundaries = [focus_boundary] if focus_boundary else sorted(by_boundary.keys())
    boundaries = [b for b in boundaries if b in by_boundary]
    if not boundaries:
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    ax_ref, ax_jump, ax_back, ax_flux = axes.ravel()

    for boundary in boundaries:
        entries = sorted(by_boundary[boundary], key=lambda item: int(item["step"]))
        steps = _to_array(entries, "step")
        ax_ref.plot(steps, _to_array(entries, "reflection_coeff"), label=boundary)
        ax_jump.plot(steps, _to_array(entries, "p_jump_rel"), label=boundary)
        ax_back.plot(steps, _to_array(entries, "backflow_frac"), label=boundary)
        ax_flux.plot(steps, _to_array(entries, "incoming_flux_density"), label=boundary)

    ax_ref.set_title("Reflection Coefficient")
    ax_ref.set_ylabel("char_in_rms / char_out_rms")
    ax_jump.set_title("Relative Pressure Jump")
    ax_jump.set_ylabel("p_jump_rel")
    ax_back.set_title("Backflow Fraction")
    ax_back.set_xlabel("Step")
    ax_back.set_ylabel("fraction")
    ax_flux.set_title("Incoming Flux Density")
    ax_flux.set_xlabel("Step")
    ax_flux.set_ylabel("rho*un (incoming)")

    for ax in (ax_ref, ax_jump, ax_back, ax_flux):
        ax.grid(True, alpha=0.3)
    if len(boundaries) > 1:
        ax_ref.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(png_path, dpi=170)
    plt.close(fig)


def plot_penetration(pen_rows: list[dict[str, float | str]], png_path: str, focus_boundary: str | None) -> None:
    by_boundary: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in pen_rows:
        by_boundary[str(row["boundary"])].append(row)

    boundaries = [focus_boundary] if focus_boundary else sorted(by_boundary.keys())
    boundaries = [b for b in boundaries if b in by_boundary]
    if not boundaries:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax_jump, ax_refl = axes

    for boundary in boundaries:
        entries = by_boundary[boundary]
        depth_groups: dict[int, list[dict[str, float | str]]] = defaultdict(list)
        for item in entries:
            depth_groups[int(float(item["depth_cells"]))].append(item)
        depths = sorted(depth_groups.keys())
        mean_jump = []
        p95_jump = []
        mean_refl = []
        p95_refl = []
        for depth in depths:
            group = depth_groups[depth]
            p_jump = np.array([float(g["p_jump_rel"]) for g in group], dtype=np.float64)
            refl = np.array([float(g["reflection_coeff"]) for g in group], dtype=np.float64)
            mean_jump.append(float(np.mean(p_jump)))
            p95_jump.append(float(np.percentile(p_jump, 95)))
            mean_refl.append(float(np.mean(refl)))
            p95_refl.append(float(np.percentile(refl, 95)))
        ax_jump.plot(depths, mean_jump, marker="o", label=f"{boundary} mean")
        ax_jump.plot(depths, p95_jump, marker="x", linestyle="--", label=f"{boundary} p95")
        ax_refl.plot(depths, mean_refl, marker="o", label=f"{boundary} mean")
        ax_refl.plot(depths, p95_refl, marker="x", linestyle="--", label=f"{boundary} p95")

    ax_jump.set_title("Penetration of p_jump_rel")
    ax_jump.set_xlabel("Depth from boundary [cells]")
    ax_jump.set_ylabel("p_jump_rel")
    ax_jump.grid(True, alpha=0.3)
    ax_jump.legend(loc="upper right", fontsize=8)

    ax_refl.set_title("Penetration of reflection coeff")
    ax_refl.set_xlabel("Depth from boundary [cells]")
    ax_refl.set_ylabel("char_in_rms / char_out_rms")
    ax_refl.grid(True, alpha=0.3)
    ax_refl.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(png_path, dpi=170)
    plt.close(fig)


def print_problem_windows(rows: list[dict[str, float | str]], boundary: str | None) -> None:
    by_boundary: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        by_boundary[str(row["boundary"])].append(row)

    targets = [boundary] if boundary else sorted(by_boundary.keys())
    for bnd in targets:
        if bnd not in by_boundary:
            continue
        entries = by_boundary[bnd]
        worst_ref = sorted(entries, key=lambda r: float(r["reflection_coeff"]), reverse=True)[:3]
        worst_jump = sorted(entries, key=lambda r: float(r["p_jump_rel"]), reverse=True)[:3]
        worst_back = sorted(entries, key=lambda r: float(r["backflow_frac"]), reverse=True)[:3]
        print(f"Problem windows for {bnd}")
        print("  top reflection steps:", ", ".join(str(int(r["step"])) for r in worst_ref))
        print("  top p_jump_rel steps:", ", ".join(str(int(r["step"])) for r in worst_jump))
        print("  top backflow steps:", ", ".join(str(int(r["step"])) for r in worst_back))


def main() -> None:
    args = parse_args()
    depths = parse_depths(args.penetration_depths)
    selected = collect_plotfiles(args.plotdir, args.step_min, args.step_max, args.stride)
    if not selected:
        raise SystemExit(f"No plotfiles found in {args.plotdir} after applying filters.")

    rows: list[dict[str, float | str]] = []
    penetration_rows: list[dict[str, float | str]] = []
    print(f"Processing {len(selected)} plotfiles from {args.plotdir}")
    for index, (step, plotfile) in enumerate(selected, start=1):
        print(f"  [{index:4d}/{len(selected):4d}] step={step}")
        fields = load_fields(plotfile)
        for boundary in args.boundaries:
            metrics = diagnose_boundary(fields, boundary, args.strip_cells, args.ref_gap_cells, args.gamma)
            metrics["step"] = step
            metrics["time"] = float(fields["time"])
            rows.append(metrics)
            for depth in depths:
                try:
                    pmetrics = diagnose_penetration(
                        fields,
                        boundary,
                        args.strip_cells,
                        args.ref_gap_cells,
                        args.gamma,
                        depth,
                    )
                except ValueError:
                    continue
                pmetrics["step"] = step
                pmetrics["time"] = float(fields["time"])
                penetration_rows.append(pmetrics)

    write_csv(rows, args.output)
    pen_output = args.penetration_output or f"{output_stem(args.output)}_penetration.csv"
    write_penetration_csv(penetration_rows, pen_output)
    if not args.no_plots:
        ts_png = f"{output_stem(args.output)}_timeseries.png"
        pen_png = f"{output_stem(args.output)}_penetration.png"
        plot_timeseries(rows, ts_png, args.focus_boundary)
        plot_penetration(penetration_rows, pen_png, args.focus_boundary)
        print(f"Wrote {ts_png}")
        print(f"Wrote {pen_png}")
    print_problem_windows(rows, args.focus_boundary)
    print_summary(rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {pen_output}")


if __name__ == "__main__":
    main()