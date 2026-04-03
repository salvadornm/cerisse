#!/usr/bin/env python3
"""
3D Numerical Schlieren three-view visualization from AMReX plotfiles.

Reads a 3D AMReX plotfile, extracts three orthogonal slices at a given origin
(or domain center by default), computes numerical Schlieren on each 2D slice,
and generates a three-view layout: Front (YZ), Side (XZ), Top (XY).

Uses the same direct AMReX reader approach as the 2D schlieren.py — no yt
dependency.  Supports multi-level AMR compositing.

Numerical Schlieren:  S = exp(-k * |∇ρ| / max|∇ρ|)

Usage
-----
  python schlieren_3d.py plot/plt00010
  python schlieren_3d.py plot/plt00010 --origin 90 90 72
  python schlieren_3d.py plot/plt00010 --plot Mach --origin 90 90 72
  python schlieren_3d.py plot/plt00010 --stl geometry.stl --mask-solid
  python schlieren_3d.py --k 30 --cmap bone

All CLI options
---------------
    plotfile (positional, optional): path to AMReX plotfile directory
    --plotdir PATH: directory containing plt* folders (default: plot)
    --origin X Y Z: slice origin in physical coordinates (default: domain center)
    --level INT: AMR output level (-1 = finest; 0..N = specific level)
    --k FLOAT: Schlieren contrast parameter (default: 20)
    --plot SPEC: select what to plot
        <varname>            -> raw field (e.g. pressure, temperature)
        Mach                 -> derived Mach number
        grad:<varname>       -> |grad(varname)|
        loggrad:<varname>    -> log10(|grad(varname)|)
        schlieren:<varname>  -> exp(-k*|grad|/max)  [default: schlieren:Density]
    --log: apply log10 to selected quantity (except schlieren mode)
    --raw: use pre-log quantity for loggrad mode
    --cmap NAME: matplotlib colormap (default: gray_r for schlieren)
    --vmin/--vmax FLOAT: color limits
    --gamma FLOAT: specific heat ratio for Mach (default: 1.4)
    --mask-solid: overlay solid regions (sld != 0)
    --mask-color NAME: solid overlay color (default: white)
    --mask-alpha FLOAT: solid overlay opacity 0-1 (default: 0.5)
    --stl PATH: overlay STL cross-section contours
    --stl-color NAME: STL contour color (default: lime)
    --stl-lw FLOAT: STL contour linewidth (default: 1.0)
    --figsize W H: figure size in inches (default: auto)
    -o, --output PATH: output filename
"""

import argparse
import glob
import os
import re
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection


# ─── Caches ──────────────────────────────────────────────────────────────────

_HEADER_CACHE = {}
_CELL_H_CACHE = {}


# ─── 3D AMReX plotfile reader ────────────────────────────────────────────────

def _parse_cell_h_3d(level_dir):
    """Parse Cell_H for 3D boxes and FabOnDisk entries."""
    cached = _CELL_H_CACHE.get(level_dir)
    if cached is not None:
        return cached

    with open(os.path.join(level_dir, "Cell_H")) as f:
        ch = [l.rstrip() for l in f.readlines()]

    fab_line_start = next(
        i for i, l in enumerate(ch) if l.strip().startswith("((")
    )

    boxes = []
    idx = fab_line_start
    while ch[idx].strip() != ")":
        line = ch[idx].strip()
        m = re.match(
            r"\(\((\d+),(\d+),(\d+)\)\s+\((\d+),(\d+),(\d+)\)", line
        )
        if m:
            boxes.append(tuple(int(m.group(g)) for g in range(1, 7)))
        idx += 1

    fabs_info = []
    for k in range(idx, len(ch)):
        if ch[k].strip().startswith("FabOnDisk:"):
            parts = ch[k].strip().split()
            fabs_info.append((parts[1], int(parts[2])))

    assert len(fabs_info) == len(boxes), (
        f"Expected {len(boxes)} fabs, got {len(fabs_info)}"
    )
    _CELL_H_CACHE[level_dir] = (boxes, fabs_info)
    return boxes, fabs_info


def _parse_header_3d(pltdir):
    """Parse AMReX Header for a 3D plotfile (cached)."""
    cached = _HEADER_CACHE.get(pltdir)
    if cached is not None:
        return cached

    with open(os.path.join(pltdir, "Header")) as f:
        hdr = [l.rstrip() for l in f.readlines()]

    nvars = int(hdr[1])
    varnames = [hdr[2 + i].strip() for i in range(nvars)]
    base = 2 + nvars
    ndim = int(hdr[base])
    if ndim != 3:
        raise RuntimeError(f"Expected 3D plotfile, got ndim={ndim}")

    sim_time = float(hdr[base + 1])
    finest_level = int(hdr[base + 2])
    lo = list(map(float, hdr[base + 3].split()))
    hi = list(map(float, hdr[base + 4].split()))

    if finest_level > 0:
        ref_ratios = list(map(int, hdr[base + 5].split()))
    else:
        ref_ratios = []

    # Locate domain-box line
    box_line = None
    for k in range(base + 4, min(base + 20, len(hdr))):
        if hdr[k].strip().startswith("(("):
            box_line = hdr[k].strip()
            break
    if box_line is None:
        raise RuntimeError(f"Could not locate domain box line in Header: {pltdir}")

    all_boxes = re.findall(
        r"\(\((\d+),(\d+),(\d+)\)\s+\((\d+),(\d+),(\d+)\)\s+"
        r"\((\d+),(\d+),(\d+)\)\)",
        box_line,
    )
    if not all_boxes:
        raise RuntimeError(f"Failed to parse 3D domain boxes from Header: {pltdir}")

    info = {
        "nvars": nvars,
        "varnames": varnames,
        "ndim": ndim,
        "sim_time": sim_time,
        "finest_level": finest_level,
        "lo": lo,
        "hi": hi,
        "ref_ratios": ref_ratios,
        "all_boxes": all_boxes,
    }
    _HEADER_CACHE[pltdir] = info
    return info


def _upsample_3d_nearest(arr, ratio):
    """Upsample a 3D array by nearest-neighbor repetition."""
    return np.repeat(
        np.repeat(np.repeat(arr, ratio, axis=0), ratio, axis=1),
        ratio, axis=2,
    )


def read_amrex_3d(pltdir, varname, out_level=-1):
    """Read a single variable from a 3D AMReX plotfile.

    Returns (field, x_coords, y_coords, z_coords, dx, dy, dz, meta).
    field has shape (Nx, Ny, Nz) with index field[i_x, j_y, k_z].
    """
    header = _parse_header_3d(pltdir)
    nvars = header["nvars"]
    varnames = header["varnames"]
    finest_level = header["finest_level"]
    lo = header["lo"]
    hi = header["hi"]
    ref_ratios = header["ref_ratios"]
    all_boxes = header["all_boxes"]

    target_level = finest_level if out_level < 0 else min(out_level, finest_level)
    fbox = all_boxes[target_level]
    Nx = int(fbox[3]) - int(fbox[0]) + 1
    Ny = int(fbox[4]) - int(fbox[1]) + 1
    Nz = int(fbox[5]) - int(fbox[2]) + 1

    if varname not in varnames:
        raise ValueError(
            f"Variable '{varname}' not found. Available: {varnames}"
        )
    var_idx = varnames.index(varname)

    field = np.zeros((Nx, Ny, Nz), dtype=np.float64)

    for lev in range(target_level + 1):
        ratio = 1
        for r in range(lev, target_level):
            ratio *= ref_ratios[r]

        level_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, fabs_info = _parse_cell_h_3d(level_dir)
        file_handles = {}

        def _get_fh(fname_local, _fh=file_handles, _ld=level_dir):
            fh = _fh.get(fname_local)
            if fh is None:
                fh = open(os.path.join(_ld, fname_local), "rb")
                _fh[fname_local] = fh
            return fh

        try:
            for b in range(len(boxes)):
                blo_i, blo_j, blo_k, bhi_i, bhi_j, bhi_k = boxes[b]
                ni = bhi_i - blo_i + 1
                nj = bhi_j - blo_j + 1
                nk = bhi_k - blo_k + 1

                fname, offset = fabs_info[b]
                bf = _get_fh(fname)
                bf.seek(offset)
                bf.readline()  # skip FAB header
                raw = np.fromfile(bf, dtype=np.float64, count=nvars * ni * nj * nk)
                if raw.size != nvars * ni * nj * nk:
                    raise RuntimeError(
                        f"Short FAB read: expected {nvars*ni*nj*nk}, got {raw.size}"
                    )
                # AMReX Fortran order: (comp, k, j, i)
                fab = raw.reshape((nvars, nk, nj, ni))
                patch = fab[var_idx].transpose(2, 1, 0)  # → (ni, nj, nk)

                if ratio == 1:
                    field[blo_i:bhi_i+1, blo_j:bhi_j+1, blo_k:bhi_k+1] = patch
                else:
                    fine = _upsample_3d_nearest(patch, ratio)
                    fi = blo_i * ratio
                    fj = blo_j * ratio
                    fk = blo_k * ratio
                    field[fi:fi+ni*ratio, fj:fj+nj*ratio, fk:fk+nk*ratio] = fine
        finally:
            for fh in file_handles.values():
                fh.close()

    dx = (hi[0] - lo[0]) / Nx
    dy = (hi[1] - lo[1]) / Ny
    dz = (hi[2] - lo[2]) / Nz
    x_coords = lo[0] + (np.arange(Nx) + 0.5) * dx
    y_coords = lo[1] + (np.arange(Ny) + 0.5) * dy
    z_coords = lo[2] + (np.arange(Nz) + 0.5) * dz

    meta = {
        "nvars": nvars,
        "varnames": varnames,
        "ndim": 3,
        "sim_time": header["sim_time"],
        "finest_level": finest_level,
        "output_level": target_level,
        "lo": lo,
        "hi": hi,
        "Nx": Nx, "Ny": Ny, "Nz": Nz,
    }
    return field, x_coords, y_coords, z_coords, dx, dy, dz, meta


def read_amrex_3d_multi(pltdir, requested_vars, out_level=-1):
    """Read multiple variables from a 3D AMReX plotfile in a single FAB sweep.

    Returns (fields_dict, x_coords, y_coords, z_coords, dx, dy, dz, meta).
    """
    header = _parse_header_3d(pltdir)
    nvars = header["nvars"]
    varnames = header["varnames"]
    finest_level = header["finest_level"]
    lo = header["lo"]
    hi = header["hi"]
    ref_ratios = header["ref_ratios"]
    all_boxes = header["all_boxes"]

    target_level = finest_level if out_level < 0 else min(out_level, finest_level)
    fbox = all_boxes[target_level]
    Nx = int(fbox[3]) - int(fbox[0]) + 1
    Ny = int(fbox[4]) - int(fbox[1]) + 1
    Nz = int(fbox[5]) - int(fbox[2]) + 1

    var_to_idx = {}
    for name in requested_vars:
        if name not in varnames:
            raise ValueError(f"Variable '{name}' not found. Available: {varnames}")
        var_to_idx[name] = varnames.index(name)

    fields = {name: np.zeros((Nx, Ny, Nz), dtype=np.float64) for name in requested_vars}

    for lev in range(target_level + 1):
        ratio = 1
        for r in range(lev, target_level):
            ratio *= ref_ratios[r]

        level_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, fabs_info = _parse_cell_h_3d(level_dir)
        file_handles = {}

        def _get_fh(fname_local, _fh=file_handles, _ld=level_dir):
            fh = _fh.get(fname_local)
            if fh is None:
                fh = open(os.path.join(_ld, fname_local), "rb")
                _fh[fname_local] = fh
            return fh

        try:
            for b in range(len(boxes)):
                blo_i, blo_j, blo_k, bhi_i, bhi_j, bhi_k = boxes[b]
                ni = bhi_i - blo_i + 1
                nj = bhi_j - blo_j + 1
                nk = bhi_k - blo_k + 1

                fname, offset = fabs_info[b]
                bf = _get_fh(fname)
                bf.seek(offset)
                bf.readline()
                raw = np.fromfile(bf, dtype=np.float64, count=nvars * ni * nj * nk)
                if raw.size != nvars * ni * nj * nk:
                    raise RuntimeError(
                        f"Short FAB read: expected {nvars*ni*nj*nk}, got {raw.size}"
                    )
                fab = raw.reshape((nvars, nk, nj, ni))

                if ratio == 1:
                    for name, vidx in var_to_idx.items():
                        fields[name][blo_i:bhi_i+1, blo_j:bhi_j+1, blo_k:bhi_k+1] = \
                            fab[vidx].transpose(2, 1, 0)
                else:
                    fi = blo_i * ratio
                    fj = blo_j * ratio
                    fk = blo_k * ratio
                    for name, vidx in var_to_idx.items():
                        coarse = fab[vidx].transpose(2, 1, 0)
                        fine = _upsample_3d_nearest(coarse, ratio)
                        fields[name][fi:fi+ni*ratio, fj:fj+nj*ratio, fk:fk+nk*ratio] = fine
        finally:
            for fh in file_handles.values():
                fh.close()

    dx = (hi[0] - lo[0]) / Nx
    dy = (hi[1] - lo[1]) / Ny
    dz = (hi[2] - lo[2]) / Nz
    x_coords = lo[0] + (np.arange(Nx) + 0.5) * dx
    y_coords = lo[1] + (np.arange(Ny) + 0.5) * dy
    z_coords = lo[2] + (np.arange(Nz) + 0.5) * dz

    meta = {
        "nvars": nvars,
        "varnames": varnames,
        "ndim": 3,
        "sim_time": header["sim_time"],
        "finest_level": finest_level,
        "output_level": target_level,
        "lo": lo,
        "hi": hi,
        "Nx": Nx, "Ny": Ny, "Nz": Nz,
    }
    return fields, x_coords, y_coords, z_coords, dx, dy, dz, meta


# ─── Gradient / Schlieren ────────────────────────────────────────────────────

def _diff_axis0(arr, spacing):
    """Central difference along axis 0 with one-sided boundaries."""
    out = np.empty_like(arr)
    inv2h = 0.5 / spacing
    invh = 1.0 / spacing
    out[1:-1, :] = (arr[2:, :] - arr[:-2, :]) * inv2h
    out[0, :] = (arr[1, :] - arr[0, :]) * invh
    out[-1, :] = (arr[-1, :] - arr[-2, :]) * invh
    return out


def _diff_axis1(arr, spacing):
    """Central difference along axis 1 with one-sided boundaries."""
    out = np.empty_like(arr)
    inv2h = 0.5 / spacing
    invh = 1.0 / spacing
    out[:, 1:-1] = (arr[:, 2:] - arr[:, :-2]) * inv2h
    out[:, 0] = (arr[:, 1] - arr[:, 0]) * invh
    out[:, -1] = (arr[:, -1] - arr[:, -2]) * invh
    return out


def grad_magnitude_2d(field, d0, d1):
    """Compute |∇field| on a 2D array. d0 = axis-0 spacing, d1 = axis-1."""
    g0 = _diff_axis0(field, d0)
    g1 = _diff_axis1(field, d1)
    return np.hypot(g0, g1)


def schlieren_2d(field, d0, d1, k=20.0):
    """Numerical Schlieren: exp(-k * |∇field| / max|∇field|).

    Non-finite cells (NaN/inf from IBM solid) are replaced with the field
    median before computing gradients, then restored as NaN in the output.
    """
    clean = field.copy()
    bad_mask = ~np.isfinite(clean)
    if np.any(bad_mask):
        med = np.nanmedian(clean[np.isfinite(clean)])
        clean[bad_mask] = med if np.isfinite(med) else 0.0
    grad = grad_magnitude_2d(clean, d0, d1)
    gmax = np.nanmax(grad[np.isfinite(grad)]) if np.any(np.isfinite(grad)) else 0.0
    if gmax > 0:
        out = np.exp(-k * grad / gmax)
    else:
        out = np.ones_like(grad)
    out[bad_mask] = np.nan
    return out


# ─── STL cross-section (vectorized) ─────────────────────────────────────────

def _read_stl_binary(path):
    """Read binary STL using numpy structured array (fast)."""
    with open(path, "rb") as f:
        f.read(80)  # header
        ntri = struct.unpack("<I", f.read(4))[0]
        dt = np.dtype([
            ("normal", "<f4", (3,)),
            ("vertices", "<f4", (3, 3)),
            ("attr", "<u2"),
        ])
        data = np.fromfile(f, dtype=dt, count=ntri)
    return data["vertices"].astype(np.float64), data["normal"].astype(np.float64)


def stl_cross_section(stl_path, axis, value):
    """Compute plane-triangle intersection segments (vectorized).

    Parameters
    ----------
    stl_path : str
        Path to binary STL file.
    axis : int
        0=x, 1=y, 2=z — the normal direction of the slice plane.
    value : float
        Coordinate value along *axis* where the plane sits.

    Returns
    -------
    segments : ndarray, shape (M, 2, 2)
        M line segments; each row is [[u0,v0],[u1,v1]] in the two remaining
        axes (sorted ascending, e.g. axis=0 → (y,z) columns).
    """
    verts, _ = _read_stl_binary(stl_path)   # (N, 3, 3)
    d = verts[:, :, axis] - value            # (N, 3)
    axes_2d = [a for a in range(3) if a != axis]
    ax0, ax1 = axes_2d

    # Per-edge crossing: edge01 (v0→v1), edge12 (v1→v2), edge20 (v2→v0)
    cross01 = d[:, 0] * d[:, 1] < 0
    cross12 = d[:, 1] * d[:, 2] < 0
    cross20 = d[:, 2] * d[:, 0] < 0

    # Triangles with exactly 2 edge crossings produce a segment
    ncross = cross01.astype(np.int8) + cross12.astype(np.int8) + cross20.astype(np.int8)
    valid = ncross == 2

    def _intersect_edge(mask, ea, eb):
        """Compute 2D intersection points for masked triangles on edge ea→eb."""
        idx = np.where(mask)[0]
        if idx.size == 0:
            return np.empty((0, 2))
        s0 = d[idx, ea]
        s1 = d[idx, eb]
        t = (s0 / (s0 - s1))[:, np.newaxis]
        p = verts[idx, ea, :] + t * (verts[idx, eb, :] - verts[idx, ea, :])
        return np.column_stack([p[:, ax0], p[:, ax1]])

    # Three cases: which pair of edges actually crosses
    case_01_12 = valid & cross01 & cross12
    case_12_20 = valid & cross12 & cross20
    case_01_20 = valid & cross01 & cross20

    parts = []
    for mask, ea1, eb1, ea2, eb2 in [
        (case_01_12, 0, 1, 1, 2),
        (case_12_20, 1, 2, 2, 0),
        (case_01_20, 0, 1, 2, 0),
    ]:
        p1 = _intersect_edge(mask, ea1, eb1)
        p2 = _intersect_edge(mask, ea2, eb2)
        if p1.shape[0] > 0:
            segs = np.stack([p1, p2], axis=1)  # (n, 2, 2)
            parts.append(segs)

    if parts:
        return np.concatenate(parts, axis=0)
    return np.empty((0, 2, 2))


# ─── Derived field loaders ───────────────────────────────────────────────────

def _load_mach_3d(pltdir, gamma=1.4, out_level=-1):
    """Load Mach number = |u| / a from 3D plotfile."""
    header = _parse_header_3d(pltdir)
    varnames = header["varnames"]
    vset = {v.lower() for v in varnames}
    vmap = {v.lower(): v for v in varnames}

    if "density" not in vset or "pressure" not in vset:
        raise ValueError("Mach requires Density and pressure in plotfile")

    eps = 1e-30

    if "x_velocity" in vset and "y_velocity" in vset and "z_velocity" in vset:
        need = [vmap["density"], vmap["pressure"],
                vmap["x_velocity"], vmap["y_velocity"], vmap["z_velocity"]]
        fields, xc, yc, zc, ddx, ddy, ddz, meta = read_amrex_3d_multi(
            pltdir, need, out_level=out_level)
        rho = fields[vmap["density"]]
        p = fields[vmap["pressure"]]
        u = fields[vmap["x_velocity"]]
        v = fields[vmap["y_velocity"]]
        w = fields[vmap["z_velocity"]]
    elif "xmom" in vset and "ymom" in vset and "zmom" in vset:
        need = [vmap["density"], vmap["pressure"],
                vmap["xmom"], vmap["ymom"], vmap["zmom"]]
        fields, xc, yc, zc, ddx, ddy, ddz, meta = read_amrex_3d_multi(
            pltdir, need, out_level=out_level)
        rho = fields[vmap["density"]]
        p = fields[vmap["pressure"]]
        rho_safe = np.maximum(rho, eps)
        u = fields[vmap["xmom"]] / rho_safe
        v = fields[vmap["ymom"]] / rho_safe
        w = fields[vmap["zmom"]] / rho_safe
    else:
        raise ValueError(
            "Mach requires (x_velocity,y_velocity,z_velocity) or "
            "(Xmom,Ymom,Zmom) in plotfile"
        )

    speed = np.sqrt(u * u + v * v + w * w)
    a = np.sqrt(gamma * np.maximum(p, 0.0) / np.maximum(rho, eps))
    mach = speed / np.maximum(a, eps)
    return mach, xc, yc, zc, ddx, ddy, ddz, meta


def _load_plot_field_3d(pltdir, varname, gamma=1.4, out_level=-1):
    """Load a raw or derived 3D field.

    Returns (field_3d, x, y, z, dx, dy, dz, meta, label).
    """
    vlow = varname.strip().lower()
    if vlow in ("mach", "mach_number"):
        f, xc, yc, zc, ddx, ddy, ddz, meta = _load_mach_3d(
            pltdir, gamma, out_level)
        return f, xc, yc, zc, ddx, ddy, ddz, meta, "Mach"

    f, xc, yc, zc, ddx, ddy, ddz, meta = read_amrex_3d(
        pltdir, varname, out_level)
    return f, xc, yc, zc, ddx, ddy, ddz, meta, varname


# ─── Plot-spec parsing ──────────────────────────────────────────────────────

def _parse_plot_spec(spec):
    """Parse --plot selector into (mode, varname)."""
    parts = spec.strip().split(":", 1)
    if len(parts) == 1:
        return "field", parts[0].strip()
    mode = parts[0].strip().lower()
    varname = parts[1].strip()
    if mode not in ("field", "grad", "loggrad", "schlieren"):
        raise ValueError(
            f"Invalid mode '{mode}'. Use field, grad, loggrad, or schlieren."
        )
    return mode, varname


def _safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


# ─── Plotfile auto-discovery ─────────────────────────────────────────────────

def _collect_plotfiles(plotdir):
    candidates = sorted(
        [d for d in glob.glob(os.path.join(plotdir, "plt*"))
         if not d.endswith(".temp") and ".old." not in d],
        key=lambda p: int(os.path.basename(p).replace("plt", "")),
    )
    return candidates


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="3D Numerical Schlieren three-view from AMReX plotfiles."
    )
    parser.add_argument("plotfile", nargs="?", help="Path to 3D plotfile")
    parser.add_argument("--plotdir", default="plot",
                        help="Directory containing plt* folders (default: plot)")
    parser.add_argument("--origin", nargs=3, type=float, default=None,
                        metavar=("X", "Y", "Z"),
                        help="Slice origin (default: domain center)")
    parser.add_argument("--level", type=int, default=-1,
                        help="AMR output level (-1 = finest)")
    parser.add_argument("--k", type=float, default=20.0,
                        help="Schlieren contrast parameter (default: 20)")
    parser.add_argument("--plot", default="schlieren:Density",
                        help="Plot spec (default: schlieren:Density)")
    parser.add_argument("--log", action="store_true",
                        help="Apply log10 (except schlieren mode)")
    parser.add_argument("--raw", action="store_true",
                        help="Use pre-log quantity for loggrad mode")
    parser.add_argument("--cmap", default=None,
                        help="Colormap (default: gray_r for schlieren)")
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--gamma", type=float, default=1.4,
                        help="Specific heat ratio for Mach (default: 1.4)")
    parser.add_argument("--mask-solid", action="store_true",
                        help="Overlay solid regions (sld != 0)")
    parser.add_argument("--mask-color", default="white",
                        help="Solid overlay color (default: white)")
    parser.add_argument("--mask-alpha", type=float, default=0.5,
                        help="Solid overlay opacity 0-1 (default: 0.5)")
    parser.add_argument("--stl", default=None,
                        help="STL file for cross-section overlay")
    parser.add_argument("--stl-color", default="lime",
                        help="STL contour color (default: lime)")
    parser.add_argument("--stl-lw", type=float, default=1.0,
                        help="STL contour linewidth (default: 1.0)")
    parser.add_argument("--figsize", nargs=2, type=float, default=None,
                        metavar=("W", "H"),
                        help="Figure size in inches (default: auto)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output filename")
    args = parser.parse_args()

    # ── Resolve plotfile ──────────────────────────────────────────────────
    plotfile = args.plotfile
    if not plotfile:
        sel = _collect_plotfiles(args.plotdir)
        if not sel:
            print(f"No plotfiles found in {args.plotdir}/")
            sys.exit(1)
        plotfile = sel[-1]
        print(f"Auto-selected: {plotfile}")

    # ── Parse plot spec ───────────────────────────────────────────────────
    mode_spec, varname = _parse_plot_spec(args.plot)

    # ── Load 3D field ─────────────────────────────────────────────────────
    print(f"Reading {plotfile} ...", flush=True)
    try:
        field, xc, yc, zc, dx, dy, dz, meta, var_label = _load_plot_field_3d(
            plotfile, varname, args.gamma, args.level,
        )
    except ValueError as exc:
        header = _parse_header_3d(plotfile)
        avail = list(header["varnames"]) + ["Mach (derived)"]
        print(f"Error: {exc}")
        print(f"Available variables: {', '.join(avail)}")
        sys.exit(1)

    print(f"  Grid: {meta['Nx']}x{meta['Ny']}x{meta['Nz']}, "
          f"dx={dx:.5g}, dy={dy:.5g}, dz={dz:.5g}")
    print(f"  Domain: [{meta['lo'][0]},{meta['lo'][1]},{meta['lo'][2]}] -> "
          f"[{meta['hi'][0]},{meta['hi'][1]},{meta['hi'][2]}]")
    n_nan = int(np.isnan(field).sum())
    n_inf = int(np.isinf(field).sum())
    finite = field[np.isfinite(field)]
    if finite.size > 0:
        extra = ""
        if n_nan > 0:
            extra += f"  ({n_nan} NaN"
        if n_inf > 0:
            extra += f", {n_inf} inf" if extra else f"  ({n_inf} inf"
        if extra:
            extra += " cells)"
        print(f"  {var_label} range: [{finite.min():.6g}, {finite.max():.6g}]{extra}")
    else:
        print(f"  {var_label}: all cells are NaN/inf")

    # ── Resolve origin ────────────────────────────────────────────────────
    if args.origin is not None:
        ox, oy, oz = args.origin
    else:
        ox = 0.5 * (meta["lo"][0] + meta["hi"][0])
        oy = 0.5 * (meta["lo"][1] + meta["hi"][1])
        oz = 0.5 * (meta["lo"][2] + meta["hi"][2])

    ix = int(np.argmin(np.abs(xc - ox)))
    iy = int(np.argmin(np.abs(yc - oy)))
    iz = int(np.argmin(np.abs(zc - oz)))
    print(f"  Slice origin: ({ox},{oy},{oz}) -> indices ({ix},{iy},{iz})")

    # ── Extract 2D slices from the raw field ──────────────────────────────
    # front: x=const → YZ plane,  shape (Nz, Ny) — rows=Z vert, cols=Y horiz
    # side:  y=const → XZ plane,  shape (Nz, Nx) — rows=Z vert, cols=X horiz
    # top:   z=const → XY plane,  shape (Ny, Nx) — rows=Y vert, cols=X horiz
    slice_front = field[ix, :, :].T   # (Nz, Ny)
    slice_side  = field[:, iy, :].T   # (Nz, Nx)
    slice_top   = field[:, :, iz].T   # (Ny, Nx)

    # ── Build plot data for each slice ────────────────────────────────────
    def _process_slice(raw_2d, d_vert, d_horiz):
        """Apply mode (schlieren/grad/loggrad/field) to a 2D slice.

        raw_2d shape is (n_vert, n_horiz).
        d_vert = cell spacing along axis-0 (vertical), d_horiz along axis-1.
        """
        # Replace non-finite (IBM solid: NaN/inf) with median for gradient modes
        bad_mask = ~np.isfinite(raw_2d)
        has_bad = np.any(bad_mask)

        if mode_spec == "schlieren":
            data = schlieren_2d(raw_2d, d_vert, d_horiz, args.k)
            default_cmap = "gray_r"
        elif mode_spec in ("grad", "loggrad"):
            clean = raw_2d.copy()
            if has_bad:
                finite_vals = clean[np.isfinite(clean)]
                med = np.median(finite_vals) if finite_vals.size > 0 else 0.0
                clean[bad_mask] = med
            data = grad_magnitude_2d(clean, d_vert, d_horiz)
            if has_bad:
                data[bad_mask] = np.nan
            if mode_spec == "loggrad" and not args.raw:
                data = np.log10(np.maximum(data, 1e-30))
            default_cmap = "inferno"
        else:  # "field"
            data = raw_2d.copy()
            default_cmap = "viridis"

        if mode_spec != "schlieren" and args.log:
            data = np.log10(np.maximum(data, 1e-30))

        return data, default_cmap

    # front: vert=Z(dz), horiz=Y(dy)
    # side:  vert=Z(dz), horiz=X(dx)
    # top:   vert=Y(dy), horiz=X(dx)
    data_front, dcmap = _process_slice(slice_front, dz, dy)
    data_side,  _     = _process_slice(slice_side,  dz, dx)
    data_top,   _     = _process_slice(slice_top,   dy, dx)

    cmap = args.cmap if args.cmap else dcmap

    # ── Solid mask ────────────────────────────────────────────────────────
    mask_front = mask_side = mask_top = None
    if args.mask_solid:
        try:
            sld = read_amrex_3d(plotfile, "sld", args.level)[0]
            mask_front = (sld[ix, :, :].T != 0)
            mask_side  = (sld[:, iy, :].T != 0)
            mask_top   = (sld[:, :, iz].T != 0)
            n_tot = int(mask_front.sum() + mask_side.sum() + mask_top.sum())
            print(f"  Solid mask: {n_tot} cells across 3 slices")
        except ValueError:
            print("  Warning: 'sld' variable not found, skipping solid mask")

    # ── STL cross-sections ────────────────────────────────────────────────
    stl_segs = [None, None, None]
    if args.stl:
        for vi, (axis, val) in enumerate([
            (0, xc[ix]), (1, yc[iy]), (2, zc[iz])
        ]):
            segs = stl_cross_section(args.stl, axis, val)
            stl_segs[vi] = segs
            print(f"  STL axis={axis} val={val:.2f}: {len(segs)} segments")

    # ── Unified color limits ──────────────────────────────────────────────
    if mode_spec == "schlieren":
        vmin = 0.0 if args.vmin is None else args.vmin
        vmax = 1.0 if args.vmax is None else args.vmax
    else:
        all_vals = np.concatenate([
            data_front.ravel(), data_side.ravel(), data_top.ravel()
        ])
        finite = all_vals[np.isfinite(all_vals)]
        vmin = float(np.min(finite)) if (args.vmin is None and finite.size > 0) else args.vmin
        vmax = float(np.max(finite)) if (args.vmax is None and finite.size > 0) else args.vmax

    # ── Plot ──────────────────────────────────────────────────────────────
    figsize = tuple(args.figsize) if args.figsize else (18, 5.5)
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    view_configs = [
        # (ax, data, horiz_coords, vert_coords, h_label, v_label, title, mask, stl)
        (axes[0], data_front, yc, zc, "Y", "Z",
         f"Front (x={xc[ix]:.1f})", mask_front, stl_segs[0]),
        (axes[1], data_side, xc, zc, "X", "Z",
         f"Side (y={yc[iy]:.1f})", mask_side, stl_segs[1]),
        (axes[2], data_top, xc, yc, "X", "Y",
         f"Top (z={zc[iz]:.1f})", mask_top, stl_segs[2]),
    ]

    last_pcm = None
    for ax, data, hc, vc, hl, vl, title, mask, segs in view_configs:
        pcm = ax.pcolormesh(hc, vc, data, shading="auto", cmap=cmap,
                            vmin=vmin, vmax=vmax)
        last_pcm = pcm

        # Solid mask overlay
        if mask is not None and np.any(mask):
            rgba = matplotlib.colors.to_rgba(args.mask_color, args.mask_alpha)
            overlay_cmap = ListedColormap([rgba])
            overlay_cmap.set_bad(alpha=0)
            overlay = np.where(mask, 1.0, np.nan)
            ax.pcolormesh(hc, vc, overlay, shading="auto",
                          cmap=overlay_cmap, vmin=0, vmax=1)

        # STL cross-section overlay
        n_stl = 0
        if segs is not None and len(segs) > 0:
            lc = LineCollection(segs, colors=args.stl_color,
                                linewidths=args.stl_lw, zorder=3)
            ax.add_collection(lc)
            n_stl = len(segs)

        title_str = f"{title}  ({n_stl} STL segs)" if n_stl > 0 else title
        ax.set_title(title_str, fontsize=10)
        ax.set_xlabel(hl)
        ax.set_ylabel(vl)
        ax.set_aspect("equal")

    # ── Suptitle ──────────────────────────────────────────────────────────
    if mode_spec == "schlieren":
        main_label = f"Schlieren({var_label})"
    elif mode_spec == "grad":
        main_label = rf"|$\nabla${var_label}|"
    elif mode_spec == "loggrad":
        main_label = rf"log$_{{10}}$|$\nabla${var_label}|"
    else:
        main_label = var_label
    if args.log and mode_spec != "schlieren":
        main_label = rf"log$_{{10}}$({main_label})"

    sim_time = meta.get("sim_time", None)
    basename = os.path.basename(os.path.normpath(plotfile))
    if sim_time is not None:
        fig.suptitle(f"{main_label} — {basename}  (t={sim_time:.6g})",
                     fontsize=12)
    else:
        fig.suptitle(f"{main_label} — {basename}", fontsize=12)

    # ── Shared colorbar ───────────────────────────────────────────────────
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
    fig.colorbar(last_pcm, cax=cbar_ax)

    plt.tight_layout(rect=[0, 0, 0.92, 0.95])

    outname = args.output if args.output else f"schlieren_3d_{_safe_name(var_label)}.png"
    fig.savefig(outname, dpi=200, bbox_inches="tight")
    print(f"Saved: {outname}")
    plt.close()


if __name__ == "__main__":
    main()
