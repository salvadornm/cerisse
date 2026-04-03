#!/usr/bin/env python3
"""
Unified Numerical Schlieren from AMReX plotfiles (2D and 3D).

Auto-detects plotfile dimensionality from the Header.  Override with --2d/--3d.

2D mode
-------
  python schlieren.py plot/plt00100
  python schlieren.py plot/plt00100 --plot Mach --mirror r --rz r
  python schlieren.py --time-average --step-min 500 --step-max 1000

3D mode (three-view layout: Front YZ, Side XZ, Top XY)
-------
  python schlieren.py plot/plt00010
  python schlieren.py plot/plt00010 --origin 90 90 72
  python schlieren.py plot/plt00010 --stl geometry.stl --mask-solid

Plot specification (both modes):
  <varname>            raw field (Density, pressure, ...)
  schlieren:<varname>  exp(-k*|grad|/max)  [default: schlieren:Density]
  grad:<varname>       |grad(varname)|
  loggrad:<varname>    log10(|grad(varname)|)
  Mach                 derived Mach number
  Vorticity            derived 2D vorticity (2D only)
  Divergence           derived 2D divergence (2D only)
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
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Parse caches: reduce repeated Header/Cell_H parsing across multi-variable loads.
_HEADER_CACHE = {}
_CELL_H_CACHE = {}

# Separate 3D caches (different box tuple format prevents cross-contamination).
_HEADER_CACHE_3D = {}
_CELL_H_CACHE_3D = {}


def _detect_ndim(pltdir):
    """Read ndim from Header without a full parse."""
    with open(os.path.join(pltdir, "Header")) as f:
        hdr = [l.rstrip() for l in f.readlines()]
    nvars = int(hdr[1])
    return int(hdr[2 + nvars])


# ---------------------------------------------------------------------------
# Direct AMReX plotfile reader (no yt dependency — avoids yt RZ bug)
# Supports multi-level AMR: composites finest available data everywhere.
# ---------------------------------------------------------------------------
def _parse_cell_h(level_dir):
    """Parse Cell_H to get fab boxes and FabOnDisk entries."""
    cached = _CELL_H_CACHE.get(level_dir)
    if cached is not None:
        return cached

    with open(os.path.join(level_dir, "Cell_H")) as f:
        ch = [l.rstrip() for l in f.readlines()]

    fab_line_start = next(i for i, l in enumerate(ch) if l.strip().startswith("(("))

    boxes = []
    idx = fab_line_start
    while ch[idx].strip() != ")":
        line = ch[idx].strip()
        m = re.match(r"\(\((\d+),(\d+)\).*?\((\d+),(\d+)\)", line)
        if m:
            boxes.append(
                (int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)))
            )
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
    return _CELL_H_CACHE[level_dir]


def _get_header_info(pltdir):
    """Parse AMReX Header once per plotfile and cache metadata."""
    cached = _HEADER_CACHE.get(pltdir)
    if cached is not None:
        return cached

    with open(os.path.join(pltdir, "Header")) as f:
        hdr = [l.rstrip() for l in f.readlines()]

    nvars = int(hdr[1])
    varnames = [hdr[2 + i].strip() for i in range(nvars)]
    base = 2 + nvars
    ndim = int(hdr[base])
    sim_time = float(hdr[base + 1])
    finest_level = int(hdr[base + 2])
    lo = list(map(float, hdr[base + 3].split()))
    hi = list(map(float, hdr[base + 4].split()))

    if finest_level > 0:
        ref_ratios = list(map(int, hdr[base + 5].split()))
    else:
        ref_ratios = []

    box_line = None
    box_line_idx = None
    for k in range(base + 4, min(base + 16, len(hdr))):
        if hdr[k].strip().startswith("(("):
            box_line = hdr[k].strip()
            box_line_idx = k
            break

    if box_line is None:
        raise RuntimeError(f"Could not locate domain box line in Header: {pltdir}")

    all_boxes = re.findall(
        r"\(\((\d+),(\d+)\)\s+\((\d+),(\d+)\)\s+\((\d+),(\d+)\)\)",
        box_line,
    )
    coord_line_idx = box_line_idx + 2 + (finest_level + 1)
    coord_sys = int(hdr[coord_line_idx].strip())
    coord_name = _coord_name_from_id(coord_sys, ndim)

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
        "coord_sys": coord_sys,
        "coord_name": coord_name,
    }
    _HEADER_CACHE[pltdir] = info
    return info


def _coord_name_from_id(coord_sys, ndim):
    """Map AMReX coord_sys integer to readable geometry name."""
    if coord_sys == 0:
        return "xy" if ndim == 2 else "cartesian"
    if coord_sys == 1:
        return "rz"
    if coord_sys == 2:
        return "spherical"
    return f"unknown({coord_sys})"


def get_box_stats(pltdir):
    """Return per-level box statistics from an AMReX plotfile.

    Returns a list of dicts, one per level, each containing:
      level, nboxes, min_size, max_size, total_cells
    """
    header = _get_header_info(pltdir)
    finest_level = header["finest_level"]
    stats = []
    for lev in range(finest_level + 1):
        level_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, _ = _parse_cell_h(level_dir)
        sizes = []
        for blo_i, blo_j, bhi_i, bhi_j in boxes:
            ni = bhi_i - blo_i + 1
            nj = bhi_j - blo_j + 1
            sizes.append((ni, nj))
        nboxes = len(boxes)
        cells = [s[0] * s[1] for s in sizes]
        total = sum(cells)
        min_sz = min(sizes) if sizes else (0, 0)
        max_sz = max(sizes) if sizes else (0, 0)
        min_cells = min(cells) if cells else 0
        max_cells = max(cells) if cells else 0
        stats.append({
            "level": lev,
            "nboxes": nboxes,
            "min_size": min_sz,
            "max_size": max_sz,
            "min_cells": min_cells,
            "max_cells": max_cells,
            "total_cells": total,
        })
    return stats


def print_box_stats(pltdir):
    """Print box statistics for each AMR level."""
    stats = get_box_stats(pltdir)
    total_boxes = sum(s["nboxes"] for s in stats)
    total_cells = sum(s["total_cells"] for s in stats)
    print(f"  AMR boxes: {total_boxes} total, {total_cells:,} cells")
    for s in stats:
        lev = s["level"]
        nb = s["nboxes"]
        mn = s["min_size"]
        mx = s["max_size"]
        tc = s["total_cells"]
        print(
            f"    Level {lev}: {nb} boxes, "
            f"min {mn[0]}×{mn[1]}, max {mx[0]}×{mx[1]}, "
            f"{tc:,} cells"
        )


def _upsample_to_finest(arr, ratio, method):
    """Upsample a coarse patch to finest resolution for visualization."""
    if ratio == 1:
        return arr
    if method == "nearest":
        return np.repeat(np.repeat(arr, ratio, axis=0), ratio, axis=1)

    # Linear prolongation at fine-cell centers to reduce block artifacts.
    ni, nj = arr.shape
    x_coarse = np.arange(ni, dtype=np.float64)
    y_coarse = np.arange(nj, dtype=np.float64)
    x_fine = (np.arange(ni * ratio, dtype=np.float64) + 0.5) / ratio - 0.5
    y_fine = (np.arange(nj * ratio, dtype=np.float64) + 0.5) / ratio - 0.5
    x_fine = np.clip(x_fine, 0.0, ni - 1.0)
    y_fine = np.clip(y_fine, 0.0, nj - 1.0)

    tmp = np.empty((ni * ratio, nj), dtype=np.float64)
    for j in range(nj):
        tmp[:, j] = np.interp(x_fine, x_coarse, arr[:, j])

    out = np.empty((ni * ratio, nj * ratio), dtype=np.float64)
    for i in range(ni * ratio):
        out[i, :] = np.interp(y_fine, y_coarse, tmp[i, :])
    return out


def _dilate_mask_4n(mask, iters):
    """Dilate a boolean mask with 4-neighbor stencil."""
    out = mask
    for _ in range(max(0, iters)):
        p = np.pad(out, 1, mode="edge")
        out = (
            p[1:-1, 1:-1] | p[:-2, 1:-1] | p[2:, 1:-1] |
            p[1:-1, :-2] | p[1:-1, 2:]
        )
    return out


def _smooth_along_seams(field, source_id, strength=0.25, width=1, passes=1):
    """Apply local smoothing only near patch seams to reduce stitch lines."""
    seam = np.zeros_like(source_id, dtype=bool)
    seam[1:, :] |= source_id[1:, :] != source_id[:-1, :]
    seam[:-1, :] |= source_id[1:, :] != source_id[:-1, :]
    seam[:, 1:] |= source_id[:, 1:] != source_id[:, :-1]
    seam[:, :-1] |= source_id[:, 1:] != source_id[:, :-1]

    mask = _dilate_mask_4n(seam, width)
    if not np.any(mask):
        return field

    alpha = float(np.clip(strength, 0.0, 1.0))
    out = field.copy()
    for _ in range(max(1, passes)):
        p = np.pad(out, 1, mode="edge")
        avg = (
            p[1:-1, 1:-1] + p[:-2, 1:-1] + p[2:, 1:-1] +
            p[1:-1, :-2] + p[1:-1, 2:]
        ) / 5.0
        out[mask] = (1.0 - alpha) * out[mask] + alpha * avg[mask]
    return out


def _safe_name(text):
    """Convert arbitrary text to a filesystem-safe token."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def _parse_plot_spec(spec):
    """Parse --plot selector into (mode, varname)."""
    parts = spec.strip().split(":", 1)
    if len(parts) == 1:
        return "field", parts[0].strip()
    mode = parts[0].strip().lower()
    varname = parts[1].strip()
    if mode not in ("field", "grad", "loggrad", "schlieren"):
        raise ValueError(
            f"Invalid --plot mode '{mode}'. Use field, grad, loggrad, or schlieren."
        )
    if not varname:
        raise ValueError("--plot requires a variable name after ':'")
    return mode, varname


def _resolve_plot_axes(coord_name, horiz_choice_rz, horiz_choice_xy, coord0, coord1, data_ij):
    """Map internal data [coord0, coord1] to plotted (x,y) axes."""
    if coord_name == "rz":
        horiz = horiz_choice_rz
        if horiz == "z":
            x_coords = coord1
            y_coords = coord0
            data = data_ij
            x_axis_name, y_axis_name = "z", "r"
        else:
            x_coords = coord0
            y_coords = coord1
            data = data_ij.T
            x_axis_name, y_axis_name = "r", "z"
    else:
        horiz = horiz_choice_xy
        if horiz == "y":
            x_coords = coord1
            y_coords = coord0
            data = data_ij
            x_axis_name, y_axis_name = "y", "x"
        else:
            x_coords = coord0
            y_coords = coord1
            data = data_ij.T
            x_axis_name, y_axis_name = "x", "y"
    return x_coords, y_coords, data, x_axis_name, y_axis_name


def _get_level_box_overlays(pltdir, out_level=-1):
    """Return physical AMR patch boxes grouped by level."""
    header = _get_header_info(pltdir)
    lo = header["lo"]
    hi = header["hi"]
    finest_level = header["finest_level"]
    all_boxes = header["all_boxes"]

    if out_level is None or int(out_level) < 0:
        target_level = finest_level
    else:
        target_level = min(int(out_level), finest_level)

    level_boxes = []
    for lev in range(target_level + 1):
        level_box = all_boxes[lev]
        n0 = int(level_box[2]) - int(level_box[0]) + 1
        n1 = int(level_box[3]) - int(level_box[1]) + 1
        d0 = (hi[0] - lo[0]) / float(n0)
        d1 = (hi[1] - lo[1]) / float(n1)

        lev_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, _ = _parse_cell_h(lev_dir)
        lev_boxes = []
        for ilo, jlo, ihi, jhi in boxes:
            lev_boxes.append(
                {
                    "coord0_min": lo[0] + ilo * d0,
                    "coord0_max": lo[0] + (ihi + 1) * d0,
                    "coord1_min": lo[1] + jlo * d1,
                    "coord1_max": lo[1] + (jhi + 1) * d1,
                }
            )
        level_boxes.append(lev_boxes)
    return level_boxes


def _map_box_to_plot_axes(coord_name, horiz_choice_rz, horiz_choice_xy, box):
    """Map a physical box from internal coord axes to displayed plot axes."""
    if coord_name == "rz":
        if horiz_choice_rz == "z":
            return box["coord1_min"], box["coord1_max"], box["coord0_min"], box["coord0_max"]
        return box["coord0_min"], box["coord0_max"], box["coord1_min"], box["coord1_max"]

    if horiz_choice_xy == "y":
        return box["coord1_min"], box["coord1_max"], box["coord0_min"], box["coord0_max"]
    return box["coord0_min"], box["coord0_max"], box["coord1_min"], box["coord1_max"]


def _iter_mirrored_boxes(x0, x1, y0, y1, mirror_axes):
    """Yield original and mirrored display-space boxes."""
    xs = [(x0, x1)]
    ys = [(y0, y1)]
    if "x" in mirror_axes:
        xs.append((-x1, -x0))
    if "y" in mirror_axes:
        ys.append((-y1, -y0))

    for xx0, xx1 in xs:
        for yy0, yy1 in ys:
            yield xx0, xx1, yy0, yy1


def _draw_amr_boxes(ax, pltdir, meta, args, mirror_axes):
    """Overlay AMR patch boxes on the current axes."""
    if not args.boxes:
        return

    level_boxes = _get_level_box_overlays(pltdir, out_level=args.level)
    default_alpha = 0.18 if args.box_fill else 1.0
    alpha = args.box_alpha if args.box_alpha is not None else default_alpha
    level_cmap = plt.get_cmap("tab10")

    for lev, lev_boxes in enumerate(level_boxes):
        edge_color = args.box_color if args.box_color else level_cmap(lev % level_cmap.N)
        face_color = edge_color if args.box_fill else "none"
        label = f"Level {lev}"
        added_label = False

        for box in lev_boxes:
            x0, x1, y0, y1 = _map_box_to_plot_axes(meta["coord_name"], args.rz, args.xy, box)
            for xx0, xx1, yy0, yy1 in _iter_mirrored_boxes(x0, x1, y0, y1, mirror_axes):
                ax.add_patch(
                    Rectangle(
                        (xx0, yy0),
                        xx1 - xx0,
                        yy1 - yy0,
                        facecolor=face_color,
                        edgecolor=edge_color,
                        fill=args.box_fill,
                        linewidth=args.box_lw,
                        alpha=alpha,
                        label=label if not added_label else None,
                    )
                )
                added_label = True

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper right", framealpha=0.85)


def _load_solid_mask(pltdir, out_level=-1, prolong_method="nearest",
                     seam_smooth=0.0, seam_width=1, seam_passes=1):
    """Load a boolean solid mask from the 'sld' variable.

    Returns a 2D boolean array where True = solid (sld != 0), or None
    if the 'sld' variable is not found in the plotfile.
    """
    try:
        sld, _, _, _, _, _ = read_amrex_plotfile(
            pltdir, "sld", out_level=out_level,
            prolong_method=prolong_method, seam_smooth=seam_smooth,
            seam_width=seam_width, seam_passes=seam_passes,
        )
    except (ValueError, KeyError):
        return None
    return sld != 0


def _normalize_mirror_axes(tokens, x_axis_name, y_axis_name):
    """Normalize mirror tokens into plotted axes set {'x','y'}."""
    out = set()
    for tok in tokens:
        t = tok.lower()
        if t in ("x", "y"):
            out.add(t)
        elif t == x_axis_name:
            out.add("x")
        elif t == y_axis_name:
            out.add("y")
    return out


def read_amrex_plotfile(
    pltdir,
    varname,
    out_level=-1,
    prolong_method="linear",
    seam_smooth=0.0,
    seam_width=1,
    seam_passes=1,
):
    """Read a single variable from a (possibly multi-level) AMReX plotfile.

    Returns (field_2d, r_coords, z_coords, dr, dz).
    field_2d is composited at the finest-level resolution,
    with shape (Nr, Nz) and index [i_r, j_z].
    """
    # ---- Parse Header (cached) ----
    header = _get_header_info(pltdir)
    nvars = header["nvars"]
    varnames = header["varnames"]
    ndim = header["ndim"]
    sim_time = header["sim_time"]
    finest_level = header["finest_level"]
    lo = header["lo"]
    hi = header["hi"]
    ref_ratios = header["ref_ratios"]
    all_boxes = header["all_boxes"]
    if out_level is None or int(out_level) < 0:
        target_level = finest_level
    else:
        target_level = min(int(out_level), finest_level)

    fine_box = all_boxes[target_level]
    Nr = int(fine_box[2]) - int(fine_box[0]) + 1
    Nz = int(fine_box[3]) - int(fine_box[1]) + 1

    coord_sys = header["coord_sys"]
    coord_name = header["coord_name"]

    var_idx = varnames.index(varname)

    # ---- Composite field at finest-level resolution ----
    field = np.zeros((Nr, Nz), dtype=np.float64)
    source_id = np.full((Nr, Nz), -1, dtype=np.int32)
    patch_id = 0

    for lev in range(target_level + 1):
        # Refinement factor from this level to finest
        ratio = 1
        for r in range(lev, target_level):
            ratio *= ref_ratios[r]

        level_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, fabs_info = _parse_cell_h(level_dir)

        # Many fabs share the same file; keep descriptors open per level to
        # reduce open/close overhead during large frame sweeps.
        file_handles = {}

        def _get_fh(fname_local):
            fh = file_handles.get(fname_local)
            if fh is None:
                fh = open(os.path.join(level_dir, fname_local), "rb")
                file_handles[fname_local] = fh
            return fh

        try:
            for i in range(len(boxes)):
                blo_i, blo_j, bhi_i, bhi_j = boxes[i]
                ni = bhi_i - blo_i + 1
                nj = bhi_j - blo_j + 1

                fname, offset = fabs_info[i]
                bf = _get_fh(fname)
                bf.seek(offset)
                bf.readline()
                raw = np.fromfile(bf, dtype=np.float64, count=nvars * ni * nj)
                if raw.size != nvars * ni * nj:
                    raise RuntimeError(
                        f"Short FAB read in {pltdir}, Level_{lev}/{fname} at offset {offset}: "
                        f"expected {nvars * ni * nj}, got {raw.size}"
                    )
                fab = raw.reshape((nvars, nj, ni))  # (comp, j, i)

                if ratio == 1:
                    field[blo_i : bhi_i + 1, blo_j : bhi_j + 1] = fab[var_idx].T
                    source_id[blo_i : bhi_i + 1, blo_j : bhi_j + 1] = patch_id
                else:
                    coarse = fab[var_idx].T  # (ni, nj)
                    fine = _upsample_to_finest(coarse, ratio, prolong_method)
                    fi = blo_i * ratio
                    fj = blo_j * ratio
                    field[fi : fi + ni * ratio, fj : fj + nj * ratio] = fine
                    source_id[fi : fi + ni * ratio, fj : fj + nj * ratio] = patch_id
                patch_id += 1
        finally:
            for fh in file_handles.values():
                fh.close()

    if seam_smooth > 0.0:
        field = _smooth_along_seams(
            field,
            source_id,
            strength=seam_smooth,
            width=seam_width,
            passes=seam_passes,
        )

    dr = (hi[0] - lo[0]) / Nr
    dz = (hi[1] - lo[1]) / Nz
    r_coords = lo[0] + (np.arange(Nr) + 0.5) * dr
    z_coords = lo[1] + (np.arange(Nz) + 0.5) * dz

    meta = {
        "ndim": ndim,
        "coord_sys": coord_sys,
        "coord_name": coord_name,
        "finest_level": finest_level,
        "output_level": target_level,
        "sim_time": sim_time,
    }
    return field, r_coords, z_coords, dr, dz, meta


def read_amrex_plotfile_multi(
    pltdir,
    varnames,
    out_level=-1,
    prolong_method="linear",
    seam_smooth=0.0,
    seam_width=1,
    seam_passes=1,
):
    """Read multiple variables from one AMReX plotfile in a single FAB sweep.

    Returns (fields, r_coords, z_coords, dr, dz, meta), where fields is a
    dict keyed by requested variable names.
    """
    requested = list(dict.fromkeys(varnames))
    if not requested:
        raise ValueError("read_amrex_plotfile_multi requires at least one variable")

    header = _get_header_info(pltdir)
    nvars = header["nvars"]
    available = header["varnames"]
    ndim = header["ndim"]
    sim_time = header["sim_time"]
    finest_level = header["finest_level"]
    lo = header["lo"]
    hi = header["hi"]
    ref_ratios = header["ref_ratios"]
    all_boxes = header["all_boxes"]

    if out_level is None or int(out_level) < 0:
        target_level = finest_level
    else:
        target_level = min(int(out_level), finest_level)

    fine_box = all_boxes[target_level]
    nr = int(fine_box[2]) - int(fine_box[0]) + 1
    nz = int(fine_box[3]) - int(fine_box[1]) + 1

    coord_sys = header["coord_sys"]
    coord_name = header["coord_name"]

    var_to_idx = {}
    for name in requested:
        if name not in available:
            raise ValueError(f"Variable '{name}' not found")
        var_to_idx[name] = available.index(name)

    fields = {name: np.zeros((nr, nz), dtype=np.float64) for name in requested}
    source_id = np.full((nr, nz), -1, dtype=np.int32)
    patch_id = 0

    for lev in range(target_level + 1):
        ratio = 1
        for r in range(lev, target_level):
            ratio *= ref_ratios[r]

        level_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, fabs_info = _parse_cell_h(level_dir)
        file_handles = {}

        def _get_fh(fname_local):
            fh = file_handles.get(fname_local)
            if fh is None:
                fh = open(os.path.join(level_dir, fname_local), "rb")
                file_handles[fname_local] = fh
            return fh

        try:
            for i in range(len(boxes)):
                blo_i, blo_j, bhi_i, bhi_j = boxes[i]
                ni = bhi_i - blo_i + 1
                nj = bhi_j - blo_j + 1

                fname, offset = fabs_info[i]
                bf = _get_fh(fname)
                bf.seek(offset)
                bf.readline()
                raw = np.fromfile(bf, dtype=np.float64, count=nvars * ni * nj)
                if raw.size != nvars * ni * nj:
                    raise RuntimeError(
                        f"Short FAB read in {pltdir}, Level_{lev}/{fname} at offset {offset}: "
                        f"expected {nvars * ni * nj}, got {raw.size}"
                    )
                fab = raw.reshape((nvars, nj, ni))

                if ratio == 1:
                    for name, vidx in var_to_idx.items():
                        fields[name][blo_i : bhi_i + 1, blo_j : bhi_j + 1] = fab[vidx].T
                    source_id[blo_i : bhi_i + 1, blo_j : bhi_j + 1] = patch_id
                else:
                    fi = blo_i * ratio
                    fj = blo_j * ratio
                    for name, vidx in var_to_idx.items():
                        coarse = fab[vidx].T
                        fine = _upsample_to_finest(coarse, ratio, prolong_method)
                        fields[name][fi : fi + ni * ratio, fj : fj + nj * ratio] = fine
                    source_id[fi : fi + ni * ratio, fj : fj + nj * ratio] = patch_id
                patch_id += 1
        finally:
            for fh in file_handles.values():
                fh.close()

    if seam_smooth > 0.0:
        for name in requested:
            fields[name] = _smooth_along_seams(
                fields[name],
                source_id,
                strength=seam_smooth,
                width=seam_width,
                passes=seam_passes,
            )

    dr = (hi[0] - lo[0]) / nr
    dz = (hi[1] - lo[1]) / nz
    r_coords = lo[0] + (np.arange(nr) + 0.5) * dr
    z_coords = lo[1] + (np.arange(nz) + 0.5) * dz

    meta = {
        "ndim": ndim,
        "coord_sys": coord_sys,
        "coord_name": coord_name,
        "finest_level": finest_level,
        "output_level": target_level,
        "sim_time": sim_time,
    }
    return fields, r_coords, z_coords, dr, dz, meta


# ---------------------------------------------------------------------------
# 3D AMReX plotfile reader
# ---------------------------------------------------------------------------
def _parse_cell_h_3d(level_dir):
    """Parse Cell_H for 3D boxes and FabOnDisk entries."""
    cached = _CELL_H_CACHE_3D.get(level_dir)
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
    _CELL_H_CACHE_3D[level_dir] = (boxes, fabs_info)
    return boxes, fabs_info


def _parse_header_3d(pltdir):
    """Parse AMReX Header for a 3D plotfile (cached)."""
    cached = _HEADER_CACHE_3D.get(pltdir)
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
    _HEADER_CACHE_3D[pltdir] = info
    return info


def _get_level_box_overlays_3d(pltdir, out_level=-1):
    """Return physical 3D AMR patch boxes grouped by level.

    Each box is a dict with keys x_min, x_max, y_min, y_max, z_min, z_max.
    """
    header = _parse_header_3d(pltdir)
    lo = header["lo"]
    hi = header["hi"]
    finest_level = header["finest_level"]
    all_boxes = header["all_boxes"]
    target_level = finest_level if out_level < 0 else min(out_level, finest_level)

    level_boxes = []
    for lev in range(target_level + 1):
        fbox = all_boxes[lev]
        Nx = int(fbox[3]) - int(fbox[0]) + 1
        Ny = int(fbox[4]) - int(fbox[1]) + 1
        Nz = int(fbox[5]) - int(fbox[2]) + 1
        dx = (hi[0] - lo[0]) / Nx
        dy = (hi[1] - lo[1]) / Ny
        dz = (hi[2] - lo[2]) / Nz

        lev_dir = os.path.join(pltdir, f"Level_{lev}")
        boxes, _ = _parse_cell_h_3d(lev_dir)
        lev_boxes = []
        for blo_i, blo_j, blo_k, bhi_i, bhi_j, bhi_k in boxes:
            lev_boxes.append({
                "x_min": lo[0] + blo_i * dx,
                "x_max": lo[0] + (bhi_i + 1) * dx,
                "y_min": lo[1] + blo_j * dy,
                "y_max": lo[1] + (bhi_j + 1) * dy,
                "z_min": lo[2] + blo_k * dz,
                "z_max": lo[2] + (bhi_k + 1) * dz,
            })
        level_boxes.append(lev_boxes)
    return level_boxes


def _draw_amr_boxes_3d(ax, pltdir, args, slice_axis, slice_val,
                       horiz_key_min, horiz_key_max,
                       vert_key_min, vert_key_max,
                       normal_key_min, normal_key_max):
    """Overlay AMR patch box outlines on one 3D-slice view.

    Only boxes that contain the slice coordinate are drawn.
    """
    if not args.boxes:
        return

    level_boxes = _get_level_box_overlays_3d(pltdir, out_level=args.level)
    default_alpha = 0.18 if args.box_fill else 1.0
    alpha = args.box_alpha if args.box_alpha is not None else default_alpha
    level_cmap = plt.get_cmap("tab10")

    for lev, lev_boxes in enumerate(level_boxes):
        edge_color = args.box_color if args.box_color else level_cmap(lev % level_cmap.N)
        face_color = edge_color if args.box_fill else "none"
        label = f"Level {lev}"
        added_label = False

        for box in lev_boxes:
            # Skip boxes that don't intersect the slice plane
            if slice_val < box[normal_key_min] or slice_val > box[normal_key_max]:
                continue
            x0 = box[horiz_key_min]
            x1 = box[horiz_key_max]
            y0 = box[vert_key_min]
            y1 = box[vert_key_max]
            ax.add_patch(
                Rectangle(
                    (x0, y0), x1 - x0, y1 - y0,
                    facecolor=face_color,
                    edgecolor=edge_color,
                    fill=args.box_fill,
                    linewidth=args.box_lw,
                    alpha=alpha,
                    label=label if not added_label else None,
                )
            )
            added_label = True

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper right", framealpha=0.85,
                  fontsize=7)


def _upsample_3d(arr, ratio, method="linear"):
    """Upsample a 3D array by nearest-neighbor or trilinear interpolation."""
    if ratio == 1:
        return arr
    if method == "nearest":
        return np.repeat(
            np.repeat(np.repeat(arr, ratio, axis=0), ratio, axis=1),
            ratio, axis=2,
        )
    # Trilinear prolongation at fine-cell centres.
    ni, nj, nk = arr.shape
    def _fine_coords(n):
        c = np.arange(n, dtype=np.float64)
        f = np.clip((np.arange(n * ratio, dtype=np.float64) + 0.5) / ratio - 0.5,
                     0.0, n - 1.0)
        idx = np.clip(np.searchsorted(c, f, side="right") - 1, 0, n - 2)
        w = f - idx.astype(np.float64)
        return idx, w
    ix, wx = _fine_coords(ni)
    iy, wy = _fine_coords(nj)
    iz, wz = _fine_coords(nk)
    # axis-0
    tmp0 = arr[ix] * (1.0 - wx[:, None, None]) + arr[ix + 1] * wx[:, None, None]
    # axis-1
    tmp1 = tmp0[:, iy, :] * (1.0 - wy[None, :, None]) + tmp0[:, iy + 1, :] * wy[None, :, None]
    # axis-2
    out = tmp1[:, :, iz] * (1.0 - wz[None, None, :]) + tmp1[:, :, iz + 1] * wz[None, None, :]
    return out


def read_amrex_3d(pltdir, varname, out_level=-1, prolong_method="linear"):
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
                patch = fab[var_idx].transpose(2, 1, 0)  # -> (ni, nj, nk)

                if ratio == 1:
                    field[blo_i:bhi_i+1, blo_j:bhi_j+1, blo_k:bhi_k+1] = patch
                else:
                    fine = _upsample_3d(patch, ratio, prolong_method)
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


def read_amrex_3d_multi(pltdir, requested_vars, out_level=-1, prolong_method="linear"):
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
                        fine = _upsample_3d(coarse, ratio, prolong_method)
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


def _diff_axis0_centered(arr, spacing, out=None):
    """2nd-order central diff along axis-0 with one-sided boundaries."""
    if out is None:
        out = np.empty_like(arr)
    inv2h = 0.5 / spacing
    invh = 1.0 / spacing
    out[1:-1, :] = (arr[2:, :] - arr[:-2, :]) * inv2h
    out[0, :] = (arr[1, :] - arr[0, :]) * invh
    out[-1, :] = (arr[-1, :] - arr[-2, :]) * invh
    return out


def _diff_axis1_centered(arr, spacing, out=None):
    """2nd-order central diff along axis-1 with one-sided boundaries."""
    if out is None:
        out = np.empty_like(arr)
    inv2h = 0.5 / spacing
    invh = 1.0 / spacing
    out[:, 1:-1] = (arr[:, 2:] - arr[:, :-2]) * inv2h
    out[:, 0] = (arr[:, 1] - arr[:, 0]) * invh
    out[:, -1] = (arr[:, -1] - arr[:, -2]) * invh
    return out


def _grad_magnitude_2d(field, d0, d1):
    """Compute |grad(field)| with reusable work arrays."""
    g0 = _diff_axis0_centered(field, d0)
    g1 = _diff_axis1_centered(field, d1)
    np.hypot(g0, g1, out=g0)
    return g0


def compute_grad_rho(rho, dr, dz):
    """Compute |∇ρ| using 2nd-order finite differences."""
    return _grad_magnitude_2d(rho, dr, dz)


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
    grad = _grad_magnitude_2d(clean, d0, d1)
    gmax = np.nanmax(grad[np.isfinite(grad)]) if np.any(np.isfinite(grad)) else 0.0
    if gmax > 0:
        out = np.exp(-k * grad / gmax)
    else:
        out = np.ones_like(grad)
    out[bad_mask] = np.nan
    return out


# ---------------------------------------------------------------------------
# STL cross-section (3D only)
# ---------------------------------------------------------------------------
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
        axes (sorted ascending, e.g. axis=0 -> (y,z) columns).
    """
    verts, _ = _read_stl_binary(stl_path)   # (N, 3, 3)
    d = verts[:, :, axis] - value            # (N, 3)
    axes_2d = [a for a in range(3) if a != axis]
    ax0, ax1 = axes_2d

    cross01 = d[:, 0] * d[:, 1] < 0
    cross12 = d[:, 1] * d[:, 2] < 0
    cross20 = d[:, 2] * d[:, 0] < 0

    ncross = cross01.astype(np.int8) + cross12.astype(np.int8) + cross20.astype(np.int8)
    valid = ncross == 2

    def _intersect_edge(mask, ea, eb):
        idx = np.where(mask)[0]
        if idx.size == 0:
            return np.empty((0, 2))
        s0 = d[idx, ea]
        s1 = d[idx, eb]
        t = (s0 / (s0 - s1))[:, np.newaxis]
        p = verts[idx, ea, :] + t * (verts[idx, eb, :] - verts[idx, ea, :])
        return np.column_stack([p[:, ax0], p[:, ax1]])

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
            segs = np.stack([p1, p2], axis=1)
            parts.append(segs)

    if parts:
        return np.concatenate(parts, axis=0)
    return np.empty((0, 2, 2))


def _compute_vorticity_2d(u0, u1, d0, d1):
    """Compute 2D out-of-plane curl: omega = d(u1)/d0 - d(u0)/d1."""
    du1_d0 = _diff_axis0_centered(u1, d0)
    du0_d1 = _diff_axis1_centered(u0, d1)
    du1_d0 -= du0_d1
    return du1_d0


def _compute_divergence_2d(u0, u1, d0, d1):
    """Compute 2D divergence: div = d(u0)/d0 + d(u1)/d1."""
    du0_d0 = _diff_axis0_centered(u0, d0)
    du1_d1 = _diff_axis1_centered(u1, d1)
    du0_d0 += du1_d1
    return du0_d0


def _read_header_varnames(pltdir):
    """Read available variable names from AMReX plotfile Header."""
    return list(_get_header_info(pltdir)["varnames"])


def _compute_mach_field(
    pltdir,
    gamma,
    out_level=-1,
    prolong_method="linear",
    seam_smooth=0.0,
    seam_width=1,
    seam_passes=1,
):
    """Compute Mach = |u|/a from plotfile fields.

    Priority:
    1) x_velocity/y_velocity + pressure + Density
    2) Xmom/Ymom + pressure + Density (velocity from momentum)
    """
    varnames = _read_header_varnames(pltdir)
    vmap = {v.lower(): v for v in varnames}
    vset = set(vmap)

    if "density" not in vset or "pressure" not in vset:
        raise ValueError("Mach requires Density and pressure in plotfile")

    eps = 1e-30
    if "x_velocity" in vset and "y_velocity" in vset:
        need = [vmap["density"], vmap["pressure"], vmap["x_velocity"], vmap["y_velocity"]]
        fields, c0, c1, d0, d1, meta = read_amrex_plotfile_multi(
            pltdir,
            need,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        rho = fields[vmap["density"]]
        p = fields[vmap["pressure"]]
        u = fields[vmap["x_velocity"]]
        v = fields[vmap["y_velocity"]]
    elif "xmom" in vset and "ymom" in vset:
        need = [vmap["density"], vmap["pressure"], vmap["xmom"], vmap["ymom"]]
        fields, c0, c1, d0, d1, meta = read_amrex_plotfile_multi(
            pltdir,
            need,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        rho = fields[vmap["density"]]
        p = fields[vmap["pressure"]]
        mx = fields[vmap["xmom"]]
        my = fields[vmap["ymom"]]
        rho_safe = np.maximum(rho, eps)
        u = mx / rho_safe
        v = my / rho_safe
    else:
        raise ValueError(
            "Mach requires (x_velocity,y_velocity) or (Xmom,Ymom) in plotfile"
        )

    # Use in-place transforms to reduce peak memory on large grids.
    mach = np.hypot(u, v)
    a = np.maximum(rho, eps)
    np.divide(np.maximum(p, 0.0), a, out=a)
    a *= gamma
    np.maximum(a, eps, out=a)
    np.sqrt(a, out=a)
    np.divide(mach, a, out=mach)
    return mach, c0, c1, d0, d1, meta


def _load_plot_field(
    pltdir,
    varname,
    gamma,
    out_level=-1,
    prolong_method="linear",
    seam_smooth=0.0,
    seam_width=1,
    seam_passes=1,
):
    """Load raw or derived field requested by --plot."""
    v = varname.strip()
    vlow = v.lower()

    # Derived velocity components / magnitude to keep plotting usable even
    # when plotfiles do not store velocity as derive variables.
    if vlow in ("x_velocity", "xvel", "velx", "ux", "u0", "u"):
        u0, _, c0, c1, d0, d1, meta = _load_velocity_components(
            pltdir,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        return u0, c0, c1, d0, d1, meta, "x_velocity"

    if vlow in ("y_velocity", "yvel", "vely", "uy", "u1", "v"):
        _, u1, c0, c1, d0, d1, meta = _load_velocity_components(
            pltdir,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        return u1, c0, c1, d0, d1, meta, "y_velocity"

    if vlow in ("velocity", "speed", "velmag", "velocity_magnitude", "umag"):
        u0, u1, c0, c1, d0, d1, meta = _load_velocity_components(
            pltdir,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        umag = np.hypot(u0, u1)
        return umag, c0, c1, d0, d1, meta, "Velocity Magnitude"

    if vlow in ("mach", "mach_number"):
        field, c0, c1, d0, d1, meta = _compute_mach_field(
            pltdir,
            gamma,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        return field, c0, c1, d0, d1, meta, "Mach"

    if vlow in ("vorticity", "vort", "omega"):
        u0, u1, c0, c1, d0, d1, meta = _load_velocity_components(
            pltdir,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        vort = _compute_vorticity_2d(u0, u1, d0, d1)
        return vort, c0, c1, d0, d1, meta, "Vorticity"

    if vlow in ("divergence", "div", "dilatation"):
        u0, u1, c0, c1, d0, d1, meta = _load_velocity_components(
            pltdir,
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        div = _compute_divergence_2d(u0, u1, d0, d1)
        return div, c0, c1, d0, d1, meta, "Divergence"

    # Discrete marker fields must use nearest-neighbor prolongation
    # to avoid interpolation artifacts at AMR coarse-fine boundaries.
    _DISCRETE_FIELDS = {"sld", "ghs"}
    if vlow in _DISCRETE_FIELDS:
        prolong_method = "nearest"

    field, c0, c1, d0, d1, meta = read_amrex_plotfile(
        pltdir,
        v,
        out_level=out_level,
        prolong_method=prolong_method,
        seam_smooth=seam_smooth,
        seam_width=seam_width,
        seam_passes=seam_passes,
    )
    return field, c0, c1, d0, d1, meta, v


def _load_velocity_components(
    pltdir,
    out_level=-1,
    prolong_method="linear",
    seam_smooth=0.0,
    seam_width=1,
    seam_passes=1,
):
    """Load velocity components aligned with coord0/coord1.

    Returns u0, u1 where u0 is velocity along coord0 and u1 along coord1.
    """
    varnames = _read_header_varnames(pltdir)
    vset = {v.lower() for v in varnames}

    if "x_velocity" in varnames and "y_velocity" in varnames:
        fields, c0, c1, d0, d1, meta = read_amrex_plotfile_multi(
            pltdir,
            ["x_velocity", "y_velocity"],
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        u0 = fields["x_velocity"]
        u1 = fields["y_velocity"]
        return u0, u1, c0, c1, d0, d1, meta

    if "xmom" in vset and "ymom" in vset and "density" in vset:
        fields, c0, c1, d0, d1, meta = read_amrex_plotfile_multi(
            pltdir,
            ["Xmom", "Ymom", "Density"],
            out_level=out_level,
            prolong_method=prolong_method,
            seam_smooth=seam_smooth,
            seam_width=seam_width,
            seam_passes=seam_passes,
        )
        mx = fields["Xmom"]
        my = fields["Ymom"]
        rho = fields["Density"]
        rho_safe = np.maximum(rho, 1e-30)
        return mx / rho_safe, my / rho_safe, c0, c1, d0, d1, meta

    raise ValueError("Velocity streamlines require x_velocity/y_velocity or Xmom/Ymom+Density")


# ---------------------------------------------------------------------------
# 3D derived field loaders
# ---------------------------------------------------------------------------
def _load_mach_3d(pltdir, gamma=1.4, out_level=-1, prolong_method="linear"):
    """Load Mach number = |u| / a from a 3D plotfile."""
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
            pltdir, need, out_level=out_level, prolong_method=prolong_method)
        rho = fields[vmap["density"]]
        p = fields[vmap["pressure"]]
        u = fields[vmap["x_velocity"]]
        v = fields[vmap["y_velocity"]]
        w = fields[vmap["z_velocity"]]
    elif "xmom" in vset and "ymom" in vset and "zmom" in vset:
        need = [vmap["density"], vmap["pressure"],
                vmap["xmom"], vmap["ymom"], vmap["zmom"]]
        fields, xc, yc, zc, ddx, ddy, ddz, meta = read_amrex_3d_multi(
            pltdir, need, out_level=out_level, prolong_method=prolong_method)
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


def _load_plot_field_3d(pltdir, varname, gamma=1.4, out_level=-1,
                       prolong_method="linear"):
    """Load a raw or derived 3D field.

    Returns (field_3d, x, y, z, dx, dy, dz, meta, label).
    """
    vlow = varname.strip().lower()
    if vlow in ("mach", "mach_number"):
        f, xc, yc, zc, ddx, ddy, ddz, meta = _load_mach_3d(
            pltdir, gamma, out_level, prolong_method=prolong_method)
        return f, xc, yc, zc, ddx, ddy, ddz, meta, "Mach"

    f, xc, yc, zc, ddx, ddy, ddz, meta = read_amrex_3d(
        pltdir, varname, out_level, prolong_method=prolong_method)
    return f, xc, yc, zc, ddx, ddy, ddz, meta, varname


def _resolve_vector_axes(coord_name, horiz_choice_rz, horiz_choice_xy, u0, u1):
    """Map internal velocity [coord0,coord1] to plotted vector components."""
    if coord_name == "rz":
        if horiz_choice_rz == "z":
            # x=z uses component along coord1, y=r uses component along coord0
            return u1, u0
        # x=r, y=z
        return u0.T, u1.T

    if horiz_choice_xy == "y":
        # x=y uses component along coord1, y=x uses component along coord0
        return u1, u0
    # x=x, y=y
    return u0.T, u1.T


def _collect_plotfiles(plotdir, step_min=None, step_max=None, stride=1):
    """Collect plotfiles in plotdir with optional step-range filtering."""
    candidates = sorted(
        [d for d in glob.glob(os.path.join(plotdir, "plt*"))
         if not d.endswith(".temp") and ".old." not in d],
        key=lambda p: int(os.path.basename(p).replace("plt", "")),
    )
    if not candidates:
        return []

    steps = [int(os.path.basename(p).replace("plt", "")) for p in candidates]
    lo = step_min if step_min is not None else -np.inf
    hi = step_max if step_max is not None else np.inf
    sel = [(s, p) for s, p in zip(steps, candidates) if lo <= s <= hi]

    step_stride = max(1, int(stride))
    if step_stride > 1:
        sel = [(s, p) for s, p in sel if s % step_stride == 0]
    return sel


def _combine_bounds(b1, b2):
    """Intersect two optional [min,max] bounds."""
    def _norm(b):
        if b is None:
            return None
        lo, hi = b
        if lo is None and hi is None:
            return None
        return (lo, hi)

    b1 = _norm(b1)
    b2 = _norm(b2)
    if b1 is None:
        return b2
    if b2 is None:
        return b1

    lo_candidates = [v for v in (b1[0], b2[0]) if v is not None]
    hi_candidates = [v for v in (b1[1], b2[1]) if v is not None]
    lo = max(lo_candidates) if lo_candidates else None
    hi = min(hi_candidates) if hi_candidates else None

    if lo is not None and hi is not None and lo > hi:
        return (1.0, 0.0)  # empty interval marker
    return (lo, hi)


def _crop_1d_bounds(coords, lo, hi):
    """Return slice bounds [i0, i1) overlapping [lo,hi] for monotonic coords."""
    if coords.size == 0:
        return 0, 0

    cmin = float(np.min(coords))
    cmax = float(np.max(coords))

    if lo is None:
        lo_eff = cmin
    else:
        lo_eff = max(float(lo), cmin)
    if hi is None:
        hi_eff = cmax
    else:
        hi_eff = min(float(hi), cmax)

    if lo_eff > hi_eff:
        return 0, 0

    i0 = int(np.searchsorted(coords, lo_eff, side="left"))
    i1 = int(np.searchsorted(coords, hi_eff, side="right"))
    i0 = max(0, min(i0, coords.size))
    i1 = max(0, min(i1, coords.size))
    return i0, i1


def _build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Unified Schlieren from AMReX plotfiles (2D and 3D)."
    )
    parser.add_argument("plotfile", nargs="?", help="Path to the plotfile")

    # ── Dimension mode ────────────────────────────────────────────────────
    dim_group = parser.add_mutually_exclusive_group()
    dim_group.add_argument("--2d", dest="force_2d", action="store_true",
                           help="Force 2D mode (auto-detected by default)")
    dim_group.add_argument("--3d", dest="force_3d", action="store_true",
                           help="Force 3D three-view mode (auto-detected by default)")

    # ── 3D-specific ───────────────────────────────────────────────────────
    parser.add_argument("--origin", nargs=3, type=float, default=None,
                        metavar=("X", "Y", "Z"),
                        help="[3D] Slice origin (default: domain center)")
    parser.add_argument("--stl", default=None,
                        help="[3D] STL file for cross-section overlay")
    parser.add_argument("--stl-color", default="lime",
                        help="[3D] STL contour color (default: lime)")
    parser.add_argument("--stl-lw", type=float, default=1.0,
                        help="[3D] STL contour linewidth (default: 1.0)")
    parser.add_argument("--figsize", nargs=2, type=float, default=None,
                        metavar=("W", "H"),
                        help="[3D] Figure size in inches (default: auto)")

    # ── 2D time-window modes ──────────────────────────────────────────────
    parser.add_argument("--time-average", action="store_true",
                        help="Average selected plotfiles in time before plotting")
    parser.add_argument("--time-std", action="store_true",
                        help="Compute temporal standard deviation field over selected steps")
    parser.add_argument("--time-fluctuation", action="store_true",
                        help="Compute fluctuation field x(step)-<x> over selected steps")
    parser.add_argument("--fluct-step", type=int, default=None,
                        help="Step used for fluctuation; default is last selected step")
    parser.add_argument("--plotdir", default="plot",
                        help="Directory containing plt* folders (default: plot)")
    parser.add_argument("--step-min", type=int, default=None,
                        help="Minimum step to include for time-window stats")
    parser.add_argument("--step-max", type=int, default=None,
                        help="Maximum step to include for time-window stats")
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep every Nth step for time-window stats (default: 1)")
    parser.add_argument("--level", type=int, default=-1,
                        help="AMR output level (-1 = finest; 0..N = specific level)")
    parser.add_argument("--x-min", type=float, default=None,
                        help="Lower bound of Cartesian x for cropping")
    parser.add_argument("--x-max", type=float, default=None,
                        help="Upper bound of Cartesian x for cropping")
    parser.add_argument("--y-min", type=float, default=None,
                        help="Lower bound of Cartesian y for cropping")
    parser.add_argument("--y-max", type=float, default=None,
                        help="Upper bound of Cartesian y for cropping")
    parser.add_argument("--r-min", type=float, default=None,
                        help="Lower bound of radial r for cropping")
    parser.add_argument("--r-max", type=float, default=None,
                        help="Upper bound of radial r for cropping")
    parser.add_argument("--z-min", type=float, default=None,
                        help="Lower bound of axial z for cropping")
    parser.add_argument("--z-max", type=float, default=None,
                        help="Upper bound of axial z for cropping")
    parser.add_argument("--k", type=float, default=20.0,
                        help="Contrast parameter for exp(-k * |∇ρ|/max) (default: 20)")
    parser.add_argument("--plot", default="schlieren:Density",
                        help="<var>, grad:<var>, loggrad:<var>, schlieren:<var>, or derived fields (mach/vorticity/divergence/velocity)")
    parser.add_argument("--log", action="store_true",
                        help="Apply log10 to selected --plot quantity (except schlieren mode)")
    parser.add_argument("--raw", action="store_true",
                        help="Use pre-log quantity for selected --plot mode (except schlieren mode)")
    parser.add_argument("--cmap", default=None,
                        help="Colormap (default: gray_r for Schlieren, inferno for log/raw)")
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--rz", choices=["r", "z"], default="z",
                        help="For RZ geometry, choose horizontal axis (default: z)")
    parser.add_argument("--xy", choices=["x", "y"], default="x",
                        help="For XY geometry, choose horizontal axis (default: x)")
    parser.add_argument("--mirror", nargs="*", choices=["x", "y", "r", "z"], default=[],
                        help="Mirror axes list, e.g. --mirror x y or --mirror r")
    parser.add_argument("--prolong", choices=["linear", "nearest"], default="linear",
                        help="Coarse-to-fine visualization prolongation (default: linear)")
    parser.add_argument("--seam-smooth", type=float, default=0.0,
                        help="Local smoothing strength near box seams [0..1] (default: 0)")
    parser.add_argument("--seam-width", type=int, default=1,
                        help="Seam neighborhood width in cells (default: 1)")
    parser.add_argument("--seam-passes", type=int, default=1,
                        help="Number of seam smoothing passes (default: 1)")
    parser.add_argument("--gamma", type=float, default=1.4,
                        help="Specific heat ratio for Mach computation (default: 1.4)")
    parser.add_argument("--streamlines", action="store_true",
                        help="Overlay velocity streamlines")
    parser.add_argument("--stream-density", type=float, default=1.2,
                        help="Streamline density (default: 1.2)")
    parser.add_argument("--stream-color", default="white",
                        help="Streamline color (default: white)")
    parser.add_argument("--stream-lw", type=float, default=0.7,
                        help="Streamline line width (default: 0.7)")
    parser.add_argument("--stream-arrowsize", type=float, default=0.9,
                        help="Streamline arrow size (default: 0.9)")
    parser.add_argument("--boxes", action="store_true",
                        help="Overlay AMR patch boxes")
    parser.add_argument("--box-color", default=None,
                        help="Use one edge color for all boxes; default colors by level")
    parser.add_argument("--box-fill", action="store_true",
                        help="Fill box interiors; default is outline only")
    parser.add_argument("--box-alpha", type=float, default=None,
                        help="Box transparency (default: 0.18 for fill, 1.0 for outline)")
    parser.add_argument("--box-lw", type=float, default=0.7,
                        help="Box edge line width (default: 0.7)")
    parser.add_argument("--mask-solid", action="store_true",
                        help="Overlay solid regions (sld != 0) with a semi-transparent tint")
    parser.add_argument("--mask-color", default="white",
                        help="Overlay color for solid regions (default: white)")
    parser.add_argument("--mask-alpha", type=float, default=0.5,
                        help="Overlay opacity 0-1 (default: 0.5)")
    parser.add_argument("--flip-x", action="store_true",
                        help="Flip image horizontally (invert x-axis)")
    parser.add_argument("--flip-y", action="store_true",
                        help="Flip image vertically (invert y-axis)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output filename (default: schlieren.png)")
    return parser


def _validate_args(args):
    if args.time_std and args.time_fluctuation:
        print("Choose only one of --time-std or --time-fluctuation.")
        sys.exit(1)

    range_pairs = [
        ("x", args.x_min, args.x_max),
        ("y", args.y_min, args.y_max),
        ("r", args.r_min, args.r_max),
        ("z", args.z_min, args.z_max),
    ]
    for name, vmin, vmax in range_pairs:
        if vmin is not None and vmax is not None and vmin > vmax:
            print(f"Invalid range: --{name}-min ({vmin}) > --{name}-max ({vmax})")
            sys.exit(1)


def _load_single_field(args, varname):
    plotfile = args.plotfile
    if not plotfile:
        sel = _collect_plotfiles(args.plotdir)
        if not sel:
            print(f"No plotfiles found in {args.plotdir}/. Specify path explicitly.")
            sys.exit(1)
        plotfile = sel[-1][1]
        print(f"Auto-selected: {plotfile}")

    try:
        field, r_coords, z_coords, dr, dz, meta, var_label = _load_plot_field(
            plotfile,
            varname,
            args.gamma,
            out_level=args.level,
            prolong_method=args.prolong,
            seam_smooth=args.seam_smooth,
            seam_width=args.seam_width,
            seam_passes=args.seam_passes,
        )
    except ValueError:
        available = _read_header_varnames(plotfile)
        available.append("Mach (derived)")
        available.append("Vorticity (derived)")
        available.append("Divergence (derived)")
        available.append("x_velocity (derived)")
        available.append("y_velocity (derived)")
        available.append("velocity / speed (derived magnitude)")
        print(f"Variable '{varname}' not found.")
        print("Available variables:")
        print("  " + ", ".join(available))
        sys.exit(1)

    return {
        "field": field,
        "r_coords": r_coords,
        "z_coords": z_coords,
        "dr": dr,
        "dz": dz,
        "meta": meta,
        "var_label": var_label,
        "plotfile": plotfile,
        "time_info": None,
        "u0_mean": None,
        "u1_mean": None,
    }


def _load_time_window_field(args, varname):
    sel = _collect_plotfiles(
        args.plotdir,
        step_min=args.step_min,
        step_max=args.step_max,
        stride=args.stride,
    )
    if not sel:
        print("No plotfiles selected for time-window stats. Check --plotdir/--step-min/--step-max/--stride.")
        sys.exit(1)

    frame_steps, frame_paths = zip(*sel)
    nframes = len(frame_paths)
    print(f"Time-average over {nframes} frames: step {frame_steps[0]} -> {frame_steps[-1]}")

    if args.time_fluctuation:
        if args.fluct_step is None:
            fluct_step = frame_steps[-1]
        else:
            fluct_step = int(args.fluct_step)
            if fluct_step not in frame_steps:
                print(
                    f"--fluct-step {fluct_step} is not in selected frames. "
                    f"Selected range: {frame_steps[0]} -> {frame_steps[-1]}"
                )
                sys.exit(1)
    else:
        fluct_step = None

    field_sum = None
    field_sq_sum = None
    field_fluct_ref = None
    var_label = None
    meta = None
    r_coords = z_coords = None
    dr = dz = None

    # Performance: if streamlines are requested with mean mode, accumulate
    # velocity in the same pass to avoid re-reading all plotfiles later.
    need_vel_mean = args.streamlines and (not args.time_std) and (not args.time_fluctuation)
    u0_sum = None
    u1_sum = None
    time_first = None
    time_last = None
    fluct_time = None

    for ip, pth in enumerate(frame_paths):
        try:
            field_i, r_i, z_i, dr_i, dz_i, meta_i, var_label_i = _load_plot_field(
                pth,
                varname,
                args.gamma,
                out_level=args.level,
                prolong_method=args.prolong,
                seam_smooth=args.seam_smooth,
                seam_width=args.seam_width,
                seam_passes=args.seam_passes,
            )
        except ValueError:
            available = _read_header_varnames(pth)
            available.append("Mach (derived)")
            available.append("Vorticity (derived)")
            available.append("Divergence (derived)")
            available.append("x_velocity (derived)")
            available.append("y_velocity (derived)")
            available.append("velocity / speed (derived magnitude)")
            print(f"Variable '{varname}' not found in {pth}.")
            print("Available variables:")
            print("  " + ", ".join(available))
            sys.exit(1)

        if field_sum is None:
            field_sum = np.array(field_i, dtype=np.float64)
            if args.time_std:
                field_sq_sum = np.array(field_i * field_i, dtype=np.float64)
            r_coords, z_coords, dr, dz, meta = r_i, z_i, dr_i, dz_i, meta_i
            var_label = var_label_i
            time_first = meta_i.get("sim_time", None)
        else:
            if field_i.shape != field_sum.shape:
                print(
                    "Time-window stats require consistent grid shape across selected frames. "
                    "Try narrowing step range or set --level to a fixed AMR level (e.g. --level 0)."
                )
                print(f"Shape mismatch: {field_sum.shape} vs {field_i.shape} at {pth}")
                sys.exit(1)
            field_sum += field_i
            if args.time_std:
                field_sq_sum += field_i * field_i

        time_last = meta_i.get("sim_time", None)

        if args.time_fluctuation and frame_steps[ip] == fluct_step:
            field_fluct_ref = np.array(field_i, dtype=np.float64)
            fluct_time = meta_i.get("sim_time", None)

        if need_vel_mean:
            try:
                u0_i, u1_i, _, _, _, _, _ = _load_velocity_components(
                    pth,
                    out_level=args.level,
                    prolong_method=args.prolong,
                    seam_smooth=args.seam_smooth,
                    seam_width=args.seam_width,
                    seam_passes=args.seam_passes,
                )
            except ValueError as exc:
                print(f"Cannot plot streamlines: {exc}")
                sys.exit(1)

            if u0_sum is None:
                u0_sum = np.array(u0_i, dtype=np.float64)
                u1_sum = np.array(u1_i, dtype=np.float64)
            else:
                if u0_i.shape != u0_sum.shape or u1_i.shape != u1_sum.shape:
                    print("Cannot average streamlines: velocity grid shape mismatch across frames.")
                    sys.exit(1)
                u0_sum += u0_i
                u1_sum += u1_i

        if (ip + 1) % 20 == 0 or ip + 1 == nframes:
            print(f"  loaded {ip + 1}/{nframes} frames", flush=True)

    field_mean = field_sum / float(nframes)
    field = field_mean
    stat_mode = "mean"

    if args.time_std:
        var_field = np.maximum(field_sq_sum / float(nframes) - field_mean * field_mean, 0.0)
        field = np.sqrt(var_field)
        stat_mode = "std"
        var_label = f"std({var_label})"
    elif args.time_fluctuation:
        if field_fluct_ref is None:
            print("Internal error: fluctuation reference frame not found.")
            sys.exit(1)
        field = field_fluct_ref - field_mean
        stat_mode = "fluctuation"
        var_label = f"{var_label}'"

    time_info = {
        "nframes": nframes,
        "step_first": frame_steps[0],
        "step_last": frame_steps[-1],
        "paths": frame_paths,
        "stat_mode": stat_mode,
        "fluct_step": fluct_step,
        "time_first": time_first,
        "time_last": time_last,
        "fluct_time": fluct_time,
    }

    u0_mean = None
    u1_mean = None
    if need_vel_mean:
        u0_mean = u0_sum / float(nframes)
        u1_mean = u1_sum / float(nframes)

    return {
        "field": field,
        "r_coords": r_coords,
        "z_coords": z_coords,
        "dr": dr,
        "dz": dz,
        "meta": meta,
        "var_label": var_label,
        "plotfile": None,
        "time_info": time_info,
        "u0_mean": u0_mean,
        "u1_mean": u1_mean,
    }


def _build_base_plot_data(field, dr, dz, mode_spec, var_label, k_param):
    if mode_spec == "field":
        return field, "viridis", var_label, f"field_{_safe_name(var_label)}"

    grad_field = compute_grad_rho(field, dr, dz)
    grad_max = np.max(grad_field)
    print(f"|∇{var_label}| range: [{np.min(grad_field):.4g}, {grad_max:.4g}]")

    if mode_spec == "loggrad":
        return (
            np.log10(np.maximum(grad_field, 1e-30)),
            "inferno",
            rf"$\log_{{10}}|\nabla({var_label})|$",
            f"loggrad_{_safe_name(var_label)}",
        )
    if mode_spec == "grad":
        return (
            grad_field,
            "inferno",
            rf"$|\nabla({var_label})|$",
            f"grad_{_safe_name(var_label)}",
        )

    if grad_max > 0:
        schl = np.exp(-k_param * grad_field / grad_max)
    else:
        schl = np.ones_like(grad_field)
    return (
        schl,
        "gray_r",
        rf"Schlieren  $\exp(-{k_param:.0f}\,|\nabla({var_label})|/\max)$",
        f"schlieren_{_safe_name(var_label)}",
    )


def _apply_mode_modifiers(base_data, field, dr, dz, mode_spec, var_label, mode, args):
    plot_data = base_data
    if mode_spec == "field":
        label = var_label
    elif mode_spec == "grad":
        label = rf"$|\nabla({var_label})|$"
    elif mode_spec == "loggrad":
        label = rf"$\log_{{10}}|\nabla({var_label})|$"
    else:
        label = rf"Schlieren  $\exp(-{args.k:.0f}\,|\nabla({var_label})|/\max)$"

    if mode_spec != "schlieren":
        if args.raw and mode_spec == "loggrad":
            plot_data = compute_grad_rho(field, dr, dz)
            label = rf"$|\nabla({var_label})|$"
            mode = f"grad_{_safe_name(var_label)}"

        if args.log:
            plot_data = np.log10(np.maximum(plot_data, 1e-30))
            if mode_spec == "field":
                label = rf"$\log_{{10}}({var_label})$"
                mode = f"log_field_{_safe_name(var_label)}"
            elif mode_spec in ("grad", "loggrad"):
                label = rf"$\log_{{10}}|\nabla({var_label})|$"
                mode = f"loggrad_{_safe_name(var_label)}"
    elif args.log or args.raw:
        print("Note: --log/--raw are ignored for schlieren mode.")

    return plot_data, label, mode


def _apply_spatial_crop(args, x_coords, y_coords, data_full, x_axis_name, y_axis_name):
    x_bounds_display = (args.x_min, args.x_max)
    y_bounds_display = (args.y_min, args.y_max)
    rz_bounds = {
        "r": (args.r_min, args.r_max),
        "z": (args.z_min, args.z_max),
    }

    x_phys = rz_bounds.get(x_axis_name)
    y_phys = rz_bounds.get(y_axis_name)
    if x_axis_name == "x":
        x_phys = (args.x_min, args.x_max)
    if x_axis_name == "y":
        x_phys = (args.y_min, args.y_max)
    if y_axis_name == "x":
        y_phys = (args.x_min, args.x_max)
    if y_axis_name == "y":
        y_phys = (args.y_min, args.y_max)

    x_combined = _combine_bounds(x_bounds_display, x_phys)
    y_combined = _combine_bounds(y_bounds_display, y_phys)

    if x_combined is None and y_combined is None:
        return x_coords, y_coords, data_full, (0, x_coords.size, 0, y_coords.size)

    x_lo = x_combined[0] if x_combined is not None else None
    x_hi = x_combined[1] if x_combined is not None else None
    y_lo = y_combined[0] if y_combined is not None else None
    y_hi = y_combined[1] if y_combined is not None else None

    ix0, ix1 = _crop_1d_bounds(x_coords, x_lo, x_hi)
    iy0, iy1 = _crop_1d_bounds(y_coords, y_lo, y_hi)
    if ix1 <= ix0 or iy1 <= iy0:
        print("Requested crop window does not overlap data domain.")
        sys.exit(1)

    x_crop = x_coords[ix0:ix1]
    y_crop = y_coords[iy0:iy1]
    d_crop = data_full[iy0:iy1, ix0:ix1]
    print(
        f"Applied crop: {x_axis_name} in [{x_crop[0]:.5g}, {x_crop[-1]:.5g}], "
        f"{y_axis_name} in [{y_crop[0]:.5g}, {y_crop[-1]:.5g}]"
    )
    return x_crop, y_crop, d_crop, (ix0, ix1, iy0, iy1)


def _main_2d(args):
    """2D Schlieren visualization entry point."""
    do_time_window = args.time_average or args.time_std or args.time_fluctuation
    mode_spec, varname = _parse_plot_spec(args.plot)

    if do_time_window:
        data = _load_time_window_field(args, varname)
    else:
        data = _load_single_field(args, varname)

    field = data["field"]
    r_coords = data["r_coords"]
    z_coords = data["z_coords"]
    dr = data["dr"]
    dz = data["dz"]
    meta = data["meta"]
    var_label = data["var_label"]
    time_info = data["time_info"]
    plotfile = data["plotfile"]

    print(
        f"Geometry: {meta['coord_name']} (coord_sys={meta['coord_sys']}), "
        f"ndim={meta['ndim']}, finest_level={meta['finest_level']}, "
        f"output_level={meta.get('output_level', meta['finest_level'])}"
    )
    print(f"Grid: {field.shape[0]}×{field.shape[1]}, dr={dr:.5g}, dz={dz:.5g}")
    print(f"{var_label} range: [{field.min():.4g}, {field.max():.4g}]")
    print_box_stats(plotfile)

    base_data, default_cmap, _base_label, mode = _build_base_plot_data(
        field, dr, dz, mode_spec, var_label, args.k
    )
    plot_data, label, mode = _apply_mode_modifiers(
        base_data, field, dr, dz, mode_spec, var_label, mode, args
    )
    cmap = args.cmap if args.cmap else default_cmap

    x_coords, y_coords, data_full, x_axis_name, y_axis_name = _resolve_plot_axes(
        meta["coord_name"], args.rz, args.xy, r_coords, z_coords, plot_data
    )

    # Resolve solid mask in display orientation
    solid_mask_display = None
    if args.mask_solid:
        src_plt = plotfile if plotfile is not None else (
            time_info["paths"][-1] if time_info else None)
        if src_plt:
            solid_mask_ij = _load_solid_mask(
                src_plt,
                out_level=args.level, prolong_method=args.prolong,
                seam_smooth=args.seam_smooth,
                seam_width=args.seam_width, seam_passes=args.seam_passes,
            )
            if solid_mask_ij is not None:
                _, _, mask_tmp, _, _ = _resolve_plot_axes(
                    meta["coord_name"], args.rz, args.xy,
                    r_coords, z_coords, solid_mask_ij.astype(float),
                )
                solid_mask_display = mask_tmp > 0.5

    mirror_axes = _normalize_mirror_axes(args.mirror, x_axis_name, y_axis_name)

    if x_axis_name == "r":
        x_label = "r [m] (radial)"
    elif x_axis_name == "z":
        x_label = "z [m] (axial)"
    else:
        x_label = f"{x_axis_name} [m]"

    if y_axis_name == "r":
        y_label = "r [m] (radial)"
    elif y_axis_name == "z":
        y_label = "z [m] (axial)"
    else:
        y_label = f"{y_axis_name} [m]"

    odd_reflection_field = (var_label == "Vorticity" and mode_spec == "field")
    if "x" in mirror_axes:
        x_coords = np.concatenate([-x_coords[::-1], x_coords])
        left = data_full[:, ::-1]
        if odd_reflection_field:
            left = -left
        data_full = np.concatenate([left, data_full], axis=1)
    if "y" in mirror_axes:
        y_coords = np.concatenate([-y_coords[::-1], y_coords])
        top = data_full[::-1, :]
        if odd_reflection_field:
            top = -top
        data_full = np.concatenate([top, data_full], axis=0)

    # Mirror solid mask to match data
    if solid_mask_display is not None:
        if "x" in mirror_axes:
            solid_mask_display = np.concatenate(
                [solid_mask_display[:, ::-1], solid_mask_display], axis=1)
        if "y" in mirror_axes:
            solid_mask_display = np.concatenate(
                [solid_mask_display[::-1, :], solid_mask_display], axis=0)

    x_coords, y_coords, data_full, crop_slices = _apply_spatial_crop(
        args, x_coords, y_coords, data_full, x_axis_name, y_axis_name
    )

    # Crop solid mask to match data
    if solid_mask_display is not None:
        ix0, ix1, iy0, iy1 = crop_slices
        solid_mask_display = solid_mask_display[iy0:iy1, ix0:ix1]
        n_masked = int(np.sum(solid_mask_display))
        print(f"Solid mask: {n_masked} cells masked ({100.*n_masked/solid_mask_display.size:.1f}%)")

    # Data-level flip
    if args.flip_x:
        x_coords = -x_coords[::-1]
        data_full = data_full[:, ::-1]
        if solid_mask_display is not None:
            solid_mask_display = solid_mask_display[:, ::-1]
    if args.flip_y:
        y_coords = -y_coords[::-1]
        data_full = data_full[::-1, :]
        if solid_mask_display is not None:
            solid_mask_display = solid_mask_display[::-1, :]

    # Schlieren is physically bounded in [0, 1]. Use these defaults unless
    # user explicitly overrides vmin/vmax.
    if mode_spec == "schlieren":
        vmin_plot = 0.0 if args.vmin is None else args.vmin
        vmax_plot = 1.0 if args.vmax is None else args.vmax
    else:
        vmin_plot = args.vmin
        vmax_plot = args.vmax

    # Use nearest shading for discrete marker fields to avoid gouraud
    # interpolation artifacts that break thin (1-cell) features like GHS rings.
    _DISCRETE_FIELDS = {"sld", "ghs"}
    _is_discrete = varname.strip().lower() in _DISCRETE_FIELDS
    shading_mode = "nearest" if _is_discrete else "auto"

    fig, ax = plt.subplots(figsize=(10, 8) if mirror_axes else (10, 5))
    pcm = ax.pcolormesh(x_coords, y_coords, data_full, shading=shading_mode, cmap=cmap,
                        vmin=vmin_plot, vmax=vmax_plot)

    # Semi-transparent overlay on solid cells
    if solid_mask_display is not None:
        overlay = np.where(solid_mask_display, 1.0, np.nan)
        rgba = matplotlib.colors.to_rgba(args.mask_color, args.mask_alpha)
        overlay_cmap = ListedColormap([rgba])
        overlay_cmap.set_bad(alpha=0)
        ax.pcolormesh(x_coords, y_coords, overlay, shading=shading_mode,
                      cmap=overlay_cmap, vmin=0, vmax=1)

    # Prefer a vertical colorbar at right; fall back to horizontal only
    # when equal-aspect geometry makes the plot very narrow.
    x_span = max(float(np.ptp(x_coords)), 1.0e-30)
    y_span = max(float(np.ptp(y_coords)), 1.0e-30)
    use_horizontal_cbar = (x_span / y_span) < 0.38

    divider = make_axes_locatable(ax)
    if use_horizontal_cbar:
        cax = divider.append_axes("bottom", size="4.5%", pad=0.14)
        cbar = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    else:
        cax = divider.append_axes("right", size="3.5%", pad=0.10)
        cbar = fig.colorbar(pcm, cax=cax, orientation="vertical")
    cbar.set_label(label)

    # Keep zero tick visible when color limits span zero.
    vmin_eff, vmax_eff = pcm.get_clim()
    near_zero_floor = abs(vmin_eff) <= 1.0e-8 * max(abs(vmax_eff - vmin_eff), 1.0)
    need_zero_tick = (vmin_eff <= 0.0 <= vmax_eff) or near_zero_floor or mode_spec == "schlieren"
    if need_zero_tick:
        ticks = np.asarray(cbar.get_ticks(), dtype=float)
        span = max(abs(vmax_eff - vmin_eff), 1.0)
        tol = 1.0e-10 * span
        if ticks.size == 0 or np.min(np.abs(ticks)) > tol:
            ticks = np.sort(np.append(ticks, 0.0))
            cbar.set_ticks(ticks)

    if args.streamlines:
        if do_time_window and time_info is not None and time_info.get("stat_mode") in ("std", "fluctuation"):
            print("Streamlines are disabled for --time-std and --time-fluctuation outputs.")
            sys.exit(1)

        if do_time_window and data.get("u0_mean") is not None:
            u0 = data["u0_mean"]
            u1 = data["u1_mean"]
        else:
            src_plotfile = plotfile if plotfile is not None else time_info["paths"][-1]
            try:
                u0, u1, _, _, _, _, _ = _load_velocity_components(
                    src_plotfile,
                    out_level=args.level,
                    prolong_method=args.prolong,
                    seam_smooth=args.seam_smooth,
                    seam_width=args.seam_width,
                    seam_passes=args.seam_passes,
                )
            except ValueError as exc:
                print(f"Cannot plot streamlines: {exc}")
                sys.exit(1)

        u_plot, v_plot = _resolve_vector_axes(meta["coord_name"], args.rz, args.xy, u0, u1)
        if "x" in mirror_axes:
            u_plot = np.concatenate([-u_plot[:, ::-1], u_plot], axis=1)
            v_plot = np.concatenate([v_plot[:, ::-1], v_plot], axis=1)
        if "y" in mirror_axes:
            u_plot = np.concatenate([u_plot[::-1, :], u_plot], axis=0)
            v_plot = np.concatenate([-v_plot[::-1, :], v_plot], axis=0)

        ix0, ix1, iy0, iy1 = crop_slices
        u_plot = u_plot[iy0:iy1, ix0:ix1]
        v_plot = v_plot[iy0:iy1, ix0:ix1]

        ax.streamplot(
            x_coords,
            y_coords,
            u_plot,
            v_plot,
            density=args.stream_density,
            color=args.stream_color,
            linewidth=args.stream_lw,
            arrowsize=args.stream_arrowsize,
        )

    if args.boxes:
        src_plotfile = plotfile if plotfile is not None else time_info["paths"][-1]
        _draw_amr_boxes(ax, src_plotfile, meta, args, mirror_axes)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_aspect("equal")

    if mode_spec == "schlieren":
        title_root = f"Schlieren({var_label})"
    elif mode_spec == "grad":
        title_root = f"|grad({var_label})|"
    elif mode_spec == "loggrad":
        title_root = f"log10|grad({var_label})|"
    else:
        title_root = var_label

    if args.log and mode_spec != "schlieren":
        title_root = f"log10({title_root})"
    elif args.raw and mode_spec in ("grad", "loggrad"):
        title_root = f"|grad({var_label})|"

    if do_time_window:
        if time_info.get("stat_mode") == "std":
            ax.set_title(
                f"{title_root} — time std "
                f"[{time_info['step_first']}, {time_info['step_last']}], "
                f"N={time_info['nframes']}"
            )
        elif time_info.get("stat_mode") == "fluctuation":
            fluct_t = time_info.get("fluct_time", None)
            if fluct_t is None:
                ax.set_title(
                    f"{title_root} — fluctuation at step {time_info['fluct_step']} "
                    f"relative to mean [{time_info['step_first']}, {time_info['step_last']}]"
                )
            else:
                ax.set_title(
                    f"{title_root} — fluctuation at step {time_info['fluct_step']} (t={fluct_t:.6g}) "
                    f"relative to mean [{time_info['step_first']}, {time_info['step_last']}]"
                )
        else:
            ax.set_title(
                f"{title_root} — time average "
                f"[{time_info['step_first']}, {time_info['step_last']}], "
                f"N={time_info['nframes']}"
            )
    else:
        step = os.path.basename(os.path.normpath(plotfile)).replace("plt", "")
        sim_time = meta.get("sim_time", None)
        if sim_time is None:
            ax.set_title(f"{title_root} — step {step}")
        else:
            ax.set_title(f"{title_root} — step {step} (t={sim_time:.6g})")

    plt.tight_layout()
    outname = args.output if args.output else f"{mode}.png"
    # Higher DPI for discrete fields so 1-cell-wide features survive rasterization.
    save_dpi = 400 if _is_discrete else 200
    fig.savefig(outname, dpi=save_dpi, bbox_inches="tight")
    print(f"Saved: {outname}")
    plt.close()


# ---------------------------------------------------------------------------
# 3D main: three-view Schlieren
# ---------------------------------------------------------------------------
def _main_3d(args):
    """3D three-view Schlieren visualization entry point."""

    # ── Resolve plotfile ──────────────────────────────────────────────────
    plotfile = args.plotfile
    if not plotfile:
        sel = _collect_plotfiles(args.plotdir)
        if not sel:
            print(f"No plotfiles found in {args.plotdir}/")
            sys.exit(1)
        plotfile = sel[-1][1]
        print(f"Auto-selected: {plotfile}")

    # ── Parse plot spec ───────────────────────────────────────────────────
    mode_spec, varname = _parse_plot_spec(args.plot)

    # ── Load 3D field ─────────────────────────────────────────────────────
    print(f"Reading {plotfile} ...", flush=True)
    try:
        field, xc, yc, zc, dx, dy, dz, meta, var_label = _load_plot_field_3d(
            plotfile, varname, args.gamma, args.level,
            prolong_method=args.prolong,
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

    # ── Extract 2D slices ─────────────────────────────────────────────────
    slice_front = field[ix, :, :].T   # (Nz, Ny)
    slice_side  = field[:, iy, :].T   # (Nz, Nx)
    slice_top   = field[:, :, iz].T   # (Ny, Nx)

    # ── Build plot data for each slice ────────────────────────────────────
    def _process_slice(raw_2d, d_vert, d_horiz):
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
            data = _grad_magnitude_2d(clean, d_vert, d_horiz)
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

    data_front, dcmap = _process_slice(slice_front, dz, dy)
    data_side,  _     = _process_slice(slice_side,  dz, dx)
    data_top,   _     = _process_slice(slice_top,   dy, dx)

    cmap = args.cmap if args.cmap else dcmap

    # ── Solid mask ────────────────────────────────────────────────────────
    mask_front = mask_side = mask_top = None
    if args.mask_solid:
        try:
            sld = read_amrex_3d(plotfile, "sld", args.level,
                                prolong_method="nearest")[0]
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
        finite_v = all_vals[np.isfinite(all_vals)]
        vmin = float(np.min(finite_v)) if (args.vmin is None and finite_v.size > 0) else args.vmin
        vmax = float(np.max(finite_v)) if (args.vmax is None and finite_v.size > 0) else args.vmax

    # ── Plot ──────────────────────────────────────────────────────────────
    figsize = tuple(args.figsize) if args.figsize else (18, 5.5)
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # view_configs: (ax, data, hc, vc, hl, vl, title, mask, segs,
    #                  slice_val, horiz_min, horiz_max, vert_min, vert_max, norm_min, norm_max)
    view_configs = [
        (axes[0], data_front, yc, zc, "Y", "Z",
         f"Front (x={xc[ix]:.1f})", mask_front, stl_segs[0],
         xc[ix], "y_min", "y_max", "z_min", "z_max", "x_min", "x_max"),
        (axes[1], data_side, xc, zc, "X", "Z",
         f"Side (y={yc[iy]:.1f})", mask_side, stl_segs[1],
         yc[iy], "x_min", "x_max", "z_min", "z_max", "y_min", "y_max"),
        (axes[2], data_top, xc, yc, "X", "Y",
         f"Top (z={zc[iz]:.1f})", mask_top, stl_segs[2],
         zc[iz], "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"),
    ]

    last_pcm = None
    for (ax, data, hc, vc, hl, vl, title, mask, segs,
         sv, hkmin, hkmax, vkmin, vkmax, nkmin, nkmax) in view_configs:
        pcm = ax.pcolormesh(hc, vc, data, shading="auto", cmap=cmap,
                            vmin=vmin, vmax=vmax)
        last_pcm = pcm

        if mask is not None and np.any(mask):
            rgba = matplotlib.colors.to_rgba(args.mask_color, args.mask_alpha)
            overlay_cmap = ListedColormap([rgba])
            overlay_cmap.set_bad(alpha=0)
            overlay = np.where(mask, 1.0, np.nan)
            ax.pcolormesh(hc, vc, overlay, shading="auto",
                          cmap=overlay_cmap, vmin=0, vmax=1)

        _draw_amr_boxes_3d(ax, plotfile, args, None, sv,
                           hkmin, hkmax, vkmin, vkmax, nkmin, nkmax)

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


# ---------------------------------------------------------------------------
# Entry point: auto-detect dimension and dispatch
# ---------------------------------------------------------------------------
def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    # ── Resolve plotfile for auto-detection ────────────────────────────────
    plotfile = args.plotfile
    if not plotfile:
        do_time = getattr(args, "time_average", False) or \
                  getattr(args, "time_std", False) or \
                  getattr(args, "time_fluctuation", False)
        if not do_time:
            sel = _collect_plotfiles(args.plotdir)
            if sel:
                plotfile = sel[-1][1]

    # ── Determine mode ────────────────────────────────────────────────────
    if args.force_3d:
        ndim = 3
    elif args.force_2d:
        ndim = 2
    elif plotfile and os.path.isdir(plotfile):
        ndim = _detect_ndim(plotfile)
        print(f"Auto-detected {ndim}D plotfile")
    else:
        ndim = 2  # default for time-window modes without explicit plotfile

    # ── Validate and dispatch ─────────────────────────────────────────────
    if ndim == 3:
        _main_3d(args)
    else:
        _validate_args(args)
        _main_2d(args)


if __name__ == "__main__":
    main()
