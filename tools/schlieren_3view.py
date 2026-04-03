#!/usr/bin/env python3
"""
schlieren_3view.py — 3D numerical-Schlieren three-view plot from AMReX plotfiles.

No yt dependency. Reads AMReX binary plotfiles directly with numpy.

Usage:
  python3 schlieren_3view.py <plotfile> [options]

Examples:
  python3 schlieren_3view.py wrk/su47/plot/plt00010
  python3 schlieren_3view.py wrk/su47/plot/plt00010 --origin 90 90 72 --field Density
  python3 schlieren_3view.py wrk/su47/plot/plt00010 --stl wrk/su47/SU47_repaired.stl
"""

from __future__ import annotations
import argparse, glob, os, re, struct, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from typing import List, Optional, Tuple


# ── AMReX plotfile reader (level-0 only, single-file) ────────────────────────

def read_header(plotdir: str):
    """Parse the top-level Header file."""
    with open(os.path.join(plotdir, "Header"), "r") as f:
        lines = f.readlines()

    idx = 0
    _version = lines[idx].strip(); idx += 1
    ncomp = int(lines[idx].strip()); idx += 1
    field_names = []
    for _ in range(ncomp):
        field_names.append(lines[idx].strip()); idx += 1
    ndim = int(lines[idx].strip()); idx += 1
    _time = float(lines[idx].strip()); idx += 1
    _max_level = int(lines[idx].strip()); idx += 1

    lo = [float(x) for x in lines[idx].split()]; idx += 1
    hi = [float(x) for x in lines[idx].split()]; idx += 1

    # Skip refinement ratios (could be empty line)
    idx += 1  # blank or ratios

    # Domain box: ((lo) (hi) (type))
    box_line = lines[idx].strip(); idx += 1
    # Parse e.g. "((0,0,0) (127,127,103) (0,0,0))"
    nums = re.findall(r'-?\d+', box_line)
    domain_lo = [int(nums[i]) for i in range(ndim)]
    domain_hi = [int(nums[ndim + i]) for i in range(ndim)]
    domain_dims = [domain_hi[i] - domain_lo[i] + 1 for i in range(ndim)]

    # Skip nsteps
    idx += 1
    # Cell sizes
    dx = [float(x) for x in lines[idx].split()]; idx += 1

    return {
        "ndim": ndim,
        "ncomp": ncomp,
        "field_names": field_names,
        "lo": np.array(lo[:ndim]),
        "hi": np.array(hi[:ndim]),
        "dims": domain_dims,
        "dx": np.array(dx[:ndim]),
    }


def read_level0(plotdir: str, header: dict):
    """Read level-0 data into a single 3D array (ncomp, nz, ny, nx)."""
    ndim = header["ndim"]
    ncomp = header["ncomp"]
    dims = header["dims"]  # [nx, ny, nz]

    cell_h = os.path.join(plotdir, "Level_0", "Cell_H")
    with open(cell_h, "r") as f:
        lines = f.readlines()

    idx = 0
    _version = int(lines[idx].strip()); idx += 1
    _how = int(lines[idx].strip()); idx += 1
    _ncomp_h = int(lines[idx].strip()); idx += 1
    _nghost = int(lines[idx].strip()); idx += 1

    # Parse box array: (nboxes nghost\n boxes... )
    box_line = lines[idx].strip(); idx += 1
    nboxes = int(box_line.lstrip("(").split()[0])

    boxes = []
    for _ in range(nboxes):
        bline = lines[idx].strip(); idx += 1
        nums = [int(x) for x in re.findall(r'-?\d+', bline)]
        blo = nums[:ndim]
        bhi = nums[ndim:2*ndim]
        boxes.append((blo, bhi))
    idx += 1  # closing ")"

    # Number of FABs
    _nfabs = int(lines[idx].strip()); idx += 1

    # FabOnDisk lines
    fab_info = []
    for _ in range(nboxes):
        parts = lines[idx].strip().split(); idx += 1
        fname = parts[1]
        offset = int(parts[2])
        fab_info.append((fname, offset))

    # Allocate full grid
    if ndim == 3:
        nx, ny, nz = dims
        full = np.zeros((ncomp, nz, ny, nx), dtype=np.float64)
    else:
        nx, ny = dims[:2]
        full = np.zeros((ncomp, ny, nx), dtype=np.float64)

    # Read each FAB
    level0_dir = os.path.join(plotdir, "Level_0")
    for bi in range(nboxes):
        blo, bhi = boxes[bi]
        fname, offset = fab_info[bi]
        fpath = os.path.join(level0_dir, fname)

        with open(fpath, "rb") as f:
            f.seek(offset)
            # Read text header line
            hdr = f.readline()

            if ndim == 3:
                bx = bhi[0] - blo[0] + 1
                by = bhi[1] - blo[1] + 1
                bz = bhi[2] - blo[2] + 1
                ncells = bx * by * bz
                data = np.fromfile(f, dtype=np.float64, count=ncells * ncomp)
                # AMReX stores (comp, z, y, x) in Fortran order → read as-is
                data = data.reshape((ncomp, bz, by, bx))
                full[:, blo[2]:bhi[2]+1, blo[1]:bhi[1]+1, blo[0]:bhi[0]+1] = data
            else:
                bx = bhi[0] - blo[0] + 1
                by = bhi[1] - blo[1] + 1
                ncells = bx * by
                data = np.fromfile(f, dtype=np.float64, count=ncells * ncomp)
                data = data.reshape((ncomp, by, bx))
                full[:, blo[1]:bhi[1]+1, blo[0]:bhi[0]+1] = data

    return full


# ── Schlieren computation ─────────────────────────────────────────────────────

def schlieren_2d(field_2d: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Compute |∇field| on a 2D slice using central differences."""
    gy, gx = np.gradient(field_2d, dy, dx)
    return np.sqrt(gx**2 + gy**2)


def numerical_schlieren(grad_mag: np.ndarray, beta: float = 0.7) -> np.ndarray:
    """Convert gradient magnitude to Schlieren intensity: exp(-beta * |∇ρ| / max|∇ρ|)."""
    gmax = grad_mag.max()
    if gmax < 1e-30:
        return np.ones_like(grad_mag)
    return np.exp(-beta * grad_mag / gmax)


# ── STL geometry outline ──────────────────────────────────────────────────────

def read_stl_vertices(stl_path: str) -> np.ndarray:
    """Read all triangle vertices from a binary STL. Returns (N*3, 3) array."""
    with open(stl_path, "rb") as f:
        f.read(80)
        ntri = struct.unpack('<I', f.read(4))[0]
        dt = np.dtype([('normal','<f4',3), ('v0','<f4',3),
                       ('v1','<f4',3), ('v2','<f4',3), ('attr','<u2')])
        data = np.frombuffer(f.read(ntri * 50), dtype=dt)
    return np.vstack([data['v0'], data['v1'], data['v2']]).astype(np.float64), data


def stl_cross_section(stl_path: str, axis: int, cut_pos: float,
                      tol: float = 2.0) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Compute the cross-section of an STL at a given plane.
    Returns list of line segments as (p1[2], p2[2]) in the 2D plane.

    axis: 0=x, 1=y, 2=z
    cut_pos: position along that axis
    tol: tolerance for selecting triangles near the cutting plane
    """
    all_verts, data = read_stl_vertices(stl_path)
    v0 = data['v0'].astype(np.float64)
    v1 = data['v1'].astype(np.float64)
    v2 = data['v2'].astype(np.float64)

    # Filter triangles that straddle the cutting plane
    d0 = v0[:, axis] - cut_pos
    d1 = v1[:, axis] - cut_pos
    d2 = v2[:, axis] - cut_pos

    # Triangle crosses the plane if not all vertices are on the same side
    crosses = ~((d0 > 0) & (d1 > 0) & (d2 > 0)) & ~((d0 < 0) & (d1 < 0) & (d2 < 0))
    idx = np.where(crosses)[0]

    # Map axes: for cut along axis, the 2D coords are the other two axes
    axes_2d = [i for i in range(3) if i != axis]

    segments = []
    for i in idx:
        tri = np.array([v0[i], v1[i], v2[i]])
        d = tri[:, axis] - cut_pos
        pts = []
        for e0, e1 in [(0,1), (1,2), (2,0)]:
            if d[e0] * d[e1] < 0:  # edge crosses plane
                t = d[e0] / (d[e0] - d[e1])
                p = tri[e0] + t * (tri[e1] - tri[e0])
                pts.append(p[axes_2d])
            elif abs(d[e0]) < 1e-10:
                pts.append(tri[e0][axes_2d])
            elif abs(d[e1]) < 1e-10:
                pts.append(tri[e1][axes_2d])
        if len(pts) >= 2:
            segments.append((pts[0], pts[1]))

    return segments


def draw_stl_outline(ax, segments: list, color='lime', lw=0.3, alpha=0.8):
    """Draw STL cross-section segments on a matplotlib axes."""
    for p1, p2 in segments:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, lw=lw, alpha=alpha)


# ── Main plotting ─────────────────────────────────────────────────────────────

def get_slice_and_coords(full: np.ndarray, header: dict, comp_idx: int,
                         axis: int, slice_idx: int):
    """
    Extract a 2D slice from 3D data.
    axis: 0=x, 1=y, 2=z
    Returns: slice_2d, (h_coords, v_coords), (dh, dv)
    """
    lo = header["lo"]
    dx = header["dx"]
    dims = header["dims"]  # [nx, ny, nz]
    # full shape: (ncomp, nz, ny, nx)

    if axis == 0:  # x-normal: show (y, z)
        slc = full[comp_idx, :, :, slice_idx]  # (nz, ny)
        h = lo[1] + (np.arange(dims[1]) + 0.5) * dx[1]
        v = lo[2] + (np.arange(dims[2]) + 0.5) * dx[2]
        dh, dv = dx[1], dx[2]
    elif axis == 1:  # y-normal: show (x, z)
        slc = full[comp_idx, :, slice_idx, :]  # (nz, nx)
        h = lo[0] + (np.arange(dims[0]) + 0.5) * dx[0]
        v = lo[2] + (np.arange(dims[2]) + 0.5) * dx[2]
        dh, dv = dx[0], dx[2]
    else:  # z-normal: show (x, y)
        slc = full[comp_idx, slice_idx, :, :]  # (ny, nx)
        h = lo[0] + (np.arange(dims[0]) + 0.5) * dx[0]
        v = lo[1] + (np.arange(dims[1]) + 0.5) * dx[1]
        dh, dv = dx[0], dx[1]

    return slc, (h, v), (dh, dv)


def pos_to_index(pos: float, lo: float, dx: float, n: int) -> int:
    idx = int((pos - lo) / dx)
    return max(0, min(n - 1, idx))


def main():
    parser = argparse.ArgumentParser(
        description="3D numerical-Schlieren three-view from AMReX plotfiles (no yt).")
    parser.add_argument("plotfile", help="Path to AMReX plotfile directory")
    parser.add_argument("--origin", type=float, nargs=3, default=None,
                        metavar=("X", "Y", "Z"),
                        help="Slice origin (default: domain center)")
    parser.add_argument("--field", default="Density",
                        help="Field for Schlieren gradient (default: Density)")
    parser.add_argument("--beta", type=float, default=0.7,
                        help="Schlieren exponential coefficient (default: 0.7)")
    parser.add_argument("--stl", default=None,
                        help="STL file path (optional, draws cross-section outline)")
    parser.add_argument("--stl-color", default="lime",
                        help="STL outline color (default: lime)")
    parser.add_argument("--stl-lw", type=float, default=0.3,
                        help="STL outline linewidth (default: 0.3)")
    parser.add_argument("--cmap", default="gray_r",
                        help="Colormap (default: gray_r)")
    parser.add_argument("--mask-solid", action="store_true",
                        help="Mask solid cells (sld>0.5) with a flat color")
    parser.add_argument("--dpi", type=int, default=200,
                        help="Output DPI (default: 200)")
    parser.add_argument("--output", "-o", default=None,
                        help="Output filename (default: schlieren_3view.png)")
    args = parser.parse_args()

    # ── Load data ──
    plotdir = args.plotfile
    if not os.path.isdir(plotdir):
        # Try as a directory containing plt*
        candidates = sorted(glob.glob(os.path.join(plotdir, "plt*")),
                            key=lambda s: [int(t) if t.isdigit() else t
                                           for t in re.split(r'(\d+)', s)])
        if candidates:
            plotdir = candidates[-1]
    if not os.path.exists(os.path.join(plotdir, "Header")):
        print(f"Error: {plotdir} is not a valid AMReX plotfile.", file=sys.stderr)
        return 1

    print(f"Reading {plotdir} ...")
    header = read_header(plotdir)
    ndim = header["ndim"]
    if ndim != 3:
        print("Error: this tool is for 3D plotfiles only.", file=sys.stderr)
        return 1

    full = read_level0(plotdir, header)
    print(f"  Grid: {header['dims']}, {header['ncomp']} components")
    print(f"  Domain: {header['lo']} → {header['hi']}")
    print(f"  Fields: {header['field_names']}")

    # ── Resolve field index ──
    if args.field not in header["field_names"]:
        print(f"Error: field '{args.field}' not found. Available: {header['field_names']}",
              file=sys.stderr)
        return 1
    comp_idx = header["field_names"].index(args.field)

    # ── Resolve sld index (for masking) ──
    sld_idx = None
    if args.mask_solid and "sld" in header["field_names"]:
        sld_idx = header["field_names"].index("sld")

    # ── Slice origin ──
    lo = header["lo"]
    hi = header["hi"]
    dx = header["dx"]
    dims = header["dims"]

    if args.origin is not None:
        origin = np.array(args.origin)
    else:
        origin = 0.5 * (lo + hi)

    ix = pos_to_index(origin[0], lo[0], dx[0], dims[0])
    iy = pos_to_index(origin[1], lo[1], dx[1], dims[1])
    iz = pos_to_index(origin[2], lo[2], dx[2], dims[2])
    print(f"  Slice origin: ({origin[0]:.2f}, {origin[1]:.2f}, {origin[2]:.2f}) → indices ({ix}, {iy}, {iz})")

    # ── Three views ──
    view_cfg = [
        {"axis": 0, "idx": ix, "title": f"Front (x={origin[0]:.1f})",
         "xlabel": "Y", "ylabel": "Z"},
        {"axis": 1, "idx": iy, "title": f"Side (y={origin[1]:.1f})",
         "xlabel": "X", "ylabel": "Z"},
        {"axis": 2, "idx": iz, "title": f"Top (z={origin[2]:.1f})",
         "xlabel": "X", "ylabel": "Y"},
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), dpi=args.dpi)

    for ax, vc in zip(axes, view_cfg):
        slc, (h, v), (dh, dv) = get_slice_and_coords(
            full, header, comp_idx, vc["axis"], vc["idx"])

        # Compute Schlieren
        grad_mag = schlieren_2d(slc, dh, dv)
        schlieren = numerical_schlieren(grad_mag, beta=args.beta)

        # Optional solid mask
        if sld_idx is not None:
            sld_slc, _, _ = get_slice_and_coords(
                full, header, sld_idx, vc["axis"], vc["idx"])
            schlieren[sld_slc > 0.5] = 0.0  # black for solid

        # Plot
        extent = [h[0] - 0.5*dh, h[-1] + 0.5*dh,
                  v[0] - 0.5*dv, v[-1] + 0.5*dv]
        ax.imshow(schlieren, origin="lower", extent=extent,
                  cmap=args.cmap, vmin=0, vmax=1, aspect="equal",
                  interpolation="bilinear")

        # STL outline
        if args.stl:
            cut_pos = origin[vc["axis"]]
            segs = stl_cross_section(args.stl, vc["axis"], cut_pos)
            draw_stl_outline(ax, segs, color=args.stl_color, lw=args.stl_lw)
            ax.set_title(f"{vc['title']}  ({len(segs)} STL segs)", fontsize=10)
        else:
            ax.set_title(vc["title"], fontsize=10)

        ax.set_xlabel(vc["xlabel"])
        ax.set_ylabel(vc["ylabel"])

    fig.suptitle(f"Numerical Schlieren — {args.field}  |  {os.path.basename(plotdir)}",
                 fontsize=12, y=1.01)
    fig.tight_layout()

    outname = args.output or "schlieren_3view.png"
    fig.savefig(outname, bbox_inches="tight", dpi=args.dpi)
    plt.close(fig)
    print(f"Saved: {outname}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
