#!/usr/bin/env python3
"""
Visualise degraded ghost points from GP-Diag CSV dumps.

Usage:
  python plot_gp_degraded.py gp_degraded_s0_L0.csv [--geom diamond_wedge_2d.dat]
  python plot_gp_degraded.py gp_degraded_s*.csv --geom diamond_wedge_2d.dat --compare
  python plot_gp_degraded.py gp_degraded_s0_L0.csv --geom ConfigA1_A0.stl  (3D)

Reads the CSV produced by reportGPDiagnostics and overlays degraded GPs
on the geometry surface.  Points are coloured by severity:
  - MILD  (n_fluid > ideal/2)  : orange
  - SEVERE (n_fluid <= ideal/2) : red
0th-order GPs are marked with a star.
"""

import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── helpers ──────────────────────────────────────────────────────────────────

def load_csv(path):
    """Load a GP-Diag CSV into a structured numpy array."""
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def parse_step_level(fname):
    """Extract (step, level) from filename like gp_degraded_s2_L0.csv."""
    m = re.search(r"_s(\d+)_L(\d+)", os.path.basename(fname))
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def load_geom_2d(path):
    """Load a 2-column .dat polygon file.  Returns (N, 2) array."""
    pts = np.loadtxt(path)
    if pts.ndim == 1:
        pts = pts.reshape(-1, 2)
    # close polygon if not already closed
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])
    return pts


def load_stl_triangles(path):
    """Load an ASCII or binary STL, return (M, 3, 3) triangle vertices."""
    try:
        from stl import mesh as stl_mesh
        m = stl_mesh.Mesh.from_file(path)
        return m.vectors  # (M, 3, 3)
    except ImportError:
        print("Warning: python-stl not installed; skipping 3D geometry overlay.")
        print("  Install with:  pip install numpy-stl")
        return None


# ── plotting ─────────────────────────────────────────────────────────────────

def classify(data):
    """Return (mild_mask, severe_mask, zeroth_mask)."""
    ideal = data["ideal"]
    nf    = data["n_fluid"]
    eff   = data["eff_order"]
    mild   = nf > ideal // 2
    severe = ~mild
    zeroth = eff == 0
    return mild, severe, zeroth


def plot_2d(csv_files, geom_path=None, outfile=None, compare=False):
    """2D scatter plot of degraded GPs overlaid on geometry."""

    if compare and len(csv_files) > 1:
        # Multi-panel: one subplot per CSV
        ncols = min(len(csv_files), 4)
        nrows = (len(csv_files) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows),
                                 squeeze=False)
        axes_flat = axes.flatten()
        for i, f in enumerate(csv_files):
            ax = axes_flat[i]
            _plot_2d_single(ax, f, geom_path)
        for j in range(len(csv_files), len(axes_flat)):
            axes_flat[j].set_visible(False)
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(figsize=(10, 8))
        _plot_2d_single(ax, csv_files[0], geom_path)
        fig.tight_layout()

    out = outfile or "gp_degraded_2d.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close(fig)


def _plot_2d_single(ax, csv_file, geom_path):
    """Plot a single 2D panel."""
    data = load_csv(csv_file)
    step, lev = parse_step_level(csv_file)
    n_total = len(data)

    mild, severe, zeroth = classify(data)

    # Geometry background
    if geom_path and os.path.exists(geom_path):
        poly = load_geom_2d(geom_path)
        ax.fill(poly[:, 0], poly[:, 1], color="lightblue", alpha=0.3, label="solid")
        ax.plot(poly[:, 0], poly[:, 1], "k-", linewidth=1.0)

    # Plot IB surface points as thin dots for context
    ax.scatter(data["ib_x"], data["ib_y"], s=4, c="grey", alpha=0.4, zorder=1,
               label="IB surface pt")

    # GP→IB lines (shows GPs are near-wall)
    for i in range(n_total):
        ax.plot([data["gp_x"][i], data["ib_x"][i]],
                [data["gp_y"][i], data["ib_y"][i]],
                "k-", alpha=0.08, linewidth=0.5, zorder=1)

    # Mild degraded GPs
    if np.any(mild):
        sc1 = ax.scatter(data["gp_x"][mild], data["gp_y"][mild],
                         c=data["max_weight"][mild], cmap="YlOrRd",
                         vmin=0.25, vmax=1.0,
                         s=30, edgecolors="orange", linewidths=0.5,
                         marker="o", zorder=3, label="MILD")

    # Severe degraded GPs
    if np.any(severe):
        sc2 = ax.scatter(data["gp_x"][severe], data["gp_y"][severe],
                         c=data["max_weight"][severe], cmap="YlOrRd",
                         vmin=0.25, vmax=1.0,
                         s=60, edgecolors="red", linewidths=1.0,
                         marker="s", zorder=4, label="SEVERE")

    # 0th-order overlay
    if np.any(zeroth):
        ax.scatter(data["gp_x"][zeroth], data["gp_y"][zeroth],
                   s=120, marker="*", facecolors="none", edgecolors="magenta",
                   linewidths=1.5, zorder=5, label="0th-order")

    # Colour bar
    if np.any(mild) or np.any(severe):
        # use whichever scatter has data
        sc = sc1 if np.any(mild) else sc2
        cb = plt.colorbar(sc, ax=ax, shrink=0.6)
        cb.set_label("max single-point weight")

    # FAB boundary annotation: highlight GPs from different FABs
    unique_fabs = np.unique(data["fab"])
    if len(unique_fabs) > 1:
        for fab_id in unique_fabs:
            mask = data["fab"] == fab_id
            if np.sum(mask) > 0:
                cx = np.mean(data["gp_x"][mask])
                cy = np.mean(data["gp_y"][mask])
                ax.annotate(f"FAB{fab_id}", (cx, cy), fontsize=6, alpha=0.5,
                            ha="center")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.legend(fontsize=7, loc="upper right")

    # n_fluid breakdown text
    nf_vals, nf_counts = np.unique(data["n_fluid"], return_counts=True)
    ideal_val = data["ideal"][0] if len(data) > 0 else "?"
    txt = f"Step {step}  L{lev}  |  {n_total} degraded GPs\n"
    for v, c in zip(nf_vals, nf_counts):
        txt += f"  {v}/{ideal_val}: {c}  "
    ax.set_title(txt, fontsize=9)


def plot_3d(csv_files, geom_path=None, outfile=None):
    """3D scatter plot of degraded GPs overlaid on geometry."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    data = load_csv(csv_files[0])
    step, lev = parse_step_level(csv_files[0])
    n_total = len(data)

    mild, severe, zeroth = classify(data)

    # Geometry overlay (transparent triangles)
    if geom_path and os.path.exists(geom_path):
        tris = load_stl_triangles(geom_path)
        if tris is not None:
            # subsample for visibility if too many
            n_tris = len(tris)
            if n_tris > 5000:
                idx = np.random.default_rng(42).choice(n_tris, 5000, replace=False)
                tris = tris[idx]
            poly = Poly3DCollection(tris, alpha=0.08, facecolor="skyblue",
                                    edgecolor="grey", linewidth=0.1)
            ax.add_collection3d(poly)

    # Mild
    if np.any(mild):
        ax.scatter(data["gp_x"][mild], data["gp_y"][mild], data["gp_z"][mild],
                   c=data["max_weight"][mild], cmap="YlOrRd",
                   vmin=0.25, vmax=1.0,
                   s=20, edgecolors="orange", linewidths=0.3,
                   marker="o", label="MILD", depthshade=True)

    # Severe
    if np.any(severe):
        ax.scatter(data["gp_x"][severe], data["gp_y"][severe], data["gp_z"][severe],
                   c=data["max_weight"][severe], cmap="YlOrRd",
                   vmin=0.25, vmax=1.0,
                   s=50, edgecolors="red", linewidths=0.8,
                   marker="s", label="SEVERE", depthshade=True)

    # 0th-order
    if np.any(zeroth):
        ax.scatter(data["gp_x"][zeroth], data["gp_y"][zeroth], data["gp_z"][zeroth],
                   s=100, marker="*", facecolors="none", edgecolors="magenta",
                   linewidths=1.5, label="0th-order")

    nf_vals, nf_counts = np.unique(data["n_fluid"], return_counts=True)
    ideal_val = data["ideal"][0] if len(data) > 0 else "?"
    txt = f"Step {step}  L{lev}  |  {n_total} degraded GPs\n"
    for v, c in zip(nf_vals, nf_counts):
        txt += f"  {v}/{ideal_val}: {c}  "
    ax.set_title(txt, fontsize=9)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend(fontsize=7)

    out = outfile or "gp_degraded_3d.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close(fig)


def plot_step_evolution(csv_files, geom_path=None, outfile=None):
    """Show step-by-step evolution: how degradation changes with geometry motion."""
    # group by (step, level)
    records = []
    for f in sorted(csv_files):
        step, lev = parse_step_level(f)
        if step is None:
            continue
        data = load_csv(f)
        n_total = len(data)
        mild_m, severe_m, zeroth_m = classify(data)
        records.append({
            "step": step, "level": lev,
            "n_total": n_total,
            "n_mild": int(np.sum(mild_m)),
            "n_severe": int(np.sum(severe_m)),
            "n_zeroth": int(np.sum(zeroth_m)),
            "mean_weight": float(np.mean(data["max_weight"])) if n_total else 0,
            "max_weight": float(np.max(data["max_weight"])) if n_total else 0,
        })

    if not records:
        print("No valid CSV files found for evolution plot.")
        return

    # separate by level
    levels = sorted(set(r["level"] for r in records))

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    for lev in levels:
        recs = [r for r in records if r["level"] == lev]
        steps = [r["step"] for r in recs]
        ax = axes[0]
        ax.plot(steps, [r["n_total"] for r in recs], "o-",
                label=f"L{lev} total")
        ax.plot(steps, [r["n_mild"] for r in recs], "s--",
                label=f"L{lev} mild", alpha=0.7)
        ax.plot(steps, [r["n_severe"] for r in recs], "^--",
                label=f"L{lev} severe", alpha=0.7)

        ax2 = axes[1]
        ax2.plot(steps, [r["max_weight"] for r in recs], "o-",
                 label=f"L{lev} worst_max_wt")
        ax2.plot(steps, [r["mean_weight"] for r in recs], "s--",
                 label=f"L{lev} mean_max_wt", alpha=0.7)

    axes[0].set_ylabel("degraded GP count")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)
    axes[1].set_ylabel("max single-point weight")
    axes[1].set_xlabel("step")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)
    fig.suptitle("GP degradation evolution with geometry motion", fontsize=11)
    fig.tight_layout()

    out = outfile or "gp_degraded_evolution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Visualise degraded ghost points")
    parser.add_argument("csv", nargs="+", help="GP-Diag CSV file(s)")
    parser.add_argument("--geom", help="Geometry file (.dat for 2D, .stl for 3D)")
    parser.add_argument("--compare", action="store_true",
                        help="Multi-panel comparison of CSV files")
    parser.add_argument("--evolution", action="store_true",
                        help="Plot step-by-step evolution of degradation")
    parser.add_argument("-o", "--output", help="Output image filename")
    parser.add_argument("--dim", type=int, choices=[2, 3], default=None,
                        help="Force 2D or 3D mode (auto-detected from z coords)")
    args = parser.parse_args()

    # Expand globs (in case shell didn't)
    csv_files = []
    for pat in args.csv:
        expanded = sorted(glob.glob(pat))
        csv_files.extend(expanded if expanded else [pat])

    if not csv_files or not os.path.exists(csv_files[0]):
        print(f"Error: no CSV files found: {args.csv}")
        sys.exit(1)

    # Evolution mode
    if args.evolution:
        plot_step_evolution(csv_files, args.geom, args.output)
        return

    # Auto-detect dimension
    sample = load_csv(csv_files[0])
    is_3d = args.dim == 3 if args.dim else (np.any(np.abs(sample["gp_z"]) > 1e-12))

    if is_3d:
        plot_3d(csv_files, args.geom, args.output)
    else:
        plot_2d(csv_files, args.geom, args.output, compare=args.compare)


if __name__ == "__main__":
    main()
