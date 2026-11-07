#!/usr/bin/env python3
"""
3D FSI Moving Sphere — Three-view schlieren movie.

Renders three orthogonal slices (XY at z=0, XZ at y=0, YZ at x=0) of
the density field with solid regions (sld != 0) masked in white.

Output: sphere_fsi_3view.mp4
"""

import os, glob, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import yt
yt.set_log_level("error")

# ============================================================================
PLOT_DIR  = "./plot"
OUTPUT    = "sphere_fsi_3view.mp4"
DPI       = 180
FPS       = 6

# Density range for colormap
RHO_MIN, RHO_MAX = 0.3, 5.0

# ============================================================================
plt_dirs = sorted(glob.glob(os.path.join(PLOT_DIR, "plt[0-9]*")))
print(f"Found {len(plt_dirs)} plot files")

frames_dir = os.path.join(PLOT_DIR, "_frames_3view")
os.makedirs(frames_dir, exist_ok=True)

# Custom colormap: inferno but with white for solid
cmap_base = plt.cm.inferno.copy()

for idx, pd in enumerate(plt_dirs):
    step = int(re.search(r'(\d+)', os.path.basename(pd)).group(1))
    ds = yt.load(pd)
    time_val = float(ds.current_time)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.suptitle(f"Step {step},  t = {time_val:.4e} s", fontsize=14, fontweight='bold')

    slice_configs = [
        ("z", 0.0, "XY plane (z=0)", "x", "y"),
        ("y", 0.0, "XZ plane (y=0)", "x", "z"),
        ("x", 0.0, "YZ plane (x=0)", "y", "z"),
    ]

    for ax, (axis, coord, title, xlabel, ylabel) in zip(axes, slice_configs):
        # Create fixed-resolution buffer for both density and sld
        slc = ds.slice(axis, coord)

        # Get domain extent for the two in-plane axes
        ax_map = {'x': 0, 'y': 1, 'z': 2}
        i1 = ax_map[xlabel]
        i2 = ax_map[ylabel]
        lo = [float(ds.domain_left_edge[i1]), float(ds.domain_left_edge[i2])]
        hi = [float(ds.domain_right_edge[i1]), float(ds.domain_right_edge[i2])]

        res = 512  # pixels per side

        # Extract density and sld as fixed-resolution buffers
        frb = yt.FixedResolutionBuffer(slc, (lo[0], hi[0], lo[1], hi[1]), (res, res))
        rho = np.array(frb[("boxlib", "Density")])
        sld = np.array(frb[("boxlib", "sld")])

        # Create RGBA image: density field with solid masked white
        norm = mcolors.LogNorm(vmin=RHO_MIN, vmax=RHO_MAX)
        rho_clipped = np.clip(rho, RHO_MIN, RHO_MAX)
        rgba = cmap_base(norm(rho_clipped))

        # Mask solid regions (sld != 0) as white
        solid_mask = sld > 0.5
        rgba[solid_mask] = [1.0, 1.0, 1.0, 1.0]  # white

        ax.imshow(rgba, origin='lower', extent=[lo[0], hi[0], lo[1], hi[1]],
                  aspect='equal', interpolation='bilinear')
        ax.set_xlabel(f"{xlabel} (m)", fontsize=11)
        ax.set_ylabel(f"{ylabel} (m)", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.tick_params(labelsize=9)

    # Add a shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_base, norm=mcolors.LogNorm(vmin=RHO_MIN, vmax=RHO_MAX))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02, shrink=0.85)
    cbar.set_label("Density (kg/m³)", fontsize=11)

    plt.tight_layout(rect=[0, 0, 0.95, 0.94])
    frame_path = os.path.join(frames_dir, f"frame_{idx:04d}.png")
    fig.savefig(frame_path, dpi=DPI, facecolor='white')
    plt.close(fig)

    print(f"  [{idx+1}/{len(plt_dirs)}] Step {step}, t={time_val:.4e}")

# ============================================================================
# Assemble movie
# ============================================================================
print(f"\nAssembling movie ({FPS} fps)...")
frame_files = sorted(glob.glob(os.path.join(frames_dir, "frame_*.png")))

if frame_files:
    # Use glob pattern for ffmpeg
    os.system(
        f"ffmpeg -y -framerate {FPS} -i {frames_dir}/frame_%04d.png "
        f"-vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' "
        f"-c:v libx264 -pix_fmt yuv420p -crf 18 {OUTPUT} 2>/dev/null"
    )

    if os.path.exists(OUTPUT):
        size_mb = os.path.getsize(OUTPUT) / 1e6
        print(f"Done: {OUTPUT} ({size_mb:.1f} MB, {len(frame_files)} frames, {FPS} fps)")
        print(f"Duration: {len(frame_files)/FPS:.1f} s")
    else:
        print("ffmpeg failed!")
else:
    print("No frames generated!")
