#!/usr/bin/env python3
"""
Generate a single 2x2 figure for 3D flow data:
- Front view: center slice normal to x
- Side view : center slice normal to y
- Top view  : center slice normal to z
- 3D panel  : volume rendering (or fallback)

Examples:
  /home/qiaoj/testcerisse/cerisse/.venv/bin/python tools/four_view_3d.py wrk/3Dsphere/plot
  /home/qiaoj/testcerisse/cerisse/.venv/bin/python tools/four_view_3d.py wrk/3Dsphere/plot --field pressure --vr-field Density --output /tmp/four_view.png
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import tempfile
from typing import Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yt


yt.funcs.mylog.setLevel("ERROR")

FieldType = Tuple[str, str]


def _natural_key(text: str) -> List[object]:
    parts = re.split(r"(\d+)", os.path.basename(text))
    key: List[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part)
    return key


def _is_plotfile_dir(path: str) -> bool:
    return os.path.isdir(path) and os.path.exists(os.path.join(path, "Header"))


def resolve_plotfile(path: str) -> str:
    if _is_plotfile_dir(path):
        return path

    if os.path.isdir(path):
        candidates = [p for p in glob.glob(os.path.join(path, "plt*")) if _is_plotfile_dir(p)]
        if not candidates:
            raise FileNotFoundError(f"No valid plt* directories found in '{path}'")
        candidates.sort(key=_natural_key)
        return candidates[-1]

    raise FileNotFoundError(f"Path '{path}' is not a valid plotfile directory")


def _resolve_field(ds: yt.data_objects.static_output.Dataset, requested: Optional[str]) -> FieldType:
    fields = list(ds.field_list)

    if requested is None:
        preferred = [
            ("boxlib", "x_velocity"),
            ("boxlib", "y_velocity"),
            ("boxlib", "z_velocity"),
            ("boxlib", "pressure"),
            ("boxlib", "Density"),
        ]
        for field in preferred:
            if field in fields:
                return field
        if fields:
            return fields[0]
        raise ValueError("Dataset has no available fields.")

    if ":" in requested:
        ftype, fname = requested.split(":", 1)
        field = (ftype, fname)
        if field in fields:
            return field
        raise ValueError(f"Field '{requested}' not found.")

    direct = ("boxlib", requested)
    if direct in fields:
        return direct

    matches = [f for f in fields if f[1] == requested]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"Field '{requested}' is ambiguous: {matches}. Use type:name format.")
    raise ValueError(f"Field '{requested}' not found.")


def _add_mach_field(ds: yt.data_objects.static_output.Dataset) -> FieldType:
    """Register a Mach number derived field: |v| / c_s."""
    mach_field = ("gas", "Mach")

    def _mach(field, data):
        vx = np.asarray(data[("boxlib", "x_velocity")], dtype=np.float64)
        vy = np.asarray(data[("boxlib", "y_velocity")], dtype=np.float64)
        vz = np.asarray(data[("boxlib", "z_velocity")], dtype=np.float64)
        p = np.asarray(data[("boxlib", "pressure")], dtype=np.float64)
        rho = np.asarray(data[("boxlib", "Density")], dtype=np.float64)
        gamma = 1.4
        cs = np.sqrt(gamma * np.abs(p) / np.where(rho > 0, rho, 1.0e-30))
        speed = np.sqrt(vx**2 + vy**2 + vz**2)
        return speed / cs

    ds.add_field(mach_field, function=_mach, sampling_type="cell",
                 units="", force_override=True)
    return mach_field


def _add_sld_masked_field(ds: yt.data_objects.static_output.Dataset, base_field: FieldType) -> FieldType:
    """Mask solid cells (sld != 0) with NaN so they appear as the nan_color (white)."""
    if base_field == ("boxlib", "sld"):
        return base_field

    masked_name = f"maskedNaN_{base_field[0]}_{base_field[1]}"
    masked_field = ("gas", masked_name)

    def _masked(field, data):
        raw = np.asarray(data[base_field], dtype=np.float64)
        sld = np.asarray(data[("boxlib", "sld")], dtype=np.float64)
        out = raw.copy()
        out[~np.isclose(sld, 0.0, atol=1.0e-14)] = np.nan
        return out

    units = str(ds.field_info[base_field].units)
    ds.add_field(
        masked_field,
        function=_masked,
        sampling_type="cell",
        units=units,
        force_override=True,
    )
    return masked_field


def _save_center_slice(
    ds: yt.data_objects.static_output.Dataset,
    axis: str,
    field: FieldType,
    out_prefix: str,
    cmap: str,
    log_scale: bool,
    buff_size: int,
    vmin: Optional[float],
    vmax: Optional[float],
) -> str:
    slc = yt.SlicePlot(ds, axis, field, center="c")
    slc.set_buff_size((buff_size, buff_size))
    slc.set_cmap(field, cmap)
    slc.set_log(field, log_scale)
    if (vmin is not None) and (vmax is not None):
        slc.set_zlim(field, vmin, vmax)
    # Make NaN cells (solid bodies) white via matplotlib's "bad" color.
    slc.set_background_color(field, color="white")
    slc.annotate_title(f"Center slice ({axis})")
    out_files = slc.save(out_prefix)
    return out_files[0] if isinstance(out_files, list) else out_files


def _save_volume_render(
    ds: yt.data_objects.static_output.Dataset,
    field: FieldType,
    out_path: str,
    log_scale: bool,
) -> str:
    scene = yt.create_scene(ds, field=field)

    source = scene[0]
    source.set_field(field)
    source.set_log(log_scale)

    cam = scene.camera
    center = ds.domain_center
    width = ds.domain_width
    cam.focus = center
    cam.position = center + 1.8 * width
    cam.north_vector = [0.0, 0.0, 1.0]
    cam.width = 1.2 * width

    scene.save(out_path, sigma_clip=4.0)
    return out_path


def _save_qcriterion_isosurface(
    ds: yt.data_objects.static_output.Dataset,
    out_path: str,
    level: int,
    stride: int,
    iso_value: Optional[float],
    iso_quantile: float,
    cmap: str,
    q_vmin: Optional[float],
    q_vmax: Optional[float],
    q_bar_title: str,
    auto_tune: bool,
) -> Tuple[str, float, int, int]:
    if stride < 1:
        raise ValueError("--q-stride must be >= 1")

    max_level = int(ds.index.max_level)
    if level < 0 or level > max_level:
        raise ValueError(f"--q-level must be in [0, {max_level}]")

    left = np.array(ds.domain_left_edge.d, dtype=np.float64)
    right = np.array(ds.domain_right_edge.d, dtype=np.float64)

    base_dims = np.array(ds.domain_dimensions, dtype=np.int64)
    dims = base_dims * (2**level)

    cg = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)

    u = np.asarray(cg[("boxlib", "x_velocity")], dtype=np.float64)
    v = np.asarray(cg[("boxlib", "y_velocity")], dtype=np.float64)
    w = np.asarray(cg[("boxlib", "z_velocity")], dtype=np.float64)
    sld = np.asarray(cg[("boxlib", "sld")], dtype=np.float64)

    # Enforce non-fluid region to zero as requested.
    fluid = np.isclose(sld, 0.0, atol=1.0e-14)
    u = np.where(fluid, u, 0.0)
    v = np.where(fluid, v, 0.0)
    w = np.where(fluid, w, 0.0)

    if stride > 1:
        u = u[::stride, ::stride, ::stride]
        v = v[::stride, ::stride, ::stride]
        w = w[::stride, ::stride, ::stride]

    spacing = (right - left) / dims
    spacing = spacing * stride

    dudx, dudy, dudz = np.gradient(u, spacing[0], spacing[1], spacing[2], edge_order=1)
    dvdx, dvdy, dvdz = np.gradient(v, spacing[0], spacing[1], spacing[2], edge_order=1)
    dwdx, dwdy, dwdz = np.gradient(w, spacing[0], spacing[1], spacing[2], edge_order=1)

    s_xx = dudx
    s_yy = dvdy
    s_zz = dwdz
    s_xy = 0.5 * (dudy + dvdx)
    s_xz = 0.5 * (dudz + dwdx)
    s_yz = 0.5 * (dvdz + dwdy)

    o_xy = 0.5 * (dudy - dvdx)
    o_xz = 0.5 * (dudz - dwdx)
    o_yz = 0.5 * (dvdz - dwdy)

    s_norm2 = s_xx**2 + s_yy**2 + s_zz**2 + 2.0 * (s_xy**2 + s_xz**2 + s_yz**2)
    o_norm2 = 2.0 * (o_xy**2 + o_xz**2 + o_yz**2)
    q = 0.5 * (o_norm2 - s_norm2)

    q_positive = q[q > 0.0]
    if q_positive.size == 0:
        raise ValueError("Q-criterion has no positive region in this dataset at current sampling.")

    # Scan quantiles for diagnostics and optional auto tuning.
    q50 = float(np.quantile(q_positive, 0.50))
    q70 = float(np.quantile(q_positive, 0.70))
    q80 = float(np.quantile(q_positive, 0.80))
    q85 = float(np.quantile(q_positive, 0.85))
    q90 = float(np.quantile(q_positive, 0.90))
    q95 = float(np.quantile(q_positive, 0.95))
    q98 = float(np.quantile(q_positive, 0.98))

    if auto_tune:
        # Auto mode targets clear structures: moderate iso, compressed color range.
        if iso_value is None:
            iso = q85
        else:
            iso = float(iso_value)

        if (q_vmin is None) and (q_vmax is None):
            q_vmin = q70
            q_vmax = q98
    elif iso_value is None:
        iso = float(np.quantile(q_positive, iso_quantile))
        quantile_candidates = [
            iso_quantile,
            max(iso_quantile - 0.05, 0.50),
            max(iso_quantile - 0.10, 0.40),
            max(iso_quantile - 0.20, 0.30),
        ]
    else:
        iso = float(iso_value)
        quantile_candidates = []

    print(
        "Q+ quantiles: "
        f"q50={q50:.3e}, q70={q70:.3e}, q80={q80:.3e}, q85={q85:.3e}, "
        f"q90={q90:.3e}, q95={q95:.3e}, q98={q98:.3e}"
    )
    if auto_tune:
        print(
            f"Auto tuned : q_iso={iso:.3e}, q_vmin={q_vmin:.3e}, q_vmax={q_vmax:.3e}"
        )
        print(
            "Suggestion : reuse these values with --q-iso/--q-vmin/--q-vmax for consistent comparison."
        )

    q_shape = np.array(q.shape, dtype=np.int64)
    grid = pv.ImageData()
    grid.dimensions = tuple((q_shape + 1).tolist())
    grid.origin = tuple(left.tolist())
    grid.spacing = tuple(spacing.tolist())
    grid.cell_data["Q"] = np.ascontiguousarray(q.ravel(order="F"))

    # Contour filter in VTK works on point data; convert from cell-centered Q first.
    point_grid = grid.cell_data_to_point_data()
    contour = point_grid.contour(isosurfaces=[iso], scalars="Q")
    if contour.n_cells == 0 and iso_value is None:
        for qnt in quantile_candidates[1:]:
            iso_try = float(np.quantile(q_positive, qnt))
            contour = point_grid.contour(isosurfaces=[iso_try], scalars="Q")
            if contour.n_cells > 0:
                iso = iso_try
                break

    if contour.n_cells == 0:
        raise ValueError(
            "Q-criterion contour is empty. Try smaller --q-iso or lower --q-quantile (e.g. 0.80)."
        )

    plotter = pv.Plotter(off_screen=True, window_size=(1200, 900))
    plotter.set_background("white")
    mesh_kwargs = {
        "scalars": "Q",
        "cmap": cmap,
        "smooth_shading": True,
        "opacity": 1.0,
        "scalar_bar_args": {
            "title": q_bar_title,
            "fmt": "%.2e",
            "label_font_size": 10,
            "title_font_size": 11,
        },
    }
    if (q_vmin is not None) and (q_vmax is not None):
        mesh_kwargs["clim"] = [q_vmin, q_vmax]

    plotter.add_mesh(contour, **mesh_kwargs)
    plotter.add_axes()
    plotter.show_grid()
    plotter.camera_position = "iso"
    plotter.show(screenshot=out_path, auto_close=True)

    cells_used = int(np.prod(q_shape))
    return out_path, iso, level, cells_used


def _compose_2x2(image_paths: List[str], titles: List[str], output: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), dpi=150)

    for ax, img_path, title in zip(axes.ravel(), image_paths, titles):
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.set_title(title, fontsize=12)
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate front/side/top + 3D in one 2x2 figure.")
    parser.add_argument("input_path", help="Plotfile path (pltXXXX) or directory containing plt*")
    parser.add_argument("--field", default=None, help="Field for center slices (default: x_velocity if available)")
    parser.add_argument("--vr-field", default=None, help="Field for 3D volume render (default: Density if available)")
    parser.add_argument("--cmap", default="RdBu_r", help="Colormap for slices")
    parser.add_argument("--log", action="store_true", help="Use log scale for slices")
    parser.add_argument("--slice-res", type=int, default=1400, help="Pixel resolution per slice panel")
    parser.add_argument("--slice-vmin", type=float, default=None, help="Slice colorbar lower limit")
    parser.add_argument("--slice-vmax", type=float, default=None, help="Slice colorbar upper limit")
    parser.add_argument(
        "--panel4",
        choices=["qcriterion", "volume"],
        default="qcriterion",
        help="Fourth panel type (default: qcriterion)",
    )
    parser.add_argument("--q-level", type=int, default=-1, help="AMR level for Q computation (default: max_level-1)")
    parser.add_argument("--q-stride", type=int, default=2, help="Subsampling stride for Q grid (default: 2)")
    parser.add_argument("--q-iso", type=float, default=None, help="Explicit Q isovalue")
    parser.add_argument("--q-quantile", type=float, default=0.92, help="Quantile of positive Q if --q-iso not set")
    parser.add_argument("--q-cmap", default="plasma", help="Colormap for Q-criterion isosurface")
    parser.add_argument("--q-vmin", type=float, default=None, help="Q panel colorbar lower limit")
    parser.add_argument("--q-vmax", type=float, default=None, help="Q panel colorbar upper limit")
    parser.add_argument("--q-bar-title", default="Q-criterion", help="Q panel scalar bar title")
    parser.add_argument(
        "--q-auto-tune",
        action="store_true",
        help="Auto scan Q quantiles and choose q-iso/q-vmin/q-vmax for clearer contrast",
    )
    parser.add_argument("--vr-log", action="store_true", help="Use log scale for volume rendering")
    parser.add_argument("--output", default="four_view.png", help="Output merged figure filename")

    args = parser.parse_args()

    if (args.slice_vmin is None) ^ (args.slice_vmax is None):
        print("Error: --slice-vmin and --slice-vmax must be set together.", file=sys.stderr)
        return 3

    if (args.q_vmin is None) ^ (args.q_vmax is None):
        print("Error: --q-vmin and --q-vmax must be set together.", file=sys.stderr)
        return 4

    try:
        plotfile = resolve_plotfile(args.input_path)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Loading: {plotfile}")
    ds = yt.load(plotfile)

    try:
        # Handle Mach as a special derived field before resolve.
        if args.field and args.field.lower() == "mach":
            slice_field = _add_mach_field(ds)
        else:
            slice_field = _resolve_field(ds, args.field)
        vr_field = _resolve_field(ds, args.vr_field if args.vr_field else "Density")
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    # Apply sld!=0 -> NaN masking so solid regions appear white.
    slice_field_masked = _add_sld_masked_field(ds, slice_field)
    vr_field_masked = _add_sld_masked_field(ds, vr_field)

    print(f"Slice field: {slice_field[0]}:{slice_field[1]} (masked by sld!=0 -> NaN/white)")

    with tempfile.TemporaryDirectory(prefix="fourview_") as tmpdir:
        front_path = _save_center_slice(ds, "x", slice_field_masked, os.path.join(tmpdir, "front"), args.cmap, args.log, args.slice_res, args.slice_vmin, args.slice_vmax)
        side_path = _save_center_slice(ds, "y", slice_field_masked, os.path.join(tmpdir, "side"), args.cmap, args.log, args.slice_res, args.slice_vmin, args.slice_vmax)
        top_path = _save_center_slice(ds, "z", slice_field_masked, os.path.join(tmpdir, "top"), args.cmap, args.log, args.slice_res, args.slice_vmin, args.slice_vmax)

        panel4_path = os.path.join(tmpdir, "panel4.png")
        panel4_title = ""
        if args.panel4 == "qcriterion":
            auto_q_level = max(int(ds.index.max_level) - 1, 0)
            q_level = auto_q_level if args.q_level < 0 else args.q_level
            q_path, q_iso, used_level, used_cells = _save_qcriterion_isosurface(
                ds=ds,
                out_path=panel4_path,
                level=q_level,
                stride=args.q_stride,
                iso_value=args.q_iso,
                iso_quantile=args.q_quantile,
                cmap=args.q_cmap,
                q_vmin=args.q_vmin,
                q_vmax=args.q_vmax,
                q_bar_title=args.q_bar_title,
                auto_tune=args.q_auto_tune,
            )
            panel4_path = q_path
            panel4_title = f"Q-criterion Isosurface (Q={q_iso:.3e}, L{used_level}, N={used_cells})"
            print(f"Panel4     : Q-criterion isosurface, Q={q_iso:.6e}, level={used_level}, cells={used_cells}")
        else:
            print(f"3D field   : {vr_field[0]}:{vr_field[1]} (masked by sld!=0 -> 0)")
            _save_volume_render(ds, vr_field_masked, panel4_path, args.vr_log)
            panel4_title = f"3D Volume Render ({vr_field[1]})"

        titles = [
            "Front View (x-normal center slice)",
            "Side View (y-normal center slice)",
            "Top View (z-normal center slice)",
            panel4_title,
        ]
        _compose_2x2([front_path, side_path, top_path, panel4_path], titles, args.output)

    print(f"Saved merged figure: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
