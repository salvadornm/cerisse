#!/usr/bin/env python3
"""
Plot Cerisse surf output (.vtp) files.

Examples:
    python3 tools/plot_surf.py wrk/3Dsphere/surf/geom0/surf0_00050.vtp
    python3 tools/plot_surf.py "wrk/3Dsphere/surf/geom0/surf0_*.vtp" --field Pressure
    python3 tools/plot_surf.py wrk/3Dsphere/surf/geom0 --latest --save surf_latest.png --offscreen
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Iterable, List

import pyvista as pv


SUPPORTED_SUFFIXES = (".vtp", ".vtk", ".stl", ".ply")


def natural_key(path: str) -> list[object]:
    base = os.path.basename(path)
    key: list[object] = []
    num = ""
    for ch in base:
        if ch.isdigit():
            num += ch
        else:
            if num:
                key.append(int(num))
                num = ""
            key.append(ch)
    if num:
        key.append(int(num))
    return key


def expand_inputs(inputs: Iterable[str]) -> List[str]:
    files: list[str] = []

    for item in inputs:
        if os.path.isdir(item):
            dir_matches: list[str] = []
            for suffix in SUPPORTED_SUFFIXES:
                dir_matches.extend(glob.glob(os.path.join(item, f"*{suffix}")))
            files.extend(dir_matches)
            continue

        if any(ch in item for ch in "*?[]"):
            files.extend(glob.glob(item))
            continue

        if os.path.isfile(item):
            files.append(item)
        else:
            print(f"Warning: '{item}' does not exist, skipping.", file=sys.stderr)

    # Keep only supported surface files and remove duplicates while preserving order.
    dedup: list[str] = []
    seen = set()
    for path in files:
        if not path.lower().endswith(SUPPORTED_SUFFIXES):
            continue
        apath = os.path.abspath(path)
        if apath not in seen:
            seen.add(apath)
            dedup.append(path)

    dedup.sort(key=natural_key)
    return dedup


def get_field_names(mesh: pv.DataSet) -> list[str]:
    names = list(mesh.cell_data.keys()) + list(mesh.point_data.keys())
    return [name.strip() for name in names]


def pick_field(mesh: pv.DataSet, requested: str | None) -> str | None:
    cell_keys = [k for k in mesh.cell_data.keys()]
    point_keys = [k for k in mesh.point_data.keys()]

    if requested is not None:
        for key in cell_keys + point_keys:
            if key.strip() == requested:
                return key
        raise ValueError(
            f"Field '{requested}' not found. Available fields: {get_field_names(mesh)}"
        )

    if cell_keys:
        return cell_keys[0]
    if point_keys:
        return point_keys[0]
    return None


def read_meshes(paths: Iterable[str]) -> list[tuple[str, pv.DataSet]]:
    meshes: list[tuple[str, pv.DataSet]] = []
    for path in paths:
        try:
            mesh = pv.read(path)
            meshes.append((path, mesh))
        except Exception as exc:  # pragma: no cover
            print(f"Error reading '{path}': {exc}", file=sys.stderr)
    return meshes


def print_info(path: str, mesh: pv.DataSet) -> None:
    print(f"Loaded: {path}")
    print(f"  Points: {mesh.n_points}")
    print(f"  Cells : {mesh.n_cells}")

    fields = get_field_names(mesh)
    if fields:
        print(f"  Fields: {', '.join(fields)}")
    else:
        print("  Fields: <none>")


def plot_meshes(
    meshes: list[tuple[str, pv.DataSet]],
    field: str | None,
    cmap: str,
    show_edges: bool,
    offscreen: bool,
    screenshot: str | None,
) -> None:
    if not meshes:
        raise ValueError("No valid surface files to plot.")

    n = len(meshes)
    if n == 1:
        shape = (1, 1)
    else:
        ncols = int((n**0.5) + 0.999999)
        nrows = (n + ncols - 1) // ncols
        shape = (nrows, ncols)

    plotter = pv.Plotter(shape=shape, off_screen=offscreen)

    for i, (path, mesh) in enumerate(meshes):
        row, col = divmod(i, shape[1])
        plotter.subplot(row, col)

        chosen = pick_field(mesh, field)
        label = os.path.basename(path)
        if chosen is not None:
            label += f"\n{chosen.strip()}"

        plotter.add_text(label, font_size=11)
        plotter.add_mesh(
            mesh,
            scalars=chosen,
            cmap=cmap,
            show_edges=show_edges,
            scalar_bar_args={"title": chosen.strip() if chosen else ""},
        )
        plotter.show_grid()

    if n > 1:
        plotter.link_views()

    if screenshot:
        plotter.show(screenshot=screenshot, auto_close=True)
        print(f"Saved screenshot: {screenshot}")
    else:
        plotter.show(auto_close=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Read and plot Cerisse surf files.")
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Surface files, glob patterns, or directories containing .vtp/.vtk/.stl/.ply",
    )
    parser.add_argument("--field", default=None, help="Scalar field to plot (default: first field)")
    parser.add_argument("--latest", action="store_true", help="Plot only the latest file after sorting")
    parser.add_argument("--cmap", default="jet", help="Matplotlib colormap name (default: jet)")
    parser.add_argument("--show-edges", action="store_true", help="Show mesh edges")
    parser.add_argument("--offscreen", action="store_true", help="Render without opening a GUI window")
    parser.add_argument("--save", default=None, help="Output image path (PNG recommended)")
    parser.add_argument("--list-fields", action="store_true", help="Print mesh fields and exit")

    args = parser.parse_args()

    files = expand_inputs(args.inputs)
    if not files:
        print("No valid surface files found.", file=sys.stderr)
        return 1

    if args.latest:
        files = [files[-1]]

    meshes = read_meshes(files)
    if not meshes:
        print("No mesh could be loaded.", file=sys.stderr)
        return 2

    for path, mesh in meshes:
        print_info(path, mesh)

    if args.list_fields:
        return 0

    try:
        plot_meshes(
            meshes=meshes,
            field=args.field,
            cmap=args.cmap,
            show_edges=args.show_edges,
            offscreen=args.offscreen,
            screenshot=args.save,
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
