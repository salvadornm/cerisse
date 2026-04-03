#!/usr/bin/env python3
"""
Quick 3D flow-field slicer for Cerisse/AMReX plotfiles.

Default behavior:
- Reads one plotfile (or the latest plt* under a directory)
- Creates a slice through the geometric box center
- Uses z-normal slice and field x_velocity by default

Examples:
  /home/qiaoj/testcerisse/cerisse/.venv/bin/python tools/slice_3d_center.py wrk/3dtest/plot/plt05000
  /home/qiaoj/testcerisse/cerisse/.venv/bin/python tools/slice_3d_center.py wrk/3dtest/plot --field Density --axis y
  /home/qiaoj/testcerisse/cerisse/.venv/bin/python tools/slice_3d_center.py wrk/3dtest/plot --list-fields
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import Iterable, List, Optional, Tuple

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


def _pick_default_field(field_list: Iterable[FieldType]) -> Optional[FieldType]:
    fields = list(field_list)
    preferred = [
        ("boxlib", "x_velocity"),
        ("boxlib", "y_velocity"),
        ("boxlib", "z_velocity"),
        ("boxlib", "velocity_magnitude"),
        ("boxlib", "pressure"),
        ("boxlib", "Density"),
    ]
    for f in preferred:
        if f in fields:
            return f
    return fields[0] if fields else None


def _resolve_field(ds: yt.data_objects.static_output.Dataset, requested: Optional[str]) -> FieldType:
    fields = list(ds.field_list)

    if requested is None:
        chosen = _pick_default_field(fields)
        if chosen is None:
            raise ValueError("Dataset has no available fields.")
        return chosen

    # Allow explicit type:name format, e.g. boxlib:x_velocity
    if ":" in requested:
        ftype, fname = requested.split(":", 1)
        cand = (ftype, fname)
        if cand in fields:
            return cand
        raise ValueError(f"Field '{requested}' not found.")

    # Try exact name under boxlib first.
    boxlib_cand = ("boxlib", requested)
    if boxlib_cand in fields:
        return boxlib_cand

    # Fallback: find by field name only.
    matches = [f for f in fields if f[1] == requested]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Field name '{requested}' is ambiguous. Use type:name format. Matches: {matches}"
        )

    raise ValueError(f"Field '{requested}' not found.")


def list_fields(ds: yt.data_objects.static_output.Dataset) -> None:
    print("Available fields:")
    for ftype, fname in ds.field_list:
        print(f"  {ftype}:{fname}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Quickly load 3D flow plotfile and slice at geometric-box center."
    )
    parser.add_argument("input_path", help="Plotfile path (pltXXXX) or directory containing plt*")
    parser.add_argument("--field", default=None, help="Field name (e.g. x_velocity or boxlib:x_velocity)")
    parser.add_argument("--axis", choices=["x", "y", "z"], default="z", help="Slice normal axis")
    parser.add_argument("--log", action="store_true", help="Use log scale for the field")
    parser.add_argument("--cmap", default="RdBu_r", help="Colormap name")
    parser.add_argument("--vmin", type=float, default=None, help="Colorbar min")
    parser.add_argument("--vmax", type=float, default=None, help="Colorbar max")
    parser.add_argument("--list-fields", action="store_true", help="List available fields and exit")
    parser.add_argument("--output", default=None, help="Output image filename prefix")

    args = parser.parse_args()

    try:
        plotfile = resolve_plotfile(args.input_path)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Loading: {plotfile}")
    ds = yt.load(plotfile)

    if args.list_fields:
        list_fields(ds)
        return 0

    try:
        field = _resolve_field(ds, args.field)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        print("Hint: run with --list-fields to inspect available fields.", file=sys.stderr)
        return 2

    print(f"Slicing field: {field[0]}:{field[1]}")
    print(f"Slice axis: {args.axis}")
    print("Slice center: geometric box center (yt center='c')")

    slc = yt.SlicePlot(ds, args.axis, field, center="c")
    slc.set_log(field, args.log)
    slc.set_cmap(field, args.cmap)

    if args.vmin is not None or args.vmax is not None:
        # If one bound is missing, yt will auto-compute that side.
        current_min = args.vmin if args.vmin is not None else None
        current_max = args.vmax if args.vmax is not None else None
        slc.set_zlim(field, current_min, current_max)

    slc.annotate_title(f"Center slice ({args.axis}-normal): {field[1]}")

    if args.output:
        out = slc.save(args.output)
    else:
        out = slc.save()

    if isinstance(out, list):
        print("Saved:")
        for item in out:
            print(f"  {item}")
    else:
        print(f"Saved: {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
