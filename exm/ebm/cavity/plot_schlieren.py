#!/usr/bin/env python3
"""Plot a Fig. 6-style numerical schlieren image for the cavity case.

The script finds the highest-numbered valid ``plot/plt*`` AMReX plotfile,
loads Density and EB volume fraction with yt, computes the signed streamwise
density gradient, and writes a grayscale PNG. Coordinates are normalized by
the cavity depth D and x/D = 0 is the cavity leading edge.

Examples
--------
Run from this directory using the Python environment that provides yt::

    python plot_schlieren.py
    python plot_schlieren.py --plotfile plot/plt02000 --output run2m6.png
    python plot_schlieren.py --percentile 99.0 --dpi 300
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import yt


# Geometry from eb_geometry.cpp.
CAVITY_DEPTH = 1.128e-4  # D [m]
CAVITY_LENGTH = 2.256e-4  # L = 2D [m]
CAVITY_LEADING_EDGE = 3.9 * CAVITY_DEPTH  # x coordinate [m]
RHO_INF = 1.177  # free-stream density [kg/m^3]


def plotfile_step(path: Path) -> int:
    """Return the integer suffix of a plt directory, or -1 if absent."""
    match = re.fullmatch(r"plt(\d+)", path.name)
    return int(match.group(1)) if match else -1


def latest_plotfile(plot_dir: Path) -> Path:
    """Find the highest-numbered plotfile that contains an AMReX Header."""
    candidates = [
        path
        for path in plot_dir.glob("plt[0-9]*")
        if path.is_dir() and (path / "Header").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"No valid plt* directories found in {plot_dir}")
    return max(candidates, key=plotfile_step)


def field_key(ds: yt.data_objects.static_output.Dataset, name: str) -> tuple[str, str]:
    """Resolve an AMReX field while giving a useful error if it is absent."""
    key = ("boxlib", name)
    if key not in ds.field_list:
        available = ", ".join(sorted(field[1] for field in ds.field_list))
        raise KeyError(f"Field {name!r} is unavailable. Plotfile fields: {available}")
    return key


def level_dimensions(ds: object, level: int) -> np.ndarray:
    """Return covering-grid dimensions without refining yt's dummy z-axis."""
    dimensions = np.asarray(ds.domain_dimensions, dtype=int).copy()
    dimensions[: int(ds.dimensionality)] *= int(ds.refine_by) ** level
    return dimensions


def level_coverage(ds: object, level: int, shape: tuple[int, int]) -> np.ndarray:
    """Mark cells actually covered by grids belonging to one AMR level."""
    covered = np.zeros(shape, dtype=bool)
    left = np.asarray(ds.domain_left_edge, dtype=float)
    width = np.asarray(ds.domain_right_edge - ds.domain_left_edge, dtype=float)
    spacing = width[:2] / np.asarray(shape, dtype=float)

    for grid in ds.index.grids:
        if int(grid.Level) != level:
            continue
        grid_left = np.asarray(grid.LeftEdge, dtype=float)
        grid_right = np.asarray(grid.RightEdge, dtype=float)
        lo = np.rint((grid_left[:2] - left[:2]) / spacing).astype(int)
        hi = np.rint((grid_right[:2] - left[:2]) / spacing).astype(int)
        covered[lo[0] : hi[0], lo[1] : hi[1]] = True
    return covered


def load_composite_gradient(
    plotfile: Path,
) -> tuple[object, np.ndarray, np.ndarray]:
    """Compute drho/dx per AMR level, then assemble a finest-pixel image.

    Differentiating a finest-level covering grid directly is incorrect in
    coarse-only regions: each coarse value is replicated into multiple fine
    pixels, producing a false odd-even gradient. Here derivatives are taken at
    each level's native spacing before finer grids overwrite coarser results.
    """
    yt.set_log_level("error")
    ds = yt.load(str(plotfile))
    domain_width = np.asarray(ds.domain_right_edge - ds.domain_left_edge, dtype=float)
    finest_shape = tuple(level_dimensions(ds, int(ds.max_level))[:2])
    gradient = np.zeros(finest_shape, dtype=float)
    volume_fraction = np.zeros(finest_shape, dtype=float)

    for level in range(int(ds.max_level) + 1):
        dimensions = level_dimensions(ds, level)
        grid = ds.covering_grid(
            level=level,
            left_edge=ds.domain_left_edge,
            dims=dimensions,
        )
        density_level = np.asarray(
            grid[field_key(ds, "Density")], dtype=float
        ).squeeze()
        vfrac_level = np.asarray(
            grid[field_key(ds, "vfrac")], dtype=float
        ).squeeze()
        if density_level.ndim != 2 or vfrac_level.shape != density_level.shape:
            raise ValueError(
                "Expected a 2-D Density/vfrac plotfile; got "
                f"Density shape {density_level.shape} and "
                f"vfrac shape {vfrac_level.shape}"
            )

        dx_level = domain_width[0] / density_level.shape[0]
        gradient_level = np.gradient(
            density_level, dx_level, axis=0, edge_order=2
        )
        covered_level = level_coverage(ds, level, density_level.shape)

        upscale = int(ds.refine_by) ** (int(ds.max_level) - level)
        if upscale > 1:
            gradient_level = np.repeat(
                np.repeat(gradient_level, upscale, axis=0), upscale, axis=1
            )
            vfrac_level = np.repeat(
                np.repeat(vfrac_level, upscale, axis=0), upscale, axis=1
            )
            covered_level = np.repeat(
                np.repeat(covered_level, upscale, axis=0), upscale, axis=1
            )

        gradient[covered_level] = gradient_level[covered_level]
        volume_fraction[covered_level] = vfrac_level[covered_level]

    gradient *= CAVITY_DEPTH / RHO_INF
    fluid = volume_fraction > 0.5
    # Do not use centered differences that cross an embedded solid boundary.
    fluid_x_stencil = fluid.copy()
    fluid_x_stencil[1:-1, :] &= fluid[:-2, :] & fluid[2:, :]
    gradient[fluid & ~fluid_x_stencil] = 0.0
    return ds, gradient, fluid


def render(
    ds: object,
    gradient: np.ndarray,
    fluid: np.ndarray,
    plotfile: Path,
    output: Path,
    percentile: float,
    dpi: int,
) -> float:
    """Render signed gradient in grayscale and return the applied magnitude limit."""
    finite_fluid = fluid & np.isfinite(gradient)
    values = np.abs(gradient[finite_fluid])
    values = values[values > 0.0]
    if values.size == 0:
        raise ValueError("No nonzero finite density gradients were found in fluid cells")
    limit = float(np.percentile(values, percentile))
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError(f"Invalid gradient display limit: {limit}")

    image = np.clip(gradient, -limit, limit).T
    solid = (~fluid).T

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    extent = (
        (left[0] - CAVITY_LEADING_EDGE) / CAVITY_DEPTH,
        (right[0] - CAVITY_LEADING_EDGE) / CAVITY_DEPTH,
        left[1] / CAVITY_DEPTH,
        right[1] / CAVITY_DEPTH,
    )

    cmap = plt.colormaps["gray"].copy()
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    rgba = cmap(norm(image))
    rgba[solid] = (0.0, 0.0, 0.0, 1.0)

    fig, ax = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
    ax.imshow(
        rgba,
        origin="lower",
        extent=extent,
        interpolation="bilinear",
        aspect="equal",
    )
    ax.set_xlabel(r"$x/D$ (cavity leading edge at 0)")
    ax.set_ylabel(r"$y/D$")
    ax.set_title(
        rf"Run 2M6, $M_\infty=0.6$: $D\,\partial\rho/(\rho_\infty\partial x)$"
        "\n"
        rf"step {plotfile_step(plotfile)}, $t={float(ds.current_time):.4e}$ s"
    )
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, facecolor="white")
    plt.close(fig)
    return limit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "plot",
        help="directory containing plt* plotfiles (default: cavity/plot)",
    )
    parser.add_argument(
        "--plotfile",
        type=Path,
        help="specific plotfile to use instead of the latest numbered plotfile",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "schlieren_drhodx.png",
        help="output PNG path (default: cavity/schlieren_drhodx.png)",
    )
    parser.add_argument(
        "--percentile",
        type=float,
        default=99.5,
        help="symmetric clipping percentile of |D/rho_inf drho/dx| (default: 99.5)",
    )
    parser.add_argument("--dpi", type=int, default=240, help="output resolution")
    args = parser.parse_args()
    if not 0.0 < args.percentile <= 100.0:
        parser.error("--percentile must be in (0, 100]")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    return args


def main() -> None:
    args = parse_args()
    plotfile = args.plotfile or latest_plotfile(args.plot_dir)
    plotfile = plotfile.resolve()
    if not (plotfile / "Header").is_file():
        raise FileNotFoundError(f"Not a valid AMReX plotfile: {plotfile}")

    ds, gradient, fluid = load_composite_gradient(plotfile)
    limit = render(
        ds,
        gradient,
        fluid,
        plotfile,
        args.output.resolve(),
        args.percentile,
        args.dpi,
    )
    print(f"Loaded: {plotfile}")
    print(f"Wrote:  {args.output.resolve()}")
    print(f"Symmetric gradient limit: +/-{limit:.6g}")


if __name__ == "__main__":
    main()
