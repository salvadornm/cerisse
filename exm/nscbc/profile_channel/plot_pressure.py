#!/usr/bin/env python3
"""Check the streamwise pressure gradient in the profiled NSCBC channel.

The script extracts pressure along the channel centreline from an AMReX
plotfile, compares it with fully developed plane-Poiseuille flow, and checks
whether dp/dx is constant away from the inlet and outlet.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import yt
except ImportError as exc:  # pragma: no cover - depends on the user environment
    missing = exc.name or "a Python dependency"
    raise SystemExit(
        f"Missing {missing!r}; install the analysis dependencies with "
        "python3 -m pip install numpy matplotlib yt"
    ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "plotfile",
        nargs="?",
        type=Path,
        help="AMReX plotfile (default: latest plot/plt*)",
    )
    parser.add_argument("--output", type=Path, default=Path("pressure_over_x.png"))
    parser.add_argument("--y", type=float, default=0.0, help="lineout y coordinate [m]")
    parser.add_argument("--temperature", type=float, default=300.0, help="reference T [K]")
    parser.add_argument("--umax", type=float, default=None, help="target centreline velocity [m/s]")
    parser.add_argument("--half-height", type=float, default=1.0e-3, help="channel h [m]")
    parser.add_argument(
        "--trim", type=float, default=0.1,
        help="fraction removed from each x end for the constancy check",
    )
    parser.add_argument("--min-r2", type=float, default=0.99)
    parser.add_argument(
        "--max-gradient-cv", type=float, default=0.10,
        help="maximum std(dp/dx)/abs(mean(dp/dx))",
    )
    parser.add_argument(
        "--max-slope-error", type=float, default=0.20,
        help="maximum relative error against the ideal pressure slope",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="return a nonzero status when any check fails",
    )
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def latest_plotfile() -> Path:
    candidates = [p for p in Path("plot").glob("plt*") if (p / "Header").is_file()]
    if not candidates:
        raise SystemExit("No AMReX plotfiles found under plot/plt*")

    def step(path: Path) -> int:
        match = re.search(r"(\d+)$", path.name)
        return int(match.group(1)) if match else -1

    return max(candidates, key=step)


def as_numpy(values) -> np.ndarray:
    return np.asarray(values.to_value() if hasattr(values, "to_value") else values, dtype=float)


def centreline_pressure(plotfile: Path, y: float) -> tuple[np.ndarray, np.ndarray]:
    ds = yt.load(str(plotfile))
    # yt represents a 2-D AMReX dataset with a singleton z direction.  For an
    # x-directed ray the fixed coordinates are therefore (y, z).
    ray = ds.ortho_ray(0, (y, 0.0))
    x = as_numpy(ray["x"])
    pressure = as_numpy(ray["pressure"])
    order = np.argsort(x)
    return x[order], pressure[order]


def sutherland_viscosity(temperature: float) -> float:
    # Must match visc_suth_t in src/cls/TransportCeris.h.
    return 1.458e-6 * temperature**1.5 / (110.4 + temperature)


def main() -> int:
    args = parse_args()
    if not 0.0 <= args.trim < 0.5:
        raise SystemExit("--trim must satisfy 0 <= trim < 0.5")

    plotfile = args.plotfile or latest_plotfile()
    x, pressure = centreline_pressure(plotfile, args.y)
    if x.size < 5:
        raise SystemExit(f"Only {x.size} centreline samples were found")

    gamma = 1.4
    gas_constant = 287.0
    umax = args.umax
    if umax is None:
        umax = 0.1*np.sqrt(gamma*gas_constant*args.temperature)

    mu = sutherland_viscosity(args.temperature)
    # For u(y)=umax[1-(y/h)^2], the x-momentum balance gives
    # dp/dx = mu*d2u/dy2 = -2*mu*umax/h^2.
    ideal_slope = -2.0*mu*umax/args.half_height**2

    length = x[-1] - x[0]
    mask = (x >= x[0] + args.trim*length) & (x <= x[-1] - args.trim*length)
    if np.count_nonzero(mask) < 3:
        raise SystemExit("Too few samples remain after applying --trim")

    x_check = x[mask]
    p_check = pressure[mask]
    measured_slope, intercept = np.polyfit(x_check, p_check, 1)
    p_fit = measured_slope*x_check + intercept
    residual = p_check - p_fit
    ss_res = np.sum(residual**2)
    ss_tot = np.sum((p_check - np.mean(p_check))**2)
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 0.0 else float("nan")

    # np.gradient uses one-sided stencils at the two ends of the selected
    # interval.  Those estimates are dominated by the inlet/outlet transition
    # and are not representative of the developed channel.  Keep only the
    # centred interior derivatives for the constancy check and plot.
    x_gradient = x_check[1:-1]
    dpdx = np.gradient(p_check, x_check, edge_order=2)[1:-1]
    mean_dpdx = float(np.mean(dpdx))
    gradient_cv = float(np.std(dpdx)/abs(mean_dpdx)) if mean_dpdx != 0.0 else float("inf")
    slope_error = abs((measured_slope - ideal_slope)/ideal_slope)

    # Anchor the ideal pressure at the downstream end; only its gradient is
    # prescribed by the analytical solution.
    ideal_pressure = pressure[-1] + ideal_slope*(x - x[-1])

    checks = {
        f"R^2 >= {args.min_r2:g}": bool(np.isfinite(r2) and r2 >= args.min_r2),
        f"gradient CV <= {args.max_gradient_cv:g}": gradient_cv <= args.max_gradient_cv,
        f"slope error <= {args.max_slope_error:g}": slope_error <= args.max_slope_error,
    }
    passed = all(checks.values())

    print(f"plotfile              : {plotfile}")
    print(f"samples/check samples  : {x.size}/{np.count_nonzero(mask)}")
    print(f"Sutherland viscosity   : {mu:.8e} Pa s")
    print(f"ideal dp/dx             : {ideal_slope:.8e} Pa/m")
    print(f"measured dp/dx          : {measured_slope:.8e} Pa/m")
    print(f"relative slope error    : {slope_error:.6g}")
    print(f"linear-fit R^2          : {r2:.8f}")
    print(f"local-gradient CV       : {gradient_cv:.6g}")
    for label, ok in checks.items():
        print(f"{'PASS' if ok else 'FAIL'}: {label}")
    print(f"overall                 : {'PASS' if passed else 'FAIL'}")

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)
    axes[0].plot(x, pressure, "o", ms=3, label="Cerisse centreline")
    axes[0].plot(x, ideal_pressure, "-", lw=1.8, label="Ideal Poiseuille")
    axes[0].set_ylabel("pressure [Pa]")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].plot(x_gradient, dpdx, "o-", ms=3, label="local dp/dx")
    axes[1].axhline(measured_slope, color="C1", ls="--", label="linear-fit slope")
    axes[1].axhline(ideal_slope, color="k", ls=":", label="ideal slope")
    axes[1].set_xlabel("x [m]")
    axes[1].set_ylabel("dp/dx [Pa/m]")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.suptitle(f"NSCBC channel pressure check: {'PASS' if passed else 'FAIL'}")
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)
    print(f"figure                  : {args.output}")
    if args.show:
        plt.show()
    plt.close(fig)

    return int(args.strict and not passed)


if __name__ == "__main__":
    sys.exit(main())
