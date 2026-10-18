#!/usr/bin/env python3
"""Compare matched small/wider probes (Eq. 97) and optionally plotfiles (Eq. 91)."""
import argparse
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TREF = 2.8e-4

def read_probes(path):
    # CERISSE headers contain unquoted commas inside integer index tuples.
    # Only numeric rows are CSV; column order is fixed by cns.time_probes.
    a = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if a.shape[1] != 4 or len(a) < 1 or not np.isfinite(a).all():
        raise ValueError(f"{path}: expected finite time,L1,L2,L3 rows")
    if np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError(f"{path}: times must increase (use fresh, non-restarted runs)")
    return a


def local_errors(small, wider):
    if small.shape != wider.shape or not np.allclose(small[:, 0], wider[:, 0], rtol=0, atol=1e-14):
        raise ValueError("Probe times differ; use identical fixed dt and step count. No interpolation is applied.")
    if np.any(wider[:, 1:] <= 0):
        raise ValueError("Reference pressures must be positive")
    return np.abs(small[:, 1:] - wider[:, 1:]) / wider[:, 1:]


def plotfiles(folder):
    files = {}
    for p in Path(folder).glob("plt*"):
        m = re.fullmatch(r"plt(\d+)", p.name)
        if m and (p / "Header").is_file():
            files[int(m[1])] = p
    if not files:
        raise ValueError(f"No plotfiles in {folder}")
    return files


def global_errors(small_folder, wider_folder):
    import yt  # Optional dependency: only needed for the top panel of Fig. 15.
    small, wider = plotfiles(small_folder), plotfiles(wider_folder)
    if small.keys() != wider.keys():
        raise ValueError("Small/wider plotfile steps differ")
    result = []
    for step in sorted(small):
        a, b = yt.load(str(small[step])), yt.load(str(wider[step]))
        if a.index.max_level != 0 or b.index.max_level != 0:
            raise ValueError("This comparison requires the supplied uniform level-0 grids")
        dims = np.asarray(a.domain_dimensions, dtype=int)
        da = np.asarray(a.domain_width) / dims
        db = np.asarray(b.domain_width) / np.asarray(b.domain_dimensions)
        if not np.allclose(da, db, rtol=1e-12, atol=0):
            raise ValueError("Grid spacings differ: double the wider cell count")
        offset = (np.asarray(a.domain_left_edge) - np.asarray(b.domain_left_edge)) / db
        if not np.allclose(offset, np.round(offset), rtol=0, atol=1e-9):
            raise ValueError("Grid cell centers are not aligned")
        if np.any(a.domain_left_edge < b.domain_left_edge) or np.any(a.domain_right_edge > b.domain_right_edge):
            raise ValueError("Reference does not contain the small domain")
        t = float(a.current_time)
        if abs(t-float(b.current_time)) > 1e-14:
            raise ValueError("Plotfile times differ")
        ga = a.covering_grid(0, a.domain_left_edge, dims)
        gb = b.covering_grid(0, a.domain_left_edge, dims)
        def pressure(ds, grid):
            fields = [f for f in ds.field_list if f[1] == 'pressure']
            if len(fields) != 1:
                raise ValueError("Expected exactly one saved pressure field")
            p = np.asarray(grid[fields[0]], dtype=float)
            if not np.isfinite(p).all() or np.any(p <= 0):
                raise ValueError("Nonfinite/nonpositive pressure in plotfile")
            return p
        p, pref = pressure(a, ga), pressure(b, gb)
        result.append((t, np.linalg.norm((p-pref).ravel()) / np.linalg.norm(pref.ravel())))
    return np.asarray(result)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('small', type=Path)
    ap.add_argument('wider', type=Path)
    ap.add_argument('--small-plots', type=Path)
    ap.add_argument('--wider-plots', type=Path)
    ap.add_argument('--output', type=Path, default=Path('pressure_error.png'))
    ap.add_argument('--label', default='CERISSE NSCBC')
    args = ap.parse_args()
    if bool(args.small_plots) != bool(args.wider_plots):
        ap.error('supply both --small-plots and --wider-plots')
    a, b = read_probes(args.small), read_probes(args.wider)
    errors = local_errors(a,b)
    global_data = global_errors(args.small_plots,args.wider_plots) if args.small_plots else None
    rows = 4 if global_data is not None else 3
    fig, axes = plt.subplots(rows, 1, sharex=True, figsize=(6, 2*rows), constrained_layout=True)
    offset = 0
    if global_data is not None:
        axes[0].plot(global_data[:,0]/TREF, global_data[:,1], label=args.label)
        axes[0].set_ylabel(r'$\epsilon$ (Eq. 91)')
        np.savetxt(args.output.with_suffix('.global.csv'), global_data, delimiter=',', header='time,relative_L2_error', comments='')
        offset = 1
    for n in range(3):
        axes[n+offset].plot(a[:,0]/TREF, errors[:,n], label=args.label)
        axes[n+offset].set_ylabel(r'$\epsilon_{L%d}$' % (n+1))
    for ax in axes:
        ax.axvline(.2, color='gray', linestyle='--', linewidth=1)
        ax.grid(alpha=.25)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    axes[0].legend()
    axes[-1].set_xlabel(r'$t/t_{ref}$, $t_{ref}=2.8\times10^{-4}$ s')
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    np.savetxt(args.output.with_suffix('.csv'), np.column_stack((a[:,0],a[:,0]/TREF,errors)), delimiter=',', header='time,t_over_tref,error_L1,error_L2,error_L3', comments='')
    print(f'Saved {args.output}; maximum local errors: {errors.max(axis=0)}')
    if a[-1,0] >= .2*TREF:
        print('Caution: t/tref >= 0.2 can be contaminated by reflections from the wider boundary (paper, p. 5129).')

if __name__ == '__main__':
    main()
