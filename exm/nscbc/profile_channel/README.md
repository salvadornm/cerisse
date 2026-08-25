# Profiled NSCBC channel

This two-dimensional verification case exercises a spatially varying NSCBC
inflow target at `xlo` and a pressure-relaxed non-reflecting NSCBC outflow at
`xhi`. The top and bottom are adiabatic no-slip walls. The inlet target is the
parabolic profile

```text
u(y) = umax [1 - (y/h)^2],  v = 0,  T = 300 K,
```

with `h = 1 mm` and centreline Mach number 0.1. The setup follows the structure
of the ideal-gas two-dimensional channel verification in Appendix D of Li et
al. (2025), except that this compact Cerisse example uses its built-in
adiabatic wall condition rather than Li's isothermal wall. The outlet uses
their pressure-relaxation value `sigma = 0.3`.

## Defining a target profile

Enable `USE_MANUAL_NSBC_TARGET` in the problem GNUmakefile and define this hook
inside namespace `PROB` in `prob.h`:

```cpp
#ifdef USE_MANUAL_NSBC_TARGET
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nscbc_target(amrex::Real x, amrex::Real y, amrex::Real z,
                  int dir, int side_sign, amrex::Real& u,
                  amrex::Real& v, amrex::Real& w, amrex::Real& T);
#endif
```

`side_sign` is `+1` on a low face and `-1` on a high face. Coordinates are
physical coordinates: the normal coordinate lies exactly on the domain face,
and tangential coordinates are at the ghost-cell centre. The four references
arrive initialized from the corresponding `cns.nscbc_*target` inputs, so a hook
may change only the fields and faces it owns. The hook is evaluated only for
NSCBC type 1 (subsonic target-relaxed inflow). It must be callable on the GPU:
avoid host-only state and use device-compatible AMReX math functions.

The manual hook changes velocity/temperature targets only. A mass-flux target
still overrides the normal target velocity, so use the default velocity target
mode for a prescribed normal profile.

Build and run from this directory with the normal Cerisse make/run workflow.

## Pressure-gradient check

After generating plotfiles, analyze the newest one with:

```bash
python3 plot_pressure.py
```

The script extracts centreline pressure, compares it with the fully developed
plane-Poiseuille slope `dp/dx = -2 mu umax/h^2`, and checks linear-fit `R^2`,
local-gradient variation, and slope error. It writes `pressure_over_x.png` and
prints PASS/FAIL results. Pass a plotfile explicitly or use `--strict` to make a
failed check return a nonzero exit status; run `python3 plot_pressure.py --help`
for all tolerances and physical parameters.
