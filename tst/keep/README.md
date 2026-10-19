# KEEP flux regression tests

CPU tests using real AMReX arrays (with bounds checking in the default debug
build) and a composition-dependent test EOS. No PelePhysics mechanism is needed.

```sh
make -j4
./main3d.gnu.DEBUG.ex
make -j4 -W main.cpp KEEP_TEST_BOUNDARY=EB
./main3d.gnu.DEBUG.ex
make -j4 -W main.cpp KEEP_TEST_BOUNDARY=IBM
./main3d.gnu.DEBUG.ex
make -j4 -W main.cpp KEEP_TEST_BOUNDARY=IBM KEEP_TEST_SPECIES=1
./main3d.gnu.DEBUG.ex
```

`-W main.cpp` forces recompilation when changing test-only preprocessor flags.
The default species count is 3. Tests cover orders 2/4/6, equality with the
legacy `CentralKEEP_SPEC.h` flux, species-sum conservation, composition-dependent
internal energy, physical constant-state fluxes, the normal MFIter entry point,
EB cut-cell and IBM ghost-point second-order fallback, zero solid-solid flux,
and dissipation with inaccessible covered states poisoned with NaNs. Directional
symmetry checks cover every compiled spatial direction. `DIM=1` or `DIM=2` can
also be selected; adjust the executable name accordingly.

These are CPU flux tests. They do not test GPU compilation or replace a full
geometry-specific EB/IBM simulation.

## Selecting the new scheme in a problem

`RHS.h` includes `Keep.h`. The parameter contract matches `Skew.h`:

```cpp
struct keep_parm {
  static constexpr int order = 4;       // 2, 4, or 6
  static constexpr bool dissipation = false;
  static constexpr amrex::Real C2skew = 0.1;
  static constexpr amrex::Real C4skew = 0.016;
};
using ProbRHS = rhs_dt<keep_t<keep_parm, ProbClosures>,
                      no_diffusive_t, no_source_t>;
```

Keep the existing diffusive/source classes if the problem uses them. The build
selects `eflux` normally or `eflux_ibm` for `CNS_USE_EB`/`AMREX_USE_GPIBM`.
The marker layout is the same as Skew: component 0 means solid; component 1
means an IBM ghost point or EB cut cell. Solid-solid faces are explicitly
zeroed. Any marked cell in the high-order stencil triggers a two-point KEEP
flux; this also avoids reaching deeper covered cells with fourth/sixth order.
EB wall fluxes and IBM ghost-state reconstruction remain in their existing
solver paths.

The optional shock sensor uses density and pressure. Near walls it restricts
sensor sampling to available fluid/ghost neighbors and uses two-point shock
damping without the high-order damping term. Unlike Skew's right-state-only
spectral radius, Keep uses the maximum adjacent |u_n|+c. The all-fluid damping
coefficients follow Skew's order-dependent coefficients. Dissipation changes the
preservation properties; it should be disabled for nondissipative KEEP studies.
The arithmetic species average is retained from KEEP_SPEC and is not claimed
to be a formally entropy-conservative flux for arbitrary mixtures.
