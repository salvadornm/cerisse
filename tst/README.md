# Tests

This directory contains unit and regression tests, together with standard
simulation cases such as Sod's shock tube, convective vortex, and Shu–Osher.

## `keep/`: KEEP flux regression tests

Tests the `keep_t` implementation in `src/rhs/Keep.h` using real AMReX arrays
and a composition-dependent test equation of state (EOS). No PelePhysics
mechanism is required. The tests cover:

- Second-, fourth-, and sixth-order fluxes and agreement with the legacy
  `CentralKEEP_SPEC.h` implementation.
- Constant-state fluxes, species-sum conservation, and composition-dependent
  internal energy.
- Embedded-boundary (EB) and immersed-boundary (IBM) stencil fallback,
  zero solid-solid flux, and optional dissipation near covered cells.

From the repository root, build and run the default CPU debug tests with:

```sh
cd tst/keep
make -j4
./main3d.gnu.DEBUG.ex
```

See [keep/README.md](keep/README.md) for the EB/IBM test variants, species-count
options, and instructions for selecting KEEP in a problem. These are flux-level
tests, not full EB/IBM simulations or GPU tests.

## `eos/`: equation-of-state unit tests

Uses Catch to compare CERISSE's calorifically perfect-gas EOS with its
PelePhysics-backed multispecies EOS configured as a single-species gamma-law
gas. The tests check internal energy (`RYP2E`), temperature and pressure
(`RYE2TP`), sound speed (`RYE2Cs`), and consistency of the combined
`RYE2TPCs` calculation with the separate calls.

Pressure, energy, and sound speed comparisons use floating-point ULP
tolerances. Cross-EOS temperature comparisons allow an absolute difference of
0.2 K because the two implementations use different molecular weights.

The GNUmakefile requires AMReX and PelePhysics under `lib/`, and selects the
PelePhysics Null mechanism, gamma-law EOS, and constant transport. From the
repository root:

```sh
cd tst/eos
make -j4 COMP=gnu DEBUG=TRUE USE_MPI=FALSE
./main3d.gnu.DEBUG.ex
```

Neither test executable requires a simulation input file.
