# IBM internal pipe flow

This is a small 3D immersed-boundary example for flow inside a pipe represented
by a closed STL. The pipe interior is fluid and the exterior is solid. The
inlet and outlet conditions are applied at the pipe caps by the IBM wall model.
The input uses an O2/N2 mixture, an isothermal 300 K inlet, and a 70 m/s inlet
speed. The current `inputs` file is a two-step smoke test; increase `max_step`
for longer runs.

## Build and run

From this directory, build with the Cerisse dependencies configured under the
repository's `lib` directory:

```bash
make -j4
./main3d.gnu.ex inputs
```

The case uses the BVH IBM backend on CPU (`USE_MPI = FALSE`, `USE_CUDA = FALSE`).
It relies on the IBM internal-flow and interpolation changes in the branch.

## PelePhysics BlackBox shutdown patch

The `patches/` directory contains a candidate patch for a PelePhysics
`BlackBoxFunction.cpp` shutdown failure observed after the flow computation
finishes. It is included for review and is **not applied automatically** by the
case build. From the PelePhysics source directory, first check whether it
applies, then apply it if appropriate:

```bash
patch --dry-run -p1 < ../../exm/ibm/pipe_IBM/patches/pelephysics-blackbox-registration-lifetime.patch
patch -p1 < ../../exm/ibm/pipe_IBM/patches/pelephysics-blackbox-registration-lifetime.patch
```

Rebuild after applying the patch. The path above assumes the PelePhysics source
directory is `lib/PelePhysics` in the repository checkout.
