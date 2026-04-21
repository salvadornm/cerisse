# IBM Test Cases

Test suite for the Immersed Boundary Method (IBM) module, covering all three backends
(BVH-CPU, BVH-GPU, CGAL-CPU) in 2D and 3D, plus application cases.

## Directory structure

```
ibm_tests/
  prob.h              shared prob.h for 6 basic cases (Mach 4 shock over geometry)
  circle.dat          2D geometry (circle, r=0.1)
  sphere.stl          3D geometry (sphere, r=0.1)

  2d_bvh_cpu/         BVH backend, CPU, 2D
  3d_bvh_cpu/         BVH backend, CPU, 3D
  2d_bvh_gpu/         BVH backend, GPU (CUDA), 2D
  3d_bvh_gpu/         BVH backend, GPU (CUDA), 3D
  2d_cgal/            CGAL backend, CPU, 2D  (requires Boost, CGAL, GMP, MPFR)
  3d_cgal/            CGAL backend, CPU, 3D  (requires Boost, CGAL, GMP, MPFR)

  complex_geom/       Complex polygon, BVH-CPU 2D (own prob.h + geometry)
  airfoil_static/     Static diamond wedge, Ma=2, BVH-CPU 2D (own prob.h + geometry)
  multi_body/         7 cylinders, BVH CPU+GPU 2D (own prob.h + 7 geometry files)
```

Also see `exm/fsi/` for moving-geometry (FSI) cases:
- `exm/fsi/airfoil_moving/` — 2D diamond wedge with 3-DOF motion
- `exm/fsi/sphere_moving_3d/` — 3D sphere with 3-DOF motion

## Quick start

Each subdirectory contains a `GNUmakefile` and `inputs`. To build and run:

```bash
cd 2d_bvh_cpu
make -j8
mpirun -np 2 ./main2d.gnu.TPROF.MPI.ex inputs
```

For GPU cases, a CUDA-capable GPU is required:
```bash
cd 2d_bvh_gpu
make -j8
mpirun -np 1 ./main2d.gnu.TPROF.MPI.CUDA.ex inputs
```

For GPU cases that provide a separate `GNUmakefile.gpu`:
```bash
make -f GNUmakefile.gpu -j8
```

For CGAL cases, set `GMP_HOME` and `MPFR_HOME` if not in system paths:
```bash
cd 2d_cgal
make -j8 GMP_HOME=/usr MPFR_HOME=/usr
```

## Test matrix

| Case            | Backend  | DIM | CUDA | CGAL | Geometry            | Physics           |
|-----------------|----------|-----|------|------|---------------------|-------------------|
| 2d_bvh_cpu      | BVH      | 2   | no   | no   | circle              | Mach 4 shock      |
| 3d_bvh_cpu      | BVH      | 3   | no   | no   | sphere              | Mach 4 shock      |
| 2d_bvh_gpu      | BVH      | 2   | yes  | no   | circle              | Mach 4 shock      |
| 3d_bvh_gpu      | BVH      | 3   | yes  | no   | sphere              | Mach 4 shock      |
| 2d_cgal         | CGAL     | 2   | no   | yes  | circle              | Mach 4 shock      |
| 3d_cgal         | CGAL     | 3   | no   | yes  | sphere              | Mach 4 shock      |
| complex_geom    | BVH      | 2   | no   | no   | complex polygon     | Mach 4 shock      |
| airfoil_static  | BVH      | 2   | no   | no   | diamond wedge       | Ma=2 steady flow  |
| multi_body      | BVH      | 2   | both | no   | 7 cylinders         | Mach 4 shock (WENO-Z5 + NS + isothermal noslip, 1024^2) |

## Reference hardware configuration

All cases were tested and verified on the following setup:

| Component | Specification |
|-----------|--------------|
| GPU       | NVIDIA GeForce RTX 4090 Laptop GPU (sm_89, compute capability 8.9) |
| CPU       | 32 cores |
| Compiler  | GCC 13.3.0 |
| CUDA      | 12.0 |
| MPI       | OpenMPI |
| AMReX     | 25.12 |
| OS        | Ubuntu 24.04 (WSL2, kernel 5.15) |

## Adapting GPU cases to different hardware

The GPU cases are configured with `CUDA_ARCH=8.0` (or `8.9`) in their GNUmakefiles.

**To use a different GPU**, change `CUDA_ARCH` in the GNUmakefile to match
your GPU's compute capability:

| GPU Family             | Cards (examples)             | CUDA_ARCH |
|------------------------|------------------------------|-----------|
| Volta                  | V100, Titan V                | 7.0       |
| Turing                 | RTX 2080, T4                 | 7.5       |
| Ampere                 | A100, RTX 3090               | 8.0       |
| Ada Lovelace           | RTX 4090, L40                | 8.9       |
| Hopper                 | H100                         | 9.0       |
| Blackwell              | B200                         | 10.0      |

You can query your GPU's compute capability with:
```bash
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader
```

Or override at build time without editing the file:
```bash
make -j8 CUDA_ARCH=7.0   # for V100
```

## Notes

- The 6 basic cases share a common `prob.h` in the parent directory; `complex_geom`,
  `airfoil_static`, and `multi_body` each have their own `prob.h`.
- GPU cases use `amrex.the_arena_is_managed=0`; CPU cases use `=1`.
- CGAL backend is CPU-only. If `USE_CGAL=TRUE` and `USE_CUDA=TRUE` are both set,
  the build system automatically falls back to BVH with a warning.
- `multi_body` provides both `GNUmakefile` (CPU) and `GNUmakefile.gpu` (GPU).

## Internal-energy floor (CLIP_MINTEMP)

The `CLIP_MINTEMP` flag in each GNUmakefile controls the lower bound used
by `cons2prims` to clip internal energy. Strong shock-on-body interactions
(Mach ≥ 3 impacting an IBM wall) can drive WENO/Roe reconstruction into
negative-internal-energy overshoot and crash the run. Enable the floor to
absorb the undershoot physically:

- `CLIP_MINTEMP = TRUE`  → `ei_min = R * 10 K / (γ-1) ≈ 7200 J/kg`  (T_min = 10 K)
- `CLIP_MINTEMP = FALSE` → `ei_min = p_min / (γ-1) ≈ 2.5e-8 J/kg`   (effectively off)

Recommended: set `CLIP_MINTEMP = TRUE` for `multi_body`, `complex_geom`,
`2d_bvh_{cpu,gpu}`, `3d_bvh_{cpu,gpu}`, and any other Mach ≥ 3 shock/body
case. Low-Mach / smooth-flow cases (cylinder_Re40 etc.) are insensitive.

## 2D axisymmetric (RZ) cases

The 2D cases in this directory are Cartesian (`coord_sys=0`). The RZ
example `exm/ibm/esa_2d/` exercises axisymmetric IBM; RZ + IBM + shock
has a known residual instability near the r = 0 axis that is not yet
covered by the above fixes.
