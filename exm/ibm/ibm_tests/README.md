# IBM Test Cases

Test suite for the Immersed Boundary Method (IBM) module, covering all three backends
(BVH-CPU, BVH-GPU, CGAL-CPU) in 2D and 3D, plus two application cases.

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

  complex_geom/       Multiple polygons, BVH-CPU 2D (own prob.h + geometry)
  airfoil_static/     Static diamond wedge, Ma=2, BVH-CPU 2D (own prob.h + geometry)
```

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

For CGAL cases, set `GMP_HOME` and `MPFR_HOME` if not in system paths:
```bash
cd 2d_cgal
make -j8 GMP_HOME=/usr MPFR_HOME=/usr
```

## Test matrix

| Case            | Backend  | DIM | CUDA | CGAL | Geometry          | Physics           |
|-----------------|----------|-----|------|------|-------------------|-------------------|
| 2d_bvh_cpu      | BVH      | 2   | no   | no   | circle            | Mach 4 shock      |
| 3d_bvh_cpu      | BVH      | 3   | no   | no   | sphere            | Mach 4 shock      |
| 2d_bvh_gpu      | BVH      | 2   | yes  | no   | circle            | Mach 4 shock      |
| 3d_bvh_gpu      | BVH      | 3   | yes  | no   | sphere            | Mach 4 shock      |
| 2d_cgal         | CGAL     | 2   | no   | yes  | circle            | Mach 4 shock      |
| 3d_cgal         | CGAL     | 3   | no   | yes  | sphere            | Mach 4 shock      |
| complex_geom    | BVH      | 2   | no   | no   | complex polygon   | Mach 4 shock      |
| airfoil_static  | BVH      | 2   | no   | no   | diamond wedge     | Ma=2 steady flow  |

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

The GPU cases are configured with `CUDA_ARCH=8.0` in their GNUmakefiles.
This targets Ada Lovelace / Ampere GPUs (RTX 30xx/40xx, A100, etc.).

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

- The 6 basic cases share a common `prob.h` in the parent directory; `complex_geom`
  and `airfoil_static` each have their own `prob.h` tailored to their geometry/physics.
- GPU cases use `amrex.the_arena_is_managed=0`; CPU cases use `=1`.
- CGAL backend is CPU-only. If `USE_CGAL=TRUE` and `USE_CUDA=TRUE` are both set,
  the build system automatically falls back to BVH with a warning.
