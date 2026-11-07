# FSI: Moving Sphere (3D)

3D quiescent atmosphere with a sphere undergoing prescribed 3-DOF rigid-body
motion (Lissajous translation + pitch). This case validates the full 3D FSI
pipeline: STL geometry loading, BVH construction, rigid transform application,
marker recomputation, ghost-point correction, solid-cell flood-fill, and FSI
load calculation.

## Physics

- Quiescent atmosphere at sea level (p=101325 Pa, T=300 K)
- Sphere (r=0.1 m), loaded from STL
- Adiabatic no-slip wall with Sutherland viscosity
- WENO-Z5 + viscous
- SSP-RK(3,3) time integration

## Build and run

### CPU (validation)

```bash
make -j8
mpirun -np 4 ./main3d.gnu.TPROF.MPI.ex inputs_cpu
```

### GPU (production)

```bash
make -f GNUmakefile.gpu -j8
mpirun -np 1 ./main3d.gnu.TPROF.MPI.CUDA.ex inputs_gpu
```

## GPU architecture

Default `CUDA_ARCH=8.9` (RTX 4090). Override:
```bash
make -f GNUmakefile.gpu -j8 CUDA_ARCH=8.0   # A100
```

## Verified runs

- CPU: 10 steps, 4 MPI ranks, 64^3, FSI loads computed correctly
- GPU: 10 steps, 1 MPI rank, 128^3, stable
