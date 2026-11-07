# FSI: Moving Diamond Wedge Airfoil

2D supersonic flow (quiescent atmosphere) over a diamond wedge airfoil undergoing
prescribed 3-DOF rigid-body motion (translation + pitch). This case exercises the
full IBM+FSI pipeline: geometry update, BVH rebuild, marker recomputation, ghost-point
correction, solid-cell flood-fill, and FSI load calculation.

## Physics

- Quiescent atmosphere at sea level (p=101325 Pa, T=300 K)
- Diamond wedge with 5-degree half-angle, chord=0.1 m
- Adiabatic no-slip wall with Sutherland viscosity
- WENO-Z5 + viscous (no turbulence model)
- SSP-RK(3,3) time integration

## Motion kinematics (3-DOF Lissajous)

The airfoil oscillates in translation and pitch simultaneously:

```
x(t) = 0.005 * sin(2pi * 100 * t)        [m]   5mm amplitude, 100 Hz
y(t) = 0.003 * sin(2pi * 150 * t + pi/2)  [m]   3mm amplitude, 150 Hz
theta(t) = 10 * sin(2pi * 120 * t)        [deg]  +/-10 deg, 120 Hz
```

The phase offset `pi/2` between x and y creates a Lissajous trajectory.
All motion parameters are set in the `inputs` file under the `motion.*` namespace
and can be changed at runtime without recompilation.

## Directory contents

```
airfoil_moving/
  prob.h                 problem definition (closures, IC, BC, wall model, motion)
  diamond_wedge_2d.dat   2D polygon geometry (diamond wedge)
  GNUmakefile            CPU build (USE_CUDA=FALSE)
  GNUmakefile.gpu        GPU build (USE_CUDA=TRUE)
  inputs_cpu             coarse grid, quick validation (256x128, AMR L1)
  inputs_gpu             fine grid, production run (1024x512, AMR L1)
```

## Build and run

### CPU (validation)

```bash
make -j8                    # builds main2d.gnu.TPROF.MPI.ex
mpirun -np 4 ./main2d.gnu.TPROF.MPI.ex inputs_cpu
```

### GPU (production)

```bash
make -f GNUmakefile.gpu -j8   # builds main2d.gnu.TPROF.MPI.CUDA.ex
mpirun -np 1 ./main2d.gnu.TPROF.MPI.CUDA.ex inputs_gpu
```

## Path configuration

The GNUmakefiles assume this directory is located at `exm/fsi/airfoil_moving/`
relative to the solver root. The key path variable is:

```makefile
AMR_SOLVER = $(abspath ../../../)    # -> solver root (cerisse/)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex
```

If you move this case to a different location, update `AMR_SOLVER` accordingly.
For example, if placed at `wrk/my_test/airfoil/`:

```makefile
AMR_SOLVER = $(abspath ../../../)    # adjust depth to reach solver root
```

## GPU architecture

The GPU GNUmakefile is configured with `CUDA_ARCH=8.9` (RTX 4090 Laptop GPU).
To use a different GPU, change `CUDA_ARCH`:

| GPU            | CUDA_ARCH |
|----------------|-----------|
| V100           | 7.0       |
| RTX 2080/T4    | 7.5       |
| A100/RTX 3090  | 8.0       |
| RTX 4090/L40   | 8.9       |
| H100           | 9.0       |

Override at build time:
```bash
make -f GNUmakefile.gpu -j8 CUDA_ARCH=8.0
```

## Compile flags

| Flag | CPU | GPU | Purpose |
|------|-----|-----|---------|
| `USE_GPIBM=TRUE` | yes | yes | Enable IBM module (BVH backend) |
| `USE_FSI=TRUE` | yes | yes | Enable FSI (rigid-body motion, load computation) |
| `CLIP_MINTEMP=TRUE` | yes | yes | Clip temperature floor (stability) |
| `USE_CUDA=TRUE` | no | yes | Enable CUDA GPU offloading |

## Input parameters of interest

| Parameter | CPU value | GPU value | Description |
|-----------|-----------|-----------|-------------|
| `amr.n_cell` | 256 128 | 1024 512 | Base grid resolution |
| `amr.max_level` | 1 | 1 | AMR levels (1 refinement) |
| `ib.move` | 1 | 1 | Enable geometry motion |
| `motion.amp_x` | 0.005 | 0.005 | x-oscillation amplitude [m] |
| `motion.freq_x` | 100 | 100 | x-oscillation frequency [Hz] |
| `motion.amp_theta` | 10 | 10 | Pitch amplitude [deg] |
| `fsi.rho_solid` | 2700 | 2700 | Solid density [kg/m3] (aluminum) |

## Reference test configuration

| Component | Specification |
|-----------|--------------|
| GPU       | NVIDIA RTX 4090 Laptop (sm_89) |
| CPU       | 32 cores |
| Compiler  | GCC 13.3 |
| CUDA      | 12.0 |
| AMReX     | 25.12 |

### Verified runs

- **CPU**: 10 steps, 4 MPI ranks, `inputs_cpu` — FSI loads computed correctly
- **GPU**: 1000 steps, 1 MPI rank, `inputs_gpu` — stable dt ~4.8e-7 s, no NaN/crash
