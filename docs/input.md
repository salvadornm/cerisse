---
icon: laptop-arrow-down
cover: >-
  https://images.unsplash.com/photo-1577401239170-897942555fb3?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw2fHxwcm9ibGVtfGVufDB8fHx8MTczMDk3NjczN3ww&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Input

This page provides detailed info on the `input` file options. These arguments are passed to AMReX, not all arguments are required (see [examples](https://app.gitbook.com/s/lJmXG8dNtMoIE65XawEJ/examples))

## AMReX options

The AMReX options covers control of

* the problem domain definition
* time-stepping
* gridding and load balancing
* output files
* checkpoint and restarting

The reference is available on [AMReX's documentation](https://amrex-codes.github.io/amrex/docs_html/Inputs_Chapter.html).

{% hint style="info" %}
In the tables below, _DIM_ means the number of dimensions, _Int_ means integer, and _Bool_ means boolean value (0 for False and 1 for True). If the option has no default value, a value must be given by the user.
{% endhint %}

### Problem definition and time-stepping

<table><thead><tr><th width="228">Option</th><th>Type</th><th align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong><code>max_step</code></strong></td><td>Int</td><td align="center"></td><td>Maximum number of time steps to take</td></tr><tr><td><strong><code>stop_time</code></strong></td><td>Real</td><td align="center"></td><td>Maximum time to reach</td></tr><tr><td><strong><code>time_step</code></strong></td><td>Real</td><td align="center"></td><td>dt (base level), higher level time step is based on number of subcycles</td></tr><tr><td><strong><code>cfl</code></strong></td><td>Real</td><td align="center"></td><td>CFL (incompatible option with time_step)</td></tr><tr><td><strong><code>geometry.is_periodic</code></strong></td><td>DIM * Int</td><td align="center">0 0 0</td><td>1 for true, 0 for false (one value for each coordinate direction)</td></tr><tr><td><strong><code>geometry.coord_sys</code></strong></td><td>Int</td><td align="center">0</td><td>0 = Cartesian; 1 = Cylindrical; 2 = Spherical (only support Cartesian)</td></tr><tr><td><strong><code>geometry.prob_lo</code></strong></td><td>DIM * Real</td><td align="center">0 0 0</td><td>Low corner of physical domain (physical not index space)</td></tr><tr><td><strong><code>geometry.prob_hi</code></strong></td><td>DIM * Real</td><td align="center"></td><td>High corner of physical domain (physical not index space)</td></tr><tr><td><strong><code>geometry.prob_extent</code></strong></td><td>DIM * Real</td><td align="center"></td><td>Extent of physical domain, choose between this or <code>prob_hi</code></td></tr><tr><td><strong><code>amr.n_cell</code></strong></td><td>DIM * Int</td><td align="center"></td><td>Number of cells at level 0 in each coordinate direction</td></tr></tbody></table>

### Gridding and load balancing

<table><thead><tr><th width="252">Option</th><th width="160">Type</th><th width="134" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong><code>amr.max_level</code></strong></td><td>Int</td><td align="center">0</td><td>Maximum level of refinement allowed<br>(0 when single-level)</td></tr><tr><td><strong><code>amr.ref_ratio</code></strong></td><td>Int *(nlev-1)</td><td align="center"></td><td>Refeniment ratio per level. If the number of <code>ref_ratio</code> is less than the number of levels - 1, the last entry will be automatically propagated</td></tr><tr><td><strong><code>amr.regrid_int</code></strong></td><td>Int</td><td align="center">-1</td><td>How often to regrid (in number of steps). No regridding will occur if set to &#x3C; 0</td></tr><tr><td><strong><code>amr.max_grid_size</code></strong></td><td>Int</td><td align="center">32</td><td>Maximum number of cells in each grid in all directions</td></tr><tr><td><strong><code>amr.blocking_factor</code></strong></td><td>Int</td><td align="center">8</td><td>Each grid must be divisible by blocking_factor in all directions (must be 1 or power of 2)</td></tr><tr><td><strong><code>amr.refine_grid_layout</code></strong></td><td>Bool</td><td align="center">1</td><td>Split grids in half until the number of grids is no less than the number of procs</td></tr><tr><td><strong><code>amr.n_error_buf</code></strong></td><td>DIM * Int</td><td align="center">1 1 1</td><td>Buffer in added around tagged cells</td></tr><tr><td><strong><code>amr.grid_eff</code></strong></td><td>Real</td><td align="center">0.7</td><td>Target value of the percentage of tagged cells in the grids</td></tr><tr><td><strong><code>amr.loadbalance_level0_int</code></strong></td><td>Int</td><td align="center">2</td><td>How often to do load balance (in number of steps). For single level (i.e., amr.max_level=0) only</td></tr><tr><td><strong><code>amr.loadbalance_with_workestimates</code></strong></td><td>Bool</td><td align="center">0</td><td>For multi-level runs, load balance is done during regrid and thus the load balance interval is controlled by <code>regrid_int</code></td></tr><tr><td><strong><code>amr.loadbalance_max_fac</code></strong></td><td>Real</td><td align="center">1.5</td><td>This controls the change in the maximum number of boxes that can be assigned to an MPI rank in load balancing</td></tr></tbody></table>

### Outputs and Restarting

<table><thead><tr><th width="263">Option</th><th>Type</th><th width="133" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong><code>amr.plot_files_output</code></strong></td><td>Bool</td><td align="center">1</td><td>Output plotfile or not (redundent because one can set plot_int = -1 to disable output)</td></tr><tr><td><strong><code>amr.plot_file</code></strong></td><td>String</td><td align="center">./plot/plt</td><td>Prefix of plotfile output</td></tr><tr><td><strong><code>amr.plot_int</code></strong></td><td>Int</td><td align="center">-1</td><td>Frequency of plotfile output; if -1 then no plotfiles will be written</td></tr><tr><td><strong><code>amr.derive_plot_vars</code></strong></td><td>Strings</td><td align="center">NONE</td><td>List of derived variables to plot; can use "ALL" or "NONE" to select the variables.</td></tr><tr><td><strong><code>amr.checkpoint_files_output</code></strong></td><td>Bool</td><td align="center">1</td><td>Same as plot_files_output, but for checkpoint files</td></tr><tr><td><strong><code>amr.check_file</code></strong></td><td>String</td><td align="center">chk</td><td>Prefix of checkpoint file output</td></tr><tr><td><strong><code>amr.check_int</code></strong></td><td>Int</td><td align="center">-1</td><td>Frequency of checkpoint file output; if -1 then no plotfiles will be written</td></tr><tr><td><strong><code>amr.restart</code></strong></td><td>String</td><td align="center"></td><td>If present, then the name of checkpoint file to restart from</td></tr><tr><td><strong><code>amr.plotfile_on_restart</code></strong></td><td>Bool</td><td align="center">0</td><td>Write a plotfile when immediately after restart or not</td></tr></tbody></table>

### GPU-related parameters

<table><thead><tr><th width="252">Option</th><th>Type</th><th width="163" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong><code>amrex.the_arena_init_size</code></strong></td><td>Int</td><td align="center">3/4 of total device memory</td><td>GPU device memory allocated to The_Arena (in bytes)</td></tr><tr><td></td><td></td><td align="center"></td><td></td></tr></tbody></table>

### Other parameters

<table><thead><tr><th width="244">Option</th><th>Type</th><th align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong><code>amrex.omp_threads</code></strong></td><td>String or Int</td><td align="center"><code>system</code></td><td><code>nosmt</code>: avoid using threads for virtual cores (aka Hyperthreading or SMT), as is default in OpenMP; <code>system</code>: use the environment variable <code>OMP_NUM_THREADS</code>. For Integer values, <code>OMP_NUM_THREADS</code> is ignored.</td></tr><tr><td><strong><code>amrex.fpe_trap_invalid</code></strong></td><td>Bool</td><td align="center">0</td><td>Produce error when invalid floating-point arthematic is detected. Helpful for debugging</td></tr></tbody></table>

### Boundary conditions

```ini
# 0 = Interior                               3 = Symmetry
# 1 = Inflow / UserBC                        4 = SlipWall =3
# 2 = Outflow (First Order Extrapolation)    5 = NoSlipWall (adiabatic)
cns.lo_bc = 1 5 0
cns.hi_bc = 2 5 0
```

In the above example, the bc in _x_ would be inflow (at lower boundary) and outflow (at upper boundary) while no specific boudnary will be defined in _z_

| Option          | Type       | Default | Description                           |
| --------------- | ---------- | ------- | ------------------------------------- |
| **`cns.lo_bc`** | DIM \* Int |         | BC flags at lower boundaries in x,y,z |
| **`cns.hi_bc`** | DIM \* Int |         | BC flags at upper boundaries in x,y,z |

If option "0" is selected, the corresponding `geometry.is_periodic` must also be set.

If option "1" is selected, the `bcnormal` function in `prob.H` will be activated (see boundary conditions).

### NSCBC boundary conditions (optional)

Cerisse supports **Navier-Stokes Characteristic Boundary Conditions (NSCBC)** based on the **LODI (Local One-Dimensional Inviscid)** formulation. NSCBC prescribes only incoming characteristic information and lets outgoing waves leave the domain, reducing artificial reflections from subsonic open boundaries.

NSCBC is selected independently on every low and high domain face. `cns.nscbc_lo` and `cns.nscbc_hi` contain one type per coordinate direction. The face must be non-periodic. The ordinary `cns.lo_bc` and `cns.hi_bc` settings are still required to fill the remaining physical ghost cells; use an open-boundary value such as `2` on an NSCBC inlet or outlet.

#### Available NSCBC types

| Type | Description |
| ---- | ----------- |
| **0** | Disabled; use the ordinary boundary condition. |
| **1** | Subsonic non-reflecting inflow with relaxation toward velocity and temperature targets. |
| **2** | Pure non-reflecting subsonic outflow. The incoming acoustic amplitude is set to zero, so the mean pressure is free to evolve. |
| **3** | Subsonic non-reflecting outflow with controlled relaxation toward `Ptarget`. |

The following two-dimensional example applies a target-relaxed inlet at `xlo`, a pressure-relaxed outlet at `xhi`, and no NSCBC treatment in the y direction:

```ini
cns.lo_bc = 2 5
cns.hi_bc = 2 5

cns.nscbc_lo = 1 0
cns.nscbc_hi = 3 0
cns.nscbc_order = 2

# x-low inflow targets
cns.nscbc_xlo_utarget = 20.0
cns.nscbc_xlo_vtarget = 0.0
cns.nscbc_xlo_Ttarget = 300.0
cns.nscbc_xlo_eta = 1.0

# x-high pressure-relaxed outflow
cns.nscbc_xhi_Ptarget = 101325.0
cns.nscbc_xhi_sigma = 0.28
cns.nscbc_xhi_Lchar = 1.0
cns.nscbc_xhi_Mmax = 0.1
cns.nscbc_xhi_use_transverse = 1
cns.nscbc_xhi_beta_transverse = 0.1
```

{% hint style="warning" %}
NSCBC parameters are **face-specific**. Use one of `xlo`, `xhi`, `ylo`, `yhi`, `zlo`, or `zhi` between `nscbc_` and the parameter name. For example, use `cns.nscbc_xlo_utarget`, not `cns.nscbc_utarget`.
{% endhint %}

| Option | Type | Default | Description |
| ------ | ---- | ------- | ----------- |
| **`cns.nscbc_lo`** | DIM \* Int | 0 0 0 | NSCBC type on each low face. |
| **`cns.nscbc_hi`** | DIM \* Int | 0 0 0 | NSCBC type on each high face. |
| **`cns.nscbc_order`** | Int | 2 | Order of the inward one-sided normal derivative; supported values are `1` and `2`. |

#### Face-specific parameters

In the table below, `<face>` is `xlo`, `xhi`, `ylo`, `yhi`, `zlo`, or `zhi` as applicable to the build dimension.

| Option | Type | Default | Used by | Description |
| ------ | ---- | ------- | ------- | ----------- |
| **`cns.nscbc_<face>_utarget`** | Real | 0 | Type 1 | Target x velocity. |
| **`cns.nscbc_<face>_vtarget`** | Real | 0 | Type 1 | Target y velocity. |
| **`cns.nscbc_<face>_wtarget`** | Real | 0 | Type 1 | Target z velocity. |
| **`cns.nscbc_<face>_Ttarget`** | Real | 300 K | Type 1 | Target temperature. |
| **`cns.nscbc_<face>_eta`** | Real | 1.0 | Type 1 | Velocity and temperature relaxation coefficient. |
| **`cns.nscbc_<face>_inflow_target`** | Int | 0 | Type 1 | Normal inflow target: `0` for velocity or `1` for mass flux. |
| **`cns.nscbc_<face>_mass_flux_target`** | Real | 0 | Type 1 | Positive mass flux into the domain in kg/(m2 s); used when `inflow_target = 1`. |
| **`cns.nscbc_<face>_Ptarget`** | Real | 101325 Pa | Type 3 | Target outlet pressure. |
| **`cns.nscbc_<face>_sigma`** | Real | 0.28 | Type 3 | Pressure relaxation coefficient. Set to zero for a perfectly non-reflecting pressure wave. |
| **`cns.nscbc_<face>_Lchar`** | Real | 1.0 | Type 3 | Characteristic length used by pressure relaxation. |
| **`cns.nscbc_<face>_Mmax`** | Real | 0.1 | Type 3 | Maximum/reference Mach number used to scale pressure relaxation. |
| **`cns.nscbc_<face>_use_transverse`** | Bool | 0 | Types 2, 3 | Include transverse characteristic terms where the stencil is available. |
| **`cns.nscbc_<face>_beta_transverse`** | Real | 1.0 | Type 3 | Transverse relaxation coefficient, normally between 0 and 1. |

#### Spatially varying inflow targets

A problem can replace the constant type-1 velocity and temperature targets with an arbitrary profile over the boundary face. Enable the compile-time flag in the problem `GNUmakefile`:

```makefile
USE_MANUAL_NSBC_TARGET = TRUE

ifeq ($(USE_MANUAL_NSBC_TARGET),TRUE)
  DEFINES += -DUSE_MANUAL_NSBC_TARGET
endif
```

Then define the following GPU-callable function inside namespace `PROB` in `prob.h`:

```cpp
#ifdef USE_MANUAL_NSBC_TARGET
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nscbc_target(amrex::Real x, amrex::Real y, amrex::Real z,
                  int dir, int side_sign,
                  amrex::Real& u, amrex::Real& v,
                  amrex::Real& w, amrex::Real& T)
{
    // Example: parabolic profile at the x-low inflow.
    if (dir == 0 && side_sign > 0) {
        constexpr amrex::Real half_height = 1.0e-3;
        constexpr amrex::Real umax = 20.0;
        const amrex::Real yn = y / half_height;
        u = umax * amrex::max(amrex::Real(0.0),
                              amrex::Real(1.0) - yn*yn);
        v = 0.0;
        w = 0.0;
        T = 300.0;
    }
}
#endif
```

`dir` is `0`, `1`, or `2` for x, y, or z. `side_sign` is `+1` on a low face and `-1` on a high face. The normal coordinate is exactly on the physical domain face; tangential coordinates are at the ghost-cell centre. The four target references are initialized from the corresponding face-specific input values, so the function may modify only the required fields or faces.

The callback is evaluated only for type-1 NSCBC boundaries and must use GPU-compatible code. It changes velocity and temperature targets only. When `inflow_target = 1`, `mass_flux_target` overrides the callback's normal target velocity; tangential velocity and temperature targets still come from the callback.

### Time Marching

| Option              | Type | Default | Description              |
| ------------------- | ---- | ------- | ------------------------ |
| **`cns.order_rk`**  | Int  |         | Order of RK (-2/0/1/2/3) |
| **`cns.stages_rk`** | Int  |         | RK stages (1/2/3)        |

_Runge-Kutta order options_

`-2` Original RK2 Scheme`0` (for testing) returns RHS`1` Euler scheme`2` Low storage Runge-Kutta second order SSPRK(m,2) with _m_ stages`3` Third order SSP Runge-Kutta (options stages 3 or 4)

### Geometry EB options

In the **input** file, users should specify the geometry of the embedded boundary with `eb2.geom_type`, then supply the required parameters in the format of `eb2.{geom_param}`.

| eb2.geom\_type            | additional parameters required                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`all_regular`**         | no EB, no additional parameters needed                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| **`plane`**               | `plane_point` - a point where the plane intersects, `plane_normal` - the normal vector of the plane that points into the solid                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| **`sphere`**              | `sphere_center`, `sphere_radius`, `sphere_has_fluid_inside` - bool value fluid inside or outside                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| **`cylinder`**            | `cylinder_center`, `cylinder_radius`, `cylinder_height`, `cylinder_direction` - (0,1,2) for (x,y,z), and `cylinder_has_fluid_inside`                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| **`box`**                 | `box_lo` and `box_hi` - lower and upper corners of the box, and `box_has_fluid_inside`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| **`stl`**                 | `stl_file` - the STL file name, `stl_scale` - the scaling factor in all directions, and `stl_center` - center of object in relation to the cooridinate system in the file, and `stl_reverse_normal` - essentially stl\_has\_fluid\_inside                                                                                                                                                                                                                                                                                                                                    |
| **`triangles`**           | <p><code>num_tri</code> - number of triangles, up to 5 (change the value in <code>custom_geometry.cpp</code> if needed),<br>for each triangle, <code>{i}</code> from 0 to num_tri-1, <code>tri_{i}</code><em><code>point_0</code>, <code>tri</code></em><code>{i}</code><em><code>point_1</code>, <code>tri</code></em><code>{i}_point_2</code> - three points that define the triangle, give the points in anti-clockwise direction to set solid inside of the triangle. The z-coordinate isn't really needed because the triangle will be extruded in the z-direction.</p> |
| **`redistribution_type`** | `"NoRedist"`, `"FluxRedist"`, `"StateRedist"` or `"NewRedist"` The first one, indicates no flux redistribution (default). FluxRedist/StateRedist correspond to AMReX options. NewRedist does not use AMReX functions but it is a variant of FluxRedist                                                                                                                                                                                                                                                                                                                       |

{% hint style="warning" %}
You can only choose **one** geometry type and **one** geometry. If you want to use multiple geometries, you need to define your own geometry
{% endhint %}

Below are some examples:

```ini
eb2.geom_type = all_regular

eb2.geom_type = cylinder
eb2.cylinder_direction = 2
eb2.cylinder_radius = 0.25
eb2.cylinder_center = 1.0 2.0 0.0
eb2.cylinder_has_fluid_inside = 0

eb2.geom_type = box
eb2.box_lo = -1.0 -1.0 0.0
eb2.box_hi =  1.0  2.0  0.0
eb2.box_has_fluid_inside = 0

# This gives the same geometry as the box above
eb2.geom_type = triangles 
triangles.num_tri = 2
triangles.tri_0_point_0 = -1.0  1.0 0.0
triangles.tri_0_point_1 = -1.0 -1.0 0.0
triangles.tri_0_point_2 =  1.0 -1.0 0.0
triangles.tri_1_point_0 = -1.0  1.0 0.0
triangles.tri_1_point_1 =  1.0 -1.0 0.0
triangles.tri_1_point_2 =  1.0  1.0 0.0
```
