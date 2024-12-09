---
icon: laptop-arrow-down
cover: >-
  https://images.unsplash.com/photo-1577401239170-897942555fb3?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw2fHxwcm9ibGVtfGVufDB8fHx8MTczMDk3NjczN3ww&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Input

This page provides detailed info on the `input` file options. These arguments are passed to AMReX, not all arguents are reuired, see examples

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

<table><thead><tr><th width="228">Option</th><th>Type</th><th align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>max_step</strong></td><td>Int</td><td align="center"></td><td>Maximum number of time steps to take</td></tr><tr><td><strong>stop_time</strong></td><td>Real</td><td align="center"></td><td>Maximum time to reach</td></tr><tr><td><strong>time_step</strong></td><td>Real</td><td align="center"></td><td>dt (base level), higher level time step is based on number of subcycles</td></tr><tr><td><strong>cfl</strong></td><td>Real</td><td align="center"></td><td>CFL (incompatible option with time_step)</td></tr><tr><td><strong>geometry.is_periodic</strong></td><td>DIM * Int</td><td align="center">0 0 0</td><td>1 for true, 0 for false (one value for each coordinate direction)</td></tr><tr><td><strong>geometry.coord_sys</strong></td><td>Int</td><td align="center">0</td><td>0 = Cartesian; 1 = Cylindrical; 2 = Spherical (only support Cartesian)</td></tr><tr><td><strong>geometry.prob_lo</strong></td><td>DIM * Real</td><td align="center">0 0 0</td><td>Low corner of physical domain (physical not index space)</td></tr><tr><td><strong>geometry.prob_hi</strong></td><td>DIM * Real</td><td align="center"></td><td>High corner of physical domain (physical not index space)</td></tr><tr><td><strong>geometry.prob_extent</strong></td><td>DIM * Real</td><td align="center"></td><td>Extent of physical domain, choose between this or <code>prob_hi</code></td></tr><tr><td><strong>amr.n_cell</strong></td><td>DIM * Int</td><td align="center"></td><td>Number of cells at level 0 in each coordinate direction</td></tr></tbody></table>

### Gridding and load balancing

<table><thead><tr><th width="252">Option</th><th width="160">Type</th><th width="134" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amr.max_level</strong></td><td>Int</td><td align="center">0</td><td>Maximum level of refinement allowed<br>(0 when single-level)</td></tr><tr><td><strong>amr.ref_ratio</strong></td><td>Int *(nlev-1)</td><td align="center"></td><td>Refeniment ratio per level. If the number of <code>ref_ratio</code> is less than the number of levels - 1, the last entry will be automatically propagated</td></tr><tr><td><strong>amr.regrid_int</strong></td><td>Int</td><td align="center">-1</td><td>How often to regrid (in number of steps). No regridding will occur if set to &#x3C; 0</td></tr><tr><td><strong>amr.max_grid_size</strong></td><td>Int</td><td align="center">32</td><td>Maximum number of cells in each grid in all directions</td></tr><tr><td><strong>amr.blocking_factor</strong></td><td>Int</td><td align="center">8</td><td>Each grid must be divisible by blocking_factor in all directions (must be 1 or power of 2)</td></tr><tr><td><strong>amr.refine_grid_layou</strong>t</td><td>Bool</td><td align="center">1</td><td>Split grids in half until the number of grids is no less than the number of procs</td></tr><tr><td><strong>amr.n_error_buf</strong></td><td>DIM * Int</td><td align="center">1 1 1</td><td>Buffer in added around tagged cells</td></tr><tr><td><strong>amr.grid_eff</strong></td><td>Real</td><td align="center">0.7</td><td>Target value of the percentage of tagged cells in the grids</td></tr><tr><td><strong>amr.loadbalance_level0_int</strong></td><td>Int</td><td align="center">2</td><td>How often to do load balance (in number of steps). For single level (i.e., amr.max_level=0) only</td></tr><tr><td><strong>amr.loadbalance_with_workestimates</strong></td><td>Bool</td><td align="center">0</td><td>For multi-level runs, load balance is done during regrid and thus the load balance interval is controlled by <code>regrid_int</code></td></tr><tr><td><strong>amr.loadbalance_max_fac</strong></td><td>Real</td><td align="center">1.5</td><td>This controls the change in the maximum number of boxes that can be assigned to an MPI rank in load balancing</td></tr></tbody></table>

### Outputs and Restarting

<table><thead><tr><th width="247">Option</th><th>Type</th><th width="133" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amr.plot_files_output</strong></td><td>Bool</td><td align="center">1</td><td>Output plotfile or not (redundent because one can set plot_int = -1 to disable output)</td></tr><tr><td><strong>amr.plot_file</strong></td><td>String</td><td align="center">./plot/plt</td><td>Prefix of plotfile output</td></tr><tr><td><strong>amr.plot_int</strong></td><td>Int</td><td align="center">-1</td><td>Frequency of plotfile output; if -1 then no plotfiles will be written</td></tr><tr><td><strong>amr.derive_plot_vars</strong></td><td>Strings</td><td align="center">NONE</td><td>List of derived variables to plot; can use "ALL" or "NONE" to select all or none of the variables. See the full list of derived variables available in Cerisse below.</td></tr><tr><td><strong>amr.checkpoint_files_output</strong></td><td>Bool</td><td align="center">1</td><td>Same as plot_files_output, but for checkpoint files</td></tr><tr><td><strong>amr.check_file</strong></td><td>String</td><td align="center">chk</td><td>Prefix of checkpoint file output</td></tr><tr><td><strong>amr.check_int</strong></td><td>Int</td><td align="center">-1</td><td>Frequency of checkpoint file output; if -1 then no plotfiles will be written</td></tr><tr><td><strong>amr.restart</strong></td><td>String</td><td align="center"></td><td>If present, then the name of checkpoint file to restart from</td></tr><tr><td><strong>amr.plotfile_on_restar</strong>t</td><td>Bool</td><td align="center">0</td><td>Write a plotfile when immediately after restart or not</td></tr></tbody></table>

### GPU-related parameters

<table><thead><tr><th width="252">Option</th><th>Type</th><th width="163" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amrex.the_arena_init_size</strong></td><td>Int</td><td align="center">3/4 of total device memory</td><td>GPU device memory allocated to The_Arena (in bytes)</td></tr></tbody></table>

### Other parameters

<table><thead><tr><th width="244">Option</th><th>Type</th><th align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amrex.omp_threads</strong></td><td>String or Int</td><td align="center"><code>system</code></td><td><code>nosmt</code>: avoid using threads for virtual cores (aka Hyperthreading or SMT), as is default in OpenMP; <code>system</code>: use the environment variable <code>OMP_NUM_THREADS</code>. For Integer values, <code>OMP_NUM_THREADS</code> is ignored.</td></tr><tr><td><strong>amrex.fpe_trap_invalid</strong></td><td>Bool</td><td align="center">0</td><td>Produce error when invalid floating-point arthematic is detected. Helpful for debugging</td></tr></tbody></table>

### Boundary conditions

```ini
# 0 = Interior                               3 = Symmetry
# 1 = Inflow / UserBC                        4 = SlipWall =3
# 2 = Outflow (First Order Extrapolation)    5 = NoSlipWall (adiabatic)
cns.lo_bc = 1 5 0
cns.hi_bc = 2 5 0
```

In the above example, the bc in _x_ would be inflow (at lower boundary) and outflow (at upper boundary) while no specific boudnary will be defined in _z_

| Option         | Type       | Default | Description                           |
| -------------- | ---------- | ------- | ------------------------------------- |
| **cns.lo\_bc** | DIM \* Int |         | BC flags at lower boundaries in x,y,z |
| **cns.hi\_bc** | DIM \* Int |         | BC flags at upper boundaries in x,y,z |

If option "0" is selected, the corresponding `geometry.is_periodic` must also be set.

If option "1" is selected, the `bcnormal` function in `prob.H` will be activated (see boundary conditions).

### Geometry EB options

In the **input** file, users should specify the geometry of the embedded boundary with `eb2.geom_type`, then supply the required parameters in the format of `eb2.{geom_param}`.

| eb2.geom\_type      | additional parameters required                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `all_regular`       | no EB, no additional parameters needed                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `plane`             | `plane_point` - a point where the plane intersects, `plane_normal` - the normal vector of the plane that points into the solid                                                                                                                                                                                                                                                                                                                                                                             |
| `sphere`            | `sphere_center`, `sphere_radius`, `sphere_has_fluid_inside` - bool value fluid inside or outside                                                                                                                                                                                                                                                                                                                                                                                                           |
| `cylinder`          | `cylinder_center`, `cylinder_radius`, `cylinder_height`, `cylinder_direction` - (0,1,2) for (x,y,z), and `cylinder_has_fluid_inside`                                                                                                                                                                                                                                                                                                                                                                       |
| `box`               | `box_lo` and `box_hi` - lower and upper corners of the box, and `box_has_fluid_inside`                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `stl`               | `stl_file` - the STL file name, `stl_scale` - the scaling factor in all directions, and `stl_center` - center of object in relation to the cooridinate system in the file, and `stl_reverse_normal` - essentially stl\_has\_fluid\_inside                                                                                                                                                                                                                                                                  |
| `triangles`         | <p><code>num_tri</code> - number of triangles, up to 5 (change the value in <code>custom_geometry.cpp</code> if needed),<br>for each triangle, <code>{i}</code> from 0 to num_tri-1, <code>tri_{i}_point_0</code>, <code>tri_{i}_point_1</code>, <code>tri_{i}_point_2</code> - three points that define the triangle, give the points in anti-clockwise direction to set solid inside of the triangle. The z-coordinate isn't really needed because the triangle will be extruded in the z-direction.</p> |
| `combustor`         | `far_wall_loc`, `ramp_plane1_point`, `ramp_plane2_point`, `ramp_plane2_normal`, `ramp_plane3_point`, `pipe_lo`, `pipe_hi`                                                                                                                                                                                                                                                                                                                                                                                  |
| `converging-nozzle` | `d_inlet`, `l_inlet`, `d_exit`, `l_nozzle`                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

NOTE: You can only choose one geometry type and one geometry. If you want to use multiple geometries, you need to define your own geometry

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
