# Input options

This page provides detailed info on the ```input`` file options.  
These arguments are passed to AMReX

## AMReX options

The AMReX options covers control of 

* the problem domain definition
* time-stepping
* gridding & load balancing
* compute backends
* output files
* checkpoint & restarting

The reference is available on [AMReX's documentation](https://amrex-codes.github.io/amrex/docs_html/Inputs_Chapter.html).

NOTE: In the tables below, DIM means the number of dimensions, Int means integer, and Bool means boolean value (0 for False and 1 for True). If the option has no default value, a value must be given by the user.


### Problem definition & time-stepping

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| **max_step**                    | Int           |         | Maximum number of time steps to take                         |
| **stop_time**                   | Real          |         | Maximum time to reach                                        |
| **time_step**                   | Real          |         | dt (base level), higher level time step is based on number of subcycles          |
| **cfl**                         | Real          |         | CFL  (incompatible option with time_step)  |
| **geometry.is_periodic**        | DIM * Int     | 0 0 0   | 1 for true, 0 for false (one value for each coordinate direction) |
| **geometry.coord_sys**          | Int           | 0       | 0 = Cartesian; 1 = Cylindrical; 2 = Spherical (only support Cartesian) |
| **geometry.prob_lo**            | DIM * Real    | 0 0 0   | Low corner of physical domain (physical not index space)     |
| **geometry.prob_hi**            | DIM * Real    |         | High corner of physical domain (physical not index space)    |
| **geometry.prob_extent**        | DIM * Real    |         | Extent of physical domain, choose between this or `prob_hi` |
| **amr.n_cell**                  | DIM * Int     |         | Number of cells at level 0 in each coordinate direction      |

### Gridding & load balancing

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| **amr.max_level**           | Int           | 0       | Maximum level of refinement allowed (0 when single-level)    |
| **amr.ref_ratio**           | Int *(nlev-1) |         | Refeniment ratio per level.  If the number of `ref_ratio` is less than the number of levels - 1, the last entry will be automatically propagated |
| **amr.regrid_int**          | Int           | -1      | How often to regrid (in number of steps). No regridding will occur if set to < 0 |
| amr.max_grid_size[^1]       | Int           | 32      | Maximum number of cells in each grid in all directions       |
| amr.blocking_factor[^1]     | Int           | 8       | Each grid must be divisible by blocking_factor in all directions (must be 1 or power of 2) |
| amr.refine_grid_layout[^1]  | Bool          | 1       | Split grids in half until the number of grids is no less than the number of procs |
| amr.n_error_buf[^1]         | DIM * Int     | 1 1 1   | Buffer in added around tagged cells[^2]                      |
| amr.grid_eff                | Real          | 0.7     | Taget value of the percentage of tagged cells in the grids   |
| amr.loadbalance_level0_int  | Int           | 2       | How often to do load balance (in number of steps). For single level (i.e., amr.max_level=0) only |
| amr.loadbalance_with_workestimates | Bool   | 0       | For multi-level runs, load balance is done during regrid and thus the load balance interval is controlled by `regrid_int` |
| amr.loadbalance_max_fac     | Real          | 1.5     | This controls the change in the maximum number of boxes that can be assigned to an MPI rank in load balancing |

### Outputs & restarting

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| amr.plot_files_output       | Bool          | 1       | Output plotfile or not (redundent because one can set plot_int = -1 to   disable output) |
| amr.plot_file               | String        | plt     | Prefix of plotfile output                                    |
| amr.plot_int                | Int           | -1      | Frequency of plotfile output; if -1 then no plotfiles will be written |
| amr.derive_plot_vars        | Strings       | NONE    | List of derived variables to plot; can use "ALL" or "NONE" to select all or none of the variables. See the full list of derived variables available in Cerisse below. |
| amr.checkpoint_files_output | Bool          | 1       | Same as plot_files_output, but for checkpoint files          |
| amr.check_file              | String        | plt     | Prefix of checkpoint file output                             |
| amr.check_int               | Int           | -1      | Frequency of checkpoint file output; if -1 then no plotfiles will be   written |
| amr.restart                 | String        |         | If present, then the name of checkpoint file to restart from |
| amr.plotfile_on_restart     | Bool          | 0       | Write a plotfile when immediately after restart or not       |


### GPU-related parameters
| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| amrex.the_arena_init_size   | Int           | 3/4 of total device memory | GPU device memory allocated to The_Arena() (in bytes) |

### Other parameters

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| amrex.omp_threads           | String or Int | `system`| `nosmt`: avoid using threads for virtual cores (aka Hyperthreading or SMT), as is default in OpenMP; `system`: use the environment variable `OMP_NUM_THREADS`. For Integer values, `OMP_NUM_THREADS` is ignored. |
| amrex.fpe_trap_invalid      | Bool          | 0       | Produce error when invalid floating-point arthematic is detected. Helpful for debugging |

[^1]: Each direction can be controlled separately, for example, use `blocking_factor_x`, `blocking_factor_y` to set different blocking factors in x and y directions
[^2]: For example, if `n_error_buf` = 3, all cells in the 7x7x7 box from the lower corner (i-3,j-3,k-3) to (i+3,j+3,k+3) will be tagged.

## Cerisse options

The Cerisse-specific options are usually (but not always) preceded by "cns":

### Boundary conditions

```ini
# 0 = Interior                               3 = Symmetry
# 1 = Inflow / UserBC                        4 = SlipWall =3
# 2 = Outflow (First Order Extrapolation)    5 = NoSlipWall (adiabatic)
cns.lo_bc = 1 5 0
cns.hi_bc = 2 5 0
```

In the above example, the bc in x would be inflow (at lower boundary) and outflow (at upper boundary) 

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- | ------- | ------------------------------------------------------------ |
| cns.lo_bc                   | DIM * Int     |         | BC flags at lower boundaries in x,y,z
| cns.hi_bc                   | DIM * Int     |         | BC flags at upper boundaries in x,y,z

If option "0" is selected, the corresponding `geometry.is_periodic` must also be set.

If option "1" is selected, the `bcnormal` function in `prob.H` will be activated.

There are also 3 options to control the BC for embedded boundaries, `cns.eb_no_slip` (default true), `cns.eb_isothermal` (default false), and a real number `cns.eb_wall_temp` if isothermal wall is selected.


### Numerics

```ini
cns.cfl = 0.3
cns.dt_cutoff = 5.e-20
```

Option             | Type             | Default   | Meaning
:------------------|:-----------------|:---------:|:--------------------------------
cns.cfl            | Real             | 0.3       | Acoustic Courant number
cns.dt_cutoff      | Real             | 5.e-20    | Minimum timestep size allowed
cns.fixed_dt       | Real             | 0         | Run with constant dt

