# Input options

This page provides detailed info on the input file options.  

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
| max_step                    | Int           |         | Maximum number of time steps to take                         |
| stop_time                   | Real          |         | Maximum time to reach                                        |
| geometry.is_periodic        | DIM * Int     | 0 0 0   | 1 for true, 0 for false (one value for each coordinate direction) |
| geometry.coord_sys          | Int           | 0       | 0 = Cartesian; 1 = Cylindrical; 2 = Spherical (only support Cartesian) |
| geometry.prob_lo            | DIM * Real    | 0 0 0   | Low corner of physical domain (physical not index space)     |
| geometry.prob_hi            | DIM * Real    |         | High corner of physical domain (physical not index space)    |
| geometry.prob_extent        | DIM * Real    |         | Extent (length) of physical domain, choose between this or `prob_hi` |
| amr.n_cell                  | DIM * Int     |         | Number of cells at level 0 in each coordinate direction      |

### Gridding & load balancing

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| amr.max_level               | Int           | 0       | Maximum level of refinement allowed (0 when single-level)    |
| amr.ref_ratio               | Ints          |         | Refinment ratio in each level, number of levels - 1 values can be used. If the number of `ref_ratio` supplied is less than the number of levels - 1, the last entry will be automatically propagated |
| amr.regrid_int              | Int           | -1      | How often to regrid (in number of steps). No regridding will occur if set to < 0 |
| amr.max_grid_size[^1]       | Int           | 32      | Maximum number of cells in each grid in all directions       |
| amr.blocking_factor[^1]     | Int           | 8       | Each grid must be divisible by blocking_factor in all directions (must be 1 or power of 2) |
| amr.refine_grid_layout[^1]  | Bool          | 1       | Split grids in half until the number of grids is no less than the number of procs |
| amr.n_error_buf[^1]         | DIM * Int     | 1 1 1   | Buffer in added around tagged cells[^2]                      |
| amr.grid_eff                | Real          | 0.7     | Taget value of the percentage of tagged cells in the grids   |
| amr.loadbalance_level0_int  | Int           | 2       | How often to do load balance (in number of steps). For single level (i.e., amr.max_level=0) only |
| amr.loadbalance_with_workestimates | Bool   | 0       | For multi-level runs, load balance is done during regrid and thus the load balance interval is controlled by `regrid_int` |
| amr.loadbalance_max_fac     | Real          | 1.5     | A parameter to control the load balancing                    |

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

```
# 0 = Interior                               3 = Symmetry
# 1 = Inflow / UserBC                        4 = SlipWall =3
# 2 = Outflow (First Order Extrapolation)    5 = NoSlipWall (adiabatic)
cns.lo_bc = 1 5 0
cns.hi_bc = 2 5 0
```

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- | ------- | ------------------------------------------------------------ |
| cns.lo_bc                   | DIM * Int     |         | BC flags at lower boundaries in x,y,z
| cns.hi_bc                   | DIM * Int     |         | BC flags at upper boundaries in x,y,z

If option "0" is selected, the corresponding `geometry.is_periodic` must also be set.

If option "1" is selected, the `bcnormal` function in `prob.H` will be activated.

There are also 3 options to control the BC for embedded boundaries, `cns.eb_no_slip` (default true), `cns.eb_isothermal` (default false), and a real number `cns.eb_wall_temp` if isothermal wall is selected.

### Physics

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- | ------- | ------------------------------------------------------------ |
| cns.do_visc                 | Bool          | 1       | Compute viscous fluxes or not (note that no-slip walls and isothermal walls require `do_visc`)
| cns.do_ext_src              | Bool          | 1       | Add the external source term specified in `prob.cpp`
| cns.do_react                | Bool          | 0       | Compute chemical reaction source term

See also the "Reaction options" section below.

### Numerics

```
cns.cfl = 0.3
cns.dt_cutoff = 5.e-20
cns.recon_scheme = 6  # 1: Godunov,  2: MUSCL,   3: WENO-Z3, 
                      # 4: WENO-JS5, 5: WENO-Z5, 6: TENO-5
...
```
Godunov is first-order (piecewise constant) reconstruction, MUSCL is second-order, while the number represents the formal order of the reconstruction scheme.

All reconstructions are carried out in characteristic space, so there is a choice of what variable to use in the Jacobian, speed of sound (0) or gamma (1). We also implemented a trick to reconstruct the characteristic system variable as well, which can significantly reduce numerical dissipation, but this is not guaranteed to be stable in all cases. Users are advised to try `cns.recon_char_var = 1` first, then switch it off when the case is unstable.

Option             | Type             | Default   | Meaning
:------------------|:-----------------|:---------:|:--------------------------------
cns.cfl            | Real             | 0.3       | Acoustic Courant number
cns.dt_cutoff      | Real             | 5.e-20    | Minimum timestep size allowed
cns.fixed_dt       | Real             | 0         | Run with constant dt
cns.recon_scheme   | Int (1 to 6)     | 5         | Reconstruction scheme 
cns.char_sys       | Int (0 or 1)     | 0 (sos)   | System for characteristic variable conversion
cns.recon_char_var | Int (0 or 1)     | 0         | Reconstruct characteristic system variable or not
cns.limiter_theta  | Real             | 2.0       | Parameter in MUSCL limiter, between 1 and 2; 1: minmod, 2: van Leer's MC (higher sharper) 
cns.rk_order       | Int (1 or 2)     | 1         | Switch between forward Euler and 2nd-order Runge-Kutta (note that WENO is not stable with forward Euler)
cns.clip_temp      | Real             | [C++ epsilon](https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon) | The `enforce_consistent_state()` routine adds energy to cells below this temperature to enforce non-negative temperature |

### LES models

Two types of LES models are available, traditional eddy-viscosity type and stochastic fields. (In practice, they can be used simultaneously, but why would someone want to do that?)

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| cns.do_les                  | Bool          | 0       | Use eddy viscosity type subgrid-scale model or not           |
| cns.do_pasr                 | Bool          | 0       | Use partially stirred reactor in reaction or not             |
| cns.les_model_name          | String        |         | Avaliable models are "Smagorinsky" or "WALE"                 |
| cns.C_s                     | Real          | 0.1     | Smagorinsky constant                                         |
| cns.C_I                     | Real          | 0.0066  | Yoshizawa constant (not used currently)                      |
| cns.Pr_T                    | Real          | 0.7     | Turbulent Prandtl number = mu_T/(cp*kappa_T)                 |
| cns.Sc_T                    | Real          | 0.7     | Turbulent Schmidt number = mu_T/(rho*D_T)                    |
| cns.Cm                      | Real          | 0.1     | PaSR mixing timescale coefficient                            |

For using stochastic fields, the number of fields must be decided at compile time using the `NUM_FIELD` option in `GNUmakefile`. The following options are available to control the runtime behaviour:

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| cns.do_psgs                 | Bool          | 0       | Apply the pressure correction term
| cns.do_pd_model             | Bool          | 0       | Apply the pressure-dilatation model
| cns.do_vpdf                 | Bool          | 0       | Apply the simplified Langevin model
| cns.do_spdf                 | Bool          | 0       | Apply the IEM model
| cns.do_restart_fields       | Bool          | 0       | To copy the mean field to all fields during restart, useful for restarting from a single-field run

### Derived plot variables

| Name                               | Meaning                                                      |
| ---------------------------------- | ------------------------------------------------------------ |
| temp                               | Temperature                                                  |
| pressure                           | Pressure                                                     |
| eint                               | Internal energy                                              |
| x_velocity, y_velocity, z_velocity | Velocities                                                   |
| MachNumber                         | Mach Number                                                  |
| magvort                            | Vorticity magnitude                                          |
| divu                               | Velocity divergence                                          |
| divrho                             | Density divengence                                           |
| shock_sensor                       | Pressure-based shock sensor                                  |
| cp                                 | Specific heat at constant pressure                           |
| cv                                 | Specific heat at constant volume                             |
| transport_coef                     | Transport coefficients. Note that species diffusivities are multipled by density (rho*D) |
| massfrac                           | Mass fractions                                               |
| molefrac                           | Mole fractions                                               |
| turb_viscosity                     | LES model turbulent viscosity (must run with `cns.do_les = 1`) |
| reynolds_stress                    | Subgrid velcoity variance (must compile with `NUM_FIELD > 0`)|
| var_Y                              | Subgrid species mass fraction variance (must compile with `NUM_FIELD > 0`) |
| var_p                              | Subgrid pressure variance (must compile with `NUM_FIELD > 0`)|

WARNING: By default, Cerisse plots all the conservative variables and reaction source terms, including all fields. This generates a huge plotfile. To disable this behaviour, set `cns.plot_fields = 0`.

Users can toggle the reaction source terms outputs using the `cns.plot_rho_omega` and `cns.update_heat_release` options. Similarly for the cost estimate used for load balancing, use `cns.plot_cost`. One may also want to set `cns.plot_rhoy = 0` to avoid plotting species' partial mass when `massfrac` or `molefrac` are already plotted.

## Geometry EB options

In **input** file

```
##---- GEOMETRY -----
#eb2.geom_type = all_regular

eb2.geom_type = cylinder
eb2.cylinder_direction = 2
eb2.cylinder_radius = 0.25
eb2.cylinder_center = 1.0 2.0 0.0
eb2.cylinder_has_fluid_inside = 0

#eb2.sphere_radius = 0.5
#eb2.sphere_center = 2.0 2.0  2.0
#eb2.sphere_has_fluid_inside = 0

#eb2.geom_type = box
#eb2.box_lo = 0.75  1.75  0.0
#eb2.box_hi = 1.25  2.25  0.0
#eb2.box_has_fluid_inside = 0
```

There are several basic geometries that are available in AMReX that can be easily specified in the input file:
- Box
- Plane
- Cylinder
- Sphere

Check [PeleC](https://amrex-combustion.github.io/PeleC/geometry/EB.html) 
page about EB geometry 



## Reaction options

Four chemical integrators are supplied by PelePhysics:

* "Reactor_Null" - does absolutely nothing
* "Reactor_RK64" - an explicit 6-step 4th-order Runge-Kutta integrator
* "Reactor_Cvode" - an adaptive implicit integrator with dense or matrix-free options
    * "direct_dense"
    * "direct_AJdense"
    * "GMERES"
    * "preGMRES"
* "Reactor_Arkode" - an adaptive explicit Runge-Kutta integrator

| Option                      | Type          | Default | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| cns.chem_integrator         | String        | "Reactor_Null" | Select the chemical integrator

std::string CNS::chem_integrator = ;
pele::physics::transport::TransportParams<pele::physics::PhysicsType::transport_type>
  CNS::trans_parms;
amrex::Vector<std::string> CNS::spec_names;
pele::physics::turbinflow::TurbInflow CNS::turb_inflow;

bool CNS::use_typical_vals_chem = false; // tell chem_integrator typical value of Temp
int CNS::reset_typical_vals_int = 10; // interval to reset the typical value
int CNS::rk_reaction_iter = 0;        // iterate to tightly couple chemistry
Real CNS::min_react_temp =
  300.0; // turn off reactor below this temperature (only work with CVODE)