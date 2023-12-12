# Input options

This page provides detailed info on the input file options.  

## AMReX options

The AMReX options include "max_step", "stop_time", and those preceded by "amrex", "amr", or "geometry". They cover control of 

* the problem domain definition
* time-stepping
* gridding & load balancing
* compute backends
* output files
* checkpoint & restarting

The reference is available on [AMReX's documentation](https://amrex-codes.github.io/amrex/docs_html/Inputs_Chapter.html).

## Cerisse options

The Cerisse-specific options are usually (but not always) preceded by "cns":

### Boundary conditions

```
# ================ BOUNDARY CONDITIONS ================= #
# 0 = Interior             3 = Symmetry
# 1 = Inflow / UserBC      4 = SlipWall =3
# 2 = Outflow (FOExtrap)   5 = NoSlipWall (adiabatic)
cns.lo_bc = 0 0 0
cns.hi_bc = 0 0 0
```

Option            | Type             | Default   | Meaning
:-----------------|:-----------------|:----------|:--------------------------------
cns.lo_bc         | DIM * Integer    |           | BC flags lower boundaries in x,y,z
cns.hi_bc         | DIM * Integer    |           | BC flags upper boundaries in x,y,z

If option "0" is selected, the corresponding `geometry.is_periodic` must also be set.

If option "1" is selected, the `bcnormal` in `prob.H` routine will be activated.

### Numerics

```
# ==================== WHAT NUMERICS =================== #
cns.cfl = 0.3
cns.dt_cutoff = 5.e-20
cns.recon_scheme = 6  # 1: Godunov,  2: MUSCL,   3: WENO-Z3, 
                      # 4: WENO-JS5, 5: WENO-Z5, 6: TENO-5
...
```
Godunov is first-order (piecewise constant) reconstruction, MUSCL is second-order, while the number represents the formal order of the reconstruction scheme.

All reconstructions are carried out in characteristic space, so there is a choice of what variable to use in the Jacobian (speed of sound or gamma). We also implemented a trick to reconstruct the characteristic system variable as well, which can significantly reduce numerical dissipation, but this is not guaranteed to be stable in all cases. Users are advised to try `cns.recon_char_var = 1` first, then switch it off when the case is unstable.

Option             | Type             | Default   | Meaning
:------------------|:-----------------|:---------:|:--------------------------------
cns.cfl            | Real             | 0.3       | Acoustic Courant number
cns.dt_cutoff      | Real             | 5.e-20    | Minimum timestep size allowed
cns.recon_scheme   | Integer (1-6)    | 5         | Reconstruction scheme 
cns.char_sys       | 0 or 1           | 0 (sos)   | System for characteristic variable conversion
cns.recon_char_var | 0 or 1           | 0         | Reconstruct characteristic system variable or not

### Geometry EB

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



## Options




