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

