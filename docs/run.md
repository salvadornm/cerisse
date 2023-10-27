# Running

This page explains how to run the code


## GNUmakefile 

General option, usually modified once.

```
DIM = 1
COMP = gnu
PRECISION = DOUBLE
```

## Make.package

Add problems files


## Simulation set-up

Check examples

### input file

Frequent changes


```
max_step = 100
stop_time = 0.2

geometry.is_periodic = 0 0 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     = 0.    0.     0.
geometry.prob_hi     = 1.    0.25   0.25
amr.n_cell           = 32    8      8
```

Input file options

variable | type | meaning
:----------- |:-------------:| -----------:
max_step       | integer        | number of steps of simulations
stop_time       | real        | physical time to stop simulation
geometry.is_periodic        | 3*integer       | indicates periodic direction
geometry.coord_sys  | integer |   cartesian cylindrical spherical
geometry.prob_lo   |  3*real    | domain lowersize x, y, z


### prob.h

Header file to define problem

```
#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "prob_parm.H"

/**
 * \brief Initialise state data.
 *
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 */
```

### prob.cpp


Start of file

```
#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
}
}
```

### prob_param.H

An additional parameter file can be used, for example


```
#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p_l = 1.0;   // left pressure
  amrex::Real u_l = 0.0;   // left velocity
  amrex::Real rho_l = 1.0; // left density
  amrex::Real rhoe_l;
  amrex::Real T_l;

  amrex::Real p_r = 0.1;     // right pressure
  amrex::Real u_r = 0.0;     // right velocity
  amrex::Real rho_r = 0.125; // right density
  amrex::Real rhoe_r;
  amrex::Real T_r;
};

#endif
```

## Geometry

## Options




