#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <cmath>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "PelePhysics.H"
#include "mechanism.H"

#include "CNS.H"
#include "index_macros.H"
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
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata (int i, int j, int k, 
              amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata,
              Parm const& /*parm*/,
              ProbParm const& prob_parm)
{
    // Set states
    state(i, j, k, URHO) = prob_parm.rho;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
    state(i, j, k, UEDEN) = prob_parm.rho * prob_parm.eint;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, UFS+n) = prob_parm.rho * prob_parm.massfrac[n];
    }
}

/**
 * \brief Fill external boundary conditions for ghost cells.
 * 
 * @param x         ghost cell cooridinates.
 * @param s_int     flow state inside of the domain.
 * @param s_ext     flow state to be filled.
 * @param idir      direction (0: x, 1: y, 2: z).
 * @param sgn       high or low boundary (1: low, -1: high).
 * @param time      time.
 * @param geomdata  domain geometry data.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal (const amrex::Real* /*x[AMREX_SPACEDIM]*/,
          const amrex::Real* /*s_int[LEN_STATE]*/,
          amrex::Real s_ext[LEN_STATE],
          const int idir,
          const int sgn, 
          const amrex::Real /*time*/,
          amrex::GeometryData const& /*geomdata*/,
          ProbParm const& prob_parm)
{
}

/**
 * \brief Modify state data and/or add turbulence to fields after restart.
 * 
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_prob_post_restart (int i, int j, int k, 
                       amrex::Array4<amrex::Real> const& state,
                       amrex::GeometryData const& geomdata,
                       Parm const& parm,
                       ProbParm const& prob_parm)
{
}

#endif