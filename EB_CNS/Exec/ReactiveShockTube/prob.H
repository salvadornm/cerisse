#ifndef CNS_PROB_H
#define CNS_PROB_H

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
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
prob_initdata (int i, int j, int k, 
               amrex::Array4<amrex::Real> const& state,
               amrex::GeometryData const& geomdata,
               Parm const& /*parm*/,
               ProbParm const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];

  // Set the states
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    if (x <= 6.0) {      
      state(i, j, k, nf*NVAR + URHO) = prob_parm.rho_l;
      state(i, j, k, nf*NVAR + UMX) = prob_parm.rho_l * prob_parm.u_l;
      state(i, j, k, nf*NVAR + UMY) = 0.0;
      state(i, j, k, nf*NVAR + UMZ) = 0.0;
      state(i, j, k, nf*NVAR + UEDEN) = prob_parm.rho_l * prob_parm.e_l + 0.5 * prob_parm.rho_l * prob_parm.u_l * prob_parm.u_l;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf*NVAR + UFS + n) = prob_parm.rho_l * prob_parm.massfrac_l[n];
      }        
    } else {
      state(i, j, k, nf*NVAR + URHO) = prob_parm.rho_r;
      state(i, j, k, nf*NVAR + UMX) = prob_parm.rho_r * prob_parm.u_r;
      state(i, j, k, nf*NVAR + UMY) = 0.0;
      state(i, j, k, nf*NVAR + UMZ) = 0.0;
      state(i, j, k, nf*NVAR + UEDEN) = prob_parm.rho_r * prob_parm.e_r + 0.5 * prob_parm.rho_r * prob_parm.u_r * prob_parm.u_r;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf*NVAR + UFS + n) = prob_parm.rho_r * prob_parm.massfrac_r[n];
      }    
    }
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
 * @sa CnsFillExtDir
 * @sa CnsFillExtDir::operator()
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal (const amrex::Real x[AMREX_SPACEDIM],
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
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_prob_post_restart (int i, int j, int k, 
                       amrex::Array4<amrex::Real> const& state,
                       amrex::GeometryData const& geomdata,
                       Parm const& parm,
                       ProbParm const& prob_parm)
{
}

#endif