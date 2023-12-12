#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "prob_parm.H"


AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& /*parm*/,
              ProbParm const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];

  // Free-stream  
  state(i, j, k, URHO)  = prob_parm.rho_l;
  state(i, j, k, UMX)   = prob_parm.rho_l * prob_parm.u_l;
  state(i, j, k, UMY)   = 0.0;
  state(i, j, k, UMZ)   = 0.0;
  state(i, j, k, UEDEN) =
    prob_parm.rhoe_l + 0.5 * prob_parm.rho_l * prob_parm.u_l * prob_parm.u_l;
  state(i, j, k, UFS) = prob_parm.rho_l;
}

// Inflow boundary condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], const amrex::Real s_int[LEN_STATE],
         const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/, ProbParm const& prob_parm)
{
  if (idir == 0 && sgn == 1) { 
    s_ext[URHO] = prob_parm.rho_l;
    s_ext[UMX] = prob_parm.rho_l * prob_parm.u_l;
    s_ext[UMY] = 0.0;
    s_ext[UMZ] = 0.0;
    s_ext[UEDEN] =
      prob_parm.rhoe_l + 0.5 * prob_parm.rho_l * prob_parm.u_l * prob_parm.u_l;
    s_ext[UFS] = prob_parm.rho_l;
  }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_restart(
  int i, int j, int k, amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_timestep(
  int i, int j, int k, const amrex::Real /*curtime*/, const amrex::Real dtlev,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<const amrex::Real> const& reactions,
  amrex::GeometryData const& /*geomdata*/, Parm const& /*parm*/, ProbParm const& pp)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_tag_error(
  int i, int j, int k, amrex::Array4<char> const& tagarr,
  amrex::Array4<const amrex::Real> const& /*sarr*/, int level, char tagval,
  const amrex::Real /*time*/, amrex::GeometryData const& geomdata,
  Parm const& /*parm*/, ProbParm const& /*pp*/)
{
}

#endif