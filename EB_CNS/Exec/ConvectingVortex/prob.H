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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata (
  int i, int j, int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  Parm const& /*parm*/,
  ProbParm const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Vortex functions (VI1)
  const amrex::Real xc = 5.0; const amrex::Real yc = 5.0; // vortex initial pos
  const amrex::Real R = 0.5; // vortex radius

  const amrex::Real rsq = (x - xc)*(x - xc) + (y - yc)*(y - yc);

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real cv, cp;
  eos.TY2Cv(CNS::h_prob_parm->T0, massfrac, cv);
  eos.TY2Cp(CNS::h_prob_parm->T0, massfrac, cp);

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    amrex::Real u[3], T;
    // if (rsq < R*R) {
      u[0] = prob_parm.v0*(1 - prob_parm.beta*(y - yc)/R*exp(-0.5*rsq/R/R));
      u[1] = prob_parm.v0*prob_parm.beta*(x - xc)/R*exp(-0.5*rsq/R/R);
      u[2] = 0.0;      
      T = prob_parm.T0 - 0.5*(prob_parm.v0*prob_parm.beta*prob_parm.v0*prob_parm.beta)/cp*exp(-rsq/R/R);      
    // } else {
    //   u[0] = prob_parm.v0;
    //   u[1] = 0.0;
    //   u[2] = 0.0;
    //   T = prob_parm.T0;
    // }
    
    amrex::Real rho; //= prob_parm.rho0 * pow(T/prob_parm.T0, 1/(eos.gamma-1));
    amrex::Real eint; //= cv * T;    
    eos.PYT2RE(prob_parm.p0, massfrac, T, rho, eint);

    // Set the state  
    state(i, j, k, nf*NVAR + URHO) = rho;
    state(i, j, k, nf*NVAR + UMX) = rho * u[0];
    state(i, j, k, nf*NVAR + UMY) = rho * u[1];
    state(i, j, k, nf*NVAR + UMZ) = rho * u[2];
    state(i, j, k, nf*NVAR + UEDEN) = rho * (eint + 0.5 * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
    state(i, j, k, nf*NVAR + UFS) = rho;
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal (
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  const amrex::Real* /*s_refl[NVAR]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_restart (
  int i, int j, int k, 
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  Parm const& parm,
  ProbParm const& prob_parm)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_timestep (
  int i, int j, int k, 
  const amrex::Real /*curtime*/,
  const amrex::Real dtlev,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<const amrex::Real> const& reactions,
  amrex::GeometryData const& /*geomdata*/,
  Parm const& /*parm*/,
  ProbParm const& pp)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error (int i, int j, int k, amrex::Array4<char> const& tagarr, 
                amrex::Array4<const amrex::Real> const& /*sarr*/, 
                int level, char tagval, const amrex::Real /*time*/,
                amrex::GeometryData const& geomdata,
                Parm const& /*parm*/, ProbParm const& /*pp*/)
{
}

#endif
