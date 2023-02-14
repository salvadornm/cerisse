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

  // TG functions
  amrex::Real u[3] = {0.0};
  u[0] = prob_parm.v0 * sin(prob_parm.omega_x * x / prob_parm.L) *
         cos(prob_parm.omega_y * y / prob_parm.L) *
         cos(prob_parm.omega_z * z / prob_parm.L);
  u[1] = -prob_parm.v0 * cos(prob_parm.omega_x * x / prob_parm.L) *
         sin(prob_parm.omega_y * y / prob_parm.L) *
         cos(prob_parm.omega_z * z / prob_parm.L);
  if (prob_parm.convecting) {
    u[0] += prob_parm.v0;
    u[1] += prob_parm.v0;
  }
  const amrex::Real p =
    prob_parm.p0 + prob_parm.rho0 * prob_parm.v0 * prob_parm.v0 / 16.0 *
                     (cos(2.0 * prob_parm.omega_x * x / prob_parm.L) +
                      cos(2.0 * prob_parm.omega_y * y / prob_parm.L)) *
                     (cos(2.0 * prob_parm.omega_z * z / prob_parm.L) + 2.0);
  amrex::Real rho, eint;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, prob_parm.T0, rho, eint);

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u[0];
  state(i, j, k, UMY) = rho * u[1];
  state(i, j, k, UMZ) = rho * u[2];
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
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
cns_prob_post_restart (
  int i, int j, int k, 
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  Parm const& parm,
  ProbParm const& prob_parm)
{
}


#endif
