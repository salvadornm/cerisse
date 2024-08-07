#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <cmath>

#include "CNS.H"
#include "PelePhysics.H"
#include "index_macros.H"
#include "mechanism.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, ProbParm const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
#if AMREX_SPACEDIM == 3
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
#else
  const amrex::Real z = 1.0;
#endif

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
  amrex::Real d = 0.05;
  amrex::Real Y1 = 0.5 * (tanh((x + 0.5) / d) - tanh((x - 0.5) / d));
  // amrex::Real Y1 = 0.5 * cos(prob_parm.omega_x * x / prob_parm.L) + 0.5;
  amrex::Real massfrac[2] = {Y1, 1.0 - Y1};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(p, massfrac, prob_parm.T0, rho, eint);

  // Set the state
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    state(i, j, k, nf * NVAR + URHO) = rho;
    state(i, j, k, nf * NVAR + UMX) = rho * u[0];
    state(i, j, k, nf * NVAR + UMY) = rho * u[1];
    state(i, j, k, nf * NVAR + UMZ) = rho * u[2];
    state(i, j, k, nf * NVAR + UEDEN) =
      rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, nf * NVAR + UFS + n) = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/, const amrex::Real* /*s_int[NVAR]*/,
  const amrex::Real* /*s_refl[NVAR]*/, amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/, const int /*sgn*/, const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/, ProbParm const& /*prob_parm*/)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_restart(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata, ProbParm const& prob_parm)
{
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_timestep(
//   int i, int j, int k, const amrex::Real /*curtime*/, const amrex::Real dtlev,
//   amrex::Array4<amrex::Real> const& state,
//   amrex::Array4<const amrex::Real> const& reactions,
//   amrex::GeometryData const& /*geomdata*/,  ProbParm const& pp)
// {
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
               amrex::Array4<const amrex::Real> const& /*sarr*/, int level,
               char tagval, const amrex::Real /*time*/,
               amrex::GeometryData const& geomdata, ProbParm const& /*pp*/)
{
}

#endif
