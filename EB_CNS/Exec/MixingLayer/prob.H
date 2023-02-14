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
  ProbParm const& pp)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    if (y > 0.0) {
      state(i, j, k, nf*NVAR + URHO) = pp.rho1;
      state(i, j, k, nf*NVAR + UMX) = pp.rho1 * pp.u1;
      state(i, j, k, nf*NVAR + UMY) = 0.0;
      state(i, j, k, nf*NVAR + UMZ) = 0.0;
      state(i, j, k, nf*NVAR + UEDEN) = pp.rho1 * pp.ei1 + 0.5 * pp.rho1 * pp.u1 * pp.u1;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf*NVAR + UFS+n) = pp.rho1 * pp.massfrac1[n];
      }
    } else {
      state(i, j, k, nf*NVAR + URHO) = pp.rho2;
      state(i, j, k, nf*NVAR + UMX) = pp.rho2 * pp.u2;
      state(i, j, k, nf*NVAR + UMY) = 0.0;
      state(i, j, k, nf*NVAR + UMZ) = 0.0;
      state(i, j, k, nf*NVAR + UEDEN) = pp.rho2 * pp.ei2 + 0.5 * pp.rho2 * pp.u2 * pp.u2;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf*NVAR + UFS+n) = pp.rho2 * pp.massfrac2[n];
      }
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal (
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real* /*s_int[NVAR]*/,
  const amrex::Real* /*s_refl[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& pp)
{
  amrex::Real y = x[1];
  amrex::Real u = 0.5*(pp.u1 + pp.u2) + 0.5*(pp.u1 - pp.u2)*tanh(2*y/pp.theta_w);

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    if (y > 0.0) {
      s_ext[nf*NVAR + URHO] = pp.rho1;
      s_ext[nf*NVAR + UMX] = pp.rho1 * u;
      s_ext[nf*NVAR + UMY] = 0.0;
      s_ext[nf*NVAR + UMZ] = 0.0;
      s_ext[nf*NVAR + UEDEN] = pp.rho1 * pp.ei1 + 0.5 * pp.rho1 * u * u;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf*NVAR + UFS+n] = pp.rho1 * pp.massfrac1[n];
      }
    } else {
      s_ext[nf*NVAR + URHO] = pp.rho2;
      s_ext[nf*NVAR + UMX] = pp.rho2 * u;
      s_ext[nf*NVAR + UMY] = 0.0;
      s_ext[nf*NVAR + UMZ] = 0.0;
      s_ext[nf*NVAR + UEDEN] = pp.rho2 * pp.ei2 + 0.5 * pp.rho2 * u * u;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf*NVAR + UFS+n] = pp.rho2 * pp.massfrac2[n];
      }
    }
  }
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
