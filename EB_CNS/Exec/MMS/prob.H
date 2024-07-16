#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <masa.h>

#include <cmath>

#include "CNS.H"
#include "PelePhysics.H"
#include "index_macros.H"
#include "mechanism.H"
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
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Get state from MASA
  const amrex::Real rho = masa_eval_3d_exact_rho(x, y, z);
  const amrex::Real u[AMREX_SPACEDIM] = {masa_eval_3d_exact_u(x, y, z),
                                         masa_eval_3d_exact_v(x, y, z),
                                         masa_eval_3d_exact_w(x, y, z)};
  const amrex::Real p = masa_eval_3d_exact_p(x, y, z);
  // amrex::Real rho = 0.0, u[AMREX_SPACEDIM] = {0.0}, p = 0.0;
  // std::vector<amrex::Real> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  // std::vector<amrex::Real> points = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
  // for (int i1 = 0; i1 < points.size(); ++i1) {
  // for (int i2 = 0; i2 < points.size(); ++i2) {
  //   amrex::Real n = points[i1];
  //   amrex::Real m = points[i2];
  //   rho += 0.25 * weights[i1] * weights[i2] * masa_eval_3d_exact_rho(x + 0.5 * n * dx[0], y + 0.5 * m * dx[1], z);
  //   u[0] += 0.25 * weights[i1] * weights[i2] * masa_eval_3d_exact_u(x + 0.5 * n * dx[0], y + 0.5 * m * dx[1], z);
  //   u[1] += 0.25 * weights[i1] * weights[i2] * masa_eval_3d_exact_v(x + 0.5 * n * dx[0], y + 0.5 * m * dx[1], z);
  //   u[2] += 0.25 * weights[i1] * weights[i2] * masa_eval_3d_exact_w(x + 0.5 * n * dx[0], y + 0.5 * m * dx[1], z);
  //   p += 0.25 * weights[i1] * weights[i2] * masa_eval_3d_exact_p(x + 0.5 * n * dx[0], y + 0.5 * m * dx[1], z);
  // }
  // }
  amrex::Real massfrac[NUM_SPECIES] = {1.0};

  amrex::Real eint;
  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(rho, massfrac, p, eint);
  
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_restart(
  int i, int j, int k, amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_timestep(
//   int i, int j, int k, const amrex::Real /*curtime*/, const amrex::Real dtlev,
//   amrex::Array4<amrex::Real> const& state,
//   amrex::Array4<const amrex::Real> const& reactions,
//   amrex::GeometryData const& /*geomdata*/, Parm const& /*parm*/, ProbParm const& pp)
// {
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_tag_error(
  int i, int j, int k, amrex::Array4<char> const& tagarr,
  amrex::Array4<const amrex::Real> const& /*sarr*/, int level, char tagval,
  const amrex::Real /*time*/, amrex::GeometryData const& geomdata,
  Parm const& /*parm*/, ProbParm const& /*pp*/)
{
}

#endif
