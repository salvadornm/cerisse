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
#include "bc_util.H"
#include "index_macros.H"
#include "mechanism.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
  int i, int j, int k, amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata, Parm const& /*parm*/, ProbParm const& pp)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Mean flow
  const amrex::Real u_mean = 1.5 * pp.u_b * (1 - y * y / pp.H / pp.H);

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real p, rho, ei;
  // eos.RTY2P(pp.rho_w, pp.T_w, (amrex::Real*) pp.massfrac.begin(), p);
  // amrex::Real T = pp.T_w + 0.3*pp.T_w*(1 - y*y/pp.H/pp.H);
  // eos.PYT2RE(p, (amrex::Real*) pp.massfrac.begin(), T, rho, ei);
  rho = pp.rho_w;
  eos.RTY2E(rho, pp.T_w, (amrex::Real*)pp.massfrac.begin(), ei);

  // Turbulent fluctuations (Hamzehloo et al)
  amrex::Real M_turb = 0.1;
  const amrex::Real Lx = 12 * pp.H, Ly = 2 * pp.H, Lz = 4 * pp.H;
  amrex::Real A =
    M_turb / pp.M_b * 1.5 *
    pp.u_b; // 0.1 * (2.5*log(pp.Re_tau * sqrt(pp.rho_b / pp.rho_w)) + 5.5);
  amrex::Real cx = cos(4 * M_PI * x / Lx);
  amrex::Real cy = 1 + cos(2 * M_PI * y / Ly);
  amrex::Real cz = cos(2 * M_PI * z / Lz);
  amrex::Real sx = sin(4 * M_PI * x / Lx);
  amrex::Real sy = sin(2 * M_PI * y / Ly);
  amrex::Real sz = sin(2 * M_PI * z / Lz);
  amrex::Real u_turb = A * cx * sy * sz * (Lx / 2);
  amrex::Real v_turb = -A * sx * cy * sz * (Ly / 2);
  amrex::Real w_turb = -A * sx * sy * cz * (Lz / 2);

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    state(i, j, k, nf * NVAR + URHO) = pp.rho_w;
    state(i, j, k, nf * NVAR + UMX) = pp.rho_w * (u_mean + u_turb);
    state(i, j, k, nf * NVAR + UMY) = pp.rho_w * v_turb;
    state(i, j, k, nf * NVAR + UMZ) = pp.rho_w * w_turb;
    state(i, j, k, nf * NVAR + UEDEN) =
      pp.rho_w * (ei + 0.5 * ((u_mean + u_turb) * (u_mean + u_turb) +
                              v_turb * v_turb + w_turb * w_turb));
    for (int n = 0; n < NUM_SPECIES; ++n) {
      state(i, j, k, nf * NVAR + UFS + n) = pp.rho_w * pp.massfrac[n];
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real x[AMREX_SPACEDIM], const amrex::Real* /*s_int[LEN_STATE]*/,
  const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE],
  const int /*idir*/, const int /*sgn*/, const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/, ProbParm const& pp)
{
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    bc_isothermal_wall(&s_ext[nf * NVAR], &s_refl[nf * NVAR], pp.T_w, 1, true);
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
  amrex::Array4<const amrex::Real> const& sarr, int /*level*/, char tagval,
  const amrex::Real /*time*/, amrex::GeometryData const& geomdata,
  Parm const& /*parm*/, ProbParm const& pp)
{
}

#endif
