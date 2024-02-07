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
#include "PMF.H"
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
 * @param pmf_data  Premixed flame data, read from pmf dat in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& /*parm*/,
              ProbParm const& prob_parm, pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  using amrex::Real;
  auto eos = pele::physics::PhysicsType::eos();

  // Geometry
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (i + 0.5) * dx[0];

  // PMF data
  const Real x1 = x - prob_parm.x_offset - 0.5 * dx[0];
  const Real x2 = x - prob_parm.x_offset + 0.5 * dx[0];
  amrex::GpuArray<Real, NUM_SPECIES + 4> pmf_vals = {0.0};  
  pele::physics::PMF::pmf(pmf_data, x1, x2, pmf_vals);

  Real T = pmf_vals[0];
  Real u = pmf_vals[1] + prob_parm.u_offset;
  Real X[NUM_SPECIES], Y[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) X[n] = pmf_vals[3 + n];
  eos.X2Y(X, Y);
  Real rho, ei;
  eos.PYT2RE(prob_parm.p0, Y, T, rho, ei);

  // Set state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEDEN) = rho * ei + 0.5 * rho * u * u;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, UFS + n) = rho * Y[n];
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
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 * @param pmf_data  Premixed flame data, read from pmf dat in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real x[AMREX_SPACEDIM], const amrex::Real* /*s_int[LEN_STATE]*/,
  const amrex::Real* /*s_refl[LEN_STATE]*/, amrex::Real s_ext[LEN_STATE], const int idir,
  const int sgn, const amrex::Real time, amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm, pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  using amrex::Real;
  auto eos = pele::physics::PhysicsType::eos();

  // PMF data
  amrex::GpuArray<Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  pele::physics::PMF::pmf(pmf_data, 0.0, 0.0, pmf_vals);

  Real T = pmf_vals[0];
  Real u = pmf_vals[1] + prob_parm.u_offset;
  Real X[NUM_SPECIES], Y[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) X[n] = pmf_vals[3 + n];
  eos.X2Y(X, Y);
  Real rho, ei;
  eos.PYT2RE(prob_parm.p0, Y, T, rho, ei);

  // Set state
  s_ext[URHO] = rho;
  s_ext[UMX] = rho * u;
  s_ext[UMY] = 0.0;
  s_ext[UMZ] = 0.0;
  s_ext[UEDEN] = rho * ei + 0.5 * rho * u * u;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    s_ext[UFS + n] = rho * Y[n];
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