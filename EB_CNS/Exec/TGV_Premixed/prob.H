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
              amrex::GeometryData const& geomdata, Parm const& /*parm*/,
              ProbParm const& prob_parm)
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

  amrex::Real u[3] = {0.0}, sd = 0.0;
  amrex::Real rho, eint;
  auto eos = pele::physics::PhysicsType::eos();
  if (prob_parm.mode == 0) {
    // TG functions
    u[0] = prob_parm.v0 * sin(2.0 * M_PI * x / prob_parm.L) *
           cos(2.0 * M_PI * y / prob_parm.L) * cos(2.0 * M_PI * z / prob_parm.L);
    u[1] = -prob_parm.v0 * cos(2.0 * M_PI * x / prob_parm.L) *
           sin(2.0 * M_PI * y / prob_parm.L) * cos(2.0 * M_PI * z / prob_parm.L);
    const amrex::Real p =
      prob_parm.p0 + prob_parm.rho0 * prob_parm.v0 * prob_parm.v0 / 16.0 *
                       (cos(2.0 * 2.0 * M_PI * x / prob_parm.L) +
                        cos(2.0 * 2.0 * M_PI * y / prob_parm.L)) *
                       (cos(2.0 * 2.0 * M_PI * z / prob_parm.L) + 2.0);    
    eos.PYT2RE(p, prob_parm.Y.begin(), prob_parm.T0, rho, eint);
  } else {
    // HIT functions
    const int index = i + prob_parm.inres * (j + prob_parm.inres * k);
    u[0] = prob_parm.d_uinput[index];
    u[1] = prob_parm.d_vinput[index];
    u[2] = prob_parm.d_winput[index];
    sd = std::sqrt(Real(2.0/3.0) * prob_parm.d_kinput[index]); // spatially varying
    // sd = std::sqrt(Real(2.0/3.0) * prob_parm.ksgs_avg); // spatially constant
    rho = prob_parm.rho0;
    eos.PYT2RE(prob_parm.p0, prob_parm.Y.begin(), prob_parm.T0, rho, eint);
  }

  // Set the state
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    state(i, j, k, nf * NVAR + URHO) = rho;
    state(i, j, k, nf * NVAR + UMX) = rho * u[0];
    state(i, j, k, nf * NVAR + UMY) = rho * u[1];
    state(i, j, k, nf * NVAR + UMZ) = rho * u[2];
    state(i, j, k, nf * NVAR + UEDEN) =
      rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k, nf * NVAR + UFS + n) = rho * prob_parm.Y[n];
  }

#if NUM_FIELD > 0
  if (prob_parm.mode == 1) {
    // Add subgrid-scale fluctuations

    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      // amrex::Real usgs = amrex::RandomNormal(0.0, sd);
      // amrex::Real vsgs = amrex::RandomNormal(0.0, sd);
      // amrex::Real wsgs = amrex::RandomNormal(0.0, sd);
      // amrex::Real usgs = sd * (2*(nf % 2)-1); // zero mean, unit variance, zero covariance
      // amrex::Real vsgs = sd * (2*(nf <= 4)-1);
      // amrex::Real wsgs = sd * (2*(nf % 4 < 2)-1);
      amrex::Real Tsgs = sd*sd*2.6 * (2*(nf % 2)-1);
      
      // amrex::Real ei = eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) -
      //     0.5 * ((u[0] + usgs) * (u[0] + usgs) + (u[1] + vsgs) * (u[1] + vsgs) +
      //            (u[2] + wsgs) * (u[2] + wsgs));
      // amrex::Real T = prob_parm.T0; // initial guess
      // eos.EY2T(eint, prob_parm.Y.begin(), T);
      // amrex::Real cv;
      // eos.TY2Cv(T, prob_parm.Y.begin(), cv);
      // eos.PYT2R(prob_parm.p0, prob_parm.Y.begin(), T, rho);

      state(i, j, k, nf * NVAR + URHO) = rho;
      // state(i, j, k, nf * NVAR + UMX) = rho * (u[0] + usgs);
      // state(i, j, k, nf * NVAR + UMY) = rho * (u[1] + vsgs);
      // state(i, j, k, nf * NVAR + UMZ) = rho * (u[2] + wsgs);
      state(i, j, k, nf * NVAR + UMX) = rho * (u[0]);
      state(i, j, k, nf * NVAR + UMY) = rho * (u[1]);
      state(i, j, k, nf * NVAR + UMZ) = rho * (u[2]);
      // Strong Reynolds analogy: total termperature fluctuation ~ 0
      state(i, j, k, nf * NVAR + UEDEN) =
        // rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
        rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) + Tsgs);
      for (int n = 0; n < NUM_SPECIES; n++)
        state(i, j, k, nf * NVAR + UFS + n) = rho * prob_parm.Y[n];

      // Add extra temperature fluctuation
      // amrex::Real cv;
      // eos.TY2Cv(T, prob_parm.Y.begin(), cv);
      // amrex::Real Tsgs = amrex::RandomNormal(0.0, 10.0); // 20K 99% confidence
      // ei += cv * Tsgs;
      // eos.EY2T(ei, prob_parm.Y.begin(), T);
      // eos.PYT2R(prob_parm.p0, prob_parm.Y.begin(), T, rho);

      // state(i, j, k, nf * NVAR + URHO) = rho;
      // state(i, j, k, nf * NVAR + UMX) = rho * (u[0] + usgs);
      // state(i, j, k, nf * NVAR + UMY) = rho * (u[1] + vsgs);
      // state(i, j, k, nf * NVAR + UMZ) = rho * (u[2] + wsgs);
      // state(i, j, k, nf * NVAR + UEDEN) =
      //   rho * (eint + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) /*+ cv * Tsgs*/);
      // for (int n = 0; n < NUM_SPECIES; n++)
      //   state(i, j, k, nf * NVAR + UFS + n) = rho * prob_parm.Y[n];
    }
  }
#endif
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
//   amrex::GeometryData const& /*geomdata*/, Parm const& /*parm*/, ProbParm const&
//   pp)
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
