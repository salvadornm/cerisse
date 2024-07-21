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
              amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    if (y > 0.0) {
      state(i, j, k, nf * NVAR + URHO) = pp.rho1;
      state(i, j, k, nf * NVAR + UMX) = pp.rho1 * pp.u1;
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) =
        pp.rho1 * pp.ei1 + 0.5 * pp.rho1 * pp.u1 * pp.u1;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf * NVAR + UFS + n) = pp.rho1 * pp.massfrac1[n];
      }
    } else {
      state(i, j, k, nf * NVAR + URHO) = pp.rho2;
      state(i, j, k, nf * NVAR + UMX) = pp.rho2 * pp.u2;
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) =
        pp.rho2 * pp.ei2 + 0.5 * pp.rho2 * pp.u2 * pp.u2;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf * NVAR + UFS + n) = pp.rho2 * pp.massfrac2[n];
      }
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], const amrex::Real* /*s_int[NVAR]*/,
         const amrex::Real* /*s_refl[NVAR]*/, amrex::Real s_ext[NVAR],
         const int /*idir*/, const int /*sgn*/, const amrex::Real time,
         amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  amrex::Real y = x[1];
  amrex::Real dy = geomdata.CellSize(1);
  // amrex::Real meanu = 0.5*(ua+ub);
  amrex::Real meanu =
    0.5 * (pp.u1 + pp.u2) + 0.5 * (pp.u1 - pp.u2) * tanh(2 * y / pp.theta_w);
  // amrex::Real ua = 0.5*(pp.u1 + pp.u2) + 0.5*(pp.u1 - pp.u2)*tanh(2*(y -
  // dy/2)/pp.theta_w); amrex::Real ub = 0.5*(pp.u1 + pp.u2) + 0.5*(pp.u1 -
  // pp.u2)*tanh(2*(y + dy/2)/pp.theta_w); amrex::Real stdu =
  // sqrt((ua-ub)*(ua-ub)/12.0); amrex::Real randnum; amrex::Real cf = (NUM_FIELD % 2
  // == 0) ? 1.0 : sqrt((amrex::Real)(NUM_FIELD)/(amrex::Real)(NUM_FIELD - 1)); //
  // correct variance if odd num of fields amrex::Real rndseed = 1.0; //(((int)(time
  // * 1e7) % 10) % 2 == 0) ? 1.0 : -1.0; // get the 10th decimal place, give 1 if
  // divisible by 2 so that the inflow is pesudo-random (change in about every 2
  // coarsest steps)

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    // if (nf == 0) {
    //   randnum = 0.0;
    // } else if (nf == NUM_FIELD) {
    //   randnum = (NUM_FIELD % 2 == 0) ? -1.0 : 0.0;
    // } else {
    //   // randnum = 2*(nf % 2) - 1;
    //   randnum = (2*(nf % 2) - 1) * cf;
    // }
    // } else {
    //   randnum = 2.0417 * (2*(nf-1)/(NUM_FIELD-1) - 1.0); //scaling for 8 fields
    // }

    if (y > 0.0) {
      s_ext[nf * NVAR + URHO] = pp.rho1;
      s_ext[nf * NVAR + UMX] = pp.rho1 * (meanu /*+ stdu * rndseed * randnum*/);
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] =
        pp.rho1 * pp.ei1 + 0.5 *
                             /* pp.rho1 * meanu * meanu;*/ s_ext[nf * NVAR + UMX] *
                             s_ext[nf * NVAR + UMX] / s_ext[nf * NVAR + URHO];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf * NVAR + UFS + n] = pp.rho1 * pp.massfrac1[n];
      }
      s_ext[nf * NVAR + UTEMP] = pp.T1;
    } else {
      s_ext[nf * NVAR + URHO] = pp.rho2;
      s_ext[nf * NVAR + UMX] = pp.rho2 * (meanu /*+ stdu * rndseed * randnum*/);
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] =
        pp.rho2 * pp.ei2 + 0.5 *
                             /*pp.rho2 * meanu * meanu;*/ s_ext[nf * NVAR + UMX] *
                             s_ext[nf * NVAR + UMX] / s_ext[nf * NVAR + URHO];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf * NVAR + UFS + n] = pp.rho2 * pp.massfrac2[n];
      }
      s_ext[nf * NVAR + UTEMP] = pp.T2;
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_restart(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  // Clear time avg
  if (pp.clean_aux_on_restart) {
    for (int n = 0; n < NUM_AUX; ++n) { state(i, j, k, UFA + n) = 0.0; }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_timestep(int i, int j, int k, const amrex::Real /*curtime*/,
                   const amrex::Real dtlev, amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<const amrex::Real> const& reactions,
                   amrex::GeometryData const& /*geomdata*/, ProbParm const& pp)
{
#if NUM_AUX > 0
  if (pp.record_statistics) {
    Real rho = state(i, j, k, URHO);
    Real rhoinv = Real(1.0) / rho;
    Real vz = 0.0;
    AMREX_D_TERM(Real vx = state(i, j, k, UMX) * rhoinv;
                 , Real vy = state(i, j, k, UMY) * rhoinv;
                 , vz = state(i, j, k, UMZ) * rhoinv;);
    Real ei = state(i, j, k, UEDEN) * rhoinv -
              0.5 * (AMREX_D_TERM(vx * vx, +vy * vy, +vz * vz));
    Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = state(i, j, k, UFS + n) * rhoinv;
    }
    Real T;
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, massfrac, T);

    state(i, j, k, UFA + 0) += state(i, j, k, URHO) * dtlev;         // RHO
    state(i, j, k, UFA + 1) += state(i, j, k, UMX) * dtlev;          // RHO*U
    state(i, j, k, UFA + 2) += state(i, j, k, UMY) * dtlev;          // RHO*V
    state(i, j, k, UFA + 3) += state(i, j, k, UMZ) * dtlev;          // RHO*W
    state(i, j, k, UFA + 4) += T * dtlev;                            // T
    state(i, j, k, UFA + 5) += state(i, j, k, UFS + H2_ID) * dtlev;  // RHO*Y_H2
    state(i, j, k, UFA + 6) += state(i, j, k, UFS + OH_ID) * dtlev;  // RHO*Y_OH
    state(i, j, k, UFA + 7) += state(i, j, k, UFS + H2O_ID) * dtlev; // RHO*Y_H2O
    state(i, j, k, UFA + 8) += state(i, j, k, UFS + HO2_ID) * dtlev; // RHO*Y_HO2
    state(i, j, k, UFA + 9) += state(i, j, k, UMX) * vx * dtlev;     // RHO*U*U
    state(i, j, k, UFA + 10) += state(i, j, k, UMY) * vy * dtlev;    // RHO*V*V
    state(i, j, k, UFA + 11) += state(i, j, k, UMZ) * vz * dtlev;    // RHO*W*W
    state(i, j, k, UFA + 12) += T * T * dtlev;                       // T*T
  }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
               amrex::Array4<const amrex::Real> const& /*sarr*/, int level,
               char tagval, const amrex::Real /*time*/,
               amrex::GeometryData const& geomdata, ProbParm const& /*pp*/)
{
}

#if NUM_AUX > 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_get_aux_name(amrex::Vector<std::string>& aux_name)
{
  aux_name.resize(13);
  aux_name[0] = "time_avg_rho";
  aux_name[1] = "time_avg_rhou";
  aux_name[2] = "time_avg_rhov";
  aux_name[3] = "time_avg_rhow";
  aux_name[4] = "time_avg_T";
  aux_name[5] = "time_avg_rhoH2";
  aux_name[6] = "time_avg_rhoOH";
  aux_name[7] = "time_avg_rhoH2O";
  aux_name[8] = "time_avg_rhoHO2";
  aux_name[9] = "time_avg_rhouu";
  aux_name[10] = "time_avg_rhovv";
  aux_name[11] = "time_avg_rhoww";
  aux_name[12] = "time_avg_TT";
}
#endif

#endif
