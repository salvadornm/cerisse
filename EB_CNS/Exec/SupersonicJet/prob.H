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
  const amrex::Real rsqr = AMREX_D_TERM(, y * y, + z * z);

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    if (rsqr > 0.1905 * 0.1905 && rsqr < 0.889 * 0.889) { // Air inlet
      const amrex::Real profile = 1; //- pow((sqrt(rsqr)-0.53975)/0.34925, 12);
      state(i, j, k, nf * NVAR + URHO) = pp.rho1;
      state(i, j, k, nf * NVAR + UMX) = pp.rho1 * pp.u1 * profile;
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) =
        pp.rho1 * pp.ei1 + 0.5 * pp.rho1 * pp.u1 * pp.u1 * profile * profile;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf * NVAR + UFS + n) = pp.rho1 * pp.Y1[n];
      }
    } else if (rsqr < 0.118 * 0.118) { // Fuel jet
      const amrex::Real profile = 1;   //- pow(rsqr/0.118/0.118, 6);
      state(i, j, k, nf * NVAR + URHO) = pp.rho2;
      state(i, j, k, nf * NVAR + UMX) = pp.rho2 * pp.u2 * profile;
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) =
        pp.rho2 * pp.ei2 + 0.5 * pp.rho2 * pp.u2 * pp.u2 * profile * profile;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf * NVAR + UFS + n) = pp.rho2 * pp.Y2[n];
      }
    } else if (rsqr > 0.889 * 0.889) { // Still air
      state(i, j, k, nf * NVAR + URHO) = pp.rho_a;
      state(i, j, k, nf * NVAR + UMX) = 0.0;
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) = pp.rho_a * pp.ei_a;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf * NVAR + UFS + n) = pp.rho_a * pp.Y_a[n];
      }
    } else { // Between two jets
      const amrex::Real profile =
        0.5 + 0.5 * tanh((rsqr - 0.15425 * 0.15425) / (0.03625 * 0.03625 * 2));
      state(i, j, k, nf * NVAR + URHO) = 0.5 * (pp.rho1 + pp.rho2);
      state(i, j, k, nf * NVAR + UMX) =
        state(i, j, k, nf * NVAR + URHO) * (profile * pp.u1 + (1 - profile) * pp.u2);
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) =
        state(i, j, k, nf * NVAR + URHO) *
          (profile * pp.ei1 + (1 - profile) * pp.ei2) +
        0.5 * state(i, j, k, nf * NVAR + UMX) * state(i, j, k, nf * NVAR + UMX) /
          state(i, j, k, nf * NVAR + URHO);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, nf * NVAR + UFS + n) =
          state(i, j, k, nf * NVAR + URHO) *
          (profile * pp.Y1[n] + (1 - profile) * pp.Y2[n]);
      }
    }
  }

#if NUM_AUX > 0
  // Aux variables for storing time averages
  for (int n = 0; n < NUM_AUX; ++n) { state(i, j, k, UFA + n) = 0.0; }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real x[AMREX_SPACEDIM], const amrex::Real s_int[NUM_FIELD * NVAR],
  const amrex::Real* /*s_refl[NUM_FIELD * NVAR]*/,
  amrex::Real s_ext[NUM_FIELD * NVAR], const int /*idir*/, const int /*sgn*/,
  const amrex::Real time, amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  BL_PROFILE("bcnormal()");

  const amrex::Real rsqr = AMREX_D_TERM(, x[1] * x[1], +x[2] * x[2]);

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    if (rsqr > 0.1905 * 0.1905 && rsqr < 0.889 * 0.889) { // Air inlet
      const amrex::Real profile = 1; //- pow((sqrt(rsqr)-0.53975)/0.34925, 12);
      s_ext[nf * NVAR + URHO] = pp.rho1;
      s_ext[nf * NVAR + UMX] = pp.rho1 * pp.u1 * profile;
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] =
        pp.rho1 * pp.ei1 + 0.5 * pp.rho1 * pp.u1 * pp.u1 * profile * profile;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf * NVAR + UFS + n] = pp.rho1 * pp.Y1[n];
      }
    } else if (rsqr < 0.118 * 0.118) { // Fuel jet
      const amrex::Real profile = 1;   //- pow(rsqr/0.118/0.118, 6);
      s_ext[nf * NVAR + URHO] = pp.rho2;
      s_ext[nf * NVAR + UMX] = pp.rho2 * pp.u2 * profile;
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] =
        pp.rho2 * pp.ei2 + 0.5 * pp.rho2 * pp.u2 * pp.u2 * profile * profile;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf * NVAR + UFS + n] = pp.rho2 * pp.Y2[n];
      }
    } else if (rsqr > 0.889 * 0.889) { // Still air
      s_ext[nf * NVAR + URHO] = pp.rho_a;
      s_ext[nf * NVAR + UMX] = 0.0;
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] = pp.rho_a * pp.ei_a;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf * NVAR + UFS + n] = pp.rho_a * pp.Y_a[n];
      }
    } else { // Between two jets (adiabatic no-slip wall)
      s_ext[nf * NVAR + URHO] = s_int[nf * NVAR + URHO];
      s_ext[nf * NVAR + UMX] =
        -s_int
          [nf * NVAR +
           UMX]; //* amrex::min(amrex::max((time-0.00017)/(0.000165-0.00017), 1.0),
                 //-1.0);
      s_ext[nf * NVAR + UMY] = -s_int[nf * NVAR + UMY];
      s_ext[nf * NVAR + UMZ] = -s_int[nf * NVAR + UMZ];
      s_ext[nf * NVAR + UEDEN] = s_int[nf * NVAR + UEDEN];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[nf * NVAR + UFS + n] = s_int[nf * NVAR + UFS + n];
      }
    }
  }

#if NUM_AUX > 0
  // Aux variables for storing time averages
  for (int n = 0; n < NUM_AUX; ++n) { s_ext[UFA + n] = 0.0; }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_restart(
  int i, int j, int k, amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& pp)
{
  if (pp.clean_aux_on_restart) {
    // Clear time avg
    for (int n = 0; n < NUM_AUX; ++n) { state(i, j, k, UFA + n) = 0.0; }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
               amrex::Array4<const amrex::Real> const& sarr, int level, char tagval,
               const amrex::Real /*time*/, amrex::GeometryData const& geomdata,
               Parm const& /*parm*/, ProbParm const& /*pp*/)
{
  BL_PROFILE("prob_tag_error()");

  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  const amrex::Real rsqr = y * y + z * z;

  // Tag cylinder
  if (level == 0 && rsqr <= 2.0 * 2.0) { tagarr(i, j, k) = tagval; }
  if (level == 1 && rsqr <= 1.7 * 1.7) { tagarr(i, j, k) = tagval; }
  if (level == 2 && rsqr <= 0.5 * 0.5) { tagarr(i, j, k) = tagval; }

  // Tag injector
  if (rsqr <= 0.3 * 0.3 && x < 0.1) { tagarr(i, j, k) = tagval; }

  // Tag flame
  amrex::Real Y_OH = sarr(i, j, k, UFS + OH_ID) / sarr(i, j, k, URHO);
  amrex::Real Y_HO2 = sarr(i, j, k, UFS + HO2_ID) / sarr(i, j, k, URHO);
  if (Y_OH >= 0.02 || Y_HO2 >= 1e-3) { tagarr(i, j, k) = tagval; }

  // Tag magvort in input file
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_timestep(
  int i, int j, int k, const amrex::Real /*curtime*/, const amrex::Real dtlev,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<const amrex::Real> const& reactions,
  amrex::GeometryData const& /*geomdata*/, Parm const& /*parm*/, ProbParm const& pp)
{
#if NUM_AUX > 0
  if (pp.record_statistics) {
    Real rho = state(i, j, k, URHO);
    Real rhoinv = Real(1.0) / rho;
    AMREX_D_TERM(Real vx = state(i, j, k, UMX) * rhoinv;
                 , Real vy = state(i, j, k, UMY) * rhoinv;
                 , Real vz = state(i, j, k, UMZ) * rhoinv;);
    Real ei = state(i, j, k, UEDEN) * rhoinv -
              0.5 * (AMREX_D_TERM(vx * vx, +vy * vy, +vz * vz));
    Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = state(i, j, k, UFS + n) * rhoinv;
    }

    Real T = 0.0;
    Real X[NUM_SPECIES];
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, massfrac, T);
    eos.Y2X(massfrac, X);

    // Mean
    AMREX_D_TERM(state(i, j, k, UFA + 0) += vx * dtlev;    // U
                 , state(i, j, k, UFA + 1) += vy * dtlev;  // V
                 , state(i, j, k, UFA + 2) += vz * dtlev;) // W
    state(i, j, k, UFA + 3) += T * dtlev;                  // T
    state(i, j, k, UFA + 4) += X[H2_ID] * dtlev;           // X(H2)
    state(i, j, k, UFA + 5) += X[OH_ID] * dtlev;           // X(OH)
    state(i, j, k, UFA + 6) += X[H2O_ID] * dtlev;          // X(H2O)
    state(i, j, k, UFA + 7) += X[N2_ID] * dtlev;           // X(N2)
    state(i, j, k, UFA + 8) += X[O2_ID] * dtlev;           // X(O2)
    state(i, j, k, UFA + 9) +=
      reactions(i, j, k, NUM_SPECIES) * dtlev; // heat release rate

    // RMS^2
    state(i, j, k, UFA + 10) += T * T * dtlev;                 // T
    state(i, j, k, UFA + 11) += X[H2_ID] * X[H2_ID] * dtlev;   // X(H2)
    state(i, j, k, UFA + 12) += X[OH_ID] * X[OH_ID] * dtlev;   // X(OH)
    state(i, j, k, UFA + 13) += X[H2O_ID] * X[H2O_ID] * dtlev; // X(H2O)
    state(i, j, k, UFA + 14) += X[N2_ID] * X[N2_ID] * dtlev;   // X(N2)
    state(i, j, k, UFA + 15) += X[O2_ID] * X[O2_ID] * dtlev;   // X(O2)
  }
#endif
}

#if NUM_AUX > 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_get_aux_name(amrex::Vector<std::string>& aux_name)
{
  aux_name.resize(NUM_AUX);

  aux_name[0] = "time_avg_u";
  aux_name[1] = "time_avg_v";
  aux_name[2] = "time_avg_w";
  aux_name[3] = "time_avg_T";
  aux_name[4] = "time_avg_X_H2";
  aux_name[5] = "time_avg_X_OH";
  aux_name[6] = "time_avg_X_H2O";
  aux_name[7] = "time_avg_X_N2";
  aux_name[8] = "time_avg_X_O2";
  aux_name[9] = "time_avg_hrr";
  aux_name[10] = "mean_squared_T";
  aux_name[11] = "mean_squared_X_H2";
  aux_name[12] = "mean_squared_X_OH";
  aux_name[13] = "mean_squared_X_H2O";
  aux_name[14] = "mean_squared_X_N2";
  aux_name[15] = "mean_squared_X_O2";
}
#endif

#endif
