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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Mean flow
  const int q = 4;
  const Real u_mean = Real(q) / Real(q - 1) * pp.u_b * (1.0 - std::pow(y / pp.H, q));

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real p, rho, ei;
  // eos.RTY2P(pp.rho_w, pp.T_w, (amrex::Real*) pp.massfrac.begin(), p);
  // amrex::Real T = pp.T_w + 0.3*pp.T_w*(1 - y*y/pp.H/pp.H);
  // eos.PYT2RE(p, (amrex::Real*) pp.massfrac.begin(), T, rho, ei);
  rho = pp.rho_b;
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
    state(i, j, k, nf * NVAR + URHO) = rho;
    state(i, j, k, nf * NVAR + UMX) = rho * (u_mean + u_turb);
    state(i, j, k, nf * NVAR + UMY) = rho * v_turb;
    state(i, j, k, nf * NVAR + UMZ) = rho * w_turb;
    state(i, j, k, nf * NVAR + UEDEN) =
      rho * (ei + 0.5 * ((u_mean + u_turb) * (u_mean + u_turb) +
                          v_turb * v_turb + w_turb * w_turb));
    for (int n = 0; n < NUM_SPECIES; ++n) {
      state(i, j, k, nf * NVAR + UFS + n) = rho * pp.massfrac[n];
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real x[AMREX_SPACEDIM], const amrex::Real* /*s_int[LEN_STATE]*/,
  const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE],
  const int /*idir*/, const int /*sgn*/, const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/, ProbParm const& pp)
{
  // Mean follows normal BC
  bc_isothermal_wall(&s_ext[0], &s_refl[0], pp.T_w, 1, true);

  // Fields extrap with the same gradient as the mean (FOExtrap of variance)
  for (int nf = 1; nf <= NUM_FIELD; ++nf) {
    int start = nf * NVAR;
    s_ext[start + URHO] = s_refl[start + URHO];
    // u_ext^(nf) = u_ext^(mean) + (u_int^(nf) - u_int^(mean))
    for (int dir = 0; dir < amrex::SpaceDim; ++dir)
      s_ext[start + UMX + dir] =
        s_ext[start + URHO] * (s_ext[UMX + dir] / s_ext[URHO] +
                               s_refl[start + UMX + dir] / s_refl[start + URHO] -
                               s_refl[UMX + dir] / s_refl[URHO]);
    for (int n = 0; n < NUM_SPECIES; ++n)
      s_ext[start + UFS + n] = s_refl[start + UFS + n];
    s_ext[start + UEDEN] = s_refl[start + UEDEN] -
                           0.5 / s_refl[start + URHO] *
                             (s_refl[start + UMX] * s_refl[start + UMX] +
                              s_refl[start + UMY] * s_refl[start + UMY] +
                              s_refl[start + UMZ] * s_refl[start + UMZ]) +
                           0.5 / s_ext[start + URHO] *
                             (s_ext[start + UMX] * s_ext[start + UMX] +
                              s_ext[start + UMY] * s_ext[start + UMY] +
                              s_ext[start + UMZ] * s_ext[start + UMZ]);
    // T_ext^(nf) = T_ext^(mean) + (T_int^(nf) - T_int^(mean))
    s_ext[start + UTEMP] = amrex::max(1e-5, s_ext[UTEMP] + (s_refl[start + UTEMP] - s_refl[UTEMP]));
  }

  // for (int nf = 0; nf <= NUM_FIELD; ++nf) {
  //   bc_isothermal_wall(&s_ext[nf * NVAR], &s_refl[nf * NVAR], pp.T_w, 1, true);
  // }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_restart(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata, ProbParm const& pp)
{ 
  //state(i, j, k, URHO) *= 34000/3000;
  //state(i, j, k, UMX) *= 34000/3000;
  //state(i, j, k, UMY) *= 34000/3000;
  //state(i, j, k, UMZ) = 0.0; //*= 34000/3000;
  //state(i, j, k, UEDEN) *= 34000/3000;
  //state(i, j, k, UFS) *= 34000/3000;

#if NUM_FIELD > 1
  // Initial turbulence
  for (int nf = 1; nf <= NUM_FIELD; ++nf) {
    amrex::Real fac = 1.0 + 0.02 * (amrex::Real(nf - 1) / amrex::Real(NUM_FIELD - 1) - 0.5);
    state(i, j, k, nf * NVAR + UMX) *= fac;
    state(i, j, k, nf * NVAR + UMZ) *= fac;
  }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_timestep(int i, int j, int k, const amrex::Real /*curtime*/,
                   const amrex::Real dtlev, amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<const amrex::Real> const& reactions,
                   amrex::GeometryData const& /*geomdata*/, ProbParm const& pp)
{
#if NUM_AUX == 10
  amrex::Real rho = state(i, j, k, URHO);
  amrex::Real rhoinv = Real(1.0) / rho;
  AMREX_D_TERM(Real u = state(i, j, k, UMX) * rhoinv;
             , Real v = state(i, j, k, UMY) * rhoinv;
             , Real w = state(i, j, k, UMZ) * rhoinv;);
  amrex::Real T = state(i, j, k, UTEMP);

  state(i, j, k, UFA + 0) += rho * dtlev;
  state(i, j, k, UFA + 1) += u * dtlev;
  state(i, j, k, UFA + 2) += v * dtlev;
  state(i, j, k, UFA + 3) += w * dtlev;
  state(i, j, k, UFA + 4) += T * dtlev;
  state(i, j, k, UFA + 5) += u * u * dtlev;
  state(i, j, k, UFA + 6) += v * v * dtlev;
  state(i, j, k, UFA + 7) += w * w * dtlev;
  state(i, j, k, UFA + 8) += u * v * dtlev;
  state(i, j, k, UFA + 9) += T * T * dtlev;
#endif
}

#if NUM_AUX == 10
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_get_aux_name(amrex::Vector<std::string>& aux_name)
{
  aux_name.resize(NUM_AUX);
  aux_name[0] = "sum_rho";
  aux_name[1] = "sum_u";
  aux_name[2] = "sum_v";
  aux_name[3] = "sum_w";
  aux_name[4] = "sum_T";
  aux_name[5] = "sum_uu";
  aux_name[6] = "sum_vv";
  aux_name[7] = "sum_ww";
  aux_name[8] = "sum_uv";
  aux_name[9] = "sum_TT";
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
               amrex::Array4<const amrex::Real> const& sarr, int /*level*/,
               char tagval, const amrex::Real /*time*/,
               amrex::GeometryData const& geomdata, ProbParm const& pp)
{
}

#endif
