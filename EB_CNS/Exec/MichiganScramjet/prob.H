#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <cmath>

#include "CNS.H"
#include "PelePhysics.H"
#include "bc_util.H"
#include "prob_parm.H"

// Declare your own geometry here and define the build method in prob.cpp
#if CNS_USE_EB
#include "custom_geometry.H"
class Scramjet : public CustomGeometry::Register<Scramjet>
{
public:
  static const std::string identifier() { return "Scramjet"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bl_profile(const amrex::Real y,
                                                    const amrex::Real z,
                                                    ProbParm const& pp,
                                                    amrex::Real& rho, amrex::Real& u,
                                                    amrex::Real& ei)
{
  // Make boundary layer
  amrex::Real d_wall = std::min(y, 2.54 - y); // dist from upper and lower walls
  d_wall = std::min(d_wall, 1.905 - std::abs(z)); // dist from left and right walls
  d_wall /= (0.1245 * 2.0); // d_wall/D
  const amrex::Real u_profile = (d_wall > 0.0) ? std::tanh(10.0 * d_wall) : 0.0; // u/u_inf
  u = pp.u * u_profile;

  // amrex::Real T_inf, p_inf;
  // const amrex::Real* Y = pp.Y.begin();
  // auto eos = pele::physics::PhysicsType::eos();
  // ei = pp.ei;
  // eos.RYET2P(pp.rho, Y, ei, T_inf, p_inf); // get T, p_inf

  // amrex::Real T_profile = 1.0 + (0.5 * pp.u * pp.u) / pp.cv_Tinf *
  //                                 (1.0 - u_profile * u_profile); // T/T_inf
  // amrex::Real cv, cv_inf = pp.cv_Tinf / T_inf;
  // eos.TY2Cv(T_inf * T_profile, Y, cv);
  // T_profile *= cv_inf / cv;

  // rho = pp.rho / T_profile;
  // eos.RYP2E(rho, Y, p_inf, ei); // rho, p_inf -> ei
}

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
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real z = 0.0;
  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , z = prob_lo[2] + (k + 0.5) * dx[2];)
  // const amrex::Real rsqr = AMREX_D_TERM(,x * x, + z * z);

  amrex::Real rho = pp.rho;
  amrex::Real u = pp.u;
  amrex::Real ei = pp.ei;
  if ((pp.do_bl || pp.make_bl_on_restart) && x > -5.0) { bl_profile(y, z, pp, rho, u, ei); }

  // Ventilated air
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    state(i, j, k, nf * NVAR + URHO) = rho;
    state(i, j, k, nf * NVAR + UMX) = rho * u * Real(y > -dx[1]);
    state(i, j, k, nf * NVAR + UMY) = 0.0;
    state(i, j, k, nf * NVAR + UMZ) = 0.0;
    state(i, j, k, nf * NVAR + UEDEN) =
      rho * ei + 0.5 * rho * u * u * Real(y > -dx[1]);
    for (int n = 0; n < NUM_SPECIES; ++n)
      state(i, j, k, nf * NVAR + UFS + n) = rho * pp.Y[n];
    state(i, j, k, nf * NVAR + UTEMP) = pp.T;
  }

  for (int na = 0; na < NUM_AUX; ++na) { state(i, j, k, UFA + na) = 0.0; }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real x[AMREX_SPACEDIM], const amrex::Real* /*s_int[LEN_STATE]*/,
  const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE], const int idir,
  const int /*sgn*/, const amrex::Real t, amrex::GeometryData const& geomdata,
  ProbParm const& pp)
{
  const bool is_short_domain = geomdata.ProbLo(0) > -35.0;
  
  if (idir == 0) {
    // Inflow / Ventilated air
    amrex::Real rho = pp.rho;
    amrex::Real u = pp.u;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;
    amrex::Real ei = pp.ei;
    amrex::Real T = pp.T;
    if (is_short_domain) {
      // Impose inflow velocity and some turbulence
      amrex::Real d_wall = std::min(x[1], 2.54 - x[1]); // dist from upper and lower walls
      d_wall = std::min(d_wall, 1.905 - std::abs(x[2])); // dist from left and right walls
      d_wall = d_wall > 0 ? d_wall : 1e-10;
      u = std::max((48000 + 4000 * std::log(d_wall)) * (1.0 - std::exp(-d_wall / 0.025)), 0.0);
      T = 996.0 * (1.0 - std::exp(-d_wall / 0.005)) + 300.0;
      auto eos = pele::physics::PhysicsType::eos();
      amrex::Real Y[NUM_SPECIES], p = 2700000;
      for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = pp.Y[n]; }
      eos.PYT2RE(p, Y, T, rho, ei);

      // turbulence in time, pseudo-random no correlation
      amrex::Real uf = 0.0, vf = 0.0, wf = 0.0;
      // constexpr amrex::Real freq[15] = {7, 11, 13, 38, 92, 5, 11, 13, 39, 92, 3, 7, 11, 38, 43};
      // constexpr amrex::Real shift[15] = {0, 0.3, 0.7, 1.1, 0.4, 0.3, 0, 1.1, 0.7, -0.4, 0.1, 0.2, 1.3, -0.4, 1.0};
      // for (int i = 0; i < 5; ++i) {
      //   uf += std::sin(freq[i] * 1e3 * t + shift[i] * M_PI);
      //   vf += std::sin(freq[i + 5] * 1e3 * t + shift[i + 5] * M_PI);
      //   wf += std::sin(freq[i + 10] * 1e3 * t + shift[i + 10] * M_PI);
      // }
      constexpr amrex::Real wy[11] = {0, 0, 1, 1, 1, 2, 2, 3, -2, 2, 3};
      constexpr amrex::Real wz[11] = {1, 2, 0, -1, 2, 0, 2, 0, 1, 3, 4};
      constexpr amrex::Real wt[11] = {8, -14, 6, 3, 40, 4, -50, 7, 31, -4, 13};
      constexpr amrex::Real mag[11] = {0.697, 0.513, 0.808, 0.256, 0.241, 0.282, 0.209, 0.303, 0.217, 0.306, 0.362};
      constexpr amrex::Real ang[11] = {0.733, 0.056, 0.487, 0.143, 0.877, 0.003, 0.24,  0.787, 0.829, 0.976, 0.061};
      for (int i = 0; i < 11; ++i) {
        uf += mag[i] * std::cos(2 * M_PI * (wy[i] * x[1] / 2.54 + wz[i] * x[2] / 3.81 + wt[i] * t / 4e-4 + ang[i]));
      }
      u += 8000.0 * uf;
      v += 5000.0 * vf;
      w += 5000.0 * wf;
    } else if (pp.do_bl) {
      bl_profile(x[1], (amrex::SpaceDim == 3 ? x[2] : Real(0.0)), pp, rho, u, ei);
      auto eos = pele::physics::PhysicsType::eos();
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = pp.Y[n]; }
      eos.REY2T(rho, ei, Y, T);
    }
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      s_ext[nf * NVAR + URHO] = rho;
      s_ext[nf * NVAR + UMX] = rho * u;
      s_ext[nf * NVAR + UMY] = rho * v;
      s_ext[nf * NVAR + UMZ] = rho * w;
      s_ext[nf * NVAR + UEDEN] = rho * ei + 0.5 * rho * (u * u + v * v + w * w);
      for (int n = 0; n < NUM_SPECIES; ++n)
        s_ext[nf * NVAR + UFS + n] = rho * pp.Y[n];
      s_ext[nf * NVAR + UTEMP] = T;      
    }
  } else if (idir == 1) {
    const amrex::Real rsqr = AMREX_D_TERM(, x[0] * x[0], +x[2] * x[2]);
    if (rsqr > 0.1245 * 0.1245) {
      // Cavity wall
      const amrex::Real rsqr_spark =
        AMREX_D_TERM(, (x[0] - 6.99) * (x[0] - 6.99), +x[2] * x[2]);
      if (pp.spark && rsqr_spark < 1.0) {
        for (int nf = 0; nf <= NUM_FIELD; ++nf)
          bc_isothermal_wall(&s_ext[nf * NVAR], &s_refl[nf * NVAR], 2000.0, idir);
      } else {
        // Adiabatic no-slip wall
        for (int nf = 0; nf <= NUM_FIELD; ++nf)
          bc_adiabatic_wall(&s_ext[nf * NVAR], &s_refl[nf * NVAR], idir);
      }
    } else {
      // Fuel
      // amrex::Real damping = time > 2e-4 ? 1.0 : 1.0 - exp(-5.0 * time / 2e-4);
      amrex::Real damping = 1.0;
      amrex::Real jet_profile = 1.0;
      // if (pp.do_bl) {
      //   amrex::Real r_over_D = sqrt(rsqr) / (0.1245 * 2.0);
      //   jet_profile = r_over_D < 0.35 ? 1.0 : 1.0 - pow(6.67*(r_over_D-0.35), 2.5);
      // }
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        s_ext[nf * NVAR + URHO] = pp.rho_j;
        s_ext[nf * NVAR + UMX] = 0.0;
        s_ext[nf * NVAR + UMY] = pp.rho_j * pp.v_j * damping * jet_profile;
        s_ext[nf * NVAR + UMZ] = 0.0;
        s_ext[nf * NVAR + UEDEN] =
          pp.rho_j * pp.ei_j +
          0.5 * s_ext[nf * NVAR + UMY] * s_ext[nf * NVAR + UMY] / pp.rho_j;
        for (int n = 0; n < NUM_SPECIES; ++n)
          s_ext[nf * NVAR + UFS + n] = pp.rho_j * pp.Y_jet[n];
        s_ext[nf * NVAR + UTEMP] = pp.T_j;
      }
    }
  }

  for (int na = 0; na < NUM_AUX; ++na) { s_ext[UFA + na] = 0.0; }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_restart(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata, ProbParm const& pp)
{
#if NUM_AUX > 0
  if (pp.clean_aux_on_restart) { // Clear time avg
    for (int n = 0; n < NUM_AUX; ++n) { state(i, j, k, UFA + n) = 0.0; }
  }
#endif

  if (pp.do_bl && pp.make_bl_on_restart) {
    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    amrex::Real z = 0.0;
    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                 , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                 , z = prob_lo[2] + (k + 0.5) * dx[2];)
    const amrex::Real rsqr = AMREX_D_TERM(, x * x, +z * z);

    if ((x < -0.5 + 0.2 * y * y) && (y > 0.0)) {
      // Make boundary layer
      amrex::Real rho, u, ei;
      bl_profile(y, z, pp, rho, u, ei);

      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        state(i, j, k, nf * NVAR + URHO) = rho;
        state(i, j, k, nf * NVAR + UMX) = rho * u;
        state(i, j, k, nf * NVAR + UMY) = 0.0;
        state(i, j, k, nf * NVAR + UMZ) = 0.0;
        state(i, j, k, nf * NVAR + UEDEN) = rho * ei + 0.5 * rho * u * u;
        for (int n = 0; n < NUM_SPECIES; ++n)
          state(i, j, k, nf * NVAR + UFS + n) = rho * pp.Y[n];
      }
      // } else if ((rsqr <= 0.1245 * 0.1245) && (y < -0.5)) {
      //   // Make jet profile
      //   amrex::Real r_over_D = sqrt(rsqr) / (0.1245 * 2.0);
      //   amrex::Real jet_profile = r_over_D < 0.35 ? 1.0 : 1.0 - pow(6.67*(r_over_D-0.35), 2.5); 
      //   for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      //     state(i, j, k, nf * NVAR + URHO) = pp.rho_j;
      //     state(i, j, k, nf * NVAR + UMX) = 0.0;
      //     state(i, j, k, nf * NVAR + UMY) = pp.rho_j * pp.v_j * jet_profile;
      //     state(i, j, k, nf * NVAR + UMZ) = 0.0;
      //     state(i, j, k, nf * NVAR + UEDEN) = pp.rho_j * pp.ei_j + 0.5 * pp.rho_j * pp.v_j * pp.v_j * jet_profile * jet_profile; 
      //     for (int n = 0; n < NUM_SPECIES; ++n)
      //       state(i, j, k, nf * NVAR + UFS + n) = pp.rho_j * pp.Y_jet[n];
      //   }
    }
  }

  if (pp.make_init_on_restart) {
    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

    if (x < pp.x_reset) {
      // Reset to init states (without BL)
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        state(i, j, k, nf * NVAR + URHO) = pp.rho;
        state(i, j, k, nf * NVAR + UMX) = pp.rho * pp.u * Real(y > -0.5 * dx[1]);
        state(i, j, k, nf * NVAR + UMY) = 0.0;
        state(i, j, k, nf * NVAR + UMZ) = 0.0;
        state(i, j, k, nf * NVAR + UEDEN) =
          pp.rho * pp.ei + 0.5 * pp.rho * pp.u * pp.u * Real(y > -0.5 * dx[1]);
        for (int n = 0; n < NUM_SPECIES; ++n)
          state(i, j, k, nf * NVAR + UFS + n) = pp.rho * pp.Y[n];
      }
    }
  }
  
  // Lower T
  // amrex::Real rho = state(i, j, k, URHO);
  // amrex::Real q2 = (state(i, j, k, UMX) * state(i, j, k, UMX) + state(i, j, k, UMY) * state(i, j, k, UMY) + state(i, j, k, UMZ) * state(i, j, k, UMZ)) / rho / rho;
  // amrex::Real ei = state(i, j, k, UEDEN) / rho - 0.5 * q2;
  // amrex::Real Y[NUM_SPECIES];
  // for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = state(i, j, k, UFS + n) / rho;
  // amrex::Real T, p;
  // auto eos = pele::physics::PhysicsType::eos();
  // eos.RYET2P(rho, Y, ei, T, p);
  // // AMREX_ALWAYS_ASSERT(T == state(i, j, k, UTEMP));
  // T += 20;
  // eos.PYT2RE(p, Y, T, rho, ei);
  // for (int n = 0; n < NUM_SPECIES; ++n) state(i, j, k, UFS + n) = rho * Y[n];
  // state(i, j, k, UMX) *= rho / state(i, j, k, URHO);
  // state(i, j, k, UMY) *= rho / state(i, j, k, URHO);
  // state(i, j, k, UMZ) *= rho / state(i, j, k, URHO);
  // state(i, j, k, URHO) = rho;
  // state(i, j, k, UTEMP) = T;
  // state(i, j, k, UEDEN) = rho * ei + 0.5 * rho * q2;

  // Remove shock train
  // const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  // if ((x < -0.5 + 0.2 * y * y) && (y > 0.0)) {
  //   // Make boundary layer
  //   amrex::Real rho = pp.rho, u = pp.u, ei = pp.ei;
  //   for (int nf = 0; nf <= NUM_FIELD; ++nf) {
  //     state(i, j, k, nf * NVAR + URHO) = rho;
  //     state(i, j, k, nf * NVAR + UMX) = rho * u;
  //     state(i, j, k, nf * NVAR + UMY) = 0.0;
  //     state(i, j, k, nf * NVAR + UMZ) = 0.0;
  //     state(i, j, k, nf * NVAR + UEDEN) = rho * ei + 0.5 * rho * u * u;
  //     for (int n = 0; n < NUM_SPECIES; ++n)
  //       state(i, j, k, nf * NVAR + UFS + n) = rho * pp.Y[n];
  //   }
  // }
  //
  // Lower p
  // const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // if ((x > -23.5) && (x < -3.9)) {
  //   amrex::Real rho = state(i, j, k, URHO);
  //   amrex::Real q2 = (state(i, j, k, UMX) * state(i, j, k, UMX) + state(i, j, k, UMY) * state(i, j, k, UMY) + state(i, j, k, UMZ) * state(i, j, k, UMZ)) / rho / rho;
  //   amrex::Real ei = state(i, j, k, UEDEN) / rho - 0.5 * q2;
  //   amrex::Real Y[NUM_SPECIES];
  //   for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = state(i, j, k, UFS + n) / rho;
  //   amrex::Real T, p;
  //   auto eos = pele::physics::PhysicsType::eos();
  //   eos.RYET2P(rho, Y, ei, T, p);
  //   p = 669339.8533007 + (2443667.4816626 - 669339.8533007) / (-3.79449173 + 23.52002364) * (x + 23.52002364);
  //   eos.PYT2RE(p, Y, T, rho, ei);
  //   for (int n = 0; n < NUM_SPECIES; ++n) state(i, j, k, UFS + n) = rho * Y[n];
  //   state(i, j, k, UMX) *= rho / state(i, j, k, URHO);
  //   state(i, j, k, UMY) *= rho / state(i, j, k, URHO);
  //   state(i, j, k, UMZ) *= rho / state(i, j, k, URHO);
  //   state(i, j, k, URHO) = rho;
  //   state(i, j, k, UTEMP) = T;
  //   state(i, j, k, UEDEN) = rho * ei + 0.5 * rho * q2;
  // }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_timestep(int i, int j, int k, const amrex::Real /*curtime*/,
                   const amrex::Real dtlev, amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<const amrex::Real> const& reactions,
                   amrex::GeometryData const& /*geomdata*/, ProbParm const& pp)
{
#if NUM_AUX > 0
  if (pp.record_statistics) { // Record time averages (mean field)
    amrex::Real rho = state(i, j, k, URHO);
    amrex::Real rhoinv = Real(1.0) / rho;
    AMREX_D_TERM(Real vx = state(i, j, k, UMX) * rhoinv;
                 , Real vy = state(i, j, k, UMY) * rhoinv;
                 , Real vz = state(i, j, k, UMZ) * rhoinv;);
    amrex::Real ke = 0.5 * (AMREX_D_TERM(vx * vx, +vy * vy, +vz * vz));
    amrex::Real ei = state(i, j, k, UEDEN) * rhoinv - ke;
    amrex::Real Y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = state(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real T = state(i, j, k, UTEMP), p;
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, Y, T);
    eos.RTY2P(rho, T, Y, p);

    state(i, j, k, UFA + 0) += rho * dtlev; // rho
    state(i, j, k, UFA + 1) += state(i, j, k, UMX) * dtlev;  // rhou
    state(i, j, k, UFA + 2) += state(i, j, k, UMY) * dtlev;  // rhov
    state(i, j, k, UFA + 3) += state(i, j, k, UMZ) * dtlev;  // rhow
    state(i, j, k, UFA + 4) += T * dtlev;  // T
    state(i, j, k, UFA + 5) += p * dtlev;  // p
    state(i, j, k, UFA + 6) += state(i, j, k, UFS + H2_ID) * dtlev;  // rhoY_H2
    state(i, j, k, UFA + 7) += state(i, j, k, UFS + OH_ID) * dtlev;  // Y_OH
    state(i, j, k, UFA + 8) += state(i, j, k, UFS + O2_ID) * dtlev;  // Y_O2
    state(i, j, k, UFA + 9) += state(i, j, k, UFS + HO2_ID) * dtlev; // Y_HO2
    state(i, j, k, UFA + 10) += reactions(i, j, k, NUM_SPECIES) * dtlev;  // heat release rate
    state(i, j, k, UFA + 11) += p * p * dtlev;  // p^2
    state(i, j, k, UFA + 12) += 2 * rho * ke * dtlev; // 2*k
  }
#endif
}

#if NUM_AUX > 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_get_aux_name(amrex::Vector<std::string>& aux_name)
{
  aux_name.resize(NUM_AUX);

  aux_name[0] = "mean_rho";
  aux_name[1] = "mean_rhou";
  aux_name[2] = "mean_rhov";
  aux_name[3] = "mean_rhow";
  aux_name[4] = "mean_T";
  aux_name[5] = "mean_p";
  aux_name[6] = "mean_rhoY_H2";
  aux_name[7] = "mean_rhoY_OH";
  aux_name[8] = "mean_rhoY_O2";
  aux_name[9] = "mean_rhoY_HO2";
  aux_name[10] = "mean_hrr";
  aux_name[11] = "mean_squared_p";
  aux_name[12] = "mean_squared_ke";
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
               amrex::Array4<const amrex::Real> const& sarr, int level,
               char tagval, const amrex::Real time,
               amrex::GeometryData const& geomdata, ProbParm const& /*pp*/)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  amrex::Real d_wall = std::min(y, 2.54 - y); // dist from upper and lower walls
  d_wall = std::min(d_wall, 1.905 - std::abs(z)); // dist from left and right walls
  if (level < 2 && d_wall > 0 && d_wall < 0.1245 && x < 0.5) { tagarr(i, j, k) = tagval; }
  d_wall = std::min(y, 2.54 - y); // dist from upper and lower walls
  if (level >= 2 && d_wall > 0 && d_wall < 0.1245 && x < -20) { tagarr(i, j, k) = tagval; }
}

#endif
