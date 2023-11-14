#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "PelePhysics.H"
#include "prob_parm.H"

#include "bc_util.H"

// Declare your own geometry here and define the build method in prob.cpp
#include "custom_geometry.H"
class Scramjet : public CustomGeometry::Register<Scramjet>
{
public:
  static const std::string identifier() { return "Scramjet"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};

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
              amrex::GeometryData const& geomdata, Parm const& /*parm*/,
              ProbParm const& pp)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
              , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
              , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];)
  const amrex::Real rsqr = AMREX_D_TERM(,x * x, + z * z);

  // if (rsqr > 0.1245 * 0.1245 || y > 0.0) {
    // Ventilated air
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      state(i, j, k, nf * NVAR + URHO) = pp.rho;
      state(i, j, k, nf * NVAR + UMX) = pp.rho * pp.u * Real(y > -dx[1]);
      state(i, j, k, nf * NVAR + UMY) = 0.0;
      state(i, j, k, nf * NVAR + UMZ) = 0.0;
      state(i, j, k, nf * NVAR + UEDEN) = pp.rho * pp.ei + 0.5 * pp.rho * pp.u * pp.u * Real(y > -dx[1]);
      for (int n = 0; n < NUM_SPECIES; ++n)
        state(i, j, k, nf * NVAR + UFS + n) = pp.rho * pp.Y[n];
    }
  // } else {
  //   // Fuel
  //   for (int nf = 0; nf <= NUM_FIELD; ++nf) {
  //     state(i, j, k, nf * NVAR + URHO) = pp.rho_j;
  //     state(i, j, k, nf * NVAR + UMX) = 0.0;
  //     state(i, j, k, nf * NVAR + UMY) = pp.rho_j * pp.v_j;
  //     state(i, j, k, nf * NVAR + UMZ) = 0.0;
  //     state(i, j, k, nf * NVAR + UEDEN) = pp.rho_j * pp.ei_j + 0.5 * pp.rho_j * pp.v_j * pp.v_j;
  //     for (int n = 0; n < NUM_SPECIES; ++n)
  //       state(i, j, k, nf * NVAR + UFS + n) = pp.rho_j * pp.Y_jet[n];
  //   }
  // }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], const amrex::Real s_int[LEN_STATE],
         const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& geomdata, ProbParm const& pp)
{
  if (idir == 0) {
    // Inflow / Ventilated air
    amrex::Real rho = pp.rho;
    amrex::Real u = pp.u;
    amrex::Real ei = pp.ei;
    if (pp.do_bl) { // TODO!
      amrex::Real y = amrex::min(x[1], 2.54-x[1]); // y on upper and lower walls
      y = amrex::min(y, 1.905-amrex::Math::abs(x[2])); // z on upper and lower walls
      y = amrex::max(y, 0.0); // prevent solver from complaining
      y /= (0.1245 * 2.0); // y/D      
      amrex::Real bl_profile = bl_profile = y > 0.68 ? 1.0 : 1.0085656 - exp(-7.0 * y);

      amrex::Real p, T = 0;
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = pp.Y[n]; }
      auto eos = pele::physics::PhysicsType::eos();
      eos.RYET2P(pp.rho, Y, ei, T, p);
      eos.RYP2E(rho, Y, p, ei);

      rho = pp.rho / (1.0 + 0.5 * (pp.u * pp.u) / pp.cv_Tinf * (1.0 - bl_profile * bl_profile));
      u = pp.u * bl_profile;
    }
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      s_ext[nf * NVAR + URHO] = rho;
      s_ext[nf * NVAR + UMX] = rho * u;
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] = rho * ei + 0.5 * rho * pp.u * pp.u;
      for (int n = 0; n < NUM_SPECIES; ++n)
        s_ext[nf * NVAR + UFS + n] = rho * pp.Y[n];
    }
  } else if (idir == 1) {
    const amrex::Real rsqr = AMREX_D_TERM(,x[0] * x[0], + x[2] * x[2]);
    if (rsqr > 0.1245 * 0.1245) {
      const amrex::Real rsqr_spark = AMREX_D_TERM(,(x[0] - 6.99) * (x[0] - 6.99), + x[2] * x[2]);
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
      if (pp.do_bl) {
        amrex::Real r_over_D = sqrt(rsqr) / (0.1245 * 2.0);
        jet_profile = r_over_D < 0.35 ? 1.0 : 1.0 - pow(6.67*(r_over_D-0.35), 2.5);
      }
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        s_ext[nf * NVAR + URHO] = pp.rho_j;
        s_ext[nf * NVAR + UMX] = 0.0;
        s_ext[nf * NVAR + UMY] = pp.rho_j * pp.v_j * damping * jet_profile;
        s_ext[nf * NVAR + UMZ] = 0.0;
        s_ext[nf * NVAR + UEDEN] = pp.rho_j * pp.ei_j + 0.5 * s_ext[nf * NVAR + UMY] * s_ext[nf * NVAR + UMY] / pp.rho_j;
        for (int n = 0; n < NUM_SPECIES; ++n)
          s_ext[nf * NVAR + UFS + n] = pp.rho_j * pp.Y_jet[n];
      }
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_restart(
  int i, int j, int k, amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& pp)
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
    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];)
    const amrex::Real rsqr = AMREX_D_TERM(,x * x, + z * z);

    if ((x < -0.5 + 0.2 * y * y) && (y > 0.0)) {
      // Make boundary layer
      amrex::Real d_wall = amrex::min(y, 2.54-y); // dist from upper and lower walls
      d_wall = amrex::min(d_wall, 1.905-amrex::Math::abs(z)); // dist from upper and lower walls
      // d_wall = amrex::max(d_wall, 0.0); // prevent solver from complaining
      d_wall /= (0.1245 * 2.0); // d_wall/D
      amrex::Real bl_profile = d_wall > 0.68 ? 1.0 : 1.0085656 - exp(-7.0 * d_wall);
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        state(i, j, k, nf * NVAR + URHO) = pp.rho / (1.0 + 0.5 * (pp.u * pp.u) / pp.cv_Tinf * (1.0 - bl_profile * bl_profile));;
        state(i, j, k, nf * NVAR + UMX) = state(i, j, k, nf * NVAR + URHO) * pp.u * bl_profile;
        state(i, j, k, nf * NVAR + UMY) = 0.0;
        state(i, j, k, nf * NVAR + UMZ) = 0.0;
        state(i, j, k, nf * NVAR + UEDEN) = pp.rho * pp.ei + 0.5 * pp.rho * pp.u * pp.u;
        for (int n = 0; n < NUM_SPECIES; ++n)
          state(i, j, k, nf * NVAR + UFS + n) = state(i, j, k, nf * NVAR + URHO) * pp.Y[n];
      }
    } else if ((rsqr <= 0.1245 * 0.1245) && (y < -0.5)) {
      // Make jet profile
      amrex::Real r_over_D = sqrt(rsqr) / (0.1245 * 2.0);
      amrex::Real jet_profile = r_over_D < 0.35 ? 1.0 : 1.0 - pow(6.67*(r_over_D-0.35), 2.5);
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        state(i, j, k, nf * NVAR + URHO) = pp.rho_j;
        state(i, j, k, nf * NVAR + UMX) = 0.0;
        state(i, j, k, nf * NVAR + UMY) = pp.rho_j * pp.v_j * jet_profile;
        state(i, j, k, nf * NVAR + UMZ) = 0.0;
        state(i, j, k, nf * NVAR + UEDEN) = pp.rho_j * pp.ei_j + 0.5 * pp.rho_j * pp.v_j * pp.v_j * jet_profile * jet_profile;
        for (int n = 0; n < NUM_SPECIES; ++n)
          state(i, j, k, nf * NVAR + UFS + n) = pp.rho_j * pp.Y_jet[n];
      }
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_timestep(
  int i, int j, int k, const amrex::Real /*curtime*/, const amrex::Real dtlev,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<const amrex::Real> const& reactions,
  amrex::GeometryData const& /*geomdata*/, Parm const& /*parm*/, ProbParm const& pp)
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
    amrex::Real T = 0.0, p = 0.0;
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, Y, T);
    eos.RTY2P(rho, T, Y, p);

    AMREX_D_TERM(state(i, j, k, UFA + 0) += vx * dtlev;    // u
                 , state(i, j, k, UFA + 1) += vy * dtlev;  // v
                 , state(i, j, k, UFA + 2) += vz * dtlev;) // w
    state(i, j, k, UFA + 3) += T * dtlev;                  // T
    state(i, j, k, UFA + 4) += p * dtlev;                  // p
    state(i, j, k, UFA + 5) += Y[H2_ID] * dtlev;           // Y_H2
    state(i, j, k, UFA + 6) += Y[OH_ID] * dtlev;           // Y_OH
    state(i, j, k, UFA + 7) += Y[O2_ID] * dtlev;          // Y_O2
    state(i, j, k, UFA + 8) += Y[HO2_ID] * dtlev;          // Y_HO2
    state(i, j, k, UFA + 9) +=
      reactions(i, j, k, NUM_SPECIES) * dtlev; // heat release rate
    state(i, j, k, UFA + 10) += p * p * dtlev; // p^2
    state(i, j, k, UFA + 11) += 2 * ke * dtlev; // 2*k
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
  aux_name[4] = "time_avg_p";
  aux_name[5] = "time_avg_Y_H2";
  aux_name[6] = "time_avg_Y_OH";
  aux_name[7] = "time_avg_Y_O2";
  aux_name[8] = "time_avg_Y_HO2";
  aux_name[9] = "time_avg_hrr";
  aux_name[10] = "mean_squared_p";
  aux_name[11] = "mean_squared_vel";
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_tag_error(
  int i, int j, int k, amrex::Array4<char> const& tagarr,
  amrex::Array4<const amrex::Real> const& sarr, int /*level*/, char tagval,
  const amrex::Real time, amrex::GeometryData const& geomdata, Parm const& /*parm*/,
  ProbParm const& /*pp*/)
{
}

#endif