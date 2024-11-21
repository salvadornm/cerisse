#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
  auto eos = pele::physics::PhysicsType::eos();

  // Calculate inflow conditions
  amrex::Real p0 = 590e4; // total pressure [Ba]
  amrex::Real M = 2.2;    // inflow Mach number
  amrex::Real T0;
  bool micka_fuel_cond = false; // use p0, T0 or p, T
  {
    ParmParse pp("prob");
    pp.get("T0", T0); // total temperature [K]
    pp.get("spark", CNS::h_prob_parm->spark); // spark plug (modelled as heating wall)
    pp.query("M", M);
    pp.query("record_statistics", CNS::h_prob_parm->record_statistics);
    pp.query("clean_aux_on_restart", CNS::h_prob_parm->clean_aux_on_restart);
    pp.query("do_bl", CNS::h_prob_parm->do_bl);
    pp.query("make_bl_on_restart", CNS::h_prob_parm->make_bl_on_restart);
    pp.query("micka_fuel_cond", micka_fuel_cond);
    pp.query("make_init_on_restart", CNS::h_prob_parm->make_init_on_restart);
    if (CNS::h_prob_parm->make_init_on_restart)
      pp.get("x_reset", CNS::h_prob_parm->x_reset);
  }
  if constexpr (NUM_AUX <= 0) {
    if (CNS::h_prob_parm->record_statistics)
      amrex::Abort("Please compile with NUM_AUX to record statistics");
  }
  if (CNS::h_prob_parm->make_bl_on_restart && CNS::h_prob_parm->do_bl) {
    amrex::Print() << "Restarting with BL...\n";
  }
  if (CNS::h_prob_parm->make_init_on_restart) {
    amrex::Print() << "Restarting and reinitialising data at x < "
                   << CNS::h_prob_parm->x_reset << "...\n";
  }

  CNS::h_prob_parm->Y[H2O_ID] = 1.068e-7 * T0 * T0 - 6.72e-5 * T0 + 2.986e-2; // closer to experiment
  CNS::h_prob_parm->Y[O2_ID] = 4.0 / 15.0 * (1.0 - CNS::h_prob_parm->Y[H2O_ID]);
  CNS::h_prob_parm->Y[N2_ID] = 11.0 / 15.0 * (1.0 - CNS::h_prob_parm->Y[H2O_ID]);

  amrex::Real gamma = 1.313, rho, T, p;

  if (M > 1.0) {
    // Iterate to find gamma
    for (int iter = 0; iter < 10; ++iter) {
      // Isentropic relations
      T = T0 / (1 + 0.5 * (gamma - 1) * M * M);
      p = p0 * std::pow(1 + 0.5 * (gamma - 1) * M * M, -gamma / (gamma - 1));
      eos.PYT2R(p, CNS::h_prob_parm->Y.begin(), T, rho);
      eos.RTY2G(rho, T, CNS::h_prob_parm->Y.begin(), gamma);
      // amrex::Print() << "Iter " << iter << " (gamma, T, p) = " << gamma << ", " <<
      // T << ", " << p << '\n';
    }
  } else {
    // Constant mass flow rate (constant area isolator)
    // First, calculate conditions at nozzle outlet (1)
    Real M1 = 2.2, p1, T1, gamma1 = 1.313;
    // Iterate to find gamma
    for (int iter = 0; iter < 10; ++iter) {
      // Isentropic relations
      T1 = T0 / (1 + 0.5 * (gamma1 - 1) * M1 * M1);
      p1 = p0 * std::pow(1 + 0.5 * (gamma1 - 1) * M1 * M1, -gamma1 / (gamma1 - 1));
      eos.PYT2R(p1, CNS::h_prob_parm->Y.begin(), T1, rho);
      eos.RTY2G(rho, T1, CNS::h_prob_parm->Y.begin(), gamma1);
      // amrex::Print() << "Iter " << iter << " (gamma1, T1, p1) = " << gamma1 << ",
      // " << T1 << ", " << p1 << '\n';
    }

    // Then, assume constant mass flow rate, knowing p and M
    {
      ParmParse pp("prob");
      pp.get("p", p);
    }
    // Iterate to find gamma
    for (int iter = 0; iter < 10; ++iter) {
      T = pow((p / p1) * (M / M1), 2) * (gamma / gamma1) * T1;
      eos.PYT2R(p, CNS::h_prob_parm->Y.begin(), T, rho);
      eos.RTY2G(rho, T, CNS::h_prob_parm->Y.begin(), gamma);
      // amrex::Print() << "Iter " << iter << " (gamma, T, rho) = " << gamma << ", "
      // << T << ", " << rho << '\n';
    }

    eos.PYT2R(p, CNS::h_prob_parm->Y.begin(), T, rho);
    eos.RTY2G(rho, T, CNS::h_prob_parm->Y.begin(), gamma);
  }
  CNS::h_prob_parm->T = T;
  amrex::Print() << "Inflow (gamma, T, p) = " << gamma << ", " << T << ", " << p
                 << '\n';
  // Real p = 55.410e4;
  amrex::Real cv;
  eos.TY2Cv(T, CNS::h_prob_parm->Y.begin(), cv);
  CNS::h_prob_parm->cv_Tinf = cv * T;

  eos.PYT2RE(p, CNS::h_prob_parm->Y.begin(), T, CNS::h_prob_parm->rho,
             CNS::h_prob_parm->ei);
  Real cs;
  eos.RTY2Cs(CNS::h_prob_parm->rho, T, CNS::h_prob_parm->Y.begin(), cs);
  // Real cs = std::sqrt(gamma * p / CNS::h_prob_parm->rho);
  CNS::h_prob_parm->u = M * cs;

  // Fuel conditions
  CNS::h_prob_parm->Y_jet[H2_ID] = 1.0;

  if (micka_fuel_cond) {
    p0 = 845.0e4 + (755.e4 - 845.e4) * (T0 - 1100.0) / 300.0; // linearly varying
    M = 1.0;
    T0 = 288.0;
    gamma = 1.405;
    for (int iter = 0; iter < 10; ++iter) {
      T = T0 / (1 + 0.5 * (gamma - 1));
      p = p0 * std::pow(1 + 0.5 * (gamma - 1), -gamma / (gamma - 1));
      eos.PYT2R(p, CNS::h_prob_parm->Y_jet.begin(), T, rho);
      eos.RTY2G(rho, T, CNS::h_prob_parm->Y_jet.begin(), gamma);
    }
  } else {
    // Rebeiro
    p = 845.0e4 + (755.e4 - 845.e4) * (T0 - 1100.0) / 300.0; // linearly varying
    M = 1.0;
    T = 288.0;
    eos.PYT2R(p, CNS::h_prob_parm->Y_jet.begin(), T, rho);
    eos.RTY2G(rho, T, CNS::h_prob_parm->Y_jet.begin(), gamma);
  }
  CNS::h_prob_parm->T_j = T;
  amrex::Print() << "Fuel (gamma, T, p) = " << gamma << ", " << T << ", " << p
                 << '\n';

  eos.PYT2RE(p, CNS::h_prob_parm->Y_jet.begin(), T, CNS::h_prob_parm->rho_j,
             CNS::h_prob_parm->ei_j);
  eos.RTY2Cs(CNS::h_prob_parm->rho_j, T, CNS::h_prob_parm->Y_jet.begin(), cs);
  // cs = std::sqrt(gamma * p / CNS::h_prob_parm->rho_j);
  CNS::h_prob_parm->v_j = M * cs;

  Real mdot_air = CNS::h_prob_parm->rho * CNS::h_prob_parm->u * 3.81 * 2.54;
  Real mdot_jet =
    CNS::h_prob_parm->rho_j * CNS::h_prob_parm->v_j * M_PI * 0.1245 * 0.1245;

  amrex::Print() << "Global eq ratio = " << mdot_jet / mdot_air * 34.0 << '\n';
  // Real phi = 0.27;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}

void CNS::fill_ext_src(int i, int j, int k, Real time, GeometryData const& geomdata,
                       Array4<const Real> const& /*state*/,
                       Array4<Real> const& ext_src, ProbParm const& pp)
{
}

#if CNS_USE_EB
void Scramjet::build(const Geometry& geom, const int max_coarsening_level)
{
  auto box = EB2::BoxIF({AMREX_D_DECL(-50., -10., -10.)}, {AMREX_D_DECL(4.45, 0.0, 10.)}, false);
  auto injector = EB2::CylinderIF(0.1190625, 2.0, 1, {AMREX_D_DECL(0.0, -0.635, 0.0)}, false);
  auto box_with_inj = EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>(box, injector);

  auto rear_wall = EB2::PlaneIF({AMREX_D_DECL(9.525, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 0.0, 0.0)});
  auto floor_wall = EB2::PlaneIF({AMREX_D_DECL(0.0, -10.0, 0.0)}, {AMREX_D_DECL(0.0, 1.0, 0.0)});
  auto inclined_wall = EB2::PlaneIF(
    {AMREX_D_DECL(9.525, 0.0, 0.0)},
    {AMREX_D_DECL(-sin(4.0 / 180.0 * M_PI), -cos(4.0 / 180.0 * M_PI), 0.0)});
  auto triangle = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(
    rear_wall, floor_wall, inclined_wall);

  auto all_objs =
    EB2::UnionIF<EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>,
                 EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>(
      box_with_inj, triangle);
  auto gshop = EB2::makeShop(all_objs);
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
#endif
