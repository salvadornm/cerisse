#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
  auto eos = pele::physics::PhysicsType::eos();

  // Calculate inflow conditions
  Real p0 = 590e4; // total pressure [Ba]
  Real M = 2.2; // inflow Mach number
  Real T0;
  {
    ParmParse pp("prob");
    pp.get("T0", T0); // total temperature [K]
    pp.get("spark", CNS::h_prob_parm->spark); // spark plug (modelled as heating wall)
    pp.query("record_statistics", CNS::h_prob_parm->record_statistics);
    pp.query("clean_aux_on_restart", CNS::h_prob_parm->clean_aux_on_restart);
  }
  CNS::h_prob_parm->Y[O2_ID] = 0.244;
  CNS::h_prob_parm->Y[N2_ID] = 0.671;
  CNS::h_prob_parm->Y[H2O_ID] = 0.085;

  // Iterate to find gamma
  Real gamma = 1.313, rho, T, p;
  for (int iter = 0; iter < 3; ++ iter) {
    T = T0 / (1 + 0.5 * (gamma - 1) * M * M);
    p = p0 * std::pow(1 + 0.5 * (gamma - 1) * M * M, -gamma / (gamma - 1));

    eos.PYT2R(p, CNS::h_prob_parm->Y.begin(), T, rho);
    eos.RTY2G(rho, T, CNS::h_prob_parm->Y.begin(), gamma);
  }
  amrex::Print() << "Inflow (gamma, T, p) = " << gamma << ", " << T << ", " << p << '\n';
  // Real p = 55.410e4;

  eos.PYT2RE(p, CNS::h_prob_parm->Y.begin(), T, CNS::h_prob_parm->rho, CNS::h_prob_parm->ei);

  Real cs = std::sqrt(gamma * p / CNS::h_prob_parm->rho);
  CNS::h_prob_parm->u = M * cs;
  
  // Fuel conditions
  p0 = 845.0e4 + (755.e4 - 845.e4) * (T0 - 1100.0) / 300.0; // linearly varying
  M = 1.0;
  T0 = 288.0;

  CNS::h_prob_parm->Y_jet[H2_ID] = 1.0;

  gamma = 1.405, rho, T, p;
  for (int iter = 0; iter < 3; ++ iter) {
    T = T0 / (1 + 0.5 * (gamma - 1));
    p = p0 * std::pow(1 + 0.5 * (gamma - 1), -gamma / (gamma - 1));

    eos.PYT2R(p, CNS::h_prob_parm->Y_jet.begin(), T, rho);
    eos.RTY2G(rho, T, CNS::h_prob_parm->Y_jet.begin(), gamma);
  }
  amrex::Print() << "Fuel (gamma, T, p) = " << gamma << ", " << T << ", " << p << '\n';

  eos.PYT2RE(p, CNS::h_prob_parm->Y_jet.begin(), T, CNS::h_prob_parm->rho_j, CNS::h_prob_parm->ei_j);

  cs = std::sqrt(gamma * p / CNS::h_prob_parm->rho_j);
  CNS::h_prob_parm->v_j = M * cs;

  Real mdot_air = CNS::h_prob_parm->rho * CNS::h_prob_parm->u * 3.81 * 2.54;
  Real mdot_jet = CNS::h_prob_parm->rho_j * CNS::h_prob_parm->v_j * M_PI * 0.1245 * 0.1245;

  amrex::Print() << "Global eq ratio = " << mdot_jet / mdot_air * 34.0 << '\n';
  // Real phi = 0.27;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}

void CNS::fill_ext_src(int i, int j, int k, Real time, GeometryData const& geomdata,
                       Array4<const Real> const& /*state*/,
                       Array4<Real> const& ext_src, Parm const& /*parm*/,
                       ProbParm const& pp)
{
}

void Scramjet::build(const Geometry& geom, const int max_coarsening_level)
{
  auto box = EB2::BoxIF({AMREX_D_DECL(-50., -10., -10.)}, {AMREX_D_DECL(4.45, 0.0, 10.)}, false);
  auto injector = EB2::CylinderIF(0.1245, 2.0, 1, {AMREX_D_DECL(0.0, -0.635, 0.0)}, false);
  auto box_with_inj = EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>(box, injector);

  auto rear_wall = EB2::PlaneIF({AMREX_D_DECL(9.525, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 0.0, 0.0)});
  auto floor_wall = EB2::PlaneIF({AMREX_D_DECL(0.0, -10.0, 0.0)}, {AMREX_D_DECL(0.0, 1.0, 0.0)});
  auto inclined_wall = EB2::PlaneIF({AMREX_D_DECL(9.53, 0.0, 0.0)}, {AMREX_D_DECL(-sin(4.8 / 180.0 * M_PI), -cos(4.8 / 180.0 * M_PI), 0.0)});
  auto triangle = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(rear_wall, floor_wall, inclined_wall);
  
  // auto all_objs = EB2::UnionIF<EB2::BoxIF, EB2::CylinderIF, 
  //                 EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>(box, injector, triangle);
  auto all_objs = EB2::UnionIF<EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>, 
                  EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>(box_with_inj, triangle);
  auto gshop = EB2::makeShop(all_objs);
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
