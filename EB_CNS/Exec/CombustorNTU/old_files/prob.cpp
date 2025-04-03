#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  ParmParse pp("prob");

  pp.query("inflow_T", CNS::h_prob_parm->inflow_T);
  pp.query("inflow_p", CNS::h_prob_parm->inflow_p);
  pp.query("inflow_mach", CNS::h_prob_parm->inflow_mach);
  pp.query("interior_T", CNS::h_prob_parm->interior_T);
  pp.query("interior_P", CNS::h_prob_parm->interior_p);

  CNS::h_prob_parm->massfrac[N2_ID] = 0.7; 
  CNS::h_prob_parm->massfrac[O2_ID] = 0.3; 

#ifdef AMREX_USE_GPU
  // Cannot use Gpu::copy because ProbParm is not trivailly copyable.
  Gpu::htod_memcpy_async(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#else
  std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#endif

  Gpu::HostVector<Real> inflow_state(NVAR);

  auto eos = pele::physics::PhysicsType::eos();
  Real rho, eint, cs;
  eos.PYT2RE(CNS::h_prob_parm->inflow_p, CNS::h_prob_parm->massfrac.begin(),
             CNS::h_prob_parm->inflow_T, rho, eint);
  eos.RTY2Cs(rho, CNS::h_prob_parm->inflow_T, CNS::h_prob_parm->massfrac.begin(),
             cs);
  Real v = CNS::h_prob_parm->inflow_mach * cs;

  std::cout << " Velocity Inflow=" << v << std::endl;

  inflow_state[URHO] = rho;
  inflow_state[UMX] = 0.0;
  inflow_state[UMY] = 0.0;
  inflow_state[UMZ] = rho * v;
  inflow_state[UEDEN] = rho * eint + 0.5 * rho * v * v; 
  inflow_state[UFS + N2_ID] = rho * CNS::h_prob_parm->massfrac[N2_ID];
  inflow_state[UFS + O2_ID] = rho * CNS::h_prob_parm->massfrac[O2_ID];

  Gpu::copyAsync(Gpu::hostToDevice, inflow_state.data(), inflow_state.data() + NVAR,
                 CNS::h_prob_parm->inflow_state);
  Gpu::streamSynchronize();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& /*state*/,
                       amrex::Array4<amrex::Real> const& ext_src, ProbParm const& pp)
{
}

#if CNS_USE_EB
void Combustor_NTNU::build(const Geometry& geom, const int max_coarsening_level)
{
  // cylindrical domain
  auto cchamber_1= EB2::CylinderIF(2.2,50, 2,{AMREX_D_DECL(2.2, 2.2, 5)}, true);

  // upstream section 
  auto solid_us = EB2::CylinderIF(2.2, 9, 2,{AMREX_D_DECL(2.2, 2.2, 0)}, false);
  auto plenum_end = EB2::CylinderIF(0.95, 9, 2, {AMREX_D_DECL(2.2, 2.2, 0)}, false);

  auto injector_tube = EB2::DifferenceIF<EB2::CylinderIF, EB2::CylinderIF>(solid_us, plenum_end);

  // central cylinder holding bluff body
  auto bb_holder = EB2::CylinderIF(0.25, 9, 2, {AMREX_D_DECL(2.2, 2.2, 0)}, false);

  // bluff body   
  auto bb_constructorplane = EB2::PlaneIF({AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(-1, 1, 0)});
  auto lathe_bb_constructorplane = EB2::lathe(bb_constructorplane);
  RealArray offset = Array<Real, AMREX_SPACEDIM>{2.2, 2.2, 3.85};
  auto trans_lathe_bb_cp = EB2::translate( lathe_bb_constructorplane, offset);
  auto bb_zlimit = EB2::PlaneIF({AMREX_D_DECL(2.2, 2.2, 4.5)}, {AMREX_D_DECL(0,0, -1)});
  auto cone = EB2::IntersectionIF<EB2::PlaneIF, EB2::TranslationIF<EB2::LatheIF<EB2::PlaneIF>>>
                    (bb_zlimit, trans_lathe_bb_cp);

  auto bb_and_holder = EB2::UnionIF<EB2::IntersectionIF<EB2::PlaneIF, EB2::TranslationIF<EB2::LatheIF<EB2::PlaneIF>>>, 
                       EB2::CylinderIF>(cone, bb_holder);
                       
  auto final_struct = EB2::UnionIF<EB2::CylinderIF, EB2::DifferenceIF<EB2::CylinderIF, EB2::CylinderIF>, 
                      EB2::UnionIF<EB2::IntersectionIF<EB2::PlaneIF, EB2::TranslationIF<EB2::LatheIF<EB2::PlaneIF>>>, 
                      EB2::CylinderIF>>
                    (cchamber_1, injector_tube, bb_and_holder);

  auto gshop = EB2::makeShop(final_struct);
  
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
#endif



