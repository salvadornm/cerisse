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
 
  CNS::h_prob_parm->massfrac[H_ID] = 0; 
  CNS::h_prob_parm->massfrac[H2_ID] = 0.013; 
  CNS::h_prob_parm->massfrac[O_ID] = 0; 
  CNS::h_prob_parm->massfrac[OH_ID] = 0; 
  CNS::h_prob_parm->massfrac[H2O_ID] = 0; 
  CNS::h_prob_parm->massfrac[O2_ID] = 0.227; 
  CNS::h_prob_parm->massfrac[HO2_ID] = 0; 
  CNS::h_prob_parm->massfrac[H2O2_ID] = 0;  
  CNS::h_prob_parm->massfrac[N2_ID] = 0.76; 
  CNS::h_prob_parm->massfrac[AR_ID] = 0; 
  CNS::h_prob_parm->massfrac[HE_ID] = 0; 
  CNS::h_prob_parm->massfrac[CO_ID] = 0; 
  CNS::h_prob_parm->massfrac[CO2_ID] = 0; 

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

  std::cout << "density " << rho << std::endl;
  std::cout << "speed of sound " << cs << std::endl;
  std::cout << " Velocity Inflow= " << v << std::endl;

  inflow_state[URHO] = rho;
  inflow_state[UMX] = 0.0;
  inflow_state[UMY] = 0.0;
  inflow_state[UMZ] = rho * v;
  inflow_state[UEDEN] = rho * eint + 0.5 * rho * v * v; 

  for (int n = 0; n < NUM_SPECIES; n++) {
    inflow_state[UFS + n] = rho * CNS::h_prob_parm->massfrac[n];
  }

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
void NTNU_combustor::build(const Geometry& geom, const int max_coarsening_level)
{
  const Real mm2cm=0.1;  // change mm to m; its acc mm to cm

  // centre of domain  
  const Real x0 =  22*mm2cm; 
  const Real y0 =  22*mm2cm;

  // big chamber
  const Real R_chamber = 22*mm2cm;
  const Real L_chamber = 90*mm2cm;
  const int  dir_chamber=2; // points to  Z
  const Real zchamber = 90*mm2cm; //
  
  // injector
  const Real R_inj = 9.5*mm2cm;
  const Real L_inj = 90*mm2cm;
  const int  dir_inj= 2; // points to Z
  const Real zinj = 30*mm2cm;

  // pipe
  const Real Rpipe = 2.5*mm2cm;
  const Real Lpipe = 90*mm2cm;
  const int  dir_pipe= 2; // points to Z
  const Real zpipe = 0.0*mm2cm;

  // aux z
  const Real zoffset = 38.5*mm2cm;
  const Real zbb = 45*mm2cm;

  auto cchamber_1 = EB2::CylinderIF(R_chamber,L_chamber, dir_chamber,{AMREX_D_DECL(x0, y0, zchamber)}, true);
  auto injector   = EB2::CylinderIF(R_inj,L_inj, dir_inj, {AMREX_D_DECL(x0, y0, zinj)}, false);

  auto cchamber_with_injector = EB2::DifferenceIF<EB2::CylinderIF, EB2::CylinderIF>(cchamber_1,injector);

   // central cylinder holding bluff body
   auto bb_holder = EB2::CylinderIF(Rpipe,Lpipe, dir_pipe, {AMREX_D_DECL(x0, y0, zpipe)}, false);

   // bluff body
   auto bb_constructorplane = EB2::PlaneIF({AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(-1, 1, 0)});
   auto lathe_bb_constructorplane = EB2::lathe(bb_constructorplane);
   RealArray offset = Array<Real, AMREX_SPACEDIM>{x0, y0, zoffset};
   auto trans_lathe_bb_cp = EB2::translate( lathe_bb_constructorplane, offset);
   auto bb_zlimit = EB2::PlaneIF({AMREX_D_DECL(x0, y0, zbb)}, {AMREX_D_DECL(0,0, -1)});
   auto cone = EB2::IntersectionIF<EB2::PlaneIF, EB2::TranslationIF<EB2::LatheIF<EB2::PlaneIF>>>
                     (bb_zlimit, trans_lathe_bb_cp);

  auto bb_and_holder = EB2::UnionIF<EB2::IntersectionIF<EB2::PlaneIF, EB2::TranslationIF<EB2::LatheIF<EB2::PlaneIF>>>,
                        EB2::CylinderIF>(cone, bb_holder);

  auto final_struct = EB2::UnionIF<
                      EB2::DifferenceIF<EB2::CylinderIF, EB2::CylinderIF>,
                      EB2::UnionIF<EB2::IntersectionIF<EB2::PlaneIF, EB2::TranslationIF<EB2::LatheIF<EB2::PlaneIF>>>,
                      EB2::CylinderIF>
                      >
                      (cchamber_with_injector, bb_and_holder);

  auto gshop = EB2::makeShop(final_struct);

  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
#endif
