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

  inflow_state[URHO] = rho;
  inflow_state[UMX] = 0.0;
  inflow_state[UMY] = 0.0;
  inflow_state[UMZ] = rho * v;
  inflow_state[UEDEN] = rho * eint + 0.5 * rho * v * v;
  inflow_state[UFS + N2_ID] = 0.0;
  inflow_state[UFS + O2_ID] = rho;

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
void Combustor_modif::build(const Geometry& geom, const int max_coarsening_level)
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

#if CNS_USE_EB
void Combustor_NTNU::build(const Geometry& geom, const int max_coarsening_level)
{
  // cylindrical domain
  auto cchamber_1= EB2::CylinderIF(2.2, 2,{AMREX_D_DECL(2.2, 2.2, 5)}, true);

  // upstream section 
  auto solid_us = EB2::CylinderIF(2.2, 4.5, 2,{AMREX_D_DECL(2.2, 2.2, 2.25)}, false);
  auto plenum_end = EB2::CylinderIF(0.95, 4.5, 2, {AMREX_D_DECL(2.2, 2.2, 2.25)}, false);

  auto injector_tube = EB2::DifferenceIF<EB2::CylinderIF, EB2::CylinderIF>(solid_us, plenum_end);

  // central cylinder holding bluff body
  auto bb_holder = EB2::CylinderIF(0.25, 4.5, 2, {AMREX_D_DECL(2.2, 2.2, 2.25)}, false);


  auto final_struct = EB2::UnionIF<EB2::CylinderIF, EB2::DifferenceIF<EB2::CylinderIF, EB2::CylinderIF>, EB2::CylinderIF>
                    (cchamber_1, injector_tube, bb_holder);

  auto gshop = EB2::makeShop(final_struct);
  
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
#endif

#if CNS_USE_EB
void test_bb::build(const Geometry& geom, const int max_coarsening_level)
{
  // auto bb_holder = EB2::CylinderIF(0.25, 4.5, 2, {AMREX_D_DECL(2.2, 2.2, 2.25)}, false);

  auto bottom_plane = EB2::PlaneIF({AMREX_D_DECL(0, 2.2, 4.5)}, {AMREX_D_DECL(0,1,0)});
  auto far_plane = EB2::PlaneIF({AMREX_D_DECL(0, 2.85, 4.5)}, {AMREX_D_DECL(0,0,-1)});
  auto inclined_plane = EB2::PlaneIF({AMREX_D_DECL(0, 2.85, 4.5)}, {AMREX_D_DECL(0, -cos(M_PI/4), cos(M_PI/4))});
  auto cut_plane_xplus = EB2::PlaneIF({AMREX_D_DECL(2.1999999, 2.85, 4.5)}, {AMREX_D_DECL(1,0,0)});
   auto cut_plane_xminus = EB2::PlaneIF({AMREX_D_DECL(2.0000001, 2.85, 4.5)}, {AMREX_D_DECL(-1,0,0)});

  auto _3d_tri = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(bottom_plane,far_plane,inclined_plane);
  auto _2d_tri_half = EB2::IntersectionIF<EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>, EB2::PlaneIF>(_3d_tri, cut_plane_xplus);
  auto _2d_tri = EB2::IntersectionIF<EB2::IntersectionIF<EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>, EB2::PlaneIF>,
                                    EB2::PlaneIF>(_2d_tri_half, cut_plane_xminus);

  auto bb_full  = EB2::lathe(cut_plane_xminus);
  // auto bb_plus_holder = EB2::UnionIF<EB2::CylinderIF, EB2::LatheIF<EB2::IntersectionIF<EB2::IntersectionIF<EB2::PlaneIF,
  //                                   EB2::PlaneIF, EB2::PlaneIF>, EB2::PlaneIF>>>(bb_holder, bb_full);




  //auto final_struct = EB2::UnionIF<EB2::UnionIF<EB2::CylinderIF, EB2::LatheIF<EB2::IntersectionIF<EB2::IntersectionIF<EB2::PlaneIF,
  //                                   EB2::PlaneIF, EB2::PlaneIF>, EB2::PlaneIF>>>, EB2::DifferenceIF<EB2::BoxIF,EB2::CylinderIF>>(bb_plus_holder,
  //                                   half_one);

  // auto final_struct = EB2::UnionIF<EB2::DifferenceIF<EB2::BoxIF,EB2::CylinderIF>, EB2::CylinderIF
  //                                  >
  //                                (half_one,bb_holder);

  auto final_struct = bb_full;

  auto gshop = EB2::makeShop(final_struct);
  
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
#endif

