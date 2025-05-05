#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include "custom_geometry.h"

// Geometry function
void Custom::build(const Geometry& geom, const int max_coarsening_level)
{
  
  const Real cm2m=0.01;  // change cm to m

  // centre of domain  
  const Real x0 =  0.0; // was 2.2 
  const Real y0 =  0.0;

  // big chamber
  const Real R_chamber = 2.2*cm2m;
  const Real L_chamber = 9.0*cm2m;
  const int  dir_chamber=2; // points to  Z
  const Real zchamber = 9*cm2m; //
  
  // injector
  const Real R_inj = 0.95*cm2m;
  const Real L_inj = 9.0*cm2m;
  const int  dir_inj= 2; // points to Z
  const Real zinj = 3.0*cm2m;

  // pipe
  const Real Rpipe = 0.25*cm2m;
  const Real Lpipe = 9.0*cm2m;
  const int  dir_pipe= 2; // points to Z
  const Real zpipe = 0.0*cm2m;

  // aux z
  const Real zoffset = 3.85*cm2m;
  const Real zbb = 4.5*cm2m;

  // cylindrical domain

  // NOTE: CylinderIF(radius, height, direction, center, has_fluid_inside);

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
