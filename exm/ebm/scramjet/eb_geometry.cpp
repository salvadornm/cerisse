#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include "custom_geometry.h"

// Geometry function
void Custom::build(const Geometry& geom, const int max_coarsening_level)
{
  auto box = EB2::BoxIF({AMREX_D_DECL(-5.0, -1.0, -1.0)}, {AMREX_D_DECL(0.0445, 0.0, 1.0)}, false);
  auto injector = EB2::CylinderIF(0.001245, 0.02, 1, {AMREX_D_DECL(0.0, -0.00635, 0.0)}, false);
  auto box_with_inj = EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>(box, injector);

  auto rear_wall = EB2::PlaneIF({AMREX_D_DECL(0.09525, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 0.0, 0.0)});
  auto floor_wall = EB2::PlaneIF({AMREX_D_DECL(0.0, -1.0, 0.0)}, {AMREX_D_DECL(0.0, 1.0, 0.0)});
  auto inclined_wall = EB2::PlaneIF(
    {AMREX_D_DECL(0.09525, 0.0, 0.0)},
    {AMREX_D_DECL(-std::sin(4.0 / 180.0 * M_PI), -std::cos(4.0 / 180.0 * M_PI), 0.0)});
  auto triangle = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(
    rear_wall, floor_wall, inclined_wall);

  auto all_objs =
    EB2::UnionIF<EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>,
                 EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>(
      box_with_inj, triangle);
  auto gshop = EB2::makeShop(all_objs);
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
