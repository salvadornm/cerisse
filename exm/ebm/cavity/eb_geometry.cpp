#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include "custom_geometry.h"

// Geometry function
void Custom::build(const Geometry& geom, const int max_coarsening_level)
{

  // cavity params
  const Real D = 1.128e-04;
  const Real L = 2.256e-04;

  // filled box assuming domain origin is (0,0,0)
  const Real box_lox = -1;
  const Real box_loy = -1;
  const Real box_loz = -1;
  const Real box_hix = 1;
  const Real box_hiy = 0;
  const Real box_hiz = 1;

  // cavity
  const Real cav_lox = 3.9 * D;
  const Real adj_cav_lox = cav_lox + 0.00000035;// + 0.0000001;
  const Real cav_loy = -D;
  const Real adj_cav_loy = cav_loy - 0.0000002;//+ 0.000000033;
  const Real cav_loz = -0.5;
  const Real cav_hix = cav_lox + L;
  const Real adj_cav_hix = cav_hix + 0.00000020;
  const Real cav_hiy = 0;
  const Real cav_hiz = 0.5;
  
  auto box1 = EB2::BoxIF({AMREX_D_DECL(box_lox,box_loy,box_loz)},{AMREX_D_DECL(box_hix,box_hiy,box_hiz)},false);
  auto box2 = EB2::BoxIF({AMREX_D_DECL(adj_cav_lox,adj_cav_loy,cav_loz)},{AMREX_D_DECL(adj_cav_hix,cav_hiy,cav_hiz)},false);

  auto box_nocav = EB2::DifferenceIF<EB2::BoxIF, EB2::BoxIF>(box1,box2);

  auto gshop = EB2::makeShop(box_nocav);

  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
