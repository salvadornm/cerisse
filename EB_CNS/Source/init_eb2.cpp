#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>

#include "custom_geometry.H"

using namespace amrex;

void initialize_EB2(const Geometry& geom, const int required_coarsening_level,
                    const int max_coarsening_level)
{
  BL_PROFILE("initializeEB2");

  ParmParse ppeb2("eb2");
  std::string geom_type = "all_regular";
  ppeb2.query("geom_type", geom_type);

  Vector<std::string> amrex_defaults(
    {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser", "stl"});
  if (std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) ==
      amrex_defaults.end()) {
    // Non-AMReX default EB types, get from CustomGeometry
    auto geometry = CustomGeometry::create(geom_type);
    geometry->build(geom, max_coarsening_level);
  } else {
    // Build with AMReX default EB2
    // EB2::Build(geom, required_coarsening_level, max_coarsening_level);
    // (..., ngrow, build_coarse_level_by_coarsening) what does that mean?
    EB2::Build(geom, required_coarsening_level, max_coarsening_level, 6, true);
    // If you have a hyperbolic system and the EB surface does not intersect AMR
    // level boundaries, you can set it to false to improve the robustness of
    // geometry generation.
  }
}
