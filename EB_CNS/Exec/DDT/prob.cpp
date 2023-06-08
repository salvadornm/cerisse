#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit (const int* /*init*/,
                     const int* /*name*/,
                     const int* /*namelen*/,
                     const Real* /*problo*/,
                     const Real* /*probhi*/)
{
}
}

void
CNS::fill_ext_src (int i, int j, int k, 
                   Real time,
                   GeometryData const& geomdata, 
                   Array4<const Real> const& /*state*/, 
                   Array4<Real> const& ext_src, 
                   Parm const& /*parm*/,
                   ProbParm const& pp)
{
}

void
GraVent::build (const Geometry& geom, const int max_coarsening_level)
{
  // Geometry configuration
  Real BR, S;
  {
    ParmParse pp("prob");
    pp.get("BR", BR); //blockage ratio
    pp.get("S", S); //obst spacing [cm], 10.0 or 30.0
  }    
  Real first_obst = 25.0, last_obst = 205.0;
  Real h = 0.5 * BR * 6.0; //obst height [cm]

  // number of obstacles
  int num_obst; // = (last_obst - first_obst) / S + 1;
  if (S == 30.0) {
    num_obst = 7;
  } else if (S == 10.0) {
    num_obst = 19;
  } else {
    amrex::Abort("S must be 10.0 or 30.0");
  }

  // Storing obstacles
  // std::array<std::unique_ptr<EB2::BoxIF>, 2*num_obst> obsts; // 2*num_obst for both upper and lower boxes
  amrex::Vector<std::unique_ptr<EB2::BoxIF>> obsts(2*num_obst);
  // we use std::unique_ptr here because EB2::BoxIF does not have a default constructor,
  // so we cannot make an array of EB2::BoxIF directly
  
  // Make obstacles
  for (int i = 0; i < num_obst; ++i) {
    Real x = first_obst + i * S;

    obsts[2*i] = std::make_unique<EB2::BoxIF>(RealArray({AMREX_D_DECL(x-0.6, -3.1, -15.1)}), 
                                              RealArray({AMREX_D_DECL(x+0.6, -3.0+h, 15.1)}), 0); //lo end, hi end with 0.1 extra
    obsts[2*i+1] = std::make_unique<EB2::BoxIF>(RealArray({AMREX_D_DECL(x-0.6, 3.0-h, -15.1)}), 
                                                RealArray({AMREX_D_DECL(x+0.6, 3.1, 15.1)}), 0);
  }

  // Combine and make geom
  if (S == 30.0) {
    auto all_obsts = EB2::makeUnion(*obsts[0],*obsts[1],*obsts[2],*obsts[3],*obsts[4],*obsts[5],*obsts[6],
                                    *obsts[7],*obsts[8],*obsts[9],*obsts[10],*obsts[11],*obsts[12],*obsts[13]);                                    
    auto gshop = EB2::makeShop(all_obsts);
    EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
  } else {
    auto all_obsts = EB2::makeUnion(*obsts[0],*obsts[1],*obsts[2],*obsts[3],*obsts[4],*obsts[5],*obsts[6],
                                    *obsts[7],*obsts[8],*obsts[9],*obsts[10],*obsts[11],*obsts[12],*obsts[13],
                                    *obsts[14],*obsts[15],*obsts[16],*obsts[17],*obsts[18],*obsts[19],*obsts[20],
                                    *obsts[21],*obsts[22],*obsts[23],*obsts[24],*obsts[25],*obsts[26],*obsts[27],
                                    *obsts[28],*obsts[29],*obsts[30],*obsts[31],*obsts[32],*obsts[33],*obsts[34],
                                    *obsts[35],*obsts[36],*obsts[37]);                                    
    auto gshop = EB2::makeShop(all_obsts);
    EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
  }

  // auto func_wrapper = [] (auto&&... args) -> auto { EB2::makeUnion(*args.get()...); }; // converts std::unique_ptr<EB2::BoxIF> to EB2::BoxIF
  // auto all_obsts = std::apply(func_wrapper, obsts);
  // auto gshop = EB2::makeShop(all_obsts);
  // EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}