#ifndef ebm_H_
#define ebm_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>


template <typename cls_t>
class ebm_t
{ 

public:	

  ebm_t() { 
  }

  ~ebm_t() {}


void init(const Geometry& geom, const int required_coarsening_level,
                    const int max_coarsening_level)
{
  BL_PROFILE("initializeEB2");

  Vector<std::string> amrex_defaults(
    {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser", "stl"});


  ParmParse ppeb2("eb2");
  std::string geom_type = "all_regular";
  ppeb2.query("geom_type", geom_type);

  if (std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) ==
      amrex_defaults.end()) {

    std::cout << " custom EB types" << std::endl;    
    //  auto geometry = CustomGeometry::create(geom_type); // SNM
    //   geometry->build(geom, max_coarsening_level);  //SNM
  } 
  else {
    std::cout << " AMReX default EB types" << std::endl;
    EB2::Build(geom, required_coarsening_level, max_coarsening_level, 6, true);
  }
}


};


#endif
