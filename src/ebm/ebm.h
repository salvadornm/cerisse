#ifndef ebm_H_
#define ebm_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>

//#include <EBMultiFab.h>

template <typename cls_t>
class ebm_t
{ 

public:	

  ebm_t() { 
  }

  ~ebm_t() {}

  const amrex::MultiFab* volfrac; // pointer  to a const MultiFab

  Vector<MultiFab*> ebmf_a;       // array of pointers to multifab 
   

///////////////////////////////////////////////////////////////////////////
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

    amrex::Print() << " custom EB types" << std::endl;    
    //  auto geometry = CustomGeometry::create(geom_type); // SNM
    amrex::Abort( " EB Type .. not ready yet (use cerisse0)"); 
    //   geometry->build(geom, max_coarsening_level);  //SNM
  } 
  else {
    amrex::Print() << " AMReX default EB types" << std::endl;
    EB2::Build(geom, required_coarsening_level, max_coarsening_level, 6, true);
  }
}
///////////////////////////////////////////////////////////////////////////
  void inline copytoRealMF(MultiFab& mf, int ibcomp, int mfcomp) {

    //amrex::Print() << " Copy vfrac data" << std::endl;


    //printf(" mfcomp=%d ibcomp = %d \n",mfcomp,ibcomp);


    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
      const Box& ebbox = mfi.validbox();  // box without ghost points

      //const Array4<Real>& vf        = (*volfrac).array(mfi);  // volfrac    
      Array4<const Real> vf = (*volfrac).const_array(mfi);
      const Array4<Real>& realfield = mf.array(mfi);       //  array data from mf

      amrex::ParallelFor(
          ebbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            realfield(i, j, k, mfcomp) = vf(i, j, k, ibcomp);
            //realfield(i, j, k, mfcomp) = 33.0_rt;   
          });
    }
  }
////////////////////////////////////////////////////////////////////////////


};


#endif
