#ifndef ebm_H_
#define ebm_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>

#include <EBMultiFab.h>


template <typename cls_t>
class ebm_t
{ 

public:	

  ebm_t() { 
  }

  ~ebm_t() {}

  const amrex::MultiFab* volfrac;        // pointer  to a const MultiFab
  const amrex::MultiCutFab* bndrycent;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

  // pointer to Amr class instance
  Amr* amr_p;

  // define a Multifab of 

  Vector<EBMultiFab<bool>*> bmf_a; //(bool multifab array)

  const amrex::Real vfracmin = 1.e-6;        // minimum to consider a cell "empty"
  const amrex::Real vfracmax = 1.0-vfracmin; // maximum to consider a cell "not-full"
  
  ///////////////////////////////////////////////////////////////////////////
  void init(Amr* pointer_amr,const Geometry& geom, const int required_level, 
                                                         const int max_level)
  {
    BL_PROFILE("initializeEB2");

    amrex::Print() << " Initialize EB at maxlevel = " << max_level << std::endl;


    Vector<std::string> amrex_defaults(
    {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser", "stl"});

    ParmParse ppeb2("eb2");
    std::string geom_type = "all_regular";
    ppeb2.query("geom_type", geom_type);

    if (std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) ==
      amrex_defaults.end()) {

      amrex::Print() << " custom EB types" << std::endl;    
      amrex::Abort( " EB Type .. not ready yet (use cerisse0)"); 
      //  auto geometry = CustomGeometry::create(geom_type); // SNM
      //   geometry->build(geom, max_coarsening_level);  //SNM
    } 
    else {
    amrex::Print() << " AMReX default EB types" << std::endl;
    EB2::Build(geom, required_level, max_level, 6, true);
    }

    // store pointer to AMR class (just in case)
    amr_p = pointer_amr;    
    
    // size of multifab
    bmf_a.resize(max_level + 1);

  }
  ////////////////////////////////////////////////////////////////////////////
  // create EBMultiFabs at a level and store pointers to it
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    bmf_a[lev] = new EBMultiFab<bool>(bxa, dm, 2, cls_t::NGHOST);    
  }
  ////////////////////////////////////////////////////////////////////////////
  void destroy_mf(int lev)
  {
    if (!bmf_a.empty()) { delete bmf_a.at(lev); }
  }
  ////////////////////////////////////////////////////////////////////////////
  void computeMarkers(int lev)
  {

    //amrex::Print( ) << " oo computeMarkers  lev=" << lev << std::endl;  

    auto& mfab = *bmf_a[lev];
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      auto& ibFab = mfab.get(mfi);
      const Box& bx = mfi.tilebox();
      const auto& ebMarkers = mfab.array(mfi); // boolean array

      Array4<const Real> vf = (*volfrac).const_array(mfi);

      // LOOP
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                
        // initialise to false
       
        const Real vfrac = vf(i,j,k);

        // SOLID MARKER
        ebMarkers(i, j, k, 0) = vfrac < vfracmin ;
        ebMarkers(i, j, k, 1) = (vfrac < vfracmax) && (vfrac > vfracmin);

      });
    }
  }
////////////////////////////////////////////////////////////////////////////
  void inline ebflux(const Geometry& geom, const MFIter& mfi,
                     const Array4<Real>& prims, const Array4<Real>& flx,
                     const Array4<Real>& rhs, const Array4<Real>& vfrac,
                     const cls_t* cls) {

    const Box& ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

    Array4<const Real> const& apx = areafrac[0]->const_array(mfi);
    Array4<const Real> const& apy = areafrac[1]->const_array(mfi);
#if (AMREX_SPACEDIM==3)    
    Array4<const Real> const& apz = areafrac[2]->const_array(mfi);
#endif    
    
    amrex::ParallelFor(
          ebbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

          // is a covered cell?

          Real vfracinv = 1.0/vfrac(i,j,k);
          for (int n = 0; n < cls_t::NCONS; n++) {
            // fluxes  from fluid faces
            Real fxp = flx(i+1,j,k,n); 
            Real fxm = flx(i,j,k,n);
            Real fyp = flx(i,j+1,k,n); 
            Real fym = flx(i,j,k,n);

#if (AMREX_SPACEDIM==2)            
            rhs(i, j, k, n) += vfracinv*(
            dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
            dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) );
#endif            
#if (AMREX_SPACEDIM==3)            
            Real fzp = flx(i,j,k+1,n); 
            Real fzm = flx (i,j,k,n);
            rhs(i, j, k, n) += vfracinv*(
            dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
            dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) +
            dxinv[2] * (apz(i, j, k + 1) * fzp - apz(i, j, k) * fzm) );
#endif      
          }      
 
          // wall fluxes ...(isothermal, adiabatic)
          // wall_flux(i,j,k,primscell,flux)

          });
  }                      
////////////////////////////////////////////////////////////////////////////
 // wall flux


};


#endif
