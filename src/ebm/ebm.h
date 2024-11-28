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


  // volfrac:: is in a single-component MultiFab. Data are in the range of  with zero representing
  //  covered cells and one for regular cells
  //
  // normbc:: Boundary normal is in a MultiCutFab with AMREX_SPACEDIM components 
  // representing the unit vector pointing toward the covered part.
  //
  // areafrac:: are returned in an Array of MultiCutFab pointers. 
  // For each direction, area fraction is for the face of that direction.
  // Data are in the range of  with zero representing a covered face and 
  // one an un-cut face.

  //** */
  const amrex::MultiCutFab* normbc;  // oboselete
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac; // obsolete

  //std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;
  //const amrex::MultiCutFab* bndrycent;

  /// ***

  // pointer to Amr class instance
  Amr* amr_p;

  // define arrays of multifab to store markers (IBM-style)
  Vector<EBMultiFab<bool>*> bmf_a;       

  // flag array
  Vector<const FabArray<EBCellFlagFab>*> ebflags_a;

  // multifab  pointers to vfrac 
  Vector<const MultiFab*> volmf_a;      

  // multicutfab pointers to area and normals to surface
  Vector<std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM>> areamcf_a;
  Vector<const MultiCutFab*> normmcf_a;

  const amrex::Real isodxerr = 1.e-8;        // relative isotropy error ||dx-dy||/dx
  const amrex::Real vfracmin = 1.e-6;        // minimum to consider a cell "empty"
  const amrex::Real vfracmax = 1.0-vfracmin; // maximum to consider a cell "not-full"
  
  ///////////////////////////////////////////////////////////////////////////
  void init(Amr* pointer_amr,const Geometry& geom, const int required_level, 
                                                         const int max_level)
  {
    BL_PROFILE("initializeEB2");

    amrex::Print() << " Initialize EB at maxlevel = " << max_level << std::endl;

   static_assert(AMREX_SPACEDIM > 1, "EB only supports 2D and 3D"); 

    Vector<std::string> amrex_defaults(
    {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser", "stl"});

    ParmParse ppeb2("eb2");
    std::string geom_type = "all_regular";
    ppeb2.query("geom_type", geom_type);

    // make sure dx == dy == dz if use EB
    if (geom_type != "all_regular") {      
      const Real* dx = geom.CellSize();
      if (AMREX_D_TERM(,  std::abs(dx[0] - dx[1]) > isodxerr * dx[0],
                       || std::abs(dx[0] - dx[2]) > isodxerr * dx[0])) {
      amrex::Abort("EB must have dx == dy == dz (for cut surface fluxes)\n");
      }
    }
 
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
    
    // size of multifab arrays
    bmf_a.resize(max_level + 1);
    volmf_a.resize(max_level + 1);

    // size of flag array
    ebflags_a.resize(max_level + 1);

    // size of multiCut Arrays
    areamcf_a.resize(max_level + 1);
    normmcf_a.resize(max_level + 1);
        
  }
  ////////////////////////////////////////////////////////////////////////////
  // create EBMultiFabs at a level and store pointers to it
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    // markers multifab 
    bmf_a[lev]   = new EBMultiFab<bool>(bxa, dm, 2, cls_t::NGHOST);       
  }
  ////////////////////////////////////////////////////////////////////////////
  void destroy_mf(int lev)
  {
    if (!bmf_a.empty())   { delete bmf_a.at(lev); }
  }
  ////////////////////////////////////////////////////////////////////////////
  void computeMarkers(int lev)
  {

    // amrex::Print( ) << " oo computeMarkers  lev=" << lev << std::endl;  

    auto& mfab = *bmf_a[lev];

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      auto& ibFab = mfab.get(mfi);
      const Box& bx = mfi.tilebox();
      const auto& ebMarkers = mfab.array(mfi); // boolean array
     
      // Array4<const Real> vf = (*volfrac).const_array(mfi); // obsolete     
      Array4<const Real> vf = (*volmf_a[lev]).const_array(mfi);      

      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                       
        const Real vfrac = vf(i,j,k);
        // solid marker
        ebMarkers(i, j, k, 0) = vfrac < vfracmin ;
        // next-to-solid marker:
        // cell is nor empty (empty cell vfrac=1) nor solid (vfrac=0)
        ebMarkers(i, j, k, 1) = (vfrac < vfracmax) && (vfrac > vfracmin);        
      });
      
      // correct for empty cells close to surface 
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                       
#if (AMREX_SPACEDIM==2)             
        ebMarkers(i, j, k, 1) =  ! ebMarkers(i, j, k, 0) && 
        ( ebMarkers(i, j, k, 1) ||
        ebMarkers(i+1, j, k, 0) || ebMarkers(i-1, j, k, 0) || 
        ebMarkers(i, j+1, k, 0) || ebMarkers(i, j-1, k, 0) );
#endif
#if (AMREX_SPACEDIM==3)             
        ebMarkers(i, j, k, 1) = ! ebMarkers(i, j, k, 0) && 
        ( ebMarkers(i, j, k, 1) ||
        ebMarkers(i+1, j, k, 0) || ebMarkers(i-1, j, k, 0) || 
        ebMarkers(i, j+1, k, 0) || ebMarkers(i, j-1, k, 0) ||
        ebMarkers(i, j, k+1, 0) || ebMarkers(i, j, k-1, 0) );
#endif
      });


    }
  }
////////////////////////////////////////////////////////////////////////////
  void inline ebflux(const Geometry& geom, const MFIter& mfi,
                     const Array4<Real>& prims, const Array4<Real>& flx,
                     const Array4<Real>& rhs, const cls_t* cls, int lev) {


   amrex::Print( ) << " oo ebflux  lev=" << lev << std::endl;  

    const Box& ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

    // extract EB arrays given a level and mfi
    Array4<const Real> vfrac = (*volmf_a[lev]).const_array(mfi);  
                     
    Array4<const Real> const& apx = areamcf_a[lev][0]->const_array(mfi);
    Array4<const Real> const& apy = areamcf_a[lev][1]->const_array(mfi);
#if (AMREX_SPACEDIM==3)    
    Array4<const Real> const& apz = areamcf_a[lev][2]->const_array(mfi);
#endif    
        
    // markers 
    const auto& ebMarkers = (*bmf_a[lev]).array(mfi);

    // just in case 
    const auto& flag = (*ebflags_a[lev])[mfi];

    amrex::ParallelFor(
          ebbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

          // only applied to covered cells (could be doen with flags)

          if (ebMarkers(i,j,k,1)){
            Real vfracinv = 1.0/vfrac(i,j,k);
            for (int n = 0; n < cls_t::NCONS; n++) {
              // fluxes  from fluid faces
              Real fxp = flx(i+1,j,k,n); 
              Real fxm = flx(i,j,k,n);
              Real fyp = flx(i,j+1,k,n); 
              Real fym = flx(i,j,k,n);

              // warning overwrite rhs in the cut-cells
#if (AMREX_SPACEDIM==2)            
              rhs(i, j, k, n) = -vfracinv*(
              dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
              dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) );
#endif            
#if (AMREX_SPACEDIM==3)            
              Real fzp = flx(i,j,k+1,n); 
              Real fzm = flx (i,j,k,n);
              rhs(i, j, k, n) = -vfracinv*(
              dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
              dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) +
              dxinv[2] * (apz(i, j, k + 1) * fzp - apz(i, j, k) * fzm) );
#endif      
            }

          // add wall fluxes
                   
          // wall fluxes ...(isothermal, adiabatic)
          // wall_flux(i,j,k,primscell,flux)

          }

          });
  }                      
////////////////////////////////////////////////////////////////////////////
 // wall flux


};


#endif
