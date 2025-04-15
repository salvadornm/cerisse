#ifndef ebm_H_
#define ebm_H_


#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_ParmParse.H>

// snm new
#include <AMReX_EBFluxRegister.H>
#include <AMReX_Geometry.H>

#include <AMReX_FArrayBox.H>


#include <algorithm>
#include <cmath>

#include <EBMultiFab.h>

#include <walltypes.h>

#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_Redistribution.H>


#include <CNSconstants.h>

// temp
#include <FluxRedistribute.h>

template <typename wallmodel,typename cls_t>
class ebm_t
{ 

public:	

  ebm_t() { 
  }

  ~ebm_t() {}


  // volfrac     :: is in a single-component MultiFab. Data are in the range of [0,1] with zero representing
  //                covered cells and one for regular cells
  //
  // normbc      :: Boundary normal is in a MultiCutFab with AMREX_SPACEDIM components 
  //                representing the unit vector pointing toward the covered part.
  //
  // areafrac    :: are returned in an Array of MultiCutFab pointers. 
  //                For each direction, area fraction is for the face of that direction.
  //                Data are in the range [0,1] of  with zero representing a covered face and 
  //                one an un-cut face.
  //
  // bndryarea  :: is a MultiCutFab with a single component representing the dimensionless boundary area. 
  //               when the cell is isotropic (i.e., ), it’s trivial to convert it to physical units. (*dx)
  //
  // bndrycent_a:: is a MultiCutFab with AMREX_SPACEDIM components, each component is [-0.5:0.5] respect
  //               to reguar centere

  
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

  // multicutfab pointers to relevant geometric numbers
  Vector<std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM>> areamcf_a;
  Vector<const MultiCutFab*> normmcf_a;
  Vector<const MultiCutFab*> bcareamcf_a;
  Vector<const MultiCutFab*> bndrycent_a;
  // mask
  Vector<iMultiFab> level_mask_a;  // object not pointers
 
  const amrex::Real isodxerr = 1.e-8;        // relative isotropy error ||dx-dy||/dx
  const amrex::Real vfracmin = 1.e-6;        // minimum to consider a cell "empty"
  const amrex::Real vfracmax = 1.0-vfracmin; // maximum to consider a cell "not-full"

  // variables for redistribution (originally declared static)
  amrex::Real eb_weight;
  std::string redistribution_type;

  // aux vars for redistribution
  const bool use_wts_in_divnc = false; //true
  const int srd_max_order = 2; // 2
  const amrex::Real target_volfrac = 0.5; //0.5
  const amrex::Real fac_for_deltaR = 1.0;

  
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
      amrex::Abort( " EB Type .. not ready yet (use cerisse1)"); 
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
    bcareamcf_a.resize(max_level + 1); 
    bndrycent_a.resize(max_level + 1); 

    // size of mask
    level_mask_a.resize(max_level + 1); 
        
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
  /**
  * @brief Update boolean markers solid and partially solid
  * @param lev current AMR level  
  **/
  void computeMarkers(int lev)
  {

    auto& mfab = *bmf_a[lev];

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      
      // make box including GHOST
      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);

      const auto& ebMarkers = mfab.array(mfi);           
      const auto& flag_arr = (*ebflags_a[lev]).const_array(mfi);

      amrex::ParallelFor(
        bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                       
        ebMarkers(i, j, k, 0) = flag_arr(i,j,k).isCovered();  
        ebMarkers(i, j, k, 1) = flag_arr(i,j,k).isSingleValued();    
      });

    }
  }
  ////////////////////////////////////////////////////////////////////////////
  /**
  * @brief Compute fluxes in the EB wall and add it to rhs
  * @param prims array of primtive variables
  * @param flxt   fluxes across faces (convective + viscous)
  **/
  void inline ebflux( const Geometry& geom, const MFIter& mfi,
                      const Array4<Real>& prims, 
                      std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                      const Array4<Real>& rhs, const cls_t* cls, int lev) {

    const Box& ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Real *dx = geom.CellSize();

    // extract EB arrays given a level and mfi
    Array4<const Real> vfrac = (*volmf_a[lev]).const_array(mfi);             // vfrac
    Array4<const Real> const& apx = areamcf_a[lev][0]->const_array(mfi);     // areas free in faces x
    Array4<const Real> const& apy = areamcf_a[lev][1]->const_array(mfi);     // faces in y 
#if (AMREX_SPACEDIM==3)    
    Array4<const Real> const& apz = areamcf_a[lev][2]->const_array(mfi);     //   faces in z
#endif    
    Array4<const Real> const& normxyz = (*normmcf_a[lev]).const_array(mfi);  // normal to surface (NDIM components)     
    Array4<const Real> const& bcarea = (*bcareamcf_a[lev]).const_array(mfi); // area wall
            
    // markers 
    const auto& ebMarkers = (*bmf_a[lev]).array(mfi);

    // just in case 
    const auto& flag = (*ebflags_a[lev])[mfi];

    auto const& flx_x = flxt[0]->array(); 
    auto const& flx_y = flxt[1]->array(); 
#if (AMREX_SPACEDIM==3)     
    auto const& flx_z = flxt[2]->array(); 
#endif

    amrex::ParallelFor(
        ebbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

          // only applied to covered cells (could be done with flags)

          if (ebMarkers(i,j,k,1)){
            Real vfracinv = 1.0/vfrac(i,j,k);

            // rebuild fluxes
            for (int n = 0; n < cls_t::NCONS; n++) {
              Real fxp = flx_x(i+1,j,k,n); 
              Real fxm = flx_x(i,j,k,n);
              Real fyp = flx_y(i,j+1,k,n); 
              Real fym = flx_y(i,j,k,n);

              //  overwrite rhs in the cut-cells

#if (AMREX_SPACEDIM==3)            
              Real fzp = flx_z(i,j,k+1,n); 
              Real fzm = flx_z(i,j,k,n);
              rhs(i, j, k, n) = -vfracinv*(
                dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
                dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) +
                dxinv[2] * (apz(i, j, k + 1) * fzp - apz(i, j, k) * fzm) );
#else
              rhs(i, j, k, n) = -vfracinv*(
                dxinv[0] * (apx(i + 1, j, k) * fxp - apx(i, j, k) * fxm) +
                dxinv[1] * (apy(i, j + 1, k) * fyp - apy(i, j, k) * fym) );            
#endif      
            }

            // build wall fluxes
            amrex::GpuArray<Real, cls_t::NCONS> flux_wall = {0.0};                             
            // primitive array at surface
            amrex::GpuArray<Real, cls_t::NPRIM> prim_wall = {0.0};
            for (int n = 0; n < cls_t::NPRIM; n++) {
              prim_wall[n] = prims(i,j,k,n);   // interpolate   ???????????   
              
              // use same tecniques weighted based on distance phi = sum w phi(node)/sum w
              // w is 1/r (only connected)


            }  
            // normal to surface 
            amrex::Real norm_wall [AMREX_SPACEDIM]= {0.0};
   
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              norm_wall[n] = -normxyz(i,j,k,n); // so the normal points towards the liquid
            }

            // calculate wall flux and add it to rhs
            wallmodel::wall_flux(geom,i,j,k,norm_wall,prim_wall,flux_wall,cls);             

            for (int n = 0; n < cls_t::NCONS; n++) {
              rhs(i,j,k,n) += flux_wall[n]*vfracinv*bcarea(i,j,k,0)*dxinv[0]; 
            }


          // snm  
          // if ((i==120) && (j==319))
          //  {
          //   printf(" in flux wall \n");
          //   for (int n = 0; n < cls_t::NCONS; n++) {
          //     printf(" n=%d flux=%f  \n",n,flux_wall[n]);
          //   }
          //   printf(" vfrac=%f bcarea=%f",vfrac(i,j,k),bcarea(i,j,k,0));
          //   for (int n = 0; n < cls_t::NPRIM; n++)
          //   {
          //     printf(" n=%d qwall=%f  \n",n,prim_wall[n]);
          //   }

          //  }
           // snm


          }
         

        });

  }                      
  ///////////////////////////////////////////////////////////////////////////
  /**
  * @brief Compute fluxes in the EB wall and add it to rhs
  * @param cons  array of conservative variables
  * @param divc  array with divergence of flux (initially rhs) includes wall fluxes
  * @param flxt  fluxes across faces (convective + viscous)
  * @param rhs  
  **/

  //(geom,mfi,cons,divc, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},
  //state, cls_d,level,dt,h_phys_bc);

  void inline redist (const Geometry& geom, const MFIter& mfi,
                      const Array4<Real>& cons, const Array4<Real>& divc, 
                      std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                      const Array4<Real>& rhs, const cls_t* cls, int lev, Real dt,
                      BCRec const* phys_bc) {

    const Box& ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Real *dx = geom.CellSize();

    // printf(" [ebm::redist] in ebm redist... (SNM temp) \n ");

                    
    // extract EB-related arrays given a level and mfi
    Array4<const Real> vfrac = (*volmf_a[lev]).const_array(mfi);             // vfrac
    Array4<const Real> const& apx = areamcf_a[lev][0]->const_array(mfi);     // areas fraction x
    Array4<const Real> const& apy = areamcf_a[lev][1]->const_array(mfi);     //                y 
#if (AMREX_SPACEDIM==3)    
    Array4<const Real> const& apz = areamcf_a[lev][2]->const_array(mfi);     //                z
#endif    
    Array4<const Real> const& bcarea  = (*bcareamcf_a[lev]).const_array(mfi); // bc area 
    Array4<const Real> const& bcent   = (*bndrycent_a[lev]).const_array(mfi); // bc centroid

    Array4<const int> const& lev_mask = level_mask_a[lev].const_array(mfi);   // level mask is an object (not a pointer)

    // flags
    //const auto& flag = (*ebflags_a[lev])[mfi];

    Array4<const EBCellFlag> const& flag = (*ebflags_a[lev]).const_array(mfi); ;

    //flags.const_array(mfi)  /<------------------------------------------
    // Array4<const EBCellFlag> const& fla

    // store weights and ebwight in two arrays to prep for redist
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);    
    // redistribution weights arrays
    FArrayBox redistwgt_fab(bxg, 1);       
    FArrayBox srd_update_scale_fab(bxg, 1);  
    auto const& redistwgt = redistwgt_fab.array();
    auto const& srd_update_scale = srd_update_scale_fab.array();
     
    amrex::ParallelFor(bxg, [=, this]  AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      redistwgt(i, j, k) = vfrac(i, j, k);
      srd_update_scale(i, j, k) = eb_weight; 
    });


    
    // create temporary array and fill first with weights (will be used as weights)
    FArrayBox tmpfab(bxg, 1, The_Async_Arena()); // cls_t::NCONS
    Array4<Real> scratch = tmpfab.array();    
    if (redistribution_type == "FluxRedist") {
      amrex::ParallelFor(bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        scratch(i, j, k) = redistwgt(i, j, k);
      });
    }
    // call MLRedistribution from AMReX
    int level_mask_not_covered = CNSConstants::level_mask_notcovered;

    // detail call                                                     AMReX notation           
    // ebbox                  :: box where we do the redist               bx
    // cls_t::NCONS           :: numbers of variables to solve            ncomp
    // rhs                    :: (output) rhs                             dUdt_out
    // divc                   :: initial rhs (at current step)            dUdt_in        
    // cons                    :: array of cons                           U_in
    // scratch                :: temp array                               scratch
    // flag                   :: flag array  (EBCellFlag)                 flag  
    // apx,apy,apz            :: cell openings (weighted by mesh)         apx,apy,apz
    // vfrac                  :: fraction of cell filled with fluid       vfrac
    // fcx,fcy,fz             :: fluxes                                   fcx,fcy,fz 
    // bcent                  :: Array4 of boundary centroid              ccc
    // physbc                 :: physical bc pointer                      d_bcrec_ptr
    // geom                   :: Geometry                                 lev_geom
    // dt                     :: time step                                dt  
    // redistribution_type    :: string StateRedist/FluxRedist/NoRedist   redistribution_type 
    // as_crse                :: int  1(true)/ 0(false) if false 
    // rr_drho_crse           :: p_drho_as_crse->array()                            rr_drho_crse
    // rr_flag_crse           :: p_rrflag_as_crse->const_array()                    rr_flag_crse   
    // as_fine                :: int  1(true)/0 (false)
    // dm_as_fine             :: Array4 Real  (not clear)  
    // lev_mask               :: Arry4 int    (result mask IntFab)                 lev_mask
    // level_mask_not_covered :: ghost cells not covered (const set as 2)           level_mask_not_covered
    // fac_for_deltaR         ::  (const set as as 1)                               fac_for_deltaR 
    // use_wts_in_divnc       :: use weights to computeRHS (constant, set as true)  use_wts_in_divnc 
    // 0                      ::  start loop of vars (set to 0)                     icomp
    // srd_max_order          ::  (const int set as 2)                              srd_max_order
    // target_volfrac         ::  (const set as 0.5)                                target_volfrac
    // srd_update_scale       :: Array of eb_weight (usually 1)                     srd_update_scale

    int as_crse = 0;
    int as_fine = 0;  // if 1 it crashes
    FArrayBox dm_as_fine(Box::TheUnitBox(), cls_t::NCONS, The_Async_Arena());
    FArrayBox fab_drho_as_crse(Box::TheUnitBox(), cls_t::NCONS, The_Async_Arena());
    IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
    // in cerisse this call is different, depends on
    const IArrayBox* p_rrflag_as_crse = &fab_rrflag_as_crse;
    FArrayBox* p_drho_as_crse = &fab_drho_as_crse;
        
    auto const& fcx = flxt[0]->array(); 
    auto const& fcy = flxt[1]->array(); 
#if (AMREX_SPACEDIM==3)     
    auto const& fcz = flxt[2]->array(); 
#endif
    
    // redistribution
    if (redistribution_type == "NewRedist")
    {
      cerisse_flux_redistribute( ebbox,rhs, divc, redistwgt, vfrac,flag,geom,cls_t::NCONS,dt);
    }
    else  
    {
      amrex::ApplyMLRedistribution(
      ebbox, cls_t::NCONS, rhs, divc, cons, scratch, flag, AMREX_D_DECL(apx, apy, apz), vfrac,
      AMREX_D_DECL(fcx, fcy, fcz), bcent, phys_bc, geom, dt, redistribution_type,
      as_crse, p_drho_as_crse->array(), p_rrflag_as_crse->const_array(), as_fine, dm_as_fine.array(), lev_mask,
      level_mask_not_covered, fac_for_deltaR, use_wts_in_divnc, 0, srd_max_order,
      target_volfrac, srd_update_scale);
    }                                
  }  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
};


#endif

// intrepolation based on distance
// temp
// #if (AMREX_SPACEDIM == 2)
//         int kk(0);
// #else
//         for (int kk = -nb; kk <= nb; kk++) {
// #endif
//         for (int jj = -nb; jj <= nb; jj++) {
//         for (int ii = -nb; ii <= nb; ii++) {

