#ifndef ebm_H_
#define ebm_H_


#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_ParmParse.H>

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

#include "custom_geometry.h"
#include <FluxRedistribute.h>

template <typename wallmodel, typename param, typename cls_t>
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
  // bndrycent_a:: (embedded boundary centroid)
  //               MultiCutFab with AMREX_SPACEDIM components, each component is [-0.5:0.5] 
  //               respect to regular centere
  // volcent_a:;   (volume centroid)
  //               MultiCutFab with AMREX_SPACEDIM components, each component is [-0.5:0.5] 
  //               respect to regular centere
  

  
  //std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;
  //const amrex::MultiCutFab* bndrycent;

  /// ***

  // pointer to Amr class instance
  Amr* amr_p;

  // define arrays of multifab to store markers (IBM-style)
  Vector<EBMultiFab<bool>*> bmf_a;       
  //Vector<IBMultiFab<uint8_t, GPDATA>*> bmf_a;  

  // flag array
  Vector<const FabArray<EBCellFlagFab>*> ebflags_a;

  // multifab  pointers to vfrac 
  Vector<const MultiFab*> volmf_a;      

  // multicutfab pointers to relevant geometric numbers
  Vector<std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM>> areamcf_a;
  Vector<const MultiCutFab*> normmcf_a;
  Vector<const MultiCutFab*> bcareamcf_a;
  Vector<const MultiCutFab*> bndrycent_a;
  Vector<const MultiCutFab*> volcent_a;
  // mask
  Vector<iMultiFab> level_mask_a;  // object not pointers
 
  // variables for redistribution (originally declared static)
  amrex::Real eb_weight;
  std::string redistribution_type;

  // EB constants and parameters --------------------------------------------------------
  static constexpr Real isodxerr = 1.e-8;        // relative isotropy error ||dx-dy||/dx
  //  for redistribution
  const bool use_wts_in_divnc = false; //true
  const int srd_max_order = 2; // 2
  const Real target_volfrac = 0.5; //0.5
  const Real fac_for_deltaR = 1.0;

  //  for interpolation
  static constexpr bool use_weighted_interp = false;  //false: unweighted average; true: inverse-distance weighted average
  static constexpr int nb = 1; // number of neighbours for interpolation     
  static constexpr Real vfracmin = 1.e-9; // minimum vfrac to consider a cell "not empty" 
  static constexpr Real vfracmax = 1.0 - vfracmin; // maximum vfrac to consider a cell "not full"

  ///////////////////////////////////////////////////////////////////////////
  void init(Amr* pointer_amr,const Geometry& geom, const int required_level, 
                                                         const int max_level)
  {
    BL_PROFILE("initializeEB2");

    amrex::Print() << " Initialize EB at level = " << required_level ;
    amrex::Print() << " out of "<< max_level << std::endl;

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
      // Non-AMReX default EB types, get from CustomGeometry
      auto geometry = CustomGeometry::create(geom_type);
      geometry->build(geom, max_level);
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
    volcent_a.resize(max_level + 1); 
    

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
  *        i,j,k,0):  if cells is covered                         :1 (true) 0 (false)
  *        i,j,k,1):  if cells is partially covered  (near solid) :1 (true) 0 (false)
  * @param lev current AMR level  
  **/
  void computeMarkers(int lev)
  {

    auto& mfab = *bmf_a[lev];

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      
      // make box including GHOST
      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);

      const auto& flag_arr = (*ebflags_a[lev]).const_array(mfi);
      const auto vfrac     = (*volmf_a[lev]).const_array(mfi);    

      // markers: comp 0 = solid mask, comp 1 = neighbor-of-solid mask
      const auto& ebMarkers = mfab.array(mfi);           


      const bool correct_cells = false; //treat small cells as solid
      // -------------------------
      // Pass 1: define "solid"
      // -------------------------
      int ncorr=0;
      amrex::Gpu::DeviceScalar<int> ncorr_d(0);
      int* p_ncorr = ncorr_d.dataPtr();

      amrex::ParallelFor( bxg, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          // Start from AMReX EB classification
          int is_solid = flag_arr(i,j,k).isCovered() ? 1 : 0;
          // If it's a cut cell with very small fluid volume, treat it as solid
          if (correct_cells)
          {
            if (flag_arr(i,j,k).isSingleValued() && (vfrac(i,j,k) < vfracmin)) {
              is_solid = 1; 
	      amrex::Gpu::Atomic::Add(p_ncorr, 1);
            }
          }  
          ebMarkers(i,j,k,0) = is_solid;
          // Initialize neighbor flag; will be rebuilt in Pass 2          
          ebMarkers(i,j,k,1) = flag_arr(i,j,k).isSingleValued(); //old      

        });
      ncorr = ncorr_d.dataValue();  // bring device counter back to host
      // ---------------------------------------------------------
      // Pass 2: rebuild "neighbor-of-solid" using updated marker0
      // ---------------------------------------------------------
      
      if (ncorr > 0)
      {
      // We need i±1/j±1/k±1 accesses; avoid OOB by shrinking one layer.
      amrex::Box bx_inner = bxg;
      bx_inner.grow(-1);

      if (bx_inner.ok()) // is box ok
      {
        amrex::ParallelFor( bx_inner, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
          // Only fluid cells can be "neighbors of solid"
          if (ebMarkers(i,j,k,0) != 0) {
            ebMarkers(i,j,k,1) = 0;
            return;
          }

          int nbr_solid = 0;
          // Face-neighbors (6-neighborhood in 3D, 4-neighborhood in 2D)
          nbr_solid |= (ebMarkers(i-1,j,k,0) != 0);
          nbr_solid |= (ebMarkers(i+1,j,k,0) != 0);
          nbr_solid |= (ebMarkers(i,j-1,k,0) != 0);
          nbr_solid |= (ebMarkers(i,j+1,k,0) != 0);
#if (AMREX_SPACEDIM == 3)
          nbr_solid |= (ebMarkers(i,j,k-1,0) != 0);
          nbr_solid |= (ebMarkers(i,j,k+1,0) != 0);
#endif
          ebMarkers(i,j,k,1) = nbr_solid ? 1 : 0;
          });
      } // endif box ok
      } // endif ncorr
    }  // end loop mfi

      // Optional: if you rely on marker(.,.,.,1) being valid in  ALL ghost cells
      // need to mfab.FillBoundary(geom[lev].periodicity()) outside the MFIter
  }
  ////////////////////////////////////////////////////////////////////////////
  /**
  * @brief check if the geometry is properly defined by printing out some geometric parameters
  * @param lev current AMR level  
  **/
  void check_geometry (int lev)
  {

    auto& mfab = *bmf_a[lev];

    amrex::Print() << " Check EB geometry at level " << lev << "\n";

    // init counters
    int empty_cutcells = 0;
    int distorted_cutcells = 0;
    int corr_cutcells = 0;

    amrex::Gpu::DeviceScalar<int> empty_d(0);
    amrex::Gpu::DeviceScalar<int> distorted_d(0);
    amrex::Gpu::DeviceScalar<int> corr_d(0);

    int* p_empty     = empty_d.dataPtr();
    int* p_distorted = distorted_d.dataPtr();
    int* p_corr      = corr_d.dataPtr();

    for (amrex::MFIter mfi(mfab, false); mfi.isValid(); ++mfi)
    {
      const Box& ebbox  = mfi.growntilebox(0);  // box without ghost points 

      const auto& flag = (*ebflags_a[lev])[mfi];
      FabType t = flag.getType(ebbox);
      const bool fab_with_eb     = (FabType::singlevalued == t);  

      if (!fab_with_eb) { continue;}

      // geometry arrays for this tile
      const auto apx       = areamcf_a[lev][0]->const_array(mfi);
      const auto apy       = areamcf_a[lev][1]->const_array(mfi);
#if (AMREX_SPACEDIM == 3)
      const auto apz       = areamcf_a[lev][2]->const_array(mfi);
#endif
      const auto vfrac     = (*volmf_a[lev]).const_array(mfi);    
      const auto flag_arr  = (*ebflags_a[lev]).const_array(mfi);           
      const auto normxyz   = (*normmcf_a[lev]).const_array(mfi);
      const auto bcarea    = (*bcareamcf_a[lev]).const_array(mfi);
      // markers
      const auto& ebMarkers = mfab.array(mfi);      

      amrex::ParallelFor(
        ebbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {

        if (ebMarkers(i,j,k,1) == 0) return; // not neighbour of solid


        // normal to surface (pointing towards the fluid)
        amrex::Real norm_wall[AMREX_SPACEDIM] = {0.0_rt};
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
          norm_wall[n] = -normxyz(i,j,k,n);
        }
        amrex::Real areaw = bcarea(i,j,k,0);
        // calculate error divergence 
        amrex::Real Err_A[AMREX_SPACEDIM] = {0.0_rt};
        Err_A[0] = apx(i,j,k)   - apx(i+1,j,k)   + areaw*norm_wall[0];
        Err_A[1] = apy(i,j,k)   - apy(i,j+1,k)   + areaw*norm_wall[1];
#if (AMREX_SPACEDIM == 3)
        Err_A[2] = apz(i,j,k)   - apz(i,j,k+1)   + areaw*norm_wall[2];
#endif
        amrex::Real sumError = 0.0_rt;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {sumError += Err_A[n];} 

        // check if vfrac is in the expected range
        if (vfrac(i,j,k) < vfracmin)  {
           //++empty_cutcells;          
	  amrex::Gpu::Atomic::Add(p_empty, 1);
        }                

        if (amrex::Math::abs(sumError) > 1.e-6_rt) {
           //++distorted_cutcells;
	  amrex::Gpu::Atomic::Add(p_distorted, 1);
        }        

        // these cells were corrected because they were empty
        if (!flag_arr(i,j,k).isSingleValued())
        {
          //++corr_cutcells;
	  amrex::Gpu::Atomic::Add(p_corr, 1);
        }


      });

    }


    empty_cutcells     = empty_d.dataValue();
    distorted_cutcells = distorted_d.dataValue();
    corr_cutcells      = corr_d.dataValue();

    // sum across MPI ranks
    amrex::ParallelDescriptor::ReduceIntSum(empty_cutcells);
    amrex::ParallelDescriptor::ReduceIntSum(distorted_cutcells);
    amrex::ParallelDescriptor::ReduceIntSum(corr_cutcells);
    
    amrex::Print() << " Total number of empty cut-cells:  " << empty_cutcells ;
    amrex::Print() << " and distorted cut-cells: " << distorted_cutcells << "\n";
    amrex::Print() << " empty cells corrected:   " << corr_cutcells << "\n";
    amrex::Print() << " ------------------------------------ \n";
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

    // avoid pointers
    auto dxinv = geom.InvCellSizeArray();
    auto dx    = geom.CellSizeArray();

    // extract EB arrays given a level and mfi
    Array4<const Real> vfrac = (*volmf_a[lev]).const_array(mfi);              // vfrac
    Array4<const Real> const& apx = areamcf_a[lev][0]->const_array(mfi);      // areas free in faces x
    Array4<const Real> const& apy = areamcf_a[lev][1]->const_array(mfi);      // faces in y 
#if (AMREX_SPACEDIM==3)    
    Array4<const Real> const& apz = areamcf_a[lev][2]->const_array(mfi);      // faces in z
#endif    
    Array4<const Real> const& normxyz = (*normmcf_a[lev]).const_array(mfi);      
    Array4<const Real> const& bcarea  = (*bcareamcf_a[lev]).const_array(mfi); 
    Array4<const Real> const& bc_centroid  = (*bndrycent_a[lev]).const_array(mfi); 
    Array4<const Real> const& vol_centroid = (*volcent_a[lev]).const_array(mfi);   

    // markers 
    const auto& ebMarkers = (*bmf_a[lev]).array(mfi);
    // fluxes
    auto const& flx_x = flxt[0]->array(); 
    auto const& flx_y = flxt[1]->array(); 
#if (AMREX_SPACEDIM==3)     
    auto const& flx_z = flxt[2]->array(); 
#endif

#ifdef USE_PELEPHYSICS    
    // transport properties
    auto const* ltransparm = trans_parms.device_parm();
    AMREX_ALWAYS_ASSERT(ltransparm != nullptr);
#else
    auto const* ltransparm = (trans_parm_t const*)nullptr;  //null pointer      
#endif    

    amrex::ParallelFor(
        ebbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

          // only solve flux in cut cells (not in regular or almost empty cells)
          bool solve_fluxwall = ebMarkers(i,j,k,1) && (vfrac(i,j,k) > vfracmin); 

          if (solve_fluxwall){            
            Real inv_hvfrac = dxinv[0]/(vfrac(i,j,k)+vfracmin); //only isotropic cells

            // rebuild fluxes in the cut cells
            for (int n = 0; n < cls_t::NCONS; n++) {
              
#if (AMREX_SPACEDIM==3)            
              rhs(i, j, k, n) = -inv_hvfrac*(
                (apx(i + 1, j, k) * flx_x(i+1,j,k,n)  - apx(i, j, k) * flx_x(i,j,k,n)) +
                (apy(i, j + 1, k) * flx_y(i,j+1,k,n)  - apy(i, j, k) * flx_y(i,j,k,n)) +
                (apz(i, j, k + 1) * flx_z(i,j,k+1,n)  - apz(i, j, k) * flx_z(i,j,k,n)) );
#else
              rhs(i, j, k, n) = -inv_hvfrac *(
                (apx(i + 1, j, k) * flx_x(i+1,j,k,n) - apx(i, j, k) * flx_x(i,j,k,n)) +
                (apy(i, j + 1, k) * flx_y(i,j+1,k,n) - apy(i, j, k) * flx_y(i,j,k,n)) );            
#endif      
            }

            // build wall fluxes
            amrex::GpuArray<Real, cls_t::NCONS> flux_wall = {0.0};                             
            // primitive array at surface
            amrex::GpuArray<Real, cls_t::NPRIM> prim_wall = {0.0};
      
            if (use_weighted_interp){
              Real sumw = 1e-30;                
              // position of centroid in cell units
              Real x_bc = i + 0.5 + bc_centroid(i,j,k,0);
              Real y_bc = j + 0.5 + bc_centroid(i,j,k,1) ;
#if (AMREX_SPACEDIM == 3)
              Real z_bc =  k  + 0.5 + bc_centroid(i,j,k,2);
              for (int kk = k - nb; kk <= k + nb; kk++) {
#else
              int kk = 0;              
#endif
                for (int jj = j - nb; jj <= j + nb; jj++) {
                  for (int ii = i - nb; ii <= i + nb; ii++) {
                  // loop over neighbour cells                          
                    if (!ebMarkers(ii,jj,kk,0)) {
                      // regular or cut-cell neighbour coordinates in cell units
                      Real xx = ii + 0.5;
                      Real yy = jj + 0.5;
#if (AMREX_SPACEDIM == 3)
                      Real zz = kk + 0.5;
#endif
                      if (ebMarkers(ii,jj,kk,1)) { // cut cell (use position of centroid)
                        xx += vol_centroid(ii,jj,kk,0);
                        yy += vol_centroid(ii,jj,kk,1);
#if (AMREX_SPACEDIM == 3)
                        zz += vol_centroid(ii,jj,kk,2);
#endif
                      }              
                      // distance from nighbour cell (ii,jj,kk) centroid  to target cell (i,j,k) boundary centroid (in cell units)
#if (AMREX_SPACEDIM == 3)
                      Real r = std::sqrt( (xx - x_bc)*(xx - x_bc) + (yy - y_bc)*(yy - y_bc) + (zz - z_bc)*(zz - z_bc) );                        
#else                                          
                      Real r = std::sqrt( (xx - x_bc)*(xx - x_bc) + (yy - y_bc)*(yy - y_bc) );
#endif                                                           
                      // interpolation weight based on distance 
                      r = amrex::max(r, 1e-6_rt);
                      Real w = 1.0_rt/r; sumw += w;
                      for (int n = 0; n < cls_t::NPRIM; n++) { prim_wall[n] += w*prims(ii,jj,kk,n);    } 
                    }  //  endif not-empty 
                  } //endfor ii
                } //endfor jj
#if (AMREX_SPACEDIM == 3)                
              } //end for kk
#endif              
              //--o  
              sumw = 1.0/sumw; //normalise weights          
              for (int n = 0; n < cls_t::NPRIM; n++) {
                prim_wall[n] = prim_wall[n]*sumw;             
              }
            }
            else {
              for (int n = 0; n < cls_t::NPRIM; n++) {
                prim_wall[n] = prims(i,j,k,n);               
              }
            }    
            // normal to surface (pointing towards the fluid)
            Real norm_wall[AMREX_SPACEDIM]= {0.0}; 
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              norm_wall[n] = -normxyz(i,j,k,n);  
            }

            // area/normal correction to ensure divergence free fluxes in cut cells 
            //====================================================================
            Real Err_A[AMREX_SPACEDIM]= {0.0};
            Real areaw = bcarea(i,j,k,0);
            Err_A[0] = apx(i,j,k) - apx(i+1,j,k) + areaw*norm_wall[0];
            Err_A[1] = apy(i,j,k) - apy(i,j+1,k) + areaw*norm_wall[1];
#if (AMREX_SPACEDIM == 3)                      
            Err_A[2] = apz(i,j,k) - apz(i,j,k+1) + areaw*norm_wall[2]; 
#endif
            // error measure (temp)
            Real sumError = 0.0_rt; Real sumnorm=0.0_rt;
            for (int n = 0; n < AMREX_SPACEDIM; n++) { sumError += Err_A[n];}
            for (int n = 0; n < AMREX_SPACEDIM; n++) { sumnorm += norm_wall[n];}

            // correct area projections 
            Real Areai[AMREX_SPACEDIM]= {0.0}; 
            for (int n = 0; n < AMREX_SPACEDIM; n++){
              Areai[n] = areaw*norm_wall[n] - Err_A[n];
            }  
            // recalculate normal based on corrected area projections
	    Real areanew = 0.0_rt;
	    for (int n = 0; n < AMREX_SPACEDIM; n++) {
              areanew += Areai[n]*Areai[n];
            }
            areanew = std::sqrt(areanew);
            Real normnew[AMREX_SPACEDIM]= {0.0};
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              normnew[n] = Areai[n]/areanew;
            }
            // update wall normal and area with corrected values
            areaw = areanew;
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              norm_wall[n] = normnew[n];
            }
            //=================================================================

	          wallmodel wm;  // create local instance
            // calculate wall flux and add it to rhs
            wm.wall_flux(geom,i,j,k,norm_wall,prim_wall,flux_wall,cls);      
                                  
            // calculate viscous walls
            if (param::solve_diffwall)
            {
              // from volume and area centroid compute distance to wall 
              // and project into normal direction
	            Real dis = 1e-6_rt; // min distance (cell units)   
              for (int n = 0; n < AMREX_SPACEDIM; n++) {                
                dis += (vol_centroid(i,j,k,n)- bc_centroid(i,j,k,n)) *norm_wall[n];       
              }    
	            dis = dis*dx[0];  // units
	            wm.wall_flux_diff(geom,i,j,k,dis,norm_wall,prims,prim_wall,flux_wall,cls,ltransparm);
            } 

            // add wall flux to rhs 
            for (int n = 0; n < cls_t::NCONS; n++) {
              rhs(i,j,k,n) += flux_wall[n]*areaw*inv_hvfrac; 
            }                           

          }  //end if partially covered       

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

  void inline redist (const Geometry& geom, const MFIter& mfi,
                      const Array4<Real>& cons, const Array4<Real>& divc, 
                      std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                      FluxReg* fr_as_crse, FluxReg* fr_as_fine, Array4<Real> const& dm_as_fine,
                      const Array4<Real>& rhs, const cls_t* cls, int lev, Real dt,
                      BCRec const* phys_bc) {

    const Box& ebbox  = mfi.growntilebox(0); 
    const Box& bxg    = mfi.growntilebox(cls_t::NGHOST);
    auto dx    = geom.CellSizeArray();

    // printf(" [ebm::redist] in ebm redist... (SNM temp) \n ");
                    
    // extract EB-related arrays/flags  given a level and mfi
    Array4<const Real> vfrac = (*volmf_a[lev]).const_array(mfi);               // vfrac
    Array4<const Real> const& apx = areamcf_a[lev][0]->const_array(mfi);       // area  fraction x
    Array4<const Real> const& apy = areamcf_a[lev][1]->const_array(mfi);       //                y 
#if (AMREX_SPACEDIM==3)    
    Array4<const Real> const& apz = areamcf_a[lev][2]->const_array(mfi);       //                z
#endif    
    Array4<const Real> const& bcarea  = (*bcareamcf_a[lev]).const_array(mfi);  // bc area 
    Array4<const Real> const& bcent   = (*bndrycent_a[lev]).const_array(mfi);  // bc centroid
    Array4<const int> const& lev_mask = level_mask_a[lev].const_array(mfi);    // level mask is an object (not a pointer)
    Array4<const EBCellFlag> const& flag = (*ebflags_a[lev]).const_array(mfi); // flags

     // create temporary arrays to store the weights
      FArrayBox redistwgt_fab(bxg, 1);        
      FArrayBox srd_update_scale_fab(bxg, 1); 
      int ncomp = cls_t::NCONS;
      if (redistribution_type == "FluxRedist") ncomp = 1;
      FArrayBox tmpfab(bxg, ncomp, The_Async_Arena());
      auto const& redistwgt = redistwgt_fab.array();
      auto const& srd_update_scale = srd_update_scale_fab.array();
      Array4<Real> scratch = tmpfab.array();
  
    // fill temporary arrays
    Real ebw = eb_weight;   // copy member to a plain scalar
    amrex::ParallelFor(bxg, [=]  AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      redistwgt(i, j, k)        = vfrac(i, j, k);
     srd_update_scale(i, j, k)  = ebw; // eb_weight;
    });

    //redist weights (in case of FluxRedist, otherwise use as scrap inside Reditribution)  // temp snm
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

    // int as_crse = 0;
    // int as_fine = 0;  // if 1 it crashes
    int as_crse = int(fr_as_crse != nullptr);
    int as_fine = int(fr_as_fine != nullptr);
    
    //  FArrayBox dm_as_fine(Box::TheUnitBox(), cls_t::NCONS, The_Async_Arena()); //snm (new 2 lines  below, now commented)
    //FArrayBox dm_as_fine(bxg, cls_t::NCONS, The_Async_Arena());
    //dm_as_fine.setVal<RunOn::Device>(0.0);

    FArrayBox fab_drho_as_crse(Box::TheUnitBox(), cls_t::NCONS, The_Async_Arena()); //snm ( new 2 lines below the lien belwo)
    IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
    //FArrayBox fab_drho_as_crse(bxg, cls_t::NCONS, The_Async_Arena());   
    //IArrayBox fab_rrflag_as_crse(bxg);   

    // in cerisse this call is different, depends on
    const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ? fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;
    FArrayBox* p_drho_as_crse = (fr_as_crse) ? fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
        
    auto const& fcx = flxt[0]->array(); 
    auto const& fcy = flxt[1]->array(); 
#if (AMREX_SPACEDIM==3)     
    auto const& fcz = flxt[2]->array(); 
#endif
    
    // redistribution temp snm
    if (redistribution_type == "NewRedist")
    {
      cerisse_flux_redistribute( ebbox,rhs, divc, redistwgt, vfrac,flag,geom,cls_t::NCONS,dt);
   }
   else  
   {
     amrex::ApplyMLRedistribution(
      ebbox, cls_t::NCONS, rhs, divc, cons, scratch, flag, AMREX_D_DECL(apx, apy, apz), vfrac,
      AMREX_D_DECL(fcx, fcy, fcz), bcent, phys_bc, geom, dt, redistribution_type,
      as_crse, p_drho_as_crse->array(), p_rrflag_as_crse->const_array(), as_fine, dm_as_fine, lev_mask,
      level_mask_not_covered, fac_for_deltaR, use_wts_in_divnc, 0, srd_max_order,
      target_volfrac, srd_update_scale);
   }                                

    // temp snm
    //amrex::Gpu::streamSynchronize();
    //AMREX_GPU_ERROR_CHECK();


  }  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
};


#endif

