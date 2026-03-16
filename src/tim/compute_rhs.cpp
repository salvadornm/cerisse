#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <eib.h>
#endif


using namespace amrex;

// Since we do not want to use expensive cudaMemCopy, we are storing all our
// data on the GPU to begin with. Concurrency on GPU using streams, parallel
// computation and data transfer, is not useful then. Therefore, we can have all
// grid point computations, per fab, in a single MFIter loop (single stream).

void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxReg* fr_as_crse, FluxReg* fr_as_fine) {
  BL_PROFILE("CNS::compute_rhs()");


  // Variables
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;     // for device access to closures (GPU)
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;    // for host  access to closures  (CPU)

  // local Indexing for convenience  
  constexpr int NCONS  = PROB::ProbClosures::NCONS;
  constexpr int NPRIM  = PROB::ProbClosures::NPRIM;
  constexpr int NGHOST = PROB::ProbClosures::NGHOST;

  // time
  const Real cur_time = state[State_Type].curTime();
  
 // PROB::ProbRHS prob_rhs;  //  local RHS object, lives only in this function
 // prob_rhs.init_coeffs();  // initialize diffusion coefficients if needed

  //...................................................................
  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& state = statemf.array(mfi);

    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(NGHOST);
#ifdef CNS_USE_EB     
    const Box& bxflux = mfi.growntilebox(NGHOST+1); // add 1 cell 
#else
    const Box& bxflux = mfi.growntilebox(NGHOST); 
#endif    
    
    // primitives and fluxes arrays
    FArrayBox primf(bxg, NPRIM, The_Async_Arena());
    Array4<Real> const& prims= primf.array();
    
#ifdef CNS_USE_EB     
    // auxiliar arrays for redistribution 
    FArrayBox divcfab(bxg, NCONS, The_Async_Arena());
    Array4<Real> const& divc= divcfab.array();    

    // store array cons 
    FArrayBox consfab(bxg, NCONS, The_Async_Arena());
    Array4<Real> const& cons= consfab.array();    
    amrex::ParallelFor(bxg, NCONS,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        cons(i,j,k,n) = state(i,j,k,n);
      });
#endif
    // flux arrays  
    std::array<FArrayBox ,AMREX_SPACEDIM> fluxt;
    for (int dir=0; dir < AMREX_SPACEDIM; ++dir)
    {
      fluxt[dir].resize(amrex::surroundingNodes(bxflux, dir),NCONS, The_Async_Arena() );
      fluxt[dir].setVal<RunOn::Device>(0.);
    }
     
    // We want to minimise function calls. So, we call prims2cons, flux and
    // source term evaluations once per fab from CPU, to be run on GPU.
    cls_h.cons2prims(mfi, state, prims); // <<<<<<<<

    // combine arrays if IBM & EBM are used together 
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
    //create auxiliary aray
    BaseFab<uint8_t> fab(bxg,2);
    Array4<uint8_t> const& geoMarkers = fab.array();    
#endif
    // extract markers
#ifdef AMREX_USE_GPIBM
    auto& ib_mf = *IBM::ib.bmf_a[level];
    IBM::ib.computeGPs(mfi, state, prims, cls_d, level);
    const auto& ibMarkers = ib_mf.array(mfi);
#endif
#ifdef CNS_USE_EB   
    // no need to update markers here
    auto& eb_mf = *EBM::eb.bmf_a[level];
    const auto& ebMarkers = eb_mf.array(mfi);
#endif

    // combine markers into one  (CAN BE DONE BETTER)
#if (AMREX_USE_GPIBM && CNS_USE_EB)    
    amrex::ParallelFor(bxg, 2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
    geoMarkers(i,j,k,n) = ebMarkers(i,j,k,n) && ibMarkers(i,j,k,n);
    });
#endif    
#if (AMREX_USE_GPIBM && !CNS_USE_EB)
    amrex::ParallelFor(bxg, 2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
    geoMarkers(i,j,k,n) = ibMarkers(i,j,k,n);
    });       
#endif    
#if (CNS_USE_EB && !AMREX_USE_GPIBM)
    amrex::ParallelFor(bxg, 2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
    geoMarkers(i,j,k,n) = ebMarkers(i,j,k,n);
    });
#endif
  
    // Euler/Diff Fluxes including boundary/discontinuity corrections
    // WARNING: state is the U array (cons)
#if (AMREX_USE_GPIBM || CNS_USE_EB )      
    prob_rhs.eflux_ibm(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d, geoMarkers);    
    prob_rhs.dflux_ibm(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d, geoMarkers);
#else
    prob_rhs.eflux(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d);
    prob_rhs.dflux(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d);
#endif

    // compute rhs as flux derivative, i.e.  rhs + = (flx[i] - flx[i+1])/dx
    // WARNING: state is now the RHS array
    // set RHS=0 (everywhere including ghost points)
    ParallelFor(bxg, NCONS,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {state(i,j,k,n) = 0.0;});
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
        auto const& flx = fluxt[dir].array();  
        ParallelFor(bx, NCONS,  // bxg or bx (should be bxg but needs correct flux)
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    state(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });
    }                  
                      
#if CNS_USE_EB    
    // internal geometry fluxes
    const Box&  ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const auto& flag = (*EBM::eb.ebflags_a[level])[mfi];
    FabType t = flag.getType(ebbox);

    const bool fab_with_eb     = (FabType::singlevalued == t);  
    const bool fab_with_fluid  =  !(t == amrex::FabType::covered);

    // EB flux     
    if (fab_with_eb) {
      EBM::eb.ebflux(geom,mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},state, cls_d,level);
    }

    // redistribution 
    // WARNING: state is  the RHS array, prims is the prims 
    // compute divc here
    amrex::ParallelFor(bxg, NCONS,  
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      divc(i,j,k,n) = state(i,j,k,n);
    });  
    
    // do redistribution only in box with EB
    FArrayBox dm_as_fine(Box::TheUnitBox(), NCONS, The_Async_Arena());
    if (fr_as_fine) {
        dm_as_fine.resize(amrex::grow(bx, 1), NCONS);
        dm_as_fine.setVal<RunOn::Device>(0.0);
    }
    if (eb_redistribution && fab_with_eb){
      EBM::eb.redist(geom, mfi, cons, divc, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},
                     fr_as_crse, fr_as_fine, dm_as_fine.array(), 
                     state, cls_d, level, dt, h_phys_bc);
    }                    
      

#endif 

    // Add source terms (only in fluid fabs)
#if CNS_USE_EB 
    if (fab_with_fluid ) {    
        prob_rhs.src(geom,mfi, prims, state, cls_d, dt, cur_time, geoMarkers);
    }
#else
    prob_rhs.src(geom,mfi, prims, state, cls_d, dt, cur_time);
#endif


    // Set solid point RHS to 0  (state hold RHS at this point)
#if AMREX_USE_GPIBM || CNS_USE_EB
    amrex::ParallelFor(bxg, NCONS,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      state(i,j,k,n) = state(i,j,k,n)*(1 - int(geoMarkers(i,j,k,0)));
    });
#endif



    // TODO: IBM::set_solid_state(mfi,state,cls_d)

    
    // Flux register accumulation (conservation across AMR levels)
    if (do_reflux && (fr_as_crse || fr_as_fine)) {
        const Real* dx = geom.CellSize();
        const int ncomp = PROB::ProbClosures::NCONS;
#ifdef CNS_USE_EB
        if (t == FabType::singlevalued) {
            const FArrayBox& volfrac = (*EBM::eb.volmf_a[level])[mfi];
            AMREX_D_TERM(const FArrayBox& areafracx = (*(EBM::eb.areamcf_a[level][0]))[mfi];,
                         const FArrayBox& areafracy = (*(EBM::eb.areamcf_a[level][1]))[mfi];,
                         const FArrayBox& areafracz = (*(EBM::eb.areamcf_a[level][2]))[mfi];)
            if (fr_as_crse) {
                fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, dx, dt,
                                    volfrac, {AMREX_D_DECL(&areafracx, &areafracy, &areafracz)},
                                    0, 0, ncomp, RunOn::Device);
            }
            if (fr_as_fine) {
                fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, dx, dt,
                                    volfrac, {AMREX_D_DECL(&areafracx, &areafracy, &areafracz)}, dm_as_fine,
                                    0, 0, ncomp, RunOn::Device);
            }
        } else 
#endif
        {
            if (fr_as_crse) {
                fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, dx, dt, 
                                    0, 0, ncomp, RunOn::Device);
            }
            if (fr_as_fine) {
                fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, dx, dt, 
                                    0, 0, ncomp, RunOn::Device);
            }
        }
    }


  } // end mfi loop
 
}

// 
// clip species state to ensure positivity and sum Y=1 (if NUM_SPECIES > 1)

#if NUM_SPECIES > 1
void CNS::clip_species_state(MultiFab& statemf) {
  BL_PROFILE("CNS:clip_speciesstate()");

  // local Indexing for convenience  
  constexpr int UFS    = PROB::ProbClosures::UFS;
  constexpr int NGHOST = PROB::ProbClosures::NGHOST;

  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& cons = statemf.array(mfi);

    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(NGHOST);
    
    // clip state
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
      // compute real rho 
      Real rhoreal = 0.0;
	    for (int n = 0; n < NUM_SPECIES; ++n) {
	      rhoreal += max( cons(i, j, k, UFS + n) ,0.0);
	    }
      // ensure sum Y[k] =-1
	    Real rhoinv = Real(1.0) / rhoreal;
	    Real Y[NUM_SPECIES]; Real sumY=0.0;
	    for (int n = 0; n < NUM_SPECIES; ++n) {
	      Y[n] = max( cons(i, j, k, UFS + n),0.0) * rhoinv; 
	      sumY += Y[n];
	    }
	    for (int n = 0; n < NUM_SPECIES; ++n) {Y[n] /= sumY;}
	     
      // readjust	      
	    for (int n = 0; n < NUM_SPECIES; ++n) {
	      cons(i, j, k, UFS + n) = rhoreal*Y[n]; 
	    }
    }); 

  } // end looop mfi  

}
#endif

