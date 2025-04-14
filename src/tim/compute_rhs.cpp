#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <eib.h>
#endif

#include <mandebug.h>

using namespace amrex;




// Since we do not want to use expensive cudaMemCopy, we are storing all our
// data on the GPU to begin with. Concurrency on GPU using streams, parallel
// computation and data transfer, is not useful then. Therefore, we can have all
// grid point computations, per fab, in a single MFIter loop (single stream).

void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxRegister* fr_as_crse, FluxRegister* fr_as_fine) {
  BL_PROFILE("CNS::compute_rhs()");

  // Variables
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;

  //...................................................................
  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& state = statemf.array(mfi);

    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls_h.NGHOST);

    // primitives and fluxes arrays
    FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
    Array4<Real> const& prims= primf.array();

    
#ifdef CNS_USE_EB     
    // auxiliar arrays for redistribution 
    FArrayBox divcfab(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& divc= divcfab.array();    

    // store array cons 
    FArrayBox consfab(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& cons= consfab.array();    
    amrex::ParallelFor(bxg, cls_h.NCONS,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        cons(i,j,k,n) = state(i,j,k,n);
      });
#else
    // FArrayBox consfab(bxg, cls_h.NCONS, The_Async_Arena());
    // Array4<Real> const& cons= consfab.array();    
    // amrex::ParallelFor(bxg, cls_h.NCONS,
    //   [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    //   {
    //     cons(i,j,k,n) = state(i,j,k,n);
    //   });          
#endif
    // fliux arrays 
    std::array<FArrayBox ,AMREX_SPACEDIM> fluxt;
    for (int dir=0; dir < AMREX_SPACEDIM; ++dir)
    {
      fluxt[dir].resize(amrex::surroundingNodes(bxg, dir),cls_h.NCONS, The_Async_Arena() );
      fluxt[dir].setVal<RunOn::Device>(0.);
    }
     
    // We want to minimise function calls. So, we call prims2cons, flux and
    // source term evaluations once per fab from CPU, to be run on GPU.
    cls_h.cons2prims(mfi, state, prims);

    // combine arrays if IBM & EBM are used together 
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
    //create auxiliary aray
    BaseFab<bool> fab(bxg,2);
    Array4<bool> const& geoMarkers = fab.array();    
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
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
        auto const& flx = fluxt[dir].array();  
        ParallelFor(bx, cls_h.NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    state(i, j, k, n) =
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });
    }                  

                        
#if CNS_USE_EB    
    // internal geometry fluxes
    const Box&  ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const auto& flag = (*EBM::eb.ebflags_a[level])[mfi];
    FabType t = flag.getType(ebbox);

    const bool fab_with_eb = (FabType::singlevalued == t);

    if (fab_with_eb) {
      EBM::eb.ebflux(geom,mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},state, cls_d,level);
    }
  
    // redistribution (WARNING at this point state has the rhs, prims has prims)
    
    // compute divc here
    amrex::ParallelFor(bxg, cls_h.NCONS,  
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      divc(i,j,k,n) = state(i,j,k,n);
    });  

    // do redistribution only in box with EB
    if (eb_redistribution && fab_with_eb){    
      EBM::eb.redist(geom,mfi,cons,divc, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},
                    state, cls_d,level,dt,h_phys_bc);
    }                    
      
#endif 

    // Source terms, including update mask (e.g inside IB)
    prob_rhs.src(mfi, prims, state, cls_d, dt);

    // Set solid point RHS to 0  (state hold RHS at this point)
#if AMREX_USE_GPIBM || CNS_USE_EB
    amrex::ParallelFor(bx, cls_h.NCONS,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      state(i,j,k,n) = state(i,j,k,n)*(1 - int(geoMarkers(i,j,k,0)));
    });
#endif
 

    // TODO: IBM::set_solid_state(mfi,state,cls_d)


    // // Flux register
    // if (do_reflux) {
    //   const auto dx = geom.CellSizeArray();
    //   if (fr_as_fine) {
    //     fr_as_fine->FineAdd(mfi,
    //                         {AMREX_D_DECL(&fluxes[0], &fluxes[1], &fluxes[2])},
    //                         dx.data(), dtsub, RunOn::Device);
    //   }
    //   if (fr_as_crse) {
    //     fr_as_crse->CrseAdd(mfi,
    //                         {AMREX_D_DECL(&fluxes[0], &fluxes[1], &fluxes[2])},
    //                         dx.data(), dtsub, RunOn::Device);
    //   }
    // }
  }
 
}
