#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <prob.h>
#include <limits>  // for std::numeric_limits (RZ divergence floor)

// IBM and EB marker paths use different array types; combined mode is not yet supported.
#if defined(AMREX_USE_GPIBM) && defined(CNS_USE_EB)
#error "compute_rhs.cpp expects exactly one of AMREX_USE_GPIBM or CNS_USE_EB to be enabled"
#endif

#ifdef AMREX_USE_GPIBM
#include <eib.h>
#endif


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

  // time
  const Real cur_time = state[State_Type].curTime();
  
  PROB::ProbRHS prob_rhs;  //  local RHS object, lives only in this function
  prob_rhs.init_coeffs();  // initialize diffusion coefficients if needed

#ifdef AMREX_USE_GPIBM
  // Pre-compute level-wide prims and apply IB ghost-cell values in a single kernel
  MultiFab prims_mf(statemf.boxArray(), statemf.DistributionMap(),
                    cls_h.NPRIM, cls_h.NGHOST,
                    MFInfo().SetArena(The_Async_Arena()));
  {
    BL_PROFILE_VAR("CNS::compute_rhs::cons2prims", prof_cons2prims);
    for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
      cls_h.cons2prims(mfi, statemf.array(mfi), prims_mf.array(mfi));
    }
    BL_PROFILE_VAR_STOP(prof_cons2prims);
  }
  {
    BL_PROFILE_VAR("IBM::computeAllGPs", prof_gp);
    IBM::ib.computeAllGPs(prims_mf, cls_d, level);
    BL_PROFILE_VAR_STOP(prof_gp);
  }
#endif

  //...................................................................
  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& state = statemf.array(mfi);

    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls_h.NGHOST);
#ifdef CNS_USE_EB     
    const Box& bxflux = mfi.growntilebox(cls_h.NGHOST+1); // add 1 cell 
#else
    const Box& bxflux = mfi.growntilebox(cls_h.NGHOST); 
#endif    
    
    // primitives and fluxes arrays
#ifdef AMREX_USE_GPIBM
    // Prims already filled + GP-corrected in pre-loop; just alias
    Array4<Real> const& prims = prims_mf.array(mfi);
#else
    FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
    Array4<Real> const& prims = primf.array();
#endif

    
#ifdef CNS_USE_EB     
    // auxiliary arrays for redistribution 
    FArrayBox divcfab(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& divc = divcfab.array();    

    // store array cons 
    FArrayBox consfab(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& cons = consfab.array();    
    amrex::ParallelFor(bxg, cls_h.NCONS,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        cons(i,j,k,n) = state(i,j,k,n);
      });
#endif

    // flux arrays  
    std::array<FArrayBox ,AMREX_SPACEDIM> fluxt;
    for (int dir=0; dir < AMREX_SPACEDIM; ++dir)
    {
      fluxt[dir].resize(amrex::surroundingNodes(bxflux, dir),cls_h.NCONS, The_Async_Arena() );
      fluxt[dir].setVal<RunOn::Device>(0.);
    }
     
    // We want to minimise function calls. So, we call prims2cons, flux and
    // source term evaluations once per fab from CPU, to be run on GPU.
#ifndef AMREX_USE_GPIBM
    cls_h.cons2prims(mfi, state, prims);
#endif

    // Geometry markers (one of IBM/EB is expected to be enabled).
    // Alias geoMarkers to the underlying marker MultiFab to avoid
    // allocating an auxiliary fab + doing a device copy.
#if (AMREX_USE_GPIBM && !CNS_USE_EB)
    auto& marker_mf = *IBM::ib.bmf_a[level];
    auto const& geoMarkers = marker_mf.array(mfi);
#elif (CNS_USE_EB && !AMREX_USE_GPIBM)
    // EB markers are maintained elsewhere; no need to update them here.
    auto& marker_mf = *EBM::eb.bmf_a[level];
    auto const& geoMarkers = marker_mf.array(mfi);
#endif
  
    // Euler/Diff Fluxes including boundary/discontinuity corrections
    // WARNING: state is the U array (cons)
    {
    BL_PROFILE_VAR("CNS::compute_rhs::eflux", prof_eflux);
#if (AMREX_USE_GPIBM || CNS_USE_EB)
    prob_rhs.eflux_ibm(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d, geoMarkers);    
#else
    prob_rhs.eflux(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d);
#endif
    BL_PROFILE_VAR_STOP(prof_eflux);
    }
    {
    BL_PROFILE_VAR("CNS::compute_rhs::dflux", prof_dflux);
#if (AMREX_USE_GPIBM || CNS_USE_EB)
    prob_rhs.dflux_ibm(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d, geoMarkers);
#else
    prob_rhs.dflux(geom, mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])}, state, cls_d);
#endif
    BL_PROFILE_VAR_STOP(prof_dflux);
    }

    // compute rhs as finite-volume flux divergence, i.e.
    //   rhs += ((F·A)_lo - (F·A)_hi) / V
    // WARNING: state is now the RHS array
    // set RHS=0 (everywhere including ghost points)
    ParallelFor(bxg, cls_h.NCONS, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {state(i,j,k,n) = 0.0;});

    // Geometry-aware finite-volume divergence computed on-the-fly in the same
    // kernel (no temporary volume/area arrays).
    // - Cartesian: uses constant areas/volume derived from dx.
    // - RZ (2D): uses cylindrical metrics in (r,z) with r = prob_lo[0] + i*dr.
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    auto const& fx = fluxt[0].array();
#if (AMREX_SPACEDIM >= 2)
    auto const& fy = fluxt[1].array();
#endif
#if (AMREX_SPACEDIM == 3)
    auto const& fz = fluxt[2].array();
#endif

    const bool is_rz = geom.IsRZ();

    // Performance note: computing geometric metrics per conserved component
    // (4D ParallelFor) is wasteful when NCONS is large (e.g. many species).
    // Use a 3D kernel and loop over n to reuse V/A metrics.
    const int ncons = cls_h.NCONS;

    if (is_rz) {

        const Real dr = dx[0];
        const Real dz = dx[1];
        const Real r0 = prob_lo[0];
        const int qpres   = PROB::ProbClosures::QPRES;
        const int umom_r  = PROB::ProbClosures::UMX;
        const Real inv_dz = Real(1.0) / dz;
        
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                    // Axisymmetric RZ: dir=0 is r, dir=1 is z.
                    // In this branch:
                    //   fx = flux density through r-faces (i±1/2)
                    //   fy = flux density through z-faces (j±1/2)
                    const Real r_lo = r0 + Real(i) * dr;
                    const Real r_hi = r_lo + dr;
                    const Real r_c  = r_lo + Real(0.5) * dr;

                    // r2diff = r_hi^2 - r_lo^2 = (r_hi + r_lo)*(r_hi - r_lo)
                    //        = (r_hi + r_lo) * dr
                    const Real r2diff = (r_hi + r_lo) * dr;
                    // With r0=0 and i>=0, r2diff is strictly positive. Keep a tiny
                    // positive floor to avoid division-by-zero if future setups
                    // violate this assumption. Use a Real-typed floor to work for
                    // both float and double builds.
                    const Real tiny = amrex::max(std::numeric_limits<Real>::min(), Real(1.0e-14) * dr * dr);
                    const Real inv_r2diff = Real(1.0) / amrex::max(r2diff, tiny);
                    const Real tiny_r = Real(1.0e-14) * dr;
                    const Real inv_r = Real(1.0) / amrex::max(r_c, tiny_r);

                    for (int n = 0; n < ncons; ++n) {
                        // FV form with metrics. After cancellation of common factors:
                        // - z term reduces to Cartesian difference / dz
                        // - r term is (2*r*F_r)|_lo - (2*r*F_r)|_hi over (r_hi^2 - r_lo^2)
                        
                        // r-direction: metric-consistent RZ divergence
                        //   -(1/r) d(rF_r)/dr ≈ 2(r_lo*F_lo - r_hi*F_hi) / (r_hi^2 - r_lo^2)
                        Real rhs_rz = (Real(2.0) * (r_lo * fx(i, j, k, n) - r_hi * fx(i + 1, j, k, n))) * inv_r2diff;
                        state(i, j, k, n) += rhs_rz;

                        // Axisymmetric Euler geometric source for radial momentum:
                        // +p/r term is not contained in -(1/r) d(r F_r)/dr - dF_z/dz
                        // when F_r uses the standard Cartesian-form momentum flux.
                        if (n == umom_r) {
                            state(i, j, k, n) += prims(i, j, k, qpres) * inv_r;
                        }

                        state(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * inv_dz;
                    }
                });
    } else {
#if (AMREX_SPACEDIM == 1)
        const Real invdx = Real(1.0) / dx[0];
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < ncons; ++n) {
                        state(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
                    }
                });
#elif (AMREX_SPACEDIM == 2)
        const Real invdx = Real(1.0) / dx[0];
        const Real invdy = Real(1.0) / dx[1];
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < ncons; ++n) {
                        state(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
                        state(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
                    }
                });
#else
        const Real invdx = Real(1.0) / dx[0];
        const Real invdy = Real(1.0) / dx[1];
        const Real invdz = Real(1.0) / dx[2];
        ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < ncons; ++n) {
                        state(i, j, k, n) += (fx(i, j, k, n) - fx(i + 1, j, k, n)) * invdx;
                        state(i, j, k, n) += (fy(i, j, k, n) - fy(i, j + 1, k, n)) * invdy;
                        state(i, j, k, n) += (fz(i, j, k, n) - fz(i, j, k + 1, n)) * invdz;
                    }
                });
#endif
    }

    // RZ-only geometric viscous source terms (e.g. hoop stress) that are not captured by the
    // metric FV divergence of the r/z face fluxes.
    if (is_rz) {
    #if (AMREX_USE_GPIBM || CNS_USE_EB)
        // TODO: implement rz_geometric_source_ibm
        // prob_rhs.rz_geometric_source_ibm(geom, mfi, prims, state, cls_d, geoMarkers);
    #else
        prob_rhs.rz_geometric_source(geom, mfi, prims, state, cls_d);
    #endif
    }
                      
#if CNS_USE_EB    
    // internal geometry fluxes
    const Box&  ebbox  = mfi.growntilebox(0);  // box without ghost points 
    const auto& flag = (*EBM::eb.ebflags_a[level])[mfi];
    FabType t = flag.getType(ebbox);

    const bool fab_with_eb = (FabType::singlevalued == t);  
    // EB flux     
    if (fab_with_eb) {
      EBM::eb.ebflux(geom,mfi, prims, {AMREX_D_DECL(&fluxt[0], &fluxt[1], &fluxt[2])},state, cls_d,level);
    }

    // redistribution 
    // WARNING: state is  the RHS array, prims is the prims 
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

    // Source terms
    prob_rhs.src(geom,mfi, prims, state, cls_d, dt, cur_time);

    // Set solid point RHS to 0  (state hold RHS at this point)
#if (AMREX_USE_GPIBM || CNS_USE_EB)
        amrex::ParallelFor(bxg,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // IBM: geoMarkers(i,j,k,0) stores geometry_index in a uint8_t.
            // EB : geoMarkers(i,j,k,0) is a bool covered-cell marker.
            // In both cases, nonzero => solid.
            if (geoMarkers(i,j,k,0) != 0) {
                for (int n = 0; n < ncons; ++n) {
                    state(i,j,k,n) = Real(0.0);
                }
            }
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
