#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <CNS_hydro_K.h>
#include <prob.h>
// #include <CentralKEEP.h>
// #include <Riemann.h>
// #include <High_resolution.h>
// #include <NLDE.h>
#ifdef AMREX_USE_GPIBM
#include <IBM.h>
#endif
using namespace amrex;

Real CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/) {
  BL_PROFILE("CNS::advance()");

  // Print() << "-------- advance start -------" << std::endl;
  // Print() << "time = " << time << std::endl;
  // Print() << "dt = " << dt << std::endl;
  // state[0].printTimeInterval(std::cout);
  // Print() << "---------------------------------" << std::endl;

  for (int i = 0; i < num_state_data_types; ++i) {
    state[i].allocOldData();
    state[i].swapTimeLevels(dt);
  }

  MultiFab& S1 = get_old_data(State_Type);
  MultiFab& S2 = get_new_data(State_Type);

  int ncons = d_prob_closures->NCONS;
  int nghost= d_prob_closures->NGHOST;
  // MultiFab dSdt(grids,dmap,ncons,0,MFInfo(),Factory());
  MultiFab Stemp(grids,dmap,ncons,nghost,MFInfo(),Factory());

  FluxRegister* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    CNS& fine_level = getLevel(level + 1);
    fr_as_crse = fine_level.flux_reg.get();
  }

  FluxRegister* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    fr_as_fine = flux_reg.get();
  }

  if (fr_as_crse) {
    fr_as_crse->setVal(Real(0.0));
  }

  if (order_rk == -2) {
    // Original time integration ///////////////////////////////////////////////
    // RK2 stage 1
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);

    // RK2 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n
    state[0].setNewTimeLevel(time + dt);
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine);
    // S_new = 0.5*(Sborder+S_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(S2, Real(0.5), S1, 0, Real(0.5), S1, 0, 0, ncons, 0);
    // S_new += 0.5*dt*dSdt
    MultiFab::Saxpy(S2, Real(0.5) * dt, Stemp, 0, 0, ncons, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*
    ////////////////////////////////////////////////////////////////////////////
  } else if (order_rk == 0) {  // returns rhs
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
    MultiFab::Copy(S2, Stemp, 0, 0, ncons, 0);
  } else if (order_rk == 1) {
    FillPatch(*this, Stemp, nghost, time, State_Type, 0,
              ncons);  // filled at t_n to evalulate f(t_n,y_n).
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);
  } else if (order_rk == 2) {
    // Low storage SSPRKm2 with m stages (C = m-1, Ceff=1-1/m). Where C is the
    // SSPRK coefficient, it also represents the max CFL over the whole
    // integration step (including m stages). From pg 84 Strong Stability
    // Preserving Runge–kutta And Multistep Time Discretizations
    int m = stages_rk;
    // Copy to S2 from S1
    MultiFab::Copy(S2, S1, 0, 0, ncons, 0);
    state[0].setOldTimeLevel(time);
    state[0].setNewTimeLevel(time);
    // first to m-1 stages
    // Print() << "-------- before RK stages -------" << std::endl;
    // Print() << "time = " << time << std::endl;
    // Print() << "dt = " << dt << std::endl;
    // state[0].printTimeInterval(std::cout);
    // Print() << "---------------------------------" << std::endl;
    for (int i = 1; i <= m - 1; i++) {
      FillPatch(*this, Stemp, nghost, time + dt * Real(i - 1) / (m - 1),
                State_Type, 0, ncons);
      compute_rhs(Stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine);
      MultiFab::Saxpy(S2, dt / Real(m - 1), Stemp, 0, 0, ncons, 0);
      state[State_Type].setNewTimeLevel(
          time + dt * Real(i) /
                     (m - 1));  // important to do this for correct fillpatch
                                // interpolations for the proceeding stages
    }
    // final stage
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
    compute_rhs(Stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine);
    MultiFab::LinComb(S2, Real(m - 1), S2, 0, dt, Stemp, 0, 0, ncons, 0);
    MultiFab::LinComb(S2, Real(1.0) / m, S1, 0, Real(1.0) / m, S2, 0, 0, ncons,
                      0);

    state[State_Type].setNewTimeLevel(
        time + dt);  // important to do this for correct fillpatch
                     // interpolations for the proceeding stages

    // Print() << "--------- after RK stages --------" << std::endl;
    // Print() << "time = " << time << std::endl;
    // Print() << "dt = " << dt << std::endl;
    // state[0].printTimeInterval(std::cout);
    // Print() << "----------------------------------" << std::endl;
  }

  else if (order_rk == 3) {
    if (stages_rk == 3) {
      state[0].setOldTimeLevel(time);
      // http://ketch.github.io/numipedia/methods/SSPRK33.html
      // state[0].setOldTimeLevel (time);
      FillPatch(*this, Stemp, nghost, time, State_Type, 0,
                ncons);  // filled at t_n to evalulate f(t_n,y_n).
      compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
      MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt);  // same time as upcoming FillPatch ensures we copy S2 to
                       // Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 4, fr_as_crse, fr_as_fine);
      MultiFab::Xpay(Stemp, dt, S2, 0, 0, ncons, 0);
      MultiFab::LinComb(S2, Real(3.0) / 4, S1, 0, Real(1.0) / 4, Stemp, 0, 0,
                        ncons, 0);

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      compute_rhs(Stemp, dt * Real(2.0) / 3, fr_as_crse, fr_as_fine);
      MultiFab::Xpay(Stemp, dt, S2, 0, 0, ncons, 0);
      MultiFab::LinComb(S2, Real(1.0) / 3, S1, 0, Real(2.0) / 3, Stemp, 0, 0,
                        ncons, 0);

      state[State_Type].setNewTimeLevel(
          time + dt);  // important to do this for correct fillpatch
                       // interpolations for the proceeding stages
    }

    else if (stages_rk == 4) {
      // http://ketch.github.io/numipedia/methods/SSPRK43.html and From pg 85
      // Strong Stability Preserving Runge–kutta And Multistep Time
      // Discretizations

      state[0].setOldTimeLevel(time);
      FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine);
      MultiFab::LinComb(S2, Real(1.0), S1, 0, dt / 2, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine);
      MultiFab::Saxpy(S2, dt / 2, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt);  // same time as upcoming FillPatch ensures we copy S2 to
                       // Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 6, fr_as_crse, fr_as_fine);
      MultiFab::LinComb(S2, Real(2.0) / 3, S1, 0, Real(1.0) / 3, S2, 0, 0,
                        ncons, 0);
      MultiFab::Saxpy(S2, dt / 6, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine);
      MultiFab::Saxpy(S2, dt / 2, Stemp, 0, 0, ncons, 0);

      state[State_Type].setNewTimeLevel(
          time + dt);  // important to do this for correct fillpatch
                       // interpolations for the proceeding stages
    }

    else {
      // Low storage SSPRKm3 with m=n^2, n>=3 stages (C=2, Ceff=0.5). From pg 85
      // Strong Stability Preserving Runge–kutta And Multistep Time
      // Discretizations
      // TODO Generally SSPRK(n^2,3) where n>2 - Ceff=1-1/n
      Print() << "SSPRK(m^2)3 not implemented yet" << std::endl;
      exit(0);
    }

  }

  // else if (order_rk == 4) {
  //   Print() << "SSPRK4 not implemented yet" << std::endl;
  //   exit(0);
  //   // TODO: SSPRK(10,4) C=6, Ceff=0.6
  // }
  return dt;
}

// Since we do not want to use expensive cudaMemCopy, we are storing all our
// data on the GPU to begin with. Concurrency on GPU using streams, parallel
// computation and data transfer, is not useful then. Therefore, we can have all
// grid point computations, per fab, in a single MFIter loop (single stream).

void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxRegister* fr_as_crse, FluxRegister* fr_as_fine) {
  BL_PROFILE("CNS::compute_rhs()");

  // Variables
  // TODO: introduce a struct for these variables?
  const PROB::ProbClosures& cls_d = *CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  const PROB::ProbParm& parms = *d_prob_parm;

  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& state = statemf.array(mfi);

    const Box& bxgnodal = mfi.grownnodaltilebox(-1, 0);  // extent is 0,N_cell+1
    const Box& bxg = mfi.growntilebox(cls_d.NGHOST);

    FArrayBox primf(bxg, cls_d.NPRIM , The_Async_Arena());
    FArrayBox tempf(bxgnodal, cls_d.NCONS, The_Async_Arena());
    Array4<Real> const& temp = tempf.array();
    Array4<Real> const& prims= primf.array();

    // We want to minimise function calls. So, we call prims2cons, flux and
    // source term evaluations once per fab from CPU, to be run on GPU.
    cls_h.cons2prims(mfi, state, prims);

    // Fluxes including boundary/discontinuity corrections
    // Note: we are over-writing state (cons) with flux derivative
    prob_rhs.eflux(geom, mfi, prims, temp, state, cls_h);
    prob_rhs.dflux(); //(prims,cons,nflx)

    // Source terms, including update mask (e.g inside IB)
    // prob_rhs.src(prims,cons,nflx,rhs)
  }
}

// Euler Fluxes ///////////////////////////////////////////////////////////////
//  if constexpr (PROB::do_euler==-1) {
//  MultiFab& basemf     = NLDE::Vbaseflow[level];
//  auto const& basefab  = basemf.array(mfi);
//  auto const& lambda   = lambdamf.array(mfi);
//  AMREX_D_TERM(auto const& nfabfx = numflxmf[0].array(mfi);,
//               auto const& nfabfy = numflxmf[1].array(mfi);,
//               auto const& nfabfz = numflxmf[2].array(mfi););

// AMREX_D_TERM(auto const& pfabfx = pntflxmf[0].array(mfi);,
//              auto const& pfabfy = pntflxmf[1].array(mfi);,
//              auto const& pfabfz = pntflxmf[2].array(mfi););

// NLDE::eflux(bxg,consfab,primfab,numflxmf);
//   Print() << "NLDE Euler fluxes not implemented yet" << std::endl;
//   exit(0);
// }
// else if (PROB::do_euler==1) {
//   CentralKEEP::FluxKEEP(consmf,primsmf,numflxmf);
// }
// else if (PROB::do_euler==2) {
//   HiRes::FluxWENO(consmf,primsmf,numflxmf);
// }
// else if (PROB::do_euler==3) {
//   Riemann::Flux(consmf,primsmf,numflxmf);
// }
// exit(0);
// }

//////////////////////////////////////////////////////////////////////////////

//   // Euler flux corrections (overwrite numflxmf) //
//   // Recompute fluxes on planes adjacent to physical boundaries (Order
//   reduction) if (flux_euler==1 && !(CentralKEEP::order_keep==2)) {
//     CentralKEEP::Flux_2nd_Order_KEEP(geom,primsmf,numflxmf);
//   }
//   // Order reduction near IBM

//   // Artificial dissipation (adding to numflxmf)
//   // JST artificial dissipation shock capturing
//   if (art_diss==1) {
//   // make multifab for spectral radius and sensor for artificial dissipation
//     MultiFab lambdamf; lambdamf.define(consmf.boxArray(),
//     consmf.DistributionMap(), AMREX_SPACEDIM, NGHOST); MultiFab sensormf;
//     sensormf.define(consmf.boxArray(), consmf.DistributionMap(),
//     AMREX_SPACEDIM, NGHOST); lambdamf = 0.0_rt; sensormf = 0.0_rt; for
//     (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//       const Box& bx      = mfi.tilebox();
//       const Box& bxnodal = mfi.grownnodaltilebox(-1,0);

//       auto const& statefab = consmf.array(mfi);
//       auto const& sensor   = sensormf.array(mfi);
//       auto const& lambda   = lambdamf.array(mfi);
//       auto const& prims    = primsmf.array(mfi);
//       AMREX_D_TERM(auto const& nfabfx = numflxmf[0].array(mfi);,
//                     auto const& nfabfy = numflxmf[1].array(mfi);,
//                     auto const& nfabfz = numflxmf[2].array(mfi););
//       amrex::ParallelFor(bx,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//           ComputeSensorLambda(i,j,k,primsfab,lambda,sensor,cls);
//         });

//         amrex::ParallelFor(bxnodal, ncons,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//         {
//           JSTflux(i,j,k,n,lambda,sensor,statefab,nfabfx,nfabfy,nfabfz,cls);
//         });
//     }
//   }
// }

//////////////////////////////////////////////////////////////////////////////

// Viscous Fluxes ////////////////////////////////////////////////////////////
// Gpu::streamSynchronize(); // ensure all rhs terms computed before assembly
// We have a separate MFIter loop here than the Euler fluxes and the source
// terms, so the work can be further parallised. As different MFIter loops can
// be in different GPU streams.

// Although conservative FD (finite difference) derivatives of viscous fluxes
// are not requried in the boundary layer, standard FD are likely sufficient.
// However, considering grid and flow discontinuities (coarse-interface
// flux-refluxing and viscous derivatives near shocks), conservative FD
// derivatives are preferred.
//   if (rhs_visc) {
//     Array<MultiFab,AMREX_SPACEDIM>& pntvflxmf = Vpntvflxmf[level];
// #if AMREX_USE_GPIBM
//     IBM::IBMultiFab& ibMultiFab = *IBM::ib.ibMFa[level];
// #endif
//     //for each fab
//     for (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//       const Box& bxpflx  = mfi.growntilebox(1);
//       auto const& prims  = primsmf.array(mfi);

//       AMREX_D_TERM(auto const& pfabfx = pntvflxmf[0].array(mfi);,
//                   auto const& pfabfy  = pntvflxmf[1].array(mfi);,
//                   auto const& pfabfz  = pntvflxmf[2].array(mfi););

//       AMREX_D_TERM(auto const& nfabfx = numflxmf[0].array(mfi);,
//                   auto const& nfabfy  = numflxmf[1].array(mfi);,
//                   auto const& nfabfz  = numflxmf[2].array(mfi););

//       // compute u,v,w,T derivatives and compute physical viscous fluxes
//       amrex::ParallelFor(bxpflx,[=] AMREX_GPU_DEVICE (int i, int j, int k)
//       noexcept {
//         ViscousFluxes(i, j, k, primsfab, pfabfx, pfabfy, pfabfz, dxinv, cls);
//       });

//       // Physical boundary viscous flux corrections
//       // (overwrite pfabfx, pfabfy, pfabfz)
//       // TODO:: generalise to wall boundary in x and z directions
//       const Box& bx  = mfi.tilebox();
//       //xlo
//       //xhi

//       //ylo
//       if(geom.Domain().smallEnd(1)==bx.smallEnd(1)) {
//         if ((*h_phys_bc).lo(1)==6) {
//           int jj = bx.smallEnd(1)-1;
//           IntVect small = {bxpflx.smallEnd(0), jj, bxpflx.smallEnd(2)};
//           IntVect big   = {bxpflx.bigEnd(0)  , jj, bxpflx.bigEnd(2)  };
//           Box bxboundary(small,big);

//           amrex::ParallelFor(bxboundary,[=] AMREX_GPU_DEVICE (int i, int j,
//           int k) noexcept {
//               ViscousWallFluxes(i, j, k, 0, primsfab, pfabfx, pfabfy, pfabfz,
//               dxinv, cls);
//           });
//         }
//       }
//       //yhi
//       if(geom.Domain().bigEnd(1)==bx.bigEnd(1)) {
//         if ((*h_phys_bc).hi(1)==6) {
//           int jj = bx.bigEnd(1) + 1;
//           IntVect small = {bxpflx.smallEnd(0), jj, bxpflx.smallEnd(2)};
//           IntVect big   = {bxpflx.bigEnd(0)  , jj, bxpflx.bigEnd(2)  };
//           Box bxboundary(small,big);

//           amrex::ParallelFor(bxboundary,[=] AMREX_GPU_DEVICE (int i, int j,
//           int k) noexcept {
//               ViscousWallFluxes(i, j, k, 1, primsfab, pfabfx, pfabfy, pfabfz,
//               dxinv, cls);
//           });
//         }
//       }
//       //zlo
//       //zhi

// #if AMREX_USE_GPIBM
//       //IBM GP viscous flux correction
//       auto const& ibFab = ibMultiFab.get(mfi);
//       auto const& markers = ibMultiFab.array(mfi);
//       auto const gp_ijk   = ibFab.gpData.gp_ijk.data();
//       amrex::ParallelFor(ibFab.gpData.ngps, [=] AMREX_GPU_DEVICE (int ii)
//       {
//         ViscousFluxGP(gp_ijk[ii](0),gp_ijk[ii](1),gp_ijk[ii](2),markers,primsfab,pfabfx,pfabfy,pfabfz,dxinv,cls);
//       });
// #endif

//       // compute numerical viscous fluxes (add to numflxmf)
//       const Box& bxnodal  = mfi.grownnodaltilebox(-1,0); // extent is 0,N+1
//       in all directions -- -1 means for all directions.
//       amrex::surroundingNodes(bx) does the same amrex::ParallelFor(bxnodal,
//       ncons,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
//         ViscousNumericalFluxes(i, j, k, n, pfabfx, pfabfy, pfabfz, nfabfx,
//         nfabfy, nfabfz);
//       });
//     }
//   }

// Re-fluxing ////////////////////////////////////////////////////////////////
// if constexpr (do_reflux) {
// if (fr_as_crse) {
//     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
//         const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ?
//         dx[0]*dx[2] : dx[0]*dx[1]); const Real scale = -dt*dA;
//         fr_as_crse->CrseInit(numflxmf[idim], idim, 0, 0, ncons, scale,
//         FluxRegister::ADD);
//     }
// }
// if (fr_as_fine) {
//     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
//         const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ?
//         dx[0]*dx[2] : dx[0]*dx[1]); const Real scale = dt*dA;
//         fr_as_fine->FineAdd(numflxmf[idim], idim, 0, 0, ncons, scale);
//     }
// }
// }
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// Add source term to RHS ////////////////////////////////////////////////////
// if (rhs_source) {
//   for (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi){
//     const Box& bx   = mfi.tilebox();
//     auto const& dsdtfab = dSdt.array(mfi);
//     auto const& statefab = consmf.array(mfi);

//     amrex::ParallelFor(bx,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//     { user_source(i,j,k,statefab,dsdtfab,lprobparm,cls,dx); });
//   }
// }
// //////////////////////////////////////////////////////////////////////////////

// Set solid point RHS to 0 //////////////////////////////////////////////////
// #if AMREX_USE_GPIBM
//   IBM::IBMultiFab& mfab = *IBM::ib.ibMFa[level];
//   for (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi){
//     const Box& bx   = mfi.tilebox();
//     auto const& dsdtfab = dSdt.array(mfi);
//     IBM::IBFab &fab = mfab.get(mfi);
//     Array4<bool> ibMarkers = fab.array();
//     amrex::ParallelFor(bx, ncons,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//     dsdtfab(i,j,k,n) = dsdtfab(i,j,k,n)*(1 - int(ibMarkers(i,j,k,0)));
//     });
//   }
// #endif
// }
//////////////////////////////////////////////////////////////////////////////
