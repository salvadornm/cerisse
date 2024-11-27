#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <CNS_hydro_K.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <eib.h>
#endif
using namespace amrex;

Real CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/) {
  BL_PROFILE("CNS::advance()");

  // Print() << " oo CNS::advance " << std::endl;
  // Print() << "time = " << time << std::endl;
  // Print() << "dt = " << dt << std::endl;
  // Print() << "stages_rk = " << stages_rk << std::endl;
  // Print() << "order_rk = " << order_rk << std::endl;
  // Print() << " num_state_data_types= " << num_state_data_types << std::endl;
  
  // state[0].printTimeInterval(std::cout);

  // for (int i = 0; i < num_state_data_types; ++i) {
  //   state[i].allocOldData();
  //   state[i].swapTimeLevels(dt);
  // }

  state[0].allocOldData();
  state[0].swapTimeLevels(dt);
  

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
    MultiFab::LinComb(S2, Real(0.5), S1, 0, Real(0.5), S2, 0, 0, ncons, 0);
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

#ifdef AMREX_USE_GPIBM
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  auto& ib_mf = *IBM::ib.bmf_a[level];
  FillPatch(*this, Stemp, cls_h.NGHOST, time + dt, State_Type, 0, ncons);

  for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(cls_h.NGHOST);
    FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
    Array4<Real> const& state_temp = Stemp.array(mfi);
    Array4<Real> const& state = S2.array(mfi);
    const auto& ibMarkers = ib_mf.array(mfi);
    Array4<Real> const& prims = primf.array();

    cls_h.cons2prims(mfi, state_temp, prims);
    IBM::ib.computeGPs(mfi, state_temp, prims, cls_d, level);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (ibMarkers(i, j, k, 1)) {
        state(i, j, k, cls_d->UMX) = prims(i, j, k, cls_d->QRHO) * prims(i, j, k, cls_d->QU);
        state(i, j, k, cls_d->UMY) = prims(i, j, k, cls_d->QRHO) * prims(i, j, k, cls_d->QV);
        state(i, j, k, cls_d->UMZ) = prims(i, j, k, cls_d->QRHO) * prims(i, j, k, cls_d->QW);
        state(i, j, k, cls_d->UET) = prims(i, j, k, cls_d->QPRES) / (cls_d->gamma - 1.0) +
            0.5 * (prims(i, j, k, cls_d->QU) * prims(i, j, k, cls_d->QU) +
                   prims(i, j, k, cls_d->QV) * prims(i, j, k, cls_d->QV) +
                   prims(i, j, k, cls_d->QW) * prims(i, j, k, cls_d->QW));
        state(i, j, k, cls_d->URHO) = prims(i, j, k, cls_d->QRHO);
      }
    });
  }
#endif

  // else if (order_rk == 4) {
  //   Print() << "SSPRK4 not implemented yet" << std::endl;
  //   exit(0);
  //   // TODO: SSPRK(10,4) C=6, Ceff=0.6
  // }
  return dt;
}

// void where_is_nan(const FArrayBox& fab, bool abort_on_nan = true) {
//   bool contains_any_nan = false;
//   for (int n = 0; n < fab.nComp(); ++n) {
//     IntVect where;
//     bool contains_nan = fab.contains_nan<RunOn::Gpu>(fab.box(), n, 1, where);    
//     if (contains_nan) {
//       amrex::Print() << "NAN found at " << where << " comp " << n << '\n';
//       contains_any_nan = true;
//     }
//   }
//   if (contains_any_nan && abort_on_nan) {
//     amrex::Abort();
//   }
// }

// Since we do not want to use expensive cudaMemCopy, we are storing all our
// data on the GPU to begin with. Concurrency on GPU using streams, parallel
// computation and data transfer, is not useful then. Therefore, we can have all
// grid point computations, per fab, in a single MFIter loop (single stream).

void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxRegister* fr_as_crse, FluxRegister* fr_as_fine) {
  BL_PROFILE("CNS::compute_rhs()");

  // Variables
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  // const PROB::ProbParm& parms = *d_prob_parm;

  for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
    Array4<Real> const& state = statemf.array(mfi);

    // const Box& bxgnodal = mfi.grownnodaltilebox(-1, 1);  // extent is 0,N_cell+1
    const Box& bxg = mfi.growntilebox(cls_h.NGHOST);

    FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
    FArrayBox tempf(bxg, cls_h.NCONS, The_Async_Arena());
    Array4<Real> const& temp = tempf.array();
    Array4<Real> const& prims= primf.array();

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
    // Note: we are over-writing state (cons) with flux derivative
#if (AMREX_USE_GPIBM || CNS_USE_EB )      
    prob_rhs.eflux_ibm(geom, mfi, prims, temp, state, cls_d, geoMarkers);
    prob_rhs.dflux_ibm(geom, mfi, prims, temp, state, cls_d, geoMarkers);
#else
    prob_rhs.eflux(geom, mfi, prims, temp, state, cls_d);
    prob_rhs.dflux(geom, mfi, prims, temp, state, cls_d);
#endif

    // Source terms, including update mask (e.g inside IB)
    prob_rhs.src(mfi, prims, state, cls_d, dt);

    // Set solid point RHS to 0  (state hold RHS at this point)
#if AMREX_USE_GPIBM || CNS_USE_EB
    const Box& bx   = mfi.tilebox();
    amrex::ParallelFor(bx, cls_h.NCONS,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      state(i,j,k,n) = state(i,j,k,n)*(1 - int(geoMarkers(i,j,k,0)));
    });
#endif
    // TODO: IBM::set_solid_state(mfi,state,cls_d)

#if CNS_USE_EB    
    // call to compute fluxes in cells close to boundary  !!
    // eb.ebflux(geom,mfi,vf, prims, stats, cls_dt)

    // Redistribution ************
    
#endif






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