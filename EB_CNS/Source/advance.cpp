#if CNS_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "CNS.H"
#include "pdf_model.H"

using namespace amrex;

/** @brief The main driver for a single level implementing the time advance.
 *  @param time the current simulation time
 *  @param dt the timestep to advance (e.g., go from time to time + dt)
 *  @param iteration where we are in the current AMR subcycle. Each
                     level will take a number of steps to reach the
                     final time of the coarser level below it. This
                     counter starts at 1
 *  @param ncycle  the number of subcycles at this level
 * */
Real
CNS::advance(Real time, Real dt, int iteration, int ncycle)
{
  BL_PROFILE("CNS::advance()");

  // enforce_consistent_state(); // Enforce rho = sum(rhoY)

  // Prepare data fabs
  for (int i = 0; i < num_state_data_types; ++i) {
    if (!(i == Reactions_Type && do_react)) { // do not swap I_R (we only use new I_R data)
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }
  }

  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& S_old = get_old_data(State_Type);
  MultiFab Sborder(grids, dmap, LEN_STATE, NUM_GROW, MFInfo(), Factory());
  MultiFab dSdt_old(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
  MultiFab dSdt_new(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
  dSdt_old.setVal(0.0);
  dSdt_new.setVal(0.0);

  MultiFab& I_R = get_new_data(Reactions_Type);

  get_new_data(Cost_Type).setVal(0.0);

  // Prepare flux register
#if CNS_USE_EB
  EBFluxRegister* fr_as_crse = nullptr;
  EBFluxRegister* fr_as_fine = nullptr;
#else
  YAFluxRegister* fr_as_crse = nullptr;
  YAFluxRegister* fr_as_fine = nullptr;
#endif
  if (do_reflux && level < parent->finestLevel()) {
    fr_as_crse = &getLevel(level + 1).flux_reg;
  }
  if (do_reflux && level > 0) {
    fr_as_fine = &flux_reg;
  }

  if (fr_as_crse) {
    fr_as_crse->reset();
  }

  // Start time-stepping
  // RK2 stage 1: U^* = U^n + dt*dUdt^n + dt*I_R^n
  if (verbose > 1) amrex::Print() << " >> AMR Cycle " << iteration << " of " << ncycle << std::endl;
  if (verbose > 0) amrex::Print() << " >> RK Step 1: Computing dSdt^{n}" << std::endl;
  FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, LEN_STATE);
  compute_dSdt(Sborder, dSdt_old, 0.5*dt, fr_as_crse, fr_as_fine, true);
  MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt_old, 0, 0, LEN_STATE, 0);
  if (do_react) {
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      MultiFab::Saxpy(S_new, dt, I_R, nf*NREACT, nf*NVAR+UFS, NUM_SPECIES, 0);
    }
  }
  // if (NUM_FIELD > 0) {
  //   do something
  // }
  enforce_consistent_state(); // Enforce rho = sum(rhoY)

  // RK2 stage 2: U^{n+1} = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1}) + dt*I_R^{n+1}
  if (verbose > 0) amrex::Print() << " >> RK Step 2: Computing dSdt^{n+1}" << std::endl;
  FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, LEN_STATE);
  compute_dSdt(Sborder, dSdt_new, 0.5*dt, fr_as_crse, fr_as_fine, (rk_reaction_iter < 1));
  MultiFab::LinComb(S_new, 0.5*dt, dSdt_new, 0, 0.5*dt, dSdt_old, 0, 0, LEN_STATE, 0);
  MultiFab::Add(S_new, S_old, 0, 0, LEN_STATE, 0);

#if (NUM_FIELD > 0) 
  computeAvg(S_new);
  FillPatch(*this, Sborder, 1, time+dt, State_Type, 0, LEN_STATE);
  compute_pdf_model(Sborder, dt, iteration);
  MultiFab::Copy(S_new, Sborder, 0, 0, LEN_STATE, 0);
  // computeAvg(S_new);
  // enforce_consistent_state(); // Enforce rho = sum(rhoY)
#endif

  if (do_react) { // Compute I_R^{n+1}(U^**) and do U^{n+1} = U^** + dt*I_R^{n+1}
    react_state(time, dt);
  }
  enforce_consistent_state(); // Enforce rho = sum(rhoY)
  
  // Iterate to couple chemistry (not tested)
  // if (do_react && (rk_reaction_iter > 1)) {
  //   for (int iter = 1; iter < rk_reaction_iter; ++iter) {
  //     if (verbose > 0) {
  //       amrex::Print() << " >> Re-computing dSdt^{n+1}" << std::endl;
  //     }
  //     FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, LEN_STATE);
  //     compute_dSdt(Sborder, dSdt_new, 0.5*dt, fr_as_crse, fr_as_fine, (iter == rk_reaction_iter));
      
  //     MultiFab::LinComb(S_new, 0.5*dt, dSdt_new, 0, 0.5*dt, dSdt_old, 0, 0, LEN_STATE, 0);
  //     MultiFab::Add(S_new, S_old, 0, 0, LEN_STATE, 0);      
  //     react_state(time, dt);
  //     enforce_consistent_state(); // Enforce rho = sum(rhoY)
  //   }
  // }

  // if (NUM_FIELD > 0) {
  //   computeAvg(S_new);
  //   FillPatch(*this, Sborder, 1, time+dt, State_Type, 0, LEN_STATE);
  //   compute_pdf_model(Sborder, dt);
  //   MultiFab::Copy(S_new, Sborder, 0, 0, LEN_STATE, 0);
  //   computeAvg(S_new);
  //   enforce_consistent_state(); // Enforce rho = sum(rhoY)
  // }

  return dt;
}

void
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt,
#if CNS_USE_EB
                   EBFluxRegister* fr_as_crse, EBFluxRegister* fr_as_fine,
#else
                   YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine,
#endif
                   bool write_to_flux_register)
{
  BL_PROFILE("CNS::compute_dSdt()");

  const Real* dx = geom.CellSize();
  const int ncomp = dSdt.nComp();

  int as_crse = (fr_as_crse != nullptr);
  int as_fine = (fr_as_fine != nullptr);

  MultiFab& cost = get_new_data(Cost_Type);

#if CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
  std::array<FArrayBox,AMREX_SPACEDIM> flux;
  FArrayBox dm_as_fine(Box::TheUnitBox(), ncomp);
  FArrayBox fab_drho_as_crse(Box::TheUnitBox(), ncomp);
  IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
  
  for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    amrex::Real wt = amrex::second();
      
#if CNS_USE_EB
    const auto& flag = flags[mfi];
    if (flag.getType(bx) == FabType::covered) {
      dSdt[mfi].setVal<RunOn::Device>(0.0, bx, 0, ncomp);
    } 
    else 
#endif
    {
      // flux is used to store centroid flux needed for reflux
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        flux[idim].resize(amrex::surroundingNodes(bx,idim), LEN_STATE);
        flux[idim].setVal<RunOn::Device>(0.);
      }

#if CNS_USE_EB
      if (flag.getType(amrex::grow(bx,1)) == FabType::regular) {
#endif
        Array4<Real const>    s_arr =    S.array(mfi);
        Array4<Real      > dsdt_arr = dSdt.array(mfi);

        compute_dSdt_box(bx, s_arr, dsdt_arr, {AMREX_D_DECL(&flux[0],&flux[1],&flux[2])});

        if (write_to_flux_register) {
          if (fr_as_crse) {
            fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt, RunOn::Device);
          }

          if (fr_as_fine) {
            fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt, RunOn::Device);
          }
        }
#if CNS_USE_EB
      } else {
        FArrayBox* p_drho_as_crse = (fr_as_crse) ?
            fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
        const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ?
            fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

        if (fr_as_fine) {
          dm_as_fine.resize(amrex::grow(bx,1),ncomp);
        }

        Array4<const Real> vf_arr = (*volfrac).array(mfi);
        Array4<const Real> bcent_arr = (*bndrycent).array(mfi);

        AMREX_D_TERM(Array4<const Real> const& apx = areafrac[0]->const_array(mfi);,
                     Array4<const Real> const& apy = areafrac[1]->const_array(mfi);,
                     Array4<const Real> const& apz = areafrac[2]->const_array(mfi));
        AMREX_D_TERM(Array4<const Real> const& fcx = facecent[0]->const_array(mfi);,
                     Array4<const Real> const& fcy = facecent[1]->const_array(mfi);,
                     Array4<const Real> const& fcz = facecent[2]->const_array(mfi));
        
        Array4<const Real> const&    s_arr =    S.array(mfi);
        Array4<      Real> const& dsdt_arr = dSdt.array(mfi);
        // AMREX_D_TERM(Array4<Real> xflx_arr = flux[0].array(); ,
        //              Array4<Real> yflx_arr = flux[1].array(); ,
        //              Array4<Real> zflx_arr = flux[2].array();)
        
        compute_dSdt_box_eb(bx, s_arr, dsdt_arr, {AMREX_D_DECL(&flux[0],&flux[1],&flux[2])}, //AMREX_D_DECL(xflx_arr, yflx_arr, zflx_arr),
                            flags.const_array(mfi), vf_arr,
                            AMREX_D_DECL(apx, apy, apz), AMREX_D_DECL(fcx, fcy, fcz), bcent_arr,
                            as_crse, p_drho_as_crse->array(), p_rrflag_as_crse->const_array(),
                            as_fine, dm_as_fine.array(), level_mask.const_array(mfi), dt);
        
        if (write_to_flux_register) {
          if (fr_as_crse) {
            fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                dx, dt, (*volfrac)[mfi],
                                {AMREX_D_DECL(&((*areafrac[0])[mfi]),
                                              &((*areafrac[1])[mfi]),
                                              &((*areafrac[2])[mfi]))},
                                RunOn::Device);
          }
          if (fr_as_fine) {
            fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                dx, dt, (*volfrac)[mfi],
                                {AMREX_D_DECL(&((*areafrac[0])[mfi]),
                                              &((*areafrac[1])[mfi]),
                                              &((*areafrac[2])[mfi]))},
                                dm_as_fine, RunOn::Device);
          }
        }
      }
#endif
    } //end EB not covered block

    // Record runtime for load balancing
    Gpu::streamSynchronize();
    wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();
    cost[mfi].plus<RunOn::Device>(wt, bx);

  } //end mfi loop
} //end omp block
}

/**
 * \brief This function enforces (1) all rhoY >= 0, (2) rho <- sum(rhoY), (3) T >= min_T, and update rhoU and rhoE accordingly.
*/
void
CNS::enforce_consistent_state ()
{
  BL_PROFILE("CNS::enforce_consistent_state()");

  // if (verbose > 0) {
  //   amrex::Print() << " >> CNS::enforce_consistent_state()" << std::endl;
  // }

  MultiFab& S = get_new_data(State_Type);
  // auto const& rho = Array4<Real>(state, URHO, 1); //this is how to slice array4!!
  // auto const& rhoU = Array4<Real>(state, UMX, AMREX_SPACEDIM);
  // auto const& rhoY = Array4<Real>(state, UFS, NUM_SPECIES);
  // auto const& rhoE = Array4<Real>(state, UEDEN, 1);

#if CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{  
  for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
      
#if CNS_USE_EB
    const auto& flag = flags[mfi];
    if (flag.getType(bx) != FabType::covered)
#endif
    {
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        Array4 s_arr = S.array(mfi, nf*NVAR);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const IntVect iv{AMREX_D_DECL(i, j, k)};
          const Real rhoOld = s_arr(iv, URHO);
          Real rhoinv = 1.0 / rhoOld;
                    
          // Clip species rhoYs and get new rho
          Real rhoNew = 0.0;
          Real rhoYOld;
          for (int n = 0; n < NUM_SPECIES; n++) {
            rhoYOld = s_arr(iv, UFS+n);
            s_arr(iv, UFS+n) = amrex::min<Real>(rhoOld, amrex::max<Real>(0.0, rhoYOld));
            rhoNew += s_arr(iv, UFS+n);
          }
          s_arr(iv, URHO) = rhoNew;

          // Keep specific energy and velocity, recompute rhoE and mom
          s_arr(iv, UEDEN) *= rhoNew * rhoinv;
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
            s_arr(iv, UMX+n) *= rhoNew * rhoinv;
          }

          // Calculate temperature (!) this adds energy to the point but without removing the excess ke
          rhoinv = 1.0 / rhoNew;
          Real Y[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; n++) {
            Y[n] = s_arr(iv, UFS+n) * rhoinv;
          }
          Real ke = 0.5 * rhoinv * (AMREX_D_TERM(s_arr(iv, UMX)*s_arr(iv, UMX) ,
                                                +s_arr(iv, UMY)*s_arr(iv, UMY) ,
                                                +s_arr(iv, UMZ)*s_arr(iv, UMZ)));
          Real ei = rhoinv * (s_arr(iv, UEDEN) - ke);
          auto eos = pele::physics::PhysicsType::eos();
          Real T;
          eos.REY2T(rhoNew, ei, Y, T);
          if (T < clip_temp) {
            std::cout << "Energy added to cell (" << i << "," << j << "," << k << "): T " << T << "->" << clip_temp << "\n";
            T = clip_temp;
            s_arr(iv, UEDEN) -= ei * rhoNew;
            eos.RTY2E(rhoNew, T, Y, ei);
            s_arr(iv, UEDEN) += ei * rhoNew;
          }
        });
      } // loop through NUM_FIELD
    } // EB !covered
  } // mfi loop
} // omp block
}