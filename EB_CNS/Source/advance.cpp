#if CNS_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "CNS.H"

using namespace amrex;

/** @brief The main driver for a single level implementing the time advance.
 *  @param time the current simulation time
 *  @param dt the timestep to advance (e.g., go from time to time + dt)
 *  @param iteration where we are in the current AMR subcycle. Each
                     level will take a number of steps to reach the
                     final time of the coarser level below it. This
                     counter starts at 1
 *  @param ncycle  the number of subcycles at this level
 */
Real CNS::advance(Real time, Real dt, int iteration, int ncycle)
{
  BL_PROFILE("CNS::advance()");

  // this must be put before we swap the new and old state (swapTimeLevels)
  // but is put at the beginning of the advance step rather than after RK step 2
  // so that state after reflux is also corrected
  // enforce_consistent_state(); // Enforce rho = sum(rhoY) and clip temp

  // Prepare data fabs
  for (int i = 0; i < num_state_data_types; ++i) {
    if (!(i == Reactions_Type &&
          do_react)) { // do not swap I_R (we only use new I_R data)
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }
  }

  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& S_old = get_old_data(State_Type);
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
  if (do_reflux && level > 0) { fr_as_fine = &flux_reg; }

  if (fr_as_crse) { fr_as_crse->reset(); }

  // Start time-stepping
  // RK2 stage 1: U^* = U^n + dt*dUdt^n + dt*I_R^n
  if (verbose > 1)
    amrex::Print() << " >> AMR Cycle " << iteration << " of " << ncycle << std::endl;
  if (verbose > 0)
    amrex::Print() << " >> RK Step 1: Computing dSdt^{n}" << std::endl;
  FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, LEN_STATE);
  compute_dSdt(Sborder, dSdt_old, 0.5 * dt, fr_as_crse, fr_as_fine, true);
  MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt_old, 0, 0, LEN_STATE, 0);

  if (rk_order == 2) {
    if (do_react) {
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        MultiFab::Saxpy(S_new, dt, I_R, nf * NREACT, nf * NVAR + UFS, NUM_SPECIES,
                        0);
      }
    }

    enforce_consistent_state(); // Enforce rho = sum(rhoY) and clip temp

    // RK2 stage 2: U^{n+1} = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1}) + dt*I_R^{n+1}
    if (verbose > 0)
      amrex::Print() << " >> RK Step 2: Computing dSdt^{n+1}" << std::endl;
    FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, LEN_STATE);
    compute_dSdt(Sborder, dSdt_new, 0.5 * dt, fr_as_crse, fr_as_fine,
                 (rk_reaction_iter < 1));
    MultiFab::LinComb(S_new, 0.5 * dt, dSdt_new, 0, 0.5 * dt, dSdt_old, 0, 0,
                      LEN_STATE, 0);
    MultiFab::Add(S_new, S_old, 0, 0, LEN_STATE, 0);
  }

#if (NUM_FIELD > 0)
  computeAvg(S_new);
  FillPatch(*this, Sborder, 1, time + dt, State_Type, 0, LEN_STATE);
  compute_pdf_model(Sborder, dt, iteration);
  MultiFab::Copy(S_new, Sborder, 0, 0, LEN_STATE, 0);
  computeAvg(S_new);
#endif

  enforce_consistent_state(); // Enforce rho = sum(rhoY)

  if (do_react) {
    // Compute I_R^{n+1}(U^**) and do U^{n+1} = U^** + dt*I_R^{n+1}
    react_state(time, dt);
  }

  enforce_consistent_state(); // Enforce rho = sum(rhoY)

  // Iterate to couple chemistry
  if (do_react && (rk_reaction_iter > 1)) {
    for (int iter = 1; iter < rk_reaction_iter; ++iter) {
      if (verbose > 0) {
        amrex::Print() << " >> Re-computing dSdt^{n+1}" << std::endl;
      }
      FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, LEN_STATE);
      compute_dSdt(Sborder, dSdt_new, 0.5 * dt, fr_as_crse, fr_as_fine,
                   (iter == rk_reaction_iter));

      MultiFab::LinComb(S_new, 0.5 * dt, dSdt_new, 0, 0.5 * dt, dSdt_old, 0, 0,
                        LEN_STATE, 0);
      MultiFab::Add(S_new, S_old, 0, 0, LEN_STATE, 0);
      react_state(time, dt);

      enforce_consistent_state(); // Enforce rho = sum(rhoY)
    }
  }

  return dt;
}

/**
 * @brief Compute RHS of Navier-Stokes equation.
 */
void CNS::compute_dSdt(const MultiFab& S, MultiFab& dSdt, Real dt,
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

#if CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    std::array<FArrayBox, AMREX_SPACEDIM> flux;

    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      amrex::Real wt = amrex::ParallelDescriptor::second();

#if CNS_USE_EB
      const auto& flag = flags[mfi];
      if (flag.getType(bx) == FabType::covered) {
        dSdt[mfi].setVal<RunOn::Device>(0.0, bx, 0, ncomp);
      } else
#endif
      {
        // flux is used to store centroid flux needed for reflux
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          flux[idim].resize(amrex::surroundingNodes(bx, idim), LEN_STATE,
                            The_Async_Arena());
          flux[idim].setVal<RunOn::Device>(0.);
        }

#if CNS_USE_EB
        if (flag.getType(amrex::grow(bx, 2)) == FabType::regular) {
#endif
          Array4<const Real> s_arr = S.array(mfi);
          Array4<Real> dsdt_arr = dSdt.array(mfi);

          compute_dSdt_box(bx, s_arr, dsdt_arr,
                           {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])});

          if (write_to_flux_register) {
            if (fr_as_crse) {
              fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                  dx, dt, RunOn::Device);
            }

            if (fr_as_fine) {
              fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                  dx, dt, RunOn::Device);
            }
          }
#if CNS_USE_EB
        } else {
          Array4<const Real> vf_arr = (*volfrac).array(mfi);
          Array4<const Real> bcent_arr = (*bndrycent).array(mfi);

          AMREX_D_TERM(
            Array4<const Real> const& apx = areafrac[0]->const_array(mfi);
            , Array4<const Real> const& apy = areafrac[1]->const_array(mfi);
            , Array4<const Real> const& apz = areafrac[2]->const_array(mfi));
          AMREX_D_TERM(
            Array4<const Real> const& fcx = facecent[0]->const_array(mfi);
            , Array4<const Real> const& fcy = facecent[1]->const_array(mfi);
            , Array4<const Real> const& fcz = facecent[2]->const_array(mfi));

          Array4<const Real> s_arr = S.array(mfi);
          Array4<Real> dsdt_arr = dSdt.array(mfi);

          compute_dSdt_box_eb(
            bx, s_arr, dsdt_arr, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
            flags.const_array(mfi), vf_arr, AMREX_D_DECL(apx, apy, apz),
            AMREX_D_DECL(fcx, fcy, fcz), bcent_arr,
            // , p_drho_as_crse->array(), p_rrflag_as_crse->const_array(),
            // as_fine, level_mask.const_array(mfi),
            dt);

          if (write_to_flux_register) {
            if (fr_as_crse) {
              fr_as_crse->CrseAdd(
                mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt,
                (*volfrac)[mfi],
                {AMREX_D_DECL(&((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
                              &((*areafrac[2])[mfi]))},
                RunOn::Device);
            }
            if (fr_as_fine) {
              // dm_as_fine is the mass flux due to redistrution, it is 0 because
              // C/F interface does *not* cross EB, same as in PeleC.
              FArrayBox dm_as_fine(amrex::grow(bx, 1), ncomp, The_Async_Arena());
              dm_as_fine.setVal<RunOn::Device>(0.0);

              fr_as_fine->FineAdd(
                mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt,
                (*volfrac)[mfi],
                {AMREX_D_DECL(&((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
                              &((*areafrac[2])[mfi]))},
                dm_as_fine, RunOn::Device);
            }
          }
        }
#endif
      } // end EB not covered block

      // Record runtime for load balancing
      Gpu::streamSynchronize();
      wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();
      get_new_data(Cost_Type)[mfi].plus<RunOn::Device>(wt, bx);

    } // end mfi loop
  }   // end omp block
}

/**
 * \brief This function enforces (1) all rhoY >= 0, (2) rho <- sum(rhoY),
 *        (3) T >= min_T, and update rhoU and rhoE accordingly. This is
 *        called before each compute_dSdt.
 */
void CNS::enforce_consistent_state()
{
  BL_PROFILE("CNS::enforce_consistent_state()");

  MultiFab& S = get_new_data(State_Type);

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
      const auto& flag_arr = flags.const_array(mfi);
      if (flag.getType(bx) != FabType::covered)
#endif
      {
        for (int nf = 0; nf <= NUM_FIELD; ++nf) {
          Array4 s_arr = S.array(mfi, nf * NVAR);
          FArrayBox Temp(bx, 1, The_Async_Arena());
          Array4 T = Temp.array();
          ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const IntVect iv{AMREX_D_DECL(i, j, k)};
            const Real rhoOld = s_arr(iv, URHO);
            if (rhoOld <= 0.0 && verbose > 1)
              std::cout << "enforce_consistent_state(): rhoOld" << iv << "="
                        << rhoOld << "!\n";
            Real rhoinv = 1.0 / rhoOld;

            // Clip species rhoYs and get new rho
            Real rhoNew = 0.0;
            Real rhoYOld;
            for (int n = 0; n < NUM_SPECIES; n++) {
              s_arr(iv, UFS + n) = amrex::max<Real>(0.0, s_arr(iv, UFS + n));
              rhoNew += s_arr(iv, UFS + n);
            }
            if (rhoNew <= 0.0) {
              if (verbose > 1)
                std::cout << "enforce_consistent_state(): rhoNew" << iv << "="
                          << rhoNew << "! Fixed \n";

              rhoNew = d_parm->smallr;
              // s_arr(iv, UFS + N2_ID) = rhoNew;
              for (int n = 0; n < NUM_SPECIES; n++) {
                s_arr(iv, UFS + n) = rhoNew / NUM_SPECIES;
              }
            }
            s_arr(iv, URHO) = rhoNew;

            // Keep specific energy and velocity, recompute rhoE and mom
            s_arr(iv, UEDEN) *= rhoNew * rhoinv;
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              s_arr(iv, UMX + n) *= rhoNew * rhoinv;
            }

            // Calculate temperature
            rhoinv = 1.0 / rhoNew;
            Real Y[NUM_SPECIES];
            for (int n = 0; n < NUM_SPECIES; n++) {
              Y[n] = s_arr(iv, UFS + n) * rhoinv;
            }
            const Real ke = 0.5 * rhoinv * rhoinv *
                      (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
                                    +s_arr(iv, UMY) * s_arr(iv, UMY),
                                    +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
            Real ei = rhoinv * s_arr(iv, UEDEN) - ke;
            auto eos = pele::physics::PhysicsType::eos();
            T(i, j, k) = 0.0;
            eos.REY2T(rhoNew, ei, Y, T(i, j, k));
          });

          ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const IntVect iv{AMREX_D_DECL(i, j, k)};

            if (T(iv) < clip_temp) {            
              // Find target temp by averaging over the neighbours
              // Real target_temp = 0.0;
              Real target_temp = clip_temp;
              
              // int n_neighbour = 0;
              int n_neighbour = 1;
              for (int ix = -1; ix <= 1; ++ix) {
              for (int iy = -1; iy <= 1; ++iy) {
              for (int iz = -1; iz <= 1; ++iz) {
              Real T_neighbour = T(i+ix, j+iy, k+iz);
              if ((T_neighbour > clip_temp) && (T_neighbour < 1e4) 
#if CNS_USE_EB
                  && flag_arr(i+ix, j+iy, k+iz).isRegular()
#endif
                ) {
                  target_temp += T_neighbour;
                  n_neighbour += 1;
                }
              }
              }
              }
              // if (n_neighbour > 0)
                target_temp /= Real(n_neighbour);
              // else
              //   target_temp = clip_temp;

              // Calculate temperature
              const Real rhoNew = s_arr(iv, URHO);
              const Real rhoinv = 1.0 / rhoNew;
              Real Y[NUM_SPECIES];
              for (int n = 0; n < NUM_SPECIES; n++) {
                Y[n] = s_arr(iv, UFS + n) * rhoinv;
              }
              const Real ke = 0.5 * rhoinv * rhoinv *
                        (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
                                      +s_arr(iv, UMY) * s_arr(iv, UMY),
                                      +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
              Real ei = rhoinv * s_arr(iv, UEDEN) - ke;
              auto eos = pele::physics::PhysicsType::eos();

              // // Keep KE, add energy to increase T
              // if (verbose > 1)
              //   std::cout << "Energy added to cell " << iv
              //             << "): T=" << T(iv) << "->" << target_temp << " rhoE "
              //             << s_arr(iv, UEDEN);
              // T(iv) = target_temp;
              // s_arr(iv, UEDEN) -= ei * rhoNew;
              // eos.RTY2E(rhoNew, T(iv), Y, ei);
              // s_arr(iv, UEDEN) += ei * rhoNew;
              // if (verbose > 1) std::cout << "->" << s_arr(iv, UEDEN) << std::endl;
              
              // // Keep total E, reduce vel
              // if (verbose > 1)
              //   std::cout << "Momentum removed from cell (" << i << "," << j << "," << k
              //             << "): T " << T << "->" << clip_temp << '\n';
              // // Correct
              // Real ei_new;
              // eos.RTY2E(rhoNew, clip_temp, Y, ei_new);
              // Real frac = (std::sqrt(1.0 - (ei_new - ei) / ke)); // this may be negative
              // AMREX_D_TERM(s_arr(iv, UMX) *= frac;,
              //              s_arr(iv, UMY) *= frac;,
              //              s_arr(iv, UMZ) *= frac;)

              // // Check
              // ke = 0.5 * rhoinv *
              //         (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
              //                       +s_arr(iv, UMY) * s_arr(iv, UMY),
              //                       +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
              // Real rE_new = rhoNew * (ei_new + ke);
              // if (rE_new != s_arr(iv, UEDEN)) {
              //   std::cout << rE_new << "<-" << s_arr(iv, UEDEN) << '\n';
              //   amrex::Abort();
              // }

              // Keep total E, reduce vel
              if (verbose > 1)
                std::cout << "Momentum removed from cell " << iv
                          << ": T=" << T(iv) << "->" << target_temp << " rhoE=" << s_arr(iv, UEDEN);
              // First take energy from KE
              Real ei_new, diff_ei;
              eos.RTY2E(rhoNew, target_temp, Y, ei_new);
              diff_ei = ei_new - ei;
              Real ctrl_parm = 1.0; // control param in [0,1] how much KE can be used to heat the cell
              Real fac = std::sqrt(1.0 - amrex::min(diff_ei, ctrl_parm * ke) / ke);
              AMREX_D_TERM(s_arr(iv, UMX) *= fac;,
                           s_arr(iv, UMY) *= fac;,
                           s_arr(iv, UMZ) *= fac;)
              // If not enough, add energy
              s_arr(iv, UEDEN) += rhoNew * amrex::max(diff_ei - ctrl_parm * ke, 0.0);
              if (verbose > 1) std::cout << "->" << s_arr(iv, UEDEN) << '\n';

              // // Check
              // Real ke_new = 0.5 * rhoinv * rhoinv * 
              //         (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
              //                       +s_arr(iv, UMY) * s_arr(iv, UMY),
              //                       +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
              // Real ke_chk = ke - amrex::min(diff_ei, ctrl_parm * ke);
              // Real ei_chk = rhoinv * s_arr(iv, UEDEN) - ke_new;
              // if (std::abs(ei_chk - ei_new) > d_parm->smallr) { //std::numeric_limits<amrex::Real>::epsilon()
              //   std::cout << "|ei_chk-ei_new|=" << std::abs(ei_chk - ei_new) << " ke " << ke << "->" << ke_new << '\n';
              //   std::cout << "|ke_chk-ke_new|=" << std::abs(ke_chk - ke_new) << " fac=" << fac << '\n';
              //   // amrex::Abort();
              // }
            }
          });
        } // loop through NUM_FIELD
      }   // EB !covered
    }     // mfi loop
  }       // omp block
}
