#if CNS_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include <cmath>

#include "CNS.H"
#include "derive.H"

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

  // Prepare data fabs
  for (int i = 0; i < num_state_data_types; ++i) {
    // do not swap I_R (we only use new I_R data)
    if (!(i == Reactions_Type && do_react)) {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }
  }
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

  if (rk_order <= 2) {
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab dSdt_old(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
    MultiFab dSdt_new(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
    MultiFab& I_R = get_new_data(Reactions_Type);

    // Predictor-corrector RK2
    // RK2 stage 1: U^* = U^n + dt*dUdt^n + dt*I_R^n
    if (verbose > 1) {
      amrex::Print() << " >> AMR Cycle " << iteration << " of " << ncycle << "\n";
    }
    if (verbose > 0) { amrex::Print() << " >> RK Step 1: Computing dSdt^{n}\n"; }
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
      enforce_consistent_state(S_new); // Enforce rho = sum(rhoY) and clip temp

      // RK2 stage 2: U^{n+1} = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1}) + dt*I_R^{n+1}
      if (verbose > 0) { amrex::Print() << " >> RK Step 2: Computing dSdt^{n+1}\n"; }
      FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, LEN_STATE);
      compute_dSdt(Sborder, dSdt_new, 0.5 * dt, fr_as_crse, fr_as_fine,
                   (rk_reaction_iter < 1));
      MultiFab::LinComb(S_new, 0.5 * dt, dSdt_new, 0, 0.5 * dt, dSdt_old, 0, 0,
                        LEN_STATE, 0);
      MultiFab::Add(S_new, S_old, 0, 0, LEN_STATE, 0);
    }
    enforce_consistent_state(S_new); // Enforce rho = sum(rhoY)

#if (NUM_FIELD > 0)
    computeAvg(S_new);
    FillPatch(*this, Sborder, 1, time + dt, State_Type, 0, LEN_STATE);
    compute_pdf_model(Sborder, dt, iteration);
    MultiFab::Copy(S_new, Sborder, 0, 0, LEN_STATE, 0);
    // enforce_consistent_state(S_new); // Enforce rho = sum(rhoY)
    // computeAvg(S_new);
#endif

    if (do_react) {
      // Compute I_R^{n+1}(U^**) and do U^{n+1} = U^** + dt*I_R^{n+1}
      react_state(time, dt);
      enforce_consistent_state(S_new); // Enforce rho = sum(rhoY)
    }

    // Iterate to couple chemistry
    if (do_react && (rk_reaction_iter > 1)) {
      for (int iter = 1; iter < rk_reaction_iter; ++iter) {
        if (verbose > 0) { amrex::Print() << " >> Re-computing dSdt^{n+1}\n"; }
        FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, LEN_STATE);
        compute_dSdt(Sborder, dSdt_new, 0.5 * dt, fr_as_crse, fr_as_fine,
                     (iter == rk_reaction_iter));

        MultiFab::LinComb(S_new, 0.5 * dt, dSdt_new, 0, 0.5 * dt, dSdt_old, 0, 0,
                          LEN_STATE, 0);
        MultiFab::Add(S_new, S_old, 0, 0, LEN_STATE, 0);
        react_state(time, dt);

        enforce_consistent_state(S_new); // Enforce rho = sum(rhoY)
      }
    }
  } else {
    // Use AMRLevel RK (no reaction, no SF for now)
    if (do_react || NUM_FIELD > 0) {
      amrex::Abort("AMRLevel RK not implemented for reaction or SF");
    }

    RK(rk_order, State_Type, time, dt, iteration, ncycle,
       [&](int stage, MultiFab& dSdt, MultiFab const& S, Real /*t*/, Real dtsub) {
         if (verbose > 0) {
           amrex::Print() << " >> AMRLevel::RK Step " << stage << '\n';
         }
         compute_dSdt(S, dSdt, dtsub, fr_as_crse, fr_as_fine);
       });

    enforce_consistent_state();
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
  const int ncomp = UFA;

#if CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    std::array<FArrayBox, AMREX_SPACEDIM> flux;
#if CNS_USE_EB
    FArrayBox dm_as_fine(Box::TheUnitBox(), ncomp, The_Async_Arena());
    FArrayBox fab_drho_as_crse(Box::TheUnitBox(), ncomp, The_Async_Arena());
    IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
#endif

    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      amrex::Real wt = amrex::ParallelDescriptor::second();

#if CNS_USE_EB
      const auto& flag = flags[mfi];
      if (flag.getType(bx) == FabType::covered) {
        dSdt[mfi].setVal<RunOn::Device>(0.0, bx, 0, LEN_STATE);
      } else
#endif
      {
        // flux is used to store centroid flux needed for reflux
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          flux[idim].resize(amrex::surroundingNodes(bx, idim), ncomp,
                            The_Async_Arena());
          flux[idim].setVal<RunOn::Device>(0.);
        }

        // shock sensor always needed in cns_riemann
        {
          Real time = state[State_Type].curTime();
          int* bcrec_dummy;
          cns_dershocksensor(amrex::grow(bx, 2), shock_sensor_mf[mfi], 0, 1, S[mfi],
                             geom, time, bcrec_dummy, level);
        }

        Array4<const Real> s_arr = S.array(mfi);
        Array4<Real> dsdt_arr = dSdt.array(mfi);
        Array4<const Real> ss_arr = shock_sensor_mf.array(mfi);

#if CNS_USE_EB
        if (flag.getType(amrex::grow(bx, 3)) == FabType::regular) {
#endif
          compute_dSdt_box(bx, s_arr, dsdt_arr,
                           {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, ss_arr);

          if (write_to_flux_register && do_reflux) {
            if (fr_as_crse) {
              fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                  dx, dt, 0, 0, ncomp, RunOn::Device);
            }
            if (fr_as_fine) {
              fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                  dx, dt, 0, 0, ncomp, RunOn::Device);
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

          FArrayBox* p_drho_as_crse =
            (fr_as_crse) ? fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
          const IArrayBox* p_rrflag_as_crse =
            (fr_as_crse) ? fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;
          if (fr_as_fine) {
            dm_as_fine.resize(amrex::grow(bx, 1), ncomp);
            dm_as_fine.setVal<RunOn::Device>(0.0);
          }

          compute_dSdt_box_eb(
            bx, s_arr, dsdt_arr, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
            flags.const_array(mfi), vf_arr, AMREX_D_DECL(apx, apy, apz),
            AMREX_D_DECL(fcx, fcy, fcz), bcent_arr, int(fr_as_crse != nullptr),
            p_drho_as_crse->array(), p_rrflag_as_crse->const_array(),
            int(fr_as_fine != nullptr), dm_as_fine.array(),
            level_mask.const_array(mfi), dt, ss_arr);

          if (write_to_flux_register && do_reflux) {
            if (fr_as_crse) {
              fr_as_crse->CrseAdd(
                mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt,
                (*volfrac)[mfi],
                {AMREX_D_DECL(&((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
                              &((*areafrac[2])[mfi]))},
                0, 0, ncomp, RunOn::Device);
            }
            if (fr_as_fine) {
              fr_as_fine->FineAdd(
                mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt,
                (*volfrac)[mfi],
                {AMREX_D_DECL(&((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
                              &((*areafrac[2])[mfi]))},
                dm_as_fine, 0, 0, ncomp, RunOn::Device);
            }
          }
        }
#endif
      } // end EB not covered block

      // Record runtime for load balancing
      if (do_load_balance) {
        Gpu::streamSynchronize();
        wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();
        get_new_data(Cost_Type)[mfi].plus<RunOn::Device>(wt, bx);
      }
    } // end mfi loop
  }   // end omp block
}

void CNS::enforce_consistent_state()
{
  MultiFab& S_new = get_new_data(State_Type);
  enforce_consistent_state(S_new);
}

/**
 * \brief This function (1) enforces all rhoY >= 0, (2) rho <- sum(rhoY),
 *        (3) enforces T >= min_T and update rhoU and rhoE accordingly,
 *        (4) compute T and stores in UTEMP.
 *        This is called before each compute_dSdt.
 */
void CNS::enforce_consistent_state(MultiFab& S)
{
  BL_PROFILE("CNS::enforce_consistent_state()");

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
#if (NUM_FIELD > 0)
        for (int nf = 1; nf <= NUM_FIELD; ++nf)
#else
        const int nf = 0;
#endif
        {
          auto s_arr = S.array(mfi, nf * NVAR);
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const IntVect iv{AMREX_D_DECL(i, j, k)};

            // Clip species rhoYs and get new rho
            Real rhoNew = 0.0;
            for (int ns = 0; ns < NUM_SPECIES; ++ns) {
              s_arr(iv, UFS + ns) = amrex::max<Real>(0.0, s_arr(iv, UFS + ns));
              rhoNew += s_arr(iv, UFS + ns);
            }
            if (rhoNew <= 0.0) {
#ifndef AMREX_USE_GPU
              if (verbose > 1) {
                std::cout << "enforce_consistent_state: rhoNew" << iv << "," << nf
                          << "=" << rhoNew << "! Fixed \n";
              }
#endif
              rhoNew = CNSConstants::smallr;
#if defined(N2_ID)
              s_arr(iv, UFS + N2_ID) = rhoNew;
#elif defined(AR_ID)
              s_arr(iv, UFS + AR_ID) = rhoNew;
#else
              for (int ns = 0; ns < NUM_SPECIES; ns++) {
                s_arr(iv, UFS + ns) = rhoNew / Real(NUM_SPECIES);
              }
#endif
            }
            s_arr(iv, URHO) = rhoNew;

            // Calculate temperature
            Real rhoinv = Real(1.0) / rhoNew;
            Real Y[NUM_SPECIES];
            for (int ns = 0; ns < NUM_SPECIES; ++ns) {
              Y[ns] = s_arr(iv, UFS + ns) * rhoinv;
            }
            const Real rke = 0.5 * rhoinv *
                             (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
                                           +s_arr(iv, UMY) * s_arr(iv, UMY),
                                           +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
            Real ei = rhoinv * (s_arr(iv, UEDEN) - rke);
            auto eos = pele::physics::PhysicsType::eos();
            Real T = s_arr(iv, UTEMP);
            eos.REY2T(rhoNew, ei, Y, T);
            s_arr(iv, UTEMP) = T;
          });

          amrex::ParallelFor(
            bx, NUM_FIELD + 1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int nf) noexcept {
              const IntVect iv{AMREX_D_DECL(i, j, k)};

              if (s_arr(iv, UTEMP) < clip_temp) {
                // Find target temp by averaging over the neighbours
                Real target_temp = clip_temp;
                int n_neighbour = 1;
                for (int ix = -1; ix <= 1; ++ix) {
                  for (int iy = -1; iy <= 1; ++iy) {
                    for (int iz = -1; iz <= 1; ++iz) {
                      Real T_neighbour = s_arr(i + ix, j + iy, k + iz, UTEMP);
                      if (std::isgreater(T_neighbour, clip_temp) &&
                          std::isless(T_neighbour, 10000)
#if CNS_USE_EB
                          && flag_arr(i + ix, j + iy, k + iz).isRegular()
#endif
                      ) {
                        target_temp += T_neighbour;
                        n_neighbour += 1;
                      }
                    } // iz
                  }   // iy
                }     // ix
                target_temp /= Real(n_neighbour);
                target_temp = amrex::min(target_temp, 200.0);

                // Calculate new temperature
                const Real rhoNew = s_arr(iv, URHO);
                const Real rhoinv = Real(1.0) / rhoNew;
                Real Y[NUM_SPECIES];
                for (int ns = 0; ns < NUM_SPECIES; ++ns) {
                  Y[ns] = s_arr(iv, UFS + ns) * rhoinv;
                }
                const Real ke = 0.5 * rhoinv * rhoinv *
                                (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
                                              +s_arr(iv, UMY) * s_arr(iv, UMY),
                                              +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
                Real ei = rhoinv * s_arr(iv, UEDEN) - ke;
                auto eos = pele::physics::PhysicsType::eos();

                // // Option 1: Keep KE, add ei to increase T
                // if (verbose > 1) {
                //   std::cout << "Energy added to cell " << iv
                //             << ": T=" << s_arr(iv, UTEMP) << "->" << target_temp
                //             << " ei=" << ei << "->";
                // }
                // s_arr(iv, UTEMP) = target_temp;
                // s_arr(iv, UEDEN) -= ei * rhoNew;
                // eos.RTY2E(rhoNew, target_temp, Y, ei);
                // s_arr(iv, UEDEN) += ei * rhoNew;
                // if (verbose > 1) { std::cout << ei << "\n"; }

                // Option 2: Keep total E, reduce vel
                // (because it usually fails at high ke/ei region)
#ifndef AMREX_USE_GPU
                if (verbose > 1) {
                  std::cout << "Momentum removed from cell " << iv
                            << ": T=" << s_arr(iv, UTEMP) << "->" << target_temp
                            << " rhoE=" << s_arr(iv, UEDEN) << "->";
                }
#endif
                s_arr(iv, UTEMP) = target_temp;
                // First take energy from KE
                Real ei_new;
                eos.RTY2E(rhoNew, target_temp, Y, ei_new);
                Real diff_ei = ei_new - ei;
                Real ctrl_parm = 1.0; // control param in [0,1] how much KE can be
                                      // used to heat the cell
                Real fac = std::sqrt(1.0 - amrex::min(diff_ei, ctrl_parm * ke) / ke);
                AMREX_D_TERM(s_arr(iv, UMX) *= fac;, s_arr(iv, UMY) *= fac;
                             , s_arr(iv, UMZ) *= fac;)
                // If not enough, add energy
                s_arr(iv, UEDEN) +=
                  rhoNew * amrex::max(diff_ei - ctrl_parm * ke, 0.0);
#ifndef AMREX_USE_GPU
                if (verbose > 1) { std::cout << s_arr(iv, UEDEN) << '\n'; }
#endif
              } // if T < clip_temp

              if (s_arr(iv, UTEMP) > 5000.0) {
                Real target_temp = 5000.0;

                // Calculate new temperature
                const Real rhoNew = s_arr(iv, URHO);
                const Real rhoinv = Real(1.0) / rhoNew;
                Real Y[NUM_SPECIES];
                for (int ns = 0; ns < NUM_SPECIES; ++ns) {
                  Y[ns] = s_arr(iv, UFS + ns) * rhoinv;
                }
                const Real ke = 0.5 * rhoinv * rhoinv *
                                (AMREX_D_TERM(s_arr(iv, UMX) * s_arr(iv, UMX),
                                              +s_arr(iv, UMY) * s_arr(iv, UMY),
                                              +s_arr(iv, UMZ) * s_arr(iv, UMZ)));
                Real ei = rhoinv * s_arr(iv, UEDEN) - ke;
                auto eos = pele::physics::PhysicsType::eos();

                // Option 1: Keep KE, add (remove) ei to increase (decrease) T
#ifndef AMREX_USE_GPU
                if (verbose > 1) {
                  std::cout << "Energy added to cell " << iv
                            << ": T=" << s_arr(iv, UTEMP) << "->" << target_temp
                            << " ei=" << ei << "->";
                }
#endif
                s_arr(iv, UTEMP) = target_temp;
                s_arr(iv, UEDEN) -= ei * rhoNew;
                eos.RTY2E(rhoNew, target_temp, Y, ei);
                s_arr(iv, UEDEN) += ei * rhoNew;
#ifndef AMREX_USE_GPU
                if (verbose > 1) { std::cout << ei << "\n"; }
#endif
              } // if T > 5000
            });
        } // nf loop
      }   // EB !covered
    }     // mfi loop
  }       // omp block
}
