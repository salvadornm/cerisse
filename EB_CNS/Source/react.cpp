#include <AMReX_FArrayBox.H>

#include "CNS.H"
#include "PaSR.H"

using namespace amrex;

/**
 * @brief TODO: Set typical values to help ODE solver.
 */
void CNS::set_typical_values_chem()
{
  // amrex::MultiFab& S_new = get_new_data(State_Type);
  // // amrex::Real minTemp = S_new.min(UTEMP);
  // // amrex::Real maxTemp = S_new.max(UTEMP);
  // amrex::Vector<amrex::Real> typical_values_chem(NUM_SPECIES + 1, 1e-10);

  // for (int n = 0; n < NUM_SPECIES; n++) {
  //   amrex::Real rhoYs_min = S_new.min(UFS + n);
  //   amrex::Real rhoYs_max = S_new.max(UFS + n);
  //   typical_values_chem[n] = amrex::max<amrex::Real>(
  //     0.5 * (rhoYs_min + rhoYs_max), 1.e-10);
  // }
  // // typical_values_chem[NUM_SPECIES] = 0.5 * (minTemp + maxTemp);
  // typical_values_chem[NUM_SPECIES] = 1500.0; // temp
  // reactor->set_typ_vals_ode(typical_values_chem);
}

/**
 * \brief Compute I_R and update S_new += dt*I_R
 */
void CNS::react_state(Real time, Real dt, bool init_react)
{
  BL_PROFILE("CNS::react_state()");
  const Real strt_time = amrex::ParallelDescriptor::second();

  if ((verbose > 0) && amrex::ParallelDescriptor::IOProcessor()) {
    if (init_react) {
      amrex::Print() << "Initialising reactions, using interval dt = " << dt
                     << std::endl;
    } else {
      amrex::Print() << " >> Computing reactions" << std::endl;
    }
  }

  // State Fabs
  MultiFab Sold(grids, dmap, LEN_STATE, 1, MFInfo(),
                Factory()); //= get_old_data(State_Type);
  if (!init_react) FillPatch(*this, Sold, 1, time, State_Type, 0, LEN_STATE);
  MultiFab& Snew = get_new_data(State_Type);
  MultiFab& I_R = get_new_data(Reactions_Type);
  I_R.setVal(0.0);

  // Fab for all Array4 in sundials box integration
  MultiFab STemp(grids, dmap, NUM_SPECIES + 2, 0);    //[rY, rEi, T]
  MultiFab SDotTemp(grids, dmap, NUM_SPECIES + 1, 0); // d[rY, rEi]/dt
  iMultiFab maskFab(grids, dmap, 1, 0); //= 1: do reaction, = -1: don't do reaction
  MultiFab fctCount(grids, dmap, 1, 0); // number of RHS evaluations in reactor

  STemp.setVal(0.0);
  SDotTemp.setVal(0.0);
  fctCount.setVal(0.0);

#if CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(Snew.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (MFIter mfi(I_R, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box bx = mfi.tilebox();
      Real wt = amrex::ParallelDescriptor::second(); // timing for each fab

#if CNS_USE_EB
      Array4<const Real> vfrac_arr = (*volfrac).const_array(mfi);
      if (flags[mfi].getType(bx) != amrex::FabType::covered)
#endif
      {
        bool do_react_fields = (NUM_FIELD > 0); // && do_spdf;
        int nf_start = do_react_fields ? 1 : 0;
        int nf_end = do_react_fields ? NUM_FIELD : 0;
        for (int nf = nf_start; nf <= nf_end; ++nf) {
          ///////////////////// Prepare for react /////////////////////
          const Array4 sold_arr =
            init_react
              ? Snew.array(mfi, nf * NVAR)
              : Sold.array(mfi, nf * NVAR); // don't have Sold at initialisation
          Array4 snew_arr = Snew.array(mfi, nf * NVAR);
          Array4 I_R_mean_arr = I_R.array(mfi, 0);
          Array4 I_R_arr = I_R.array(mfi, nf * NREACT);
          Array4 rY = STemp.array(mfi, 0);
          Array4 rEi = STemp.array(mfi, NUM_SPECIES);
          Array4 T = STemp.array(mfi, NUM_SPECIES + 1);
          Array4 rYsrc = SDotTemp.array(mfi, 0);
          Array4 rEisrc = SDotTemp.array(mfi, NUM_SPECIES);
          Array4 mask = maskFab.array(mfi);
          Array4 fc = fctCount.array(mfi);

          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // calculate [rY, rEi, T]
            const Real rho_old = sold_arr(i, j, k, URHO);
            Real Y[NUM_SPECIES];
            for (int n = 0; n < NUM_SPECIES; ++n) {
              rY(i, j, k, n) = sold_arr(i, j, k, UFS + n);
              Y[n] = sold_arr(i, j, k, UFS + n) / rho_old;
            }

            rEi(i, j, k) =
              sold_arr(i, j, k, UEDEN) -
              0.5 *
                (AMREX_D_TERM(sold_arr(i, j, k, UMX) * sold_arr(i, j, k, UMX),
                              +sold_arr(i, j, k, UMY) * sold_arr(i, j, k, UMY),
                              +sold_arr(i, j, k, UMZ) * sold_arr(i, j, k, UMZ))) /
                rho_old;

            const Real ei = rEi(i, j, k) / rho_old;
            auto eos = pele::physics::PhysicsType::eos();
            eos.REY2T(rho_old, ei, Y, T(i, j, k));

            if (T(i, j, k) <= 0.0)
              std::cout << "Before reaction T=" << T(i, j, k) << " @ " << i << ","
                        << j << "," << k << std::endl;

            // calculate d[rY, rEi]/dt
            for (int n = 0; n < NUM_SPECIES; ++n) {
              rYsrc(i, j, k, n) = (snew_arr(i, j, k, UFS + n) - rY(i, j, k, n)) / dt;
            }

            const Real rEinew =
              snew_arr(i, j, k, UEDEN) -
              0.5 *
                (AMREX_D_TERM(snew_arr(i, j, k, UMX) * snew_arr(i, j, k, UMX),
                              +snew_arr(i, j, k, UMY) * snew_arr(i, j, k, UMY),
                              +snew_arr(i, j, k, UMZ) * snew_arr(i, j, k, UMZ))) /
                snew_arr(i, j, k, URHO);
            rEisrc(i, j, k) = (rEinew - rEi(i, j, k)) / dt;

            // fill mask
            mask(i, j, k) = (T(i, j, k) > min_react_temp) ? 1 : -1;
#if CNS_USE_EB
            mask(i, j, k) = (vfrac_arr(i, j, k) > 0.0) ? mask(i, j, k) : -1;
#endif
          });

          /////////////////////////// React ///////////////////////////
          Real current_time = 0.0;
          reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time
#ifdef AMREX_USE_GPU
                         ,
                         amrex::Gpu::gpuStream()
#endif
          );
          amrex::Gpu::Device::streamSynchronize();

          //////////////////////// Unpack data ////////////////////////
          // Prepare (rho, velocities, mu) for PaSR
          FArrayBox qfab, mufab;
          if (do_pasr && !init_react) {
            const Box bxg1 = amrex::grow(bx, 1);
            qfab.resize(bxg1, 4, The_Async_Arena()); // [rho, u, v, w]
            auto const& qarr = qfab.array();
            amrex::ParallelFor(
              bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                qarr(i, j, k, QRHO) = sold_arr(i, j, k, URHO);
                AMREX_D_TERM(qarr(i, j, k, QU) =
                               sold_arr(i, j, k, UMX) / sold_arr(i, j, k, URHO);
                             , qarr(i, j, k, QV) =
                                 sold_arr(i, j, k, UMY) / sold_arr(i, j, k, URHO);
                             , qarr(i, j, k, QW) =
                                 sold_arr(i, j, k, UMZ) / sold_arr(i, j, k, URHO);)
              });

            mufab.resize(bx, 1, The_Async_Arena()); // [mu]
            auto const& mu = mufab.array();
            amrex::ParallelFor(
              bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Copy to input
                Real muloc;
                Real xiloc, lamloc, Ddiag[NUM_SPECIES]; // not used
                Real Tin = T(i, j, k);
                Real rhoin = sold_arr(i, j, k, URHO);
                Real yin[NUM_SPECIES];
                for (int n = 0; n < NUM_SPECIES; ++n) {
                  yin[n] = sold_arr(i, j, k, UFS + n) / rhoin;
                }

                constexpr bool get_xi = false;
                constexpr bool get_mu = true;
                constexpr bool get_lam = false;
                constexpr bool get_Ddiag = false;
                auto trans = pele::physics::PhysicsType::transport();
                auto const* ltransparm = trans_parms.device_trans_parm();
                trans.transport(get_xi, get_mu, get_lam, get_Ddiag, Tin, rhoin, yin,
                                Ddiag, muloc, xiloc, lamloc, ltransparm);

                // Copy to output
                mu(i, j, k) = muloc;
              });
          }

          // For unpack_pasr
          auto const& qarr = qfab.array();
          auto const& muarr = mufab.array();
          const auto dx = geom.CellSizeArray();
          const auto dxinv = geom.InvCellSizeArray();

          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (mask(i, j, k) != -1) {
              // Monitor problem cell
              bool any_rY_unbounded = false;
              for (int n = 0; n < NUM_SPECIES; ++n)
                any_rY_unbounded |=
                  (rY(i, j, k, n) < -1e-5 || rY(i, j, k, n) > 1.0 + 1e-5 ||
                   std::isnan(rY(i, j, k, n)));
              if (any_rY_unbounded) {
                std::cout << "Reaction causing rY=[ ";
                for (int n = 0; n < NUM_SPECIES; ++n)
                  std::cout << rY(i, j, k, n) << " ";
                std::cout << "] @ " << i << "," << j << "," << k << '\n';
              }
              if (T(i, j, k) < clip_temp) {
                std::cout << "Reaction causing T=" << T(i, j, k) << " @ " << i << ","
                          << j << "," << k << '\n';
              }

              // update U^{n+1} = U^** + dt*I_R^{n+1}
              if (!init_react) {
                amrex::Real new_rho = 0.0;
                if (do_pasr) {
                  unpack_pasr(i, j, k, snew_arr, new_rho, sold_arr, rY, rYsrc, qarr,
                              muarr, dx, dxinv, dt);
                } else {
                  // amrex::Real min_rY = 0.0;
                  // for (int n = 0; n < NUM_SPECIES; ++n) {
                  //   min_rY = amrex::min(min_rY, rY(i, j, k, n));
                  // }
                  for (int n = 0; n < NUM_SPECIES; ++n) {
                    snew_arr(i, j, k, UFS + n) = amrex::max(0.0, rY(i, j, k, n));
                    // snew_arr(i, j, k, UFS + n) = rY(i, j, k, n) - min_rY;
                    new_rho += snew_arr(i, j, k, UFS + n);
                  }
                }

                // Enforce conservation of rho
                // if (std::abs(snew_arr(i, j, k, URHO) - new_rho) > 1e-3 *
                // snew_arr(i, j, k, URHO))
                //     std::cout << "Reaction changing rho: diff=" << snew_arr(i, j,
                //     k, URHO) - new_rho
                //               << " @ " << i << "," << j << "," << k << '\n';  //
                //               for debug
                for (int n = 0; n < NUM_SPECIES; ++n) {
                  snew_arr(i, j, k, UFS + n) *= snew_arr(i, j, k, URHO) / new_rho;
                }
                // snew_arr(i, j, k, URHO) = new_rho;
              }

              // update drY/dt
              for (int n = 0; n < NUM_SPECIES; ++n) {
                I_R_arr(i, j, k, n) =
                  (snew_arr(i, j, k, UFS + n) - sold_arr(i, j, k, UFS + n)) / dt -
                  rYsrc(i, j, k, n);
              }

              // update heat release rate
              if (update_heat_release) {
                Real Y[NUM_SPECIES];
                for (int n = 0; n < NUM_SPECIES; n++) {
                  Y[n] = snew_arr(i, j, k, UFS + n) / snew_arr(i, j, k, URHO);
                }
                Real hi[NUM_SPECIES];
                auto eos = pele::physics::PhysicsType::eos();
                eos.RTY2Hi(snew_arr(i, j, k, URHO), T(i, j, k), Y, hi);
                I_R_arr(i, j, k, NUM_SPECIES) = 0.0;
                for (int n = 0; n < NUM_SPECIES; n++) {
                  I_R_arr(i, j, k, NUM_SPECIES) -= hi[n] * I_R_arr(i, j, k, n);
                }
              }

              // average I_R and put into the front NREACT entries
              if (NUM_FIELD > 0) {
                for (int n = 0; n < NREACT; ++n) {
                  I_R_mean_arr(i, j, k, n) +=
                    I_R_arr(i, j, k, n) / amrex::Real(NUM_FIELD);
                }
              }
              // }
            }
          });
        } // end fields loop
      }   // end EB not covered block

      // Record runtime for load balancing
      amrex::Gpu::streamSynchronize();
      wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();
      get_new_data(Cost_Type)[mfi].plus<amrex::RunOn::Device>(wt, bx);

    } // end mfi loop
  }   // end omp block

  if (Snew.nGrow() > 0) { Snew.FillBoundary(geom.periodicity()); }

  if (verbose >= 2) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    const int NProcs = amrex::ParallelDescriptor::NProcs();

    Real max_runtime = amrex::ParallelDescriptor::second() - strt_time;
    Real min_runtime = max_runtime;
    Real sum_runtime = max_runtime;

    amrex::ParallelDescriptor::ReduceRealMax(max_runtime, IOProc);
    amrex::ParallelDescriptor::ReduceRealMin(min_runtime, IOProc);
    amrex::ParallelDescriptor::ReduceRealSum(sum_runtime, IOProc);

    Real max_fc = fctCount.sum(0, true);
    Real min_fc = max_fc;
    Real sum_fc = max_fc;

    amrex::ParallelDescriptor::ReduceRealMax(max_fc, IOProc);
    amrex::ParallelDescriptor::ReduceRealMin(min_fc, IOProc);
    amrex::ParallelDescriptor::ReduceRealSum(sum_fc, IOProc);

    amrex::Print() << "    Runtime = [" << min_runtime << "..."
                   << sum_runtime / NProcs << "..." << max_runtime << "]\n";
    amrex::Print() << "    # RHS eval = [" << min_fc << "..." << sum_fc / NProcs
                   << "..." << max_fc << "]\n";
  }
}