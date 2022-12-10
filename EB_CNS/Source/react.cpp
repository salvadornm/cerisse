#include <AMReX_FArrayBox.H>

// #include "index_macros.H"
#include "CNS.H"
// #include "PelePhysics.H"

void
CNS::set_typical_values_chem ()
{
  amrex::MultiFab& S_new = get_new_data(State_Type);
  // amrex::Real minTemp = S_new.min(UTEMP);
  // amrex::Real maxTemp = S_new.max(UTEMP);
  amrex::Vector<amrex::Real> typical_values_chem(NUM_SPECIES + 1, 1e-10);

  for (int n = 0; n < NUM_SPECIES; n++) {
    amrex::Real rhoYs_min = S_new.min(UFS + n);
    amrex::Real rhoYs_max = S_new.max(UFS + n);
    typical_values_chem[n] = amrex::max<amrex::Real>(
      0.5 * (rhoYs_min + rhoYs_max), 1.e-10);
  }
  // typical_values_chem[NUM_SPECIES] = 0.5 * (minTemp + maxTemp);
  reactor->set_typ_vals_ode(typical_values_chem);
}

void
CNS::react_source (amrex::Real /*time*/,
                   amrex::Real dt,
                   bool react_init)
{}

/*
  // Update I_R

  BL_PROFILE("CNS::reaction_source()");

  const amrex::Real strt_time = amrex::ParallelDescriptor::second();

  AMREX_ASSERT(do_react == 1);

  if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
    if (react_init) {
      amrex::Print() << " ... Initialising reactions" << std::endl;
    } else {
      amrex::Print() << " ... Computing reactions" << std::endl;
    }
  }

  const amrex::MultiFab& S_new = get_new_data(State_Type);
  const amrex::MultiFab& S_old = get_old_data(State_Type);
//   const int ng = S_new.nGrow();
//   prefetchToDevice(S_new);

  // Create a MultiFab with all of the non-reacting source terms.
  amrex::MultiFab non_react_src(grids, dmap, LEN_STATE, ng, amrex::MFInfo(), Factory());
  non_react_src.setVal(0);
  if (!react_init) {
    // Only do this if we are not at the first step
    // Build non-reacting source term, and an S_new that does not include
    // reactions
    amrex::MultiFab::LinComb(non_react_src, 1.0, S_new, 0, -1.0, S_old, 0, 0, LEN_STATE, 0);
  }

  amrex::MultiFab& react_src = get_new_data(Reactions_Type);
  react_src.setVal(0.0);
  // prefetchToDevice(react_src);

  // for sundials box integration
  amrex::MultiFab STemp(grids, dmap, NUM_SPECIES + 2, 0);
  amrex::MultiFab extsrc_rY(grids, dmap, NUM_SPECIES, 0);
  amrex::MultiFab extsrc_rE(grids, dmap, 1, 0);
  amrex::iMultiFab dummyMask(grids, dmap, 1, 0);
  dummyMask.setVal(1);  
  amrex::MultiFab fctCount(grids, dmap, 1, 0);

  if (!react_init) {
    amrex::MultiFab::Copy(STemp, S_old, UFS, 0, NUM_SPECIES, 0);
    amrex::MultiFab::Copy(STemp, S_old, UTEMP, NUM_SPECIES, 1, 0);
    amrex::MultiFab::Copy(STemp, S_old, UEINT, NUM_SPECIES+1, 1, 0);
  } else {
    amrex::MultiFab::Copy(STemp, S_new, UFS, 0, NUM_SPECIES, 0);
    amrex::MultiFab::Copy(STemp, S_new, UTEMP, NUM_SPECIES, 1, 0);
    amrex::MultiFab::Copy(STemp, S_new, UEINT, NUM_SPECIES + 1, 1, 0);
  }
  amrex::MultiFab::Copy(extsrc_rY, non_react_src, UFS, 0, NUM_SPECIES, 0);

  auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_new,amrex::TilingIfNotGPU()); 
      mfi.isValid(); ++mfi) {

      const amrex::Box& bx = mfi.growntilebox(ng);

      // old state or the state at t=0
      auto const& sold_arr = react_init ? S_new.array(mfi) : S_old.array(mfi);
      // new state
      auto const& snew_arr = S_new.array(mfi);
      auto const& nonrs_arr = non_react_src.array(mfi);
      auto const& I_R = react_src.array(mfi);

      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(bx);
      if (typ == amrex::FabType::covered) {
        continue;
      } else if ( (typ == amrex::FabType::singlevalued) ||
                  (typ == amrex::FabType::regular) ) {
        amrex::Real wt = amrex::ParallelDescriptor::second(); // timing for each fab

        amrex::Real current_time = 0.0;

        auto const& rhoY    = STemp.array(mfi);
        auto const& T       = STemp.array(mfi, NUM_SPECIES);
        auto const& rhoE    = STemp.array(mfi, NUM_SPECIES+1);
        auto const& rhoYsrc = extsrc_rY.array(mfi);
        auto const& rhoEsrc = extsrc_rE.array(mfi);
        auto const& mask    = dummyMask.array(mfi);
        auto const& fc      = fctCount.array(mfi);

        // Calculate rhoEdot
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // work on old state
            amrex::Real rhou = sold_arr(i, j, k, UMX);
            amrex::Real rhov = sold_arr(i, j, k, UMY);
            amrex::Real rhow = sold_arr(i, j, k, UMZ);
            amrex::Real rhoinv = 1.0 / sold_arr(i, j, k, URHO);

            amrex::Real rhoe_old =
              sold_arr(i,j,k,UEDEN) - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)*rhoinv;
            // amrex::Real rhoe_old = sold_arr(i,j,k,UEINT); //should be the same
              
            // work on new state
            rhou = snew_arr(i, j, k, UMX);
            rhov = snew_arr(i, j, k, UMY);
            rhow = snew_arr(i, j, k, UMZ);
            rhoinv = 1.0 / snew_arr(i, j, k, URHO);

            rhoEsrc(i,j,k) = (snew_arr(i,j,k,UEDEN) - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)*rhoinv - rhoe_old) / dt;
            // (snew_arr(i,j,k,UEDEN) - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)*rhoinv - sold_arr(i,j,k,UEINT)) / dt;
              
            // for (int n = 0; n < NUM_SPECIES; ++n) {
            //   rhoYsrc(i,j,k,n) = (snew_arr(i,j,k,UFS+n) - sold_arr(i,j,k,UFS+n)) / dt;
            // } this is done above using copy from lincomb
          });

        reactor->react(
          bx, rhoY, rhoYsrc, T, rhoE, rhoEsrc, fc, mask, dt, current_time
#ifdef AMREX_USE_GPU
          , amrex::Gpu::gpuStream()
#endif
        );

        amrex::Gpu::Device::streamSynchronize();

        // unpack data
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // // work on old state
            // amrex::Real rhou = sold_arr(i, j, k, UMX);
            // amrex::Real rhov = sold_arr(i, j, k, UMY);
            // amrex::Real rhow = sold_arr(i, j, k, UMZ);
            // amrex::Real rho_old = sold_arr(i, j, k, URHO);
            // amrex::Real rhoInv = 1.0 / rho_old;

            // amrex::Real e_old =
            //   (sold_arr(i, j, k, UEDEN) // old total energy
            //    - 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv) // KE
            //   * rhoInv;

            // rhou = snew_arr(i, j, k, UMX);
            // rhov = snew_arr(i, j, k, UMY);
            // rhow = snew_arr(i, j, k, UMZ);
            // rhoInv = 1.0 / snew_arr(i, j, k, URHO);

            // amrex::Real rhoedot_ext =
            //   (snew_arr(i, j, k, UEDEN) // new total energy
            //    - 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) * rhoInv // KE
            //    - rho_old * e_old) // old internal energy
            //   / dt;

            // amrex::Real umnew =
            //   sold_arr(i, j, k, UMX) + dt * nonrs_arr(i, j, k, UMX);
            // amrex::Real vmnew =
            //   sold_arr(i, j, k, UMY) + dt * nonrs_arr(i, j, k, UMY);
            // amrex::Real wmnew =
            //   sold_arr(i, j, k, UMZ) + dt * nonrs_arr(i, j, k, UMZ);

            // // get new rho
            // amrex::Real rhonew = 0.0;

            // for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
            //   rhonew += rhoY(i, j, k, nsp);
            // }

            // omega_rhoY
            // Note: we subtract by s_old even though our starting state is s_new 
            // because I_R is integrated for a whole dt step
            for (int n = 0; n < NUM_SPECIES; ++n) {
              I_R(i,j,k,n) = (rhoY(i,j,k,n) - sold_arr(i,j,k,UFS+n))/dt 
                            - rhoYsrc(i,j,k,UFS+n);
            }

            // omega_rhoe
            I_R(i,j,k,NUM_SPECIES) = (rhoE - sold_arr(i,j,k,UEDEN) + 0.5*()/sold_arr(i,j,k,URHO));
            // I_R(i,j,k,NUM_SPECIES) =
            //   rho_old * e_old + dt * rhoedot_ext // new internal energy
            //    + 0.5 * (umnew * umnew + vmnew * vmnew + wmnew * wmnew) /
            //        rhonew                  // new KE
            //    - sold_arr(i, j, k, UEDEN)) // old total energy
            //     / dt -
            //   nonrs_arr(i, j, k, UEDEN);

            // heat release
            I_R(i, j, k, NUM_SPECIES + 1) = 0.0;
            auto eos = pele::physics::PhysicsType::eos();

            amrex::Real hi[NUM_SPECIES] = {0.0};

            amrex::Real Yspec[NUM_SPECIES] = {0.0};
            for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
              Yspec[nsp] = rhoY[nsp] / snew_arr(i, j, k, URHO);
            }
            eos.RTY2Hi(snew_arr(i, j, k, URHO), T(i, j, k), Yspec, hi);

            for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
              I_R(i, j, k, NUM_SPECIES + 1) -= hi[nsp] * I_R(i, j, k, nsp);
            }
          });

        wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();

//         // if (do_react_load_balance) {
//           const amrex::Box vbox = mfi.tilebox();
//           get_new_data(Cost_Type)[mfi].plus<amrex::RunOn::Device>(
//             wt, vbox);
//         // }

        // update heat release
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            I_R(i, j, k, NUM_SPECIES + 1) = 0.0;
            auto eos = pele::physics::PhysicsType::eos();

            amrex::Real hi[NUM_SPECIES] = {0.0};

            amrex::Real Yspec[NUM_SPECIES] = {0.0};
            for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
              Yspec[nsp] =
                snew_arr(i, j, k, UFS + nsp) / snew_arr(i, j, k, URHO);
            }
            eos.RTY2Hi(
              snew_arr(i, j, k, URHO), snew_arr(i, j, k, UTEMP), Yspec, hi);

            for (int nsp = 0; nsp < NUM_SPECIES; nsp++) {
              I_R(i, j, k, NUM_SPECIES + 1) -= hi[nsp] * I_R(i, j, k, nsp);
            }
          });
//       }
//     }
//   }

//   if (ng > 0) {
//     S_new.FillBoundary(geom.periodicity());
//   }

  if (verbose >= 2) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real max_runtime = amrex::ParallelDescriptor::second() - strt_time;
    amrex::Real min_runtime = max_runtime;

    amrex::ParallelDescriptor::ReduceRealMax(max_runtime, IOProc);
    amrex::ParallelDescriptor::ReduceRealMin(min_runtime, IOProc);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "CNS::react_state() runtime = [" << min_runtime << 
                                                   "..." << max_runtime << "]\n";
    }
  }
}

*/