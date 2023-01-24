#include <AMReX_FArrayBox.H>

// #include "index_macros.H"
#include "CNS.H"
// #include "PelePhysics.H"

void
CNS::set_typical_values_chem ()
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

// /** \brief Compute I_R and update S_new += dt*I_R
//  */
// void
// CNS::react_state (amrex::Real /*time*/,
//                   amrex::Real dt,
//                   bool init_react)
// {
//   BL_PROFILE("CNS::react_state()");

// //   const amrex::Real strt_time = amrex::ParallelDescriptor::second();

//   // if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
//   //   if (init_react) {
//   //     amrex::Print() << "Initialising reactions, using interval dt = " << dt << std::endl;
//   //   } else {
//   //     amrex::Print() << " >> Computing reactions" << std::endl;
//   //   }
//   // }

//   // Create dSdt of all non-reacting source terms (if dSdt does not exist)  
//   // amrex::MultiFab non_react_src(grids, dmap, NUM_SPECIES, 0, amrex::MFInfo(), Factory());
//   // non_react_src.setVal(0.0);
//   // if (react_init) {
//   //   dSdt = &non_react_src;
//   // }

//   const amrex::MultiFab& S_old = get_old_data(State_Type);
//   amrex::MultiFab& S_new = get_new_data(State_Type);
//   const int ng = S_new.nGrow();

//   amrex::MultiFab& I_R = get_new_data(Reactions_Type);
//   I_R.setVal(0.0);

//   // Fab for all Array4 in sundials box integration [rhoY; rhoYdit; T; rhoEi; rhoEidot]
//   amrex::MultiFab tempFab(grids, dmap, 3, 0); // this stores T, rEi, rEidot
//   amrex::MultiFab rYFab(grids, dmap, NUM_SPECIES, 0);
//   amrex::MultiFab rYdotFab(grids, dmap, NUM_SPECIES, 0);
//   amrex::iMultiFab maskFab(grids, dmap, 1, 0); //maskfab.setVal(1); // can use this if we use cvode
//   amrex::MultiFab fctCount(grids, dmap, 1, 0); // number of steps output

//   // Fill in rhoY and rhoYdot
//   // amrex::MultiFab::Copy(rYFab, S_old, UFS, 0, NUM_SPECIES, 0);
//   // amrex::MultiFab::Copy(rYdotFab, *dSdt, UFS, 0, NUM_SPECIES, 0);

// #if (AMREX_SPACEDIM > 1) //1D cannot have EB
//   auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
//   auto const& flags = fact.getMultiEBCellFlagFab();
// #endif
// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
// #endif
// {
//   for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//     // const amrex::Box& bx = mfi.growntilebox(ng);
//     const amrex::Box& bx = mfi.tilebox();
//     amrex::Real wt = amrex::ParallelDescriptor::second(); // timing for each fab
//     amrex::Real current_time = 0.0;
//     amrex::AllPrint() << "For box " << bx << std::endl;
    
// #if (AMREX_SPACEDIM > 1) //1D cannot have EB
//     const auto& flag_fab = flags[mfi];
//     amrex::FabType typ = flag_fab.getType(bx);
//     if (typ == amrex::FabType::covered) { 
//       continue; 
//     } 
//     else if ((typ == amrex::FabType::singlevalued) || (typ == amrex::FabType::regular)) 
// #endif
//     { 
//       auto const& rhoY    = rYFab.array(mfi);
//       auto const& rhoYsrc = rYdotFab.array(mfi);
//       auto const& T        = tempFab.array(mfi, 0);
//       auto const& rhoEi    = tempFab.array(mfi, 1);
//       auto const& rhoEisrc = tempFab.array(mfi, 2);
//       auto const& mask = maskFab.array(mfi);
//       auto const& fc   = fctCount.array(mfi);
//       auto const& sold_arr = init_react ? S_new.array(mfi) : S_old.array(mfi); //not yet initialised sold
//       auto const& snew_arr = S_new.array(mfi);
//       auto const& I_R_arr  = I_R.array(mfi);

//       // amrex::Print() << "after arr" << std::endl;

//       // Calculate Y, rhoYdot, T, rhoe, rhoedot
//       amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//         amrex::Real rho, AMREX_D_DECL(rhou, rhov, rhow), rhoinv;
//         amrex::Real Y[NUM_SPECIES];
//         const amrex::Real dtinv = 1.0/dt;

//         // access old states
//         rho = sold_arr(i, j, k, URHO);
//         AMREX_D_TERM(rhou = sold_arr(i, j, k, UMX); ,
//                      rhov = sold_arr(i, j, k, UMY); ,
//                      rhow = sold_arr(i, j, k, UMZ);)
//         rhoinv = 1.0 / rho;

//         // amrex::Print() << "mkpt 1 ("<< i << "," << j << "," << k << ")" << std::endl;

//         rhoEi(i,j,k) = sold_arr(i,j,k,UEDEN) - 0.5*(AMREX_D_TERM(rhou*rhou, + rhov*rhov, + rhow*rhow))*rhoinv;
                
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//           rhoY(i,j,k,n) = sold_arr(i, j, k, UFS+n);
//           Y[n] = sold_arr(i, j, k, UFS+n) * rhoinv;
//         }

//         auto eos = pele::physics::PhysicsType::eos();
//         eos.EY2T(rhoEi(i,j,k)*rhoinv, Y, T(i,j,k));        
        
//         // access new states
//         rho = snew_arr(i, j, k, URHO);
//         AMREX_D_TERM(rhou = snew_arr(i, j, k, UMX); ,
//                      rhov = snew_arr(i, j, k, UMY); ,
//                      rhow = snew_arr(i, j, k, UMZ);)
//         rhoinv = 1.0 / rho;

//         // calculate source terms
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//           rhoYsrc(i,j,k,n) = (snew_arr(i, j, k, UFS+n) - rhoY(i,j,k,n)) * dtinv;
//         }
        
//         rhoEisrc(i,j,k) = (snew_arr(i,j,k,UEDEN) - 0.5*(AMREX_D_TERM(rhou*rhou, + rhov*rhov, + rhow*rhow))*rhoinv - rhoEi(i,j,k)) * dtinv;
        
//         mask(i,j,k) = (T(i,j,k) > min_react_temp) ? 1 : -1;
//       });

//       // amrex::Print() << "before react()" << std::endl;

//       // PelePhysics react
//       reactor->react(bx, rhoY, rhoYsrc, T, rhoEi, rhoEisrc, fc, mask, dt, current_time
// #ifdef AMREX_USE_GPU
//         , amrex::Gpu::gpuStream()
// #endif
//       );
//       amrex::Gpu::Device::streamSynchronize();

//       // Unpack data
//       // Update U^{n+1} = U^** + dt*I_R^{n+1}
//       if (!init_react) {
//         amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//           if (mask(i,j,k) == 1) {
//             snew_arr(i, j, k, URHO) = 0.0;
//             for (int n = 0; n < NUM_SPECIES; ++n) {
//               snew_arr(i, j, k, URHO) += rhoY(i, j, k, n); 
//               snew_arr(i, j, k, UFS+n) = rhoY(i, j, k, n);
//             }
//           }
//         });
//       }

//       // amrex::Print() << "Before fill in I_R" << std::endl;

//       // Fill I_R
//       amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//         if (mask(i,j,k) == 1) {
//           const amrex::Real dtinv = 1.0/dt;

//           // update drhoY/dt
//           // amrex::Print() << "before drhoY/dt" << std::endl;
//           for (int n = 0; n < NUM_SPECIES; ++n) {
//             I_R_arr(i,j,k,n) = (rhoY(i,j,k,n) - sold_arr(i,j,k, UFS+n)) * dtinv - rhoYsrc(i,j,k,n);
//           }

//           // amrex::Print() << "after drhoY/dt" << std::endl;

//           // update heat release
//           I_R_arr(i, j, k, NUM_SPECIES) = 0.0;
//           if (!init_react) {
//             amrex::Real hi[NUM_SPECIES] = {0.0};
//             amrex::Real Y[NUM_SPECIES] = {0.0};
//             for (int n = 0; n < NUM_SPECIES; n++) {
//               Y[n] = snew_arr(i, j, k, UFS+n) / snew_arr(i, j, k, URHO);
//             }

//             auto eos = pele::physics::PhysicsType::eos();
//             eos.RTY2Hi(snew_arr(i, j, k, URHO), T(i, j, k), Y, hi);

//             // amrex::Print() << "after eos" << std::endl;

//             I_R_arr(i, j, k, NUM_SPECIES) = 0.0;
//             for (int n = 0; n < NUM_SPECIES; n++) {
//               I_R_arr(i, j, k, NUM_SPECIES) -= hi[n] * I_R_arr(i, j, k, n);
//             }
//           }
//         }
//       });

//     } // not covered

//     // Record runtime for load balancing
//     wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();
//     const amrex::Box vbox = mfi.tilebox();
//     get_new_data(Cost_Type)[mfi].plus<amrex::RunOn::Device>(wt, vbox);
//   } // mfi loop
// } //omp

/** \brief Compute I_R and update S_new += dt*I_R
 */
void
CNS::react_state (amrex::Real /*time*/,
                  amrex::Real dt,
                  bool init_react)
{

  BL_PROFILE("CNS::react_state()");
  const amrex::Real strt_time = amrex::ParallelDescriptor::second();

  if ((verbose > 0) && amrex::ParallelDescriptor::IOProcessor()) {
    if (init_react) {
      amrex::Print() << "Initialising reactions, using interval dt = " << dt << std::endl;
    } else {
      amrex::Print() << " >> Computing reactions" << std::endl;
    }
  }

  // State Fabs
  const amrex::MultiFab& Sold = get_old_data(State_Type);
  amrex::MultiFab& Snew = get_new_data(State_Type);
  amrex::MultiFab& I_R = get_new_data(Reactions_Type);
  I_R.setVal(0.0);

  // Fab for all Array4 in sundials box integration
  amrex::MultiFab  STemp(grids, dmap, NUM_SPECIES + 2, 0); //[rY, rEi, T]
  amrex::MultiFab  SDotTemp(grids, dmap, NUM_SPECIES + 1, 0); //d[rY, rEi]/dt
  amrex::iMultiFab maskFab(grids, dmap, 1, 0); //= 1: do reaction, = -1: don't do reaction
  amrex::MultiFab  fctCount(grids, dmap, 1, 0); //number of steps output

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
  auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
{
  for (amrex::MFIter mfi(I_R, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.tilebox();
    amrex::Real wt = amrex::ParallelDescriptor::second(); // timing for each fab
    
#if (AMREX_SPACEDIM > 1) //1D cannot have EB
    const auto& flag_fab = flags[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if ((typ == amrex::FabType::singlevalued) || (typ == amrex::FabType::regular)) 
#endif
    {
      // ============== Prepare for react ==============
      const amrex::Array4 sold_arr = init_react ? Snew.array(mfi) : Sold.array(mfi); //don't have Sold at initialisation
      amrex::Array4 snew_arr = Snew.array(mfi);
      amrex::Array4 I_R_arr = I_R.array(mfi);
      amrex::Array4 rY = STemp.array(mfi, 0);
      amrex::Array4 rEi = STemp.array(mfi, NUM_SPECIES);
      amrex::Array4 T = STemp.array(mfi, NUM_SPECIES+1);
      amrex::Array4 rYsrc = SDotTemp.array(mfi, 0);
      amrex::Array4 rEisrc = SDotTemp.array(mfi, NUM_SPECIES);
      amrex::Array4 mask = maskFab.array(mfi);
      amrex::Array4 fc = fctCount.array(mfi);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // calculate [rY, rEi, T]
        const amrex::Real rho_old = sold_arr(i,j,k,URHO);
        amrex::Real Y[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; ++n) {
          rY(i,j,k,n) = sold_arr(i,j,k,UFS+n);
          Y[n] = sold_arr(i,j,k,UFS+n) / rho_old;
        }

        rEi(i,j,k) = sold_arr(i,j,k,UEDEN) 
            - 0.5*(AMREX_D_TERM(sold_arr(i,j,k,UMX)*sold_arr(i,j,k,UMX),
                               +sold_arr(i,j,k,UMY)*sold_arr(i,j,k,UMY), 
                               +sold_arr(i,j,k,UMZ)*sold_arr(i,j,k,UMZ))
                  ) / rho_old;
        
        amrex::Real ei = rEi(i,j,k) / rho_old;
        auto eos = pele::physics::PhysicsType::eos();
        eos.EY2T(ei, Y, T(i,j,k));

        // calculate d[rY, rEi]/dt
        for (int n = 0; n < NUM_SPECIES; ++n) {
          rYsrc(i,j,k,n) = (snew_arr(i,j,k,UFS+n) - rY(i,j,k,n)) / dt;
        }

        rEisrc(i,j,k) = (snew_arr(i,j,k,UEDEN) 
            - 0.5*(AMREX_D_TERM(snew_arr(i,j,k,UMX)*snew_arr(i,j,k,UMX),
                               +snew_arr(i,j,k,UMY)*snew_arr(i,j,k,UMY), 
                               +snew_arr(i,j,k,UMZ)*snew_arr(i,j,k,UMZ)))
            - rEi(i,j,k)) / dt;
        
        // fill mask
        mask(i,j,k) = (T(i,j,k) > min_react_temp) ? 1 : -1;
      });

      // ===================== React =====================
      // // My dummy reactor
      // amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      //   if (mask(i,j,k) != -1) {
      //     amrex::Real rho = 0.0;
      //     for (int n = 0; n < NUM_SPECIES; ++n) {        
      //       rY(i,j,k,n) += dt*rYsrc(i,j,k,n);        
      //       rho += rY(i,j,k,n);
      //     }
      //     rY(i,j,k,N2_ID) += -1e3*rho*dt;
      //     rY(i,j,k,N_ID) += 1e3*rho*dt;      
      //   }
      // });
      amrex::Real current_time = 0.0;
      reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time
#ifdef AMREX_USE_GPU
      , amrex::Gpu::gpuStream()
#endif
      );
      amrex::Gpu::Device::streamSynchronize();

      // ================== Unpack data ==================
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (mask(i,j,k) != -1) {
          // update U^{n+1} = U^** + dt*I_R^{n+1}
          if (!init_react) {
            snew_arr(i,j,k,URHO) = 0.0;
            for (int n = 0; n < NUM_SPECIES; ++n) {
              snew_arr(i,j,k,URHO) += rY(i,j,k,n);
              snew_arr(i,j,k,UFS+n) = rY(i,j,k,n);
            }
          }

          // update drY/dt      
          for (int n = 0; n < NUM_SPECIES; ++n) {
            I_R_arr(i,j,k,n) = (rY(i,j,k,n)-sold_arr(i,j,k,UFS+n)) / dt
                              - rYsrc(i,j,k,n);
          }

          // update heat release rate
          if (update_heat_release) {
            amrex::Real Y[NUM_SPECIES];
            for (int n = 0; n < NUM_SPECIES; n++) {
              Y[n] = snew_arr(i,j,k,UFS+n) / snew_arr(i,j,k,URHO);
            }
            amrex::Real hi[NUM_SPECIES];
            auto eos = pele::physics::PhysicsType::eos();
            eos.RTY2Hi(snew_arr(i,j,k,URHO), T(i,j,k), Y, hi);
            I_R_arr(i,j,k,NUM_SPECIES) = 0.0;
            for (int n = 0; n < NUM_SPECIES; n++) {
              I_R_arr(i,j,k,NUM_SPECIES) -= hi[n] * I_R_arr(i,j,k,n);
            }
          }
        }
      });
    } //end EB not covered block

    // Record runtime for load balancing
    amrex::Gpu::streamSynchronize();
    wt = (amrex::ParallelDescriptor::second() - wt) / bx.d_numPts();
    get_new_data(Cost_Type)[mfi].plus<amrex::RunOn::Device>(wt, bx);

  } //end mfi loop
} //end omp block

  if (Snew.nGrow() > 0) {
    Snew.FillBoundary(geom.periodicity());
  }

  if (verbose >= 2) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real max_runtime = amrex::ParallelDescriptor::second() - strt_time;
    amrex::Real min_runtime = max_runtime;

    amrex::ParallelDescriptor::ReduceRealMax(max_runtime, IOProc);
    amrex::ParallelDescriptor::ReduceRealMin(min_runtime, IOProc);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "    Runtime = [" << min_runtime << 
                                    "..." << max_runtime << "]\n";
    }
  }
}