#include "pdf_model.H"

/** 
 * \brief Average field data ans write to S
*/
void
CNS::computeAvg (MultiFab& S)
{
  BL_PROFILE("CNS::computeAvg()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    auto sarr = S.array(mfi);

    amrex::ParallelFor(bx, NVAR,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      sarr(i, j, k, n) = 0.0;
      for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        sarr(i, j, k, n) += sarr(i, j, k, nf*NVAR + n);
      }
      sarr(i, j, k, n) /= amrex::Real(NUM_FIELD);
    });

    amrex::Gpu::synchronize();
  } // mfi loop
} // omp parallel
}

/** 
 * \brief Compute modelled terms for Velocity-PDF (Langevin model) / Species-PDF (IEM model)
*/
void
CNS::compute_pdf_model (amrex::MultiFab& S, amrex::Real dt, int iteration)
{
  BL_PROFILE("CNS::compute_pdf_model()");

  if (verbose > 0) {
    amrex::Print() << " >> Computing PDF source term:";
    if (do_psgs) { amrex::Print() << " p_sgs"; }
    if (do_vpdf) { amrex::Print() << " VPDF"; }
    if (do_spdf) { amrex::Print() << " SPDF"; } 
    amrex::Print() << std::endl;
  }

  const auto* dx = geom.CellSize();
  // amrex::Real dx0[AMREX_SPACEDIM];
  // for (int i = 0; i < AMREX_SPACEDIM; ++i) {
  //   dx0[i] = dx[i]*pow(2, amrex::max(level-2,0));
  // }

  // Prepare Wiener process for vpdf. It is constant over space.
  amrex::Real dW[AMREX_SPACEDIM][NUM_FIELD]; 
  WienerProcess.generate_new(sqrt(dt), 2, 0);
  WienerProcess.get_rand(nStep(), dW);

  /////////////////////////////////////////////////
  // TO TEST UNIQUERAND: TURN OFF ALL PDF MODELS //
  /////////////////////////////////////////////////
  // WienerProcess.generate_new(sqrt(dt), 2, 1);
  // WienerProcess.get_rand(nStep(), dW);

  // int f = 3;
  // amrex::Print() << "** Level " << level << " cycle " << iteration-1
  //                << ": dW[0][" << f << "] = " << dW[0][f] << std::endl;

  // for (amrex::MFIter mfi(S, true); mfi.isValid(); ++mfi) {
  //   const amrex::Box& bx = mfi.tilebox();
  //   amrex::Array4<Real> sarr = S.array(mfi);
  //   amrex::ParallelFor(bx, NUM_FIELD, [=] AMREX_GPU_DEVICE (int i, int j, int k, int nf) noexcept {
  //     for (int d = 0; d < AMREX_SPACEDIM; ++d) {
  //       sarr(i,j,k,(nf+1)*NVAR+UMX+d) += 10*dW[d][nf]; // s.t. var(u) = 10*t
  //     }
  //   });
  // }
  /////////////////////////////////////////////////
  /////////////////////////////////////////////////

  // This has been moved to random.H
  // if (do_vpdf) {
  //   if (iteration == 1 && this->level == 0) {
  //     WienerProcess.generate_new(sqrt(dt), this->level, 0);
  //   }
  //   WienerProcess.get(this->level, iteration-1, dW);
    
  //   if (amrex::ParallelDescriptor::IOProcessor()) {
  //     amrex::Real sqrtdt = sqrt(dt); //sqrt because it is standard deviation

  //     amrex::Real corr_fac = (NUM_FIELD % 2 == 0) ? 1.0 : sqrt((amrex::Real)(NUM_FIELD)/(amrex::Real)(NUM_FIELD - 1)); //to give exactly variance=1      
  //     for (int d = 0; d < AMREX_SPACEDIM; ++d) {
  //       for (int nf = 0; nf < NUM_FIELD-1; nf+=2) {
  //         dW[nf][d] = sqrtdt * corr_fac * (2.0 * (double)amrex::Random_int(2) - 1.0);
  //         dW[nf+1][d] = -dW[nf][d];
  //       }
  //       if (NUM_FIELD % 2 == 1) {
  //         int swap_id = amrex::Random_int(NUM_FIELD-2);
  //         dW[NUM_FIELD-1][d] = dW[swap_id][d];
  //         dW[swap_id][d] = 0.0;
  //       }
  //     }

  //     amrex::Real mean;
  //     amrex::Real stdd;
  //     for (int d = 0; d < AMREX_SPACEDIM; ++d) {
  //       mean = 0.0;
  //       stdd = 0.0;
  //       for (int nf = 0; nf < NUM_FIELD; ++nf) {
  //         dW[nf][d] = amrex::RandomNormal(0.0, sqrtdt);
  //         mean += dW[nf][d];
  //         stdd += dW[nf][d]*dW[nf][d];
  //       }
  //       mean /= NUM_FIELD;
  //       stdd = sqrt(stdd / NUM_FIELD - mean*mean);
  //       for (int nf = 0; nf < NUM_FIELD; ++nf) {
  //         dW[nf][d] = (dW[nf][d] - mean) * sqrtdt / stdd;
  //       }
  //     }
  //   }

  //   amrex::ParallelDescriptor::Bcast(&dW[0][0], NUM_FIELD*AMREX_SPACEDIM, 
  //                                    amrex::ParallelDescriptor::IOProcessorNumber());
  // }

#if CNS_USE_EB
  auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

  amrex::MultiFab& cost = get_new_data(Cost_Type);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{    
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    auto wt = amrex::second();

    const amrex::Box& bx = mfi.tilebox();

#if CNS_USE_EB
    const auto& flag = flags[mfi];
    if (flag.getType(bx) != amrex::FabType::covered) // Just to save some computation
#endif
    {
      amrex::Array4<Real> sarr = S.array(mfi);

      // The order of Langevin and IEM does not matter, as they act on different variables
      // The order of p_sgs and Langevin may matter, do not know
      if (do_psgs) {
        vpdf_psgs_model(bx, sarr, dt, dW, dx);
      }

      if (do_vpdf) {
        vpdf_langevin_model(bx, sarr, dt, dW, dx);
      }
      
      if (do_spdf) {
        spdf_iem_model(bx, sarr, dt, dW, dx);
      }
    }

    amrex::Gpu::streamSynchronize();

    wt = (amrex::second() - wt) / bx.d_numPts();
    cost[mfi].plus<amrex::RunOn::Device>(wt, bx);
  } // mfi loop
} // omp parallel
}