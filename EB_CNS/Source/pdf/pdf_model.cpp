#include "pdf_model.H"

/**
 * \brief Average field data and write to S
 */
void CNS::computeAvg(MultiFab& S)
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
                             sarr(i, j, k, n) += sarr(i, j, k, nf * NVAR + n);
                           }
                           sarr(i, j, k, n) /= amrex::Real(NUM_FIELD);
                         });

      amrex::Gpu::synchronize();
    } // mfi loop
  }   // omp parallel
}

/**
 * \brief Compute modelled terms for Velocity-PDF (Langevin model) / Species-PDF (IEM
 * model)
 */
void CNS::compute_pdf_model(amrex::MultiFab& S, amrex::Real dt, int iteration)
{
  BL_PROFILE("CNS::compute_pdf_model()");

  if (verbose > 0) {
    amrex::Print() << " >> Computing PDF source term:";
    if (do_psgs) amrex::Print() << " p_sgs";
    if (do_vpdf) amrex::Print() << " VPDF";
    if (do_spdf) amrex::Print() << " SPDF";
    amrex::Print() << std::endl;
  }

  const auto* dx = geom.CellSize();

  // Prepare Wiener process for vpdf. It is constant over space.
  amrex::Real dW[AMREX_SPACEDIM][NUM_FIELD];
  WienerProcess.generate_new(sqrt(dt), 0, 0);
  WienerProcess.get_rand(nStep(), dW);

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
      if (flag.getType(bx) !=
          amrex::FabType::covered) // Just to save some computation
#endif
      {
        amrex::Array4<Real> sarr = S.array(mfi);

        // The order of Langevin and IEM does not matter, as they act on different
        // variables The order of p_sgs and Langevin may matter, do not know
        if (do_psgs) { vpdf_psgs_model(bx, sarr, dt, dW, dx); }

        if (do_vpdf) { vpdf_langevin_model(bx, sarr, dt, dW, dx); }

        if (do_spdf) { spdf_iem_model(bx, sarr, dt, dW, dx); }
      }

      amrex::Gpu::streamSynchronize();

      wt = (amrex::second() - wt) / bx.d_numPts();
      cost[mfi].plus<amrex::RunOn::Device>(wt, bx);
    } // mfi loop
  }   // omp parallel
}