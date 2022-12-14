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
CNS::compute_pdf_model (amrex::MultiFab& S, amrex::Real dt)
{
    BL_PROFILE("CNS::compute_pdf_model()");

    const auto* dx = geom.CellSize();

    // Prepare Wiener process for vpdf. It is constant over space.
    amrex::Real dW[NUM_FIELD][AMREX_SPACEDIM]; 
    if (do_vpdf) {
        amrex::Real sqrtdt = sqrt(dt); //sqrt because it is standard deviation
        for (int nf = 0; nf < NUM_FIELD;      ++nf) {
        for (int  i = 0;  i < AMREX_SPACEDIM; ++i ) {
            dW[nf][i] = amrex::RandomNormal(0.0, sqrtdt);
        }
        }
    }

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
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

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
        const auto& flag = flags[mfi];
        if (flag.getType(bx) != amrex::FabType::covered) // Just to save some computation
#endif
        {
            amrex::Array4<Real> sarr = S.array(mfi);

            if (do_vpdf) {
                vpdf_langevin_model(bx, sarr, dt, dW, dx);
            }
            if (do_spdf) {
                spdf_iem_model(bx, sarr, dt);
            }
        }

        amrex::Gpu::streamSynchronize();

        wt = (amrex::second() - wt) / bx.d_numPts();
        cost[mfi].plus<amrex::RunOn::Device>(wt, bx);
    } // mfi loop
} // omp parallel
}