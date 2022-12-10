#ifndef CNS_PDF_MODEL_H_
#define CNS_PDF_MODEL_H_

#include "CNS.H"
#include "index_macros.H"

/** \brief Average field data ans write to S
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

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
vpdf_langevin_model (const amrex::Box& bx,
                     amrex::Array4<amrex::Real>& s,
                     const amrex::Real dt,
                     const amrex::Real dW[NUM_FIELD][AMREX_SPACEDIM],
                     const amrex::Real dx[AMREX_SPACEDIM])
{   
    amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Compute turbulent statistics
        amrex::Real mean_u, u, tke = 0.0;        
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
            for (int nf = 1; nf <= NUM_FIELD; ++nf) {
                u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
                tke += (u - mean_u)*(u - mean_u);
            }
        }
        tke /= 2.0 * amrex::Real(NUM_FIELD); // k = 0.5*(u^2 + v^2 + w^2)
        amrex::Real deltainv = 1.0 / std::pow(AMREX_D_TERM(dx[0],*dx[1],*dx[2]), 1.0/amrex::Real(AMREX_SPACEDIM));
        amrex::Real eps = 0.53 * std::pow(k, 1.5) * deltainv;

        // u = u + 2.075*tke/eps*(mean_u - u)*dt + sqrt(2.1*eps)*dW
        for (int nf = 1; nf <= NUM_FIELD; ++nf) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
                u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
                
                u += 2.075 * tke/eps * (mean_u - u) * dt + sqrt(2.1*eps)*dW[nf-1][dir];
                
                s(i, j, k, nf*NVAR + UMX+dir) = s(i, j, k, URHO) * u;
            }
        }
    });
}

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
spdf_iem_model (const amrex::Box& bx,
                amrex::Array4<amrex::Real>& s,
                const amrex::Real dt)
{   
    // amrex::ParallelFor(bx, 
    // [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    //     // Compute turbulent statistics
    //     amrex::Real k   = 1.0
    //     amrex::Real eps = 1.0;

    //     // u = u + k/eps*(mean_y - y)*dt
    //     for (int nf = 1; nf <= NUM_FIELD; ++nf) {
    //         for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    //             mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
    //             u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
                
    //             u += 2.075 * k/eps * (mean_u - u) * dt + sqrt(2.1*eps)*dW[nf-1][dir];
                
    //             s(i, j, k, nf*NVAR + UMX+dir) = s(i, j, k, URHO) * u;
    //         }
    //     }
    // });
}

/** \brief Compute modelled terms for Velocity-PDF (Langevin model) / Species-PDF (IEM model)
 * 
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

#endif