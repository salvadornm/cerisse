#ifndef CNS_K_H_
#define CNS_K_H_

#include <AMReX_FArrayBox.H>
#include <limits>
#include <cmath>

// #include "CNS.H"
#include "index_macros.H"
#include "PelePhysics.H"
// #include "prob.H"
// #include "parm.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
cns_estdt (int i, int j, int k, amrex::Array4<Real const> const& state,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           Parm const& parm) noexcept
{
    amrex::Real rho = state(i, j, k, URHO);    
    amrex::Real rhoinv = Real(1.0) / rho;
    AMREX_D_TERM(amrex::Real vx = state(i, j, k, UMX) * rhoinv; ,
                 amrex::Real vy = state(i, j, k, UMY) * rhoinv; ,
                 amrex::Real vz = state(i, j, k, UMZ) * rhoinv;);
    amrex::Real ei = state(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vx*vx, +vy*vy, +vz*vz));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = state(i, j, k, UFS + n) * rhoinv;
    }

    amrex::Real cs, T;
    auto eos = pele::physics::PhysicsType::eos();
    eos.EY2T(ei, massfrac, T);
    eos.RTY2Cs(rho, T, massfrac, cs);

    AMREX_D_TERM(
        amrex::Real dt  = dx[0]/(amrex::Math::abs(vx)+cs); 
        ,
        amrex::Real dty = dx[1]/(amrex::Math::abs(vy)+cs);
        dt = amrex::min(dt,dty);
        ,
        amrex::Real dtz = dx[2]/(amrex::Math::abs(vz)+cs);
        dt = amrex::min(dt,dtz););  
    
    return dt;
}


// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_compute_temperature (int i, int j, int k, 
//                          amrex::Array4<amrex::Real> const& S) noexcept
// {
//     // Reset Eint and compute T
    
//     amrex::Real rhoinv = amrex::Real(1.0) / S(i,j,k,URHO);
//     amrex::Real mx = S(i,j,k,UMX);
//     amrex::Real my = S(i,j,k,UMY);
//     amrex::Real mz = S(i,j,k,UMZ);

//     S(i,j,k,UEINT) = S(i,j,k,UEDEN) - amrex::Real(0.5) * rhoinv * (AMREX_D_TERM(mx*mx, + my*my, + mz*mz));

//     amrex::Real massfrac[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         massfrac[n] = S(i, j, k, UFS + n) * rhoinv;
//     }

//     auto eos = pele::physics::PhysicsType::eos();
//     eos.EY2T(rhoinv*S(i,j,k,UEINT), massfrac, S(i,j,k,UTEMP));
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
// void
// cns_compute_average (int i, int j, int k,
//                      amrex::Array4<amrex::Real> const& S) noexcept
// {
//     for (int n = 0; n < NVAR; ++n) {
//         S(i, j, k, n) = 0.0;
//         for (int nf = 1; nf <= NUM_FIELD; ++nf) {
//             S(i, j, k, n) += S(i, j, k, nf*NVAR + n);
//         }
//         S(i, j, k, n) /= amrex::Real(NUM_FIELD);
//     }
// }

// void
// CNS::computeAvg (MultiFab& State)
// {
//     BL_PROFILE("CNS::computeAvg()");

// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     auto const& sarrs = State.arrays();    
//     amrex::ParallelFor(State, 
//     [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
//         cns_compute_average(i, j, k, sarrs[nbx]);            
//     });

//     amrex::Gpu::synchronize();
// }

#endif
