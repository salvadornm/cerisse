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

/** 
 * \brief Compute next dt size.
*/
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

#endif
