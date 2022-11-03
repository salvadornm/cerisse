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
cns_estdt(int i, int j, int k, amrex::Array4<Real const> const& state,
          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
          Parm const& parm) noexcept
{
    using amrex::Real;

    Real rho = state(i,j,k,URHO);
    Real mx  = state(i,j,k,UMX);
    Real my  = state(i,j,k,UMY);
#if (AMREX_SPACEDIM == 3)
    Real mz  = state(i,j,k,UMZ);
#endif
    Real ei  = state(i,j,k,UEINT);
    Real rhoinv = Real(1.0)/amrex::max(rho,parm.smallr);
    Real vx = mx*rhoinv;
    Real vy = my*rhoinv;
#if (AMREX_SPACEDIM == 3)
    Real vz = mz*rhoinv;
#endif
    Real p = amrex::max((parm.eos_gamma-Real(1.0))*ei, parm.smallp);
    Real cs = std::sqrt(parm.eos_gamma*p*rhoinv);
    Real dtx = dx[0]/(amrex::Math::abs(vx)+cs);
    Real dty = dx[1]/(amrex::Math::abs(vy)+cs);
#if (AMREX_SPACEDIM == 3)
    Real dtz = dx[2]/(amrex::Math::abs(vz)+cs);
#endif

#if (AMREX_SPACEDIM == 2)
    return amrex::min(dtx,dty);
#else
    return amrex::min(dtx,amrex::min(dty,dtz));
#endif
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_compute_temperature(int i, int j, int k, 
                        amrex::Array4<amrex::Real> const& S,
                        Parm const& parm) noexcept
{
    // Reset Eint and compute T
    
    amrex::Real rhoinv = amrex::Real(1.0)/S(i,j,k,URHO);
    amrex::Real mx = S(i,j,k,UMX);
    amrex::Real my = S(i,j,k,UMY);
#if (AMREX_SPACEDIM == 2)
    S(i,j,k,UEINT) = S(i,j,k,UEDEN) - amrex::Real(0.5) * rhoinv * (mx*mx + my*my);
#else
    amrex::Real mz = S(i,j,k,UMZ);
    S(i,j,k,UEINT) = S(i,j,k,UEDEN) - amrex::Real(0.5) * rhoinv * (mx*mx + my*my + mz*mz);
#endif

    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = S(i, j, k, UFS + n) * rhoinv;
    }

    auto eos = pele::physics::PhysicsType::eos();
    eos.EY2T(rhoinv * S(i,j,k,UEINT), massfrac, S(i,j,k,UTEMP));

    // amrex::Real cv;
    // eos.TY2Cv(S(i,j,k,UTEMP), massfrac, cv);
    // S(i,j,k,UTEMP) = rhoinv * S(i,j,k,UEINT) * amrex::Real(1.0) / cv;
    // S(i,j,k,UTEMP) = rhoinv * S(i,j,k,UEINT) * amrex::Real(1.0) / parm.cv;

    // S(i,j,k,UTEMP) = 1500.; //for debug
}

#endif
