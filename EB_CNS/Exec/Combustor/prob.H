#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
inline
void
cns_initdata(int i, int j, int k, 
             amrex::Array4<amrex::Real> const& state,
             amrex::GeometryData const& geomdata,
             Parm const& parm,
             ProbParm const& prob_parm)
{
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real rho, eint;
    eos.PYT2RE(prob_parm.interior_p, (amrex::Real*) prob_parm.massfrac.begin(), prob_parm.interior_T,   
               rho, eint);
    
    state(i,j,k,URHO) = rho;

    state(i,j,k,UMX) = amrex::Real(0.0);
    state(i,j,k,UMY) = amrex::Real(0.0);
    state(i,j,k,UMZ) = amrex::Real(0.0);

    state(i,j,k,UTEMP) = prob_parm.interior_T;
    state(i,j,k,UEINT) = rho * eint;
    state(i,j,k,UEDEN) = rho * eint;

    for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * prob_parm.massfrac[n];
    }
}

#endif
