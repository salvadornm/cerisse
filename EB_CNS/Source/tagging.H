#ifndef CNS_TAGGING_H
#define CNS_TAGGING_H

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_tag_grad (int i, int j, int k,
              amrex::Array4<char> const& tag,
              amrex::Array4<amrex::Real const> const& s,
              amrex::Real grad_threshold, char tagval) noexcept
{
    amrex::Real ax, ay, az = 0.0;
    AMREX_D_TERM(
        ax = amrex::Math::abs(s(i+1,j,k) - s(i,j,k)); 
        ax = amrex::max(ax, amrex::Math::abs(s(i,j,k) - s(i-1,j,k))); 
        ,
        ay = amrex::Math::abs(s(i,j+1,k) - s(i,j,k)); 
        ay = amrex::max(ay, amrex::Math::abs(s(i,j,k) - s(i,j-1,k)));
        ,
        az = amrex::Math::abs(s(i,j,k+1) - s(i,j,k)); 
        az = amrex::max(az, amrex::Math::abs(s(i,j,k) - s(i,j,k-1))););

    if (amrex::max<amrex::Real>(ax, ay, az) >= grad_threshold) {
        tag(i, j, k) = tagval;
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_tag_error (int i, int j, int k,
               amrex::Array4<char> const& tag,
               amrex::Array4<amrex::Real const> const& s,
               amrex::Real err_threshold, char tagval) noexcept
{    
    if (s(i, j, k) >= err_threshold) {
        tag(i, j, k) = tagval;
    }
}

#endif
