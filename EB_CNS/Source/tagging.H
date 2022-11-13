#ifndef CNS_TAGGING_H
#define CNS_TAGGING_H

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_tag_graderror (int i, int j, int k,
                   amrex::Array4<char> const& tag,
                   amrex::Array4<amrex::Real const> const& s,
                   amrex::Real grad_threshold, char tagval) noexcept
{
    AMREX_D_TERM(
        amrex::Real ax = amrex::Math::abs(s(i+1,j,k) - s(i,j,k)); ,
        amrex::Real ay = amrex::Math::abs(s(i,j+1,k) - s(i,j,k)); ,
        amrex::Real az = amrex::Math::abs(s(i,j,k+1) - s(i,j,k)); )
    AMREX_D_TERM(
        ax = amrex::max(ax,amrex::Math::abs(s(i,j,k) - s(i-1,j,k))); ,
        ay = amrex::max(ay,amrex::Math::abs(s(i,j,k) - s(i,j-1,k))); ,
        az = amrex::max(az,amrex::Math::abs(s(i,j,k) - s(i,j,k-1))); )
#if AMREX_SPACEDIM >= 2
    if (amrex::max<amrex::Real>(AMREX_D_DECL(ax, ay, az)) >= grad_threshold)
#elif AMREX_SPACEDIM == 1
    if (ax >= grad_threshold)
#endif
    {
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
