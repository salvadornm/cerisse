#ifndef CNS_TAGGING_H
#define CNS_TAGGING_H

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>

#include <cmath>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_tag_grad(int i, int j, int k, amrex::Array4<char> const& tag,
             amrex::Array4<const amrex::Real> const& s,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
             amrex::Real grad_threshold, char tagval) noexcept
{
  amrex::Real ax, ay, az = 0.0;
  AMREX_D_TERM(ax = amrex::Math::abs(s(i + 1, j, k) - s(i - 1, j, k)) * dxinv[0];
               , ay = amrex::Math::abs(s(i, j + 1, k) - s(i, j - 1, k)) * dxinv[1];
               , az = amrex::Math::abs(s(i, j, k + 1) - s(i, j, k - 1)) * dxinv[2];)
  if (0.5 * amrex::max<amrex::Real>(ax, ay, az) >= grad_threshold) {
    tag(i, j, k) = tagval;
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_tag_error(int i, int j, int k, amrex::Array4<char> const& tag,
              amrex::Array4<amrex::Real const> const& s, amrex::Real err_threshold,
              char tagval) noexcept
{
  if (s(i, j, k) >= err_threshold) { tag(i, j, k) = tagval; }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_tag_range(int i, int j, int k, amrex::Array4<char> const& tag,
              amrex::Array4<amrex::Real const> const& s, amrex::Real lo_threshold,
              amrex::Real hi_threshold, char tagval) noexcept
{
  if ((s(i, j, k) >= lo_threshold) && (s(i, j, k) <= hi_threshold)) {
    tag(i, j, k) = tagval;
  }
}

#endif
