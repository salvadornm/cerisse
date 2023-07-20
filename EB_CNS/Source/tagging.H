#ifndef CNS_TAGGING_H
#define CNS_TAGGING_H

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>

#include <cmath>

// struct TaggingParm
// {
//   amrex::Real denerr = 1.0e10;
//   amrex::Real dengrad = 1.0e10;
//   amrex::Real denratio = 1.0e10;
//   amrex::Real presserr = 1.0e10;
//   amrex::Real pressgrad = 1.0e10;
//   amrex::Real velerr = 1.0e10;
//   amrex::Real velgrad = 1.0e10;
//   amrex::Real vorterr = 1.0e10;
//   amrex::Real temperr = 1.0e10;
//   amrex::Real lotemperr = -1.0e10;
//   amrex::Real tempgrad = 1.0e10;
//   amrex::Real ftracerr = 1.0e10;
//   amrex::Real ftracgrad = 1.0e10;
//   amrex::Real vfracerr = 1.0e10;
//   amrex::Real detag_eb_factor = 3.0;

//   std::string eb_refine_type = "static";

//   int max_denerr_lev = 10;
//   int max_dengrad_lev = 10;
//   int max_denratio_lev = 10;
//   int max_presserr_lev = 10;
//   int max_pressgrad_lev = 10;
//   int max_velerr_lev = 10;
//   int max_velgrad_lev = 10;
//   int max_vorterr_lev = 10;
//   int max_temperr_lev = 10;
//   int max_lotemperr_lev = 10;
//   int max_tempgrad_lev = 10;
//   int max_ftracerr_lev = 10;
//   int max_ftracgrad_lev = 10;
//   int max_vfracerr_lev = 10;
//   int max_eb_refine_lev = 10;
//   int min_eb_refine_lev = -1;
//   int adapt_eb_refined_lev = -1;

//   amrex::Vector<amrex::AMRErrorTag> err_tags;
// };

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

// TODO
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// cns_shock_sensor (int i, int j, int k,
//                   amrex::Array4<char> const& tag,
//                   amrex::Array4<amrex::Real const> const& s,
//                   amrex::Real lo_threshold, amrex::Real hi_threshold, char tagval)
//                   noexcept
// {
// }

#endif
