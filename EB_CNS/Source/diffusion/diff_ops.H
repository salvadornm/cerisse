#ifndef CNS_DIFF_OPS_H_
#define CNS_DIFF_OPS_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>

//// Normal derivative with variable order using central difference. ////

template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;

/// @brief 2nd order central difference.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff<2>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return (q(iv, comp) - q(iv - ivn, comp)) * dxinv[idir];
}

/// @brief 4th order central difference.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff<4>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);

  // Note that the normal stencil is [-1.5, -0.5, 0.5, 1.5]*dx
  //  i-2,  i-1  |  i,  i+1
  //             ^ we calculate dqdn at this interface
  Real w1 = -1.0 / 12.0; // only second order is "exact"
  Real w2 = 5.0 / 4.0;
  Real w3 = -5.0 / 4.0;
  Real w4 = 1.0 / 12.0;

  return (w1 * q(iv + ivn, comp) + w2 * q(iv, comp) + w3 * q(iv - ivn, comp) +
          w4 * q(iv - 2 * ivn, comp)) *
         dxinv[idir];
}

//// Tangential derivative with variable order using central difference. ////

template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff(amrex::IntVect iv, int idir, int tdir, int comp,
             amrex::Array4<amrex::Real const> const& q,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;

/// @brief 2nd order central difference. d/d(tdir) on face normal to idir.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff<2>(amrex::IntVect iv, int idir, int tdir, int comp,
                amrex::Array4<amrex::Real const> const& q,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  const auto ivt = amrex::IntVect::TheDimensionVector(tdir);
  return (q(iv + ivt, comp) - q(iv - ivt, comp) + q(iv - ivn + ivt, comp) -
          q(iv - ivn - ivt, comp)) *
         0.25 * dxinv[tdir];
}

/// @brief 4th order central difference. d/d(tdir) on face normal to idir.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff<4>(amrex::IntVect iv, int idir, int tdir, int comp,
                amrex::Array4<amrex::Real const> const& q,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  const auto ivt = amrex::IntVect::TheDimensionVector(tdir);

  // // Note that the transverse stencil is [-2, -1, 1, 2]*dx
  // // weights for derivatives
  Real w1 = -1.0 / 12.0;
  Real w2 = 8.0 / 12.0;
  Real w3 = -8.0 / 12.0;
  Real w4 = 1.0 / 12.0;

  Real d4 = w1 * q(iv + ivn + 2 * ivt, comp) + w2 * q(iv + ivn + ivt, comp) +
            w3 * q(iv + ivn - ivt, comp) + w4 * q(iv + ivn - 2 * ivt, comp);
  Real d3 = w1 * q(iv + 2 * ivt, comp) + w2 * q(iv + ivt, comp) +
            w3 * q(iv - ivt, comp) + w4 * q(iv - 2 * ivt, comp);
  Real d2 = w1 * q(iv - ivn + 2 * ivt, comp) + w2 * q(iv - ivn + ivt, comp) +
            w3 * q(iv - ivn - ivt, comp) + w4 * q(iv - ivn - 2 * ivt, comp);
  Real d1 = w1 * q(iv - 2 * ivn + 2 * ivt, comp) + w2 * q(iv - 2 * ivn + ivt, comp) +
            w3 * q(iv - 2 * ivn - ivt, comp) + w4 * q(iv - 2 * ivn - 2 * ivt, comp);

  Real i1 = -1.0 / 12.0;
  Real i2 = 7.0 / 12.0;

  return (i1 * d1 + i2 * d2 + i2 * d3 + i1 * d4) * dxinv[tdir];
}

//// Interpolation with variable order using polynomial. ////

template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp(amrex::IntVect iv, int idir, int comp,
       amrex::Array4<amrex::Real const> const& q) noexcept;

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp<2>(amrex::IntVect iv, int idir, int comp,
          amrex::Array4<amrex::Real const> const& q) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return 0.5 * (q(iv, comp) + q(iv - ivn, comp));
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp<4>(amrex::IntVect iv, int idir, int comp,
          amrex::Array4<amrex::Real const> const& q) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);

  // Note that the stencil is [-1.5, -0.5, 0.5, 1.5]*dx
  Real w1 = -1.0 / 16.0;
  Real w2 = 9.0 / 16.0;

  return (w1 * q(iv + ivn, comp) + w2 * q(iv, comp) + w2 * q(iv - ivn, comp) +
          w1 * q(iv - 2 * ivn, comp));
}

#endif