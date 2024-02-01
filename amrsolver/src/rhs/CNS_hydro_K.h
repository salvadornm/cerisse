#ifndef CNS_HYDRO_K_H_
#define CNS_HYDRO_K_H_

#include <AMReX_FArrayBox.H>

#include <cmath>

using namespace amrex;
// Viscous fluxes at cell centers
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ViscousFluxes(
    int i, int j, int k, auto const& prims, auto const& fx, auto const& fy,
    auto const& fz, GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
    PROB::ProbClosures const& cls) noexcept {
  // 2nd order accurate central difference of primitive vars /////////////////
  // x direction
  Real ux = prims(i, j, k, cls.QU);
  Real dudx =
      (prims(i + 1, j, k, cls.QU) - prims(i - 1, j, k, cls.QU)) * 0.5_rt * dxinv[0];
  Real dvdx =
      (prims(i + 1, j, k, cls.QV) - prims(i - 1, j, k, cls.QV)) * 0.5_rt * dxinv[0];
  Real dwdx =
      (prims(i + 1, j, k, cls.QW) - prims(i - 1, j, k, cls.QW)) * 0.5_rt * dxinv[0];
  Real dTdx =
      (prims(i + 1, j, k, cls.QT) - prims(i - 1, j, k, cls.QT)) * 0.5_rt * dxinv[0];

  // y direction
  Real uy = prims(i, j, k, cls.QV);
  Real dudy =
      (prims(i, j + 1, k, cls.QU) - prims(i, j - 1, k, cls.QU)) * 0.5_rt * dxinv[1];
  Real dvdy =
      (prims(i, j + 1, k, cls.QV) - prims(i, j - 1, k, cls.QV)) * 0.5_rt * dxinv[1];
  Real dwdy =
      (prims(i, j + 1, k, cls.QW) - prims(i, j - 1, k, cls.QW)) * 0.5_rt * dxinv[1];
  Real dTdy =
      (prims(i, j + 1, k, cls.QT) - prims(i, j - 1, k, cls.QT)) * 0.5_rt * dxinv[1];

  // z direction
  Real uz = prims(i, j, k, cls.QW);
  Real dudz =
      (prims(i, j, k + 1, cls.QU) - prims(i, j, k - 1, cls.QU)) * 0.5_rt * dxinv[2];
  Real dvdz =
      (prims(i, j, k + 1, cls.QV) - prims(i, j, k - 1, cls.QV)) * 0.5_rt * dxinv[2];
  Real dwdz =
      (prims(i, j, k + 1, cls.QW) - prims(i, j, k - 1, cls.QW)) * 0.5_rt * dxinv[2];
  Real dTdz =
      (prims(i, j, k + 1, cls.QT) - prims(i, j, k - 1, cls.QT)) * 0.5_rt * dxinv[2];

  // divergence
  Real div = dudx + dvdy + dwdz;

  // constants
  Real mu = cls.visc(prims(i, j, k, cls.QT));
  Real lambda = cls.cond(prims(i, j, k, cls.QT));
  Real r1_3 = Real(1.0) / Real(3.0);

  // viscous fluxes
  Real tauxx = Real(2.0) * mu * (dudx - r1_3 * div);
  Real tauxy = mu * (dudy + dvdx);
  Real tauxz = mu * (dudz + dwdx);

  // tauxy = tauyx
  Real tauyy = Real(2.0) * mu * (dvdy - r1_3 * div);
  Real tauyz = mu * (dvdz + dwdy);

  // tauzx = tauxz;
  // tauzy = tauyz;
  Real tauzz = Real(2.0) * mu * (dwdz - r1_3 * div);

  // assemble fluxes on LHS
  fx(i, j, k, cls.URHO) = Real(0.0);
  fx(i, j, k, cls.UMX) = -tauxx;
  fx(i, j, k, cls.UMY) = -tauxy;
  fx(i, j, k, cls.UMZ) = -tauxz;
  fx(i, j, k, cls.UET) = -lambda * dTdx - tauxx * ux - tauxy * uy - tauxz * uz;

  fy(i, j, k, cls.URHO) = Real(0.0);
  fy(i, j, k, cls.UMX) = -tauxy;
  fy(i, j, k, cls.UMY) = -tauyy;
  fy(i, j, k, cls.UMZ) = -tauyz;
  fy(i, j, k, cls.UET) = -lambda * dTdy - tauxy * ux - tauyy * uy - tauyz * uz;

  fz(i, j, k, cls.URHO) = Real(0.0);
  fz(i, j, k, cls.UMX) = -tauxz;
  fz(i, j, k, cls.UMY) = -tauyz;
  fz(i, j, k, cls.UMZ) = -tauzz;
  fz(i, j, k, cls.UET) = -lambda * dTdz - tauxz * ux - tauyz * uy - tauzz * uz;
}

// Viscous fluxes at cell centers
// TODO: remove if statement and generalise to other directions by passing
// one-sided derivative coefficients
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ViscousWallFluxes(
    int i, int j, int k, int loc, auto const& prims, auto const& fx,
    auto const& fy, auto const& fz, GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
    PROB::ProbClosures const& cls) noexcept {
  // 2nd order accurate central difference of primitive vars /////////////////
  // x direction
  Real ux = prims(i, j, k, cls.QU);
  Real dudx =
      (prims(i + 1, j, k, cls.QU) - prims(i - 1, j, k, cls.QU)) * 0.5_rt * dxinv[0];
  Real dvdx =
      (prims(i + 1, j, k, cls.QV) - prims(i - 1, j, k, cls.QV)) * 0.5_rt * dxinv[0];
  Real dwdx =
      (prims(i + 1, j, k, cls.QW) - prims(i - 1, j, k, cls.QW)) * 0.5_rt * dxinv[0];
  Real dTdx =
      (prims(i + 1, j, k, cls.QT) - prims(i - 1, j, k, cls.QT)) * 0.5_rt * dxinv[0];

  // y direction
  Real uy = prims(i, j, k, cls.QV);
  Real dudy, dvdy, dwdy, dTdy;
  // bottom boundary
  if (loc == 0) {
    dudy = (-prims(i, j + 2, k, cls.QU) + 4 * prims(i, j + 1, k, cls.QU) -
            3 * prims(i, j, k, cls.QU)) *
           0.5_rt * dxinv[1];
    dvdy = (-prims(i, j + 2, k, cls.QV) + 4 * prims(i, j + 1, k, cls.QV) -
            3 * prims(i, j, k, cls.QV)) *
           0.5_rt * dxinv[1];
    dwdy = (-prims(i, j + 2, k, cls.QW) + 4 * prims(i, j + 1, k, cls.QW) -
            3 * prims(i, j, k, cls.QW)) *
           0.5_rt * dxinv[1];
    dTdy = (-prims(i, j + 2, k, cls.QT) + 4 * prims(i, j + 1, k, cls.QT) -
            3 * prims(i, j, k, cls.QT)) *
           0.5_rt * dxinv[1];
  }
  // top boundary
  else if (loc == 1) {
    dudy = (3 * prims(i, j, k, cls.QU) - 4 * prims(i, j - 1, k, cls.QU) +
            prims(i, j - 2, k, cls.QU)) *
           0.5_rt * dxinv[1];
    dvdy = (3 * prims(i, j, k, cls.QV) - 4 * prims(i, j - 1, k, cls.QV) +
            prims(i, j - 2, k, cls.QV)) *
           0.5_rt * dxinv[1];
    dwdy = (3 * prims(i, j, k, cls.QW) - 4 * prims(i, j - 1, k, cls.QW) +
            prims(i, j - 2, k, cls.QW)) *
           0.5_rt * dxinv[1];
    dTdy = (3 * prims(i, j, k, cls.QT) - 4 * prims(i, j - 1, k, cls.QT) +
            prims(i, j - 2, k, cls.QT)) *
           0.5_rt * dxinv[1];
  }

  // z direction
  Real uz = prims(i, j, k, cls.QW);
  Real dudz =
      (prims(i, j, k + 1, cls.QU) - prims(i, j, k - 1, cls.QU)) * 0.5_rt * dxinv[2];
  Real dvdz =
      (prims(i, j, k + 1, cls.QV) - prims(i, j, k - 1, cls.QV)) * 0.5_rt * dxinv[2];
  Real dwdz =
      (prims(i, j, k + 1, cls.QW) - prims(i, j, k - 1, cls.QW)) * 0.5_rt * dxinv[2];
  Real dTdz =
      (prims(i, j, k + 1, cls.QT) - prims(i, j, k - 1, cls.QT)) * 0.5_rt * dxinv[2];

  // divergence
  Real div = dudx + dvdy + dwdz;

  // constants
  Real mu = cls.visc(prims(i, j, k, cls.QT));
  Real lambda = cls.cond(prims(i, j, k, cls.QT));
  Real r1_3 = Real(1.0) / Real(3.0);

  // viscous fluxes
  Real tauxx = Real(2.0) * mu * (dudx - r1_3 * div);
  Real tauxy = mu * (dudy + dvdx);
  Real tauxz = mu * (dudz + dwdx);

  // tauxy = tauyx
  Real tauyy = Real(2.0) * mu * (dvdy - r1_3 * div);
  Real tauyz = mu * (dvdz + dwdy);

  // tauzx = tauxz;
  // tauzy = tauyz;
  Real tauzz = Real(2.0) * mu * (dwdz - r1_3 * div);

  // assemble fluxes on LHS
  fx(i, j, k, cls.URHO) = Real(0.0);
  fx(i, j, k, cls.UMX) = -tauxx;
  fx(i, j, k, cls.UMY) = -tauxy;
  fx(i, j, k, cls.UMZ) = -tauxz;
  fx(i, j, k, cls.UET) = -lambda * dTdx - tauxx * ux - tauxy * uy - tauxz * uz;

  fy(i, j, k, cls.URHO) = Real(0.0);
  fy(i, j, k, cls.UMX) = -tauxy;
  fy(i, j, k, cls.UMY) = -tauyy;
  fy(i, j, k, cls.UMZ) = -tauyz;
  fy(i, j, k, cls.UET) = -lambda * dTdy - tauxy * ux - tauyy * uy - tauyz * uz;

  fz(i, j, k, cls.URHO) = Real(0.0);
  fz(i, j, k, cls.UMX) = -tauxz;
  fz(i, j, k, cls.UMY) = -tauyz;
  fz(i, j, k, cls.UMZ) = -tauzz;
  fz(i, j, k, cls.UET) = -lambda * dTdz - tauxz * ux - tauyz * uy - tauzz * uz;
}

// 2nd order accurate interpolation at i-1/2,j-1/2,k-1/2
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ViscousNumericalFluxes(
    int i, int j, int k, int n, const auto& pfx, const auto& pfy,
    const auto& pfz, const auto& nfx, const auto& nfy, const auto& nfz) {
  nfx(i, j, k, n) += Real(0.5) * (pfx(i - 1, j, k, n) + pfx(i, j, k, n));
  nfy(i, j, k, n) += Real(0.5) * (pfy(i, j - 1, k, n) + pfy(i, j, k, n));
  nfz(i, j, k, n) += Real(0.5) * (pfz(i, j, k - 1, n) + pfz(i, j, k, n));
}

///
/// \brief Viscous fluxes at ghost point
/// \param i index
/// \param markers solid
/// \sa CnsFillExtDir
///
/// c1 =
/// ```
/// {rst}
/// Assuming the reconstructed state :math:`p^\text{th}`
/// order of accuracy for the reconstructed state.
///
/// :math:`\frac{\partial f}{\partial x}\bigg|_i=`
/// ```
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ViscousFluxGP(
    int i, int j, int k, const Array4<bool>& markers, auto const& prims,
    const auto& fx, const auto& fy, const auto& fz,
    const GpuArray<Real, AMREX_SPACEDIM>& dxinv,
    const PROB::ProbClosures& cls) {
  Real ux, dudx, dvdx, dwdx, dTdx;
  Real uy, dudy, dvdy, dwdy, dTdy;
  Real uz, dudz, dvdz, dwdz, dTdz;
  Real mu = cls.visc(prims(i, j, k, cls.QT));
  Real lambda = cls.cond(prims(i, j, k, cls.QT));
  Real r1_3 = Real(1.0) / Real(3.0);

  // x direction //
  ux = prims(i, j, k, cls.QU);
  if (markers(i + 1, j, k, 0)) {  // the point to the right is solid then
                                  // one-sided derivative to left
    dudx = (1.5_rt * prims(i, j, k, cls.QU) - 2.0_rt * prims(i - 1, j, k, cls.QU) +
            0.5_rt * prims(i - 2, j, k, cls.QU)) *
           dxinv[0];
    dvdx = (1.5_rt * prims(i, j, k, cls.QV) - 2.0_rt * prims(i - 1, j, k, cls.QV) +
            0.5_rt * prims(i - 2, j, k, cls.QV)) *
           dxinv[0];
    dwdx = (1.5_rt * prims(i, j, k, cls.QW) - 2.0_rt * prims(i - 1, j, k, cls.QW) +
            0.5_rt * prims(i - 2, j, k, cls.QW)) *
           dxinv[0];
    dTdx = (1.5_rt * prims(i, j, k, cls.QT) - 2.0_rt * prims(i - 1, j, k, cls.QT) +
            0.5_rt * prims(i - 2, j, k, cls.QT)) *
           dxinv[0];
  } else if (markers(i - 1, j, k, 0)) {  // the point to the left is solid then
                                         // one-sided derivative to right
    dudx = (-1.5_rt * prims(i, j, k, cls.QU) + 2.0_rt * prims(i + 1, j, k, cls.QU) -
            0.5_rt * prims(i + 2, j, k, cls.QU)) *
           dxinv[1];
    dvdx = (-1.5_rt * prims(i, j, k, cls.QV) + 2.0_rt * prims(i + 1, j, k, cls.QV) -
            0.5_rt * prims(i + 2, j, k, cls.QV)) *
           dxinv[1];
    dwdx = (-1.5_rt * prims(i, j, k, cls.QW) + 2.0_rt * prims(i + 1, j, k, cls.QW) -
            0.5_rt * prims(i + 2, j, k, cls.QW)) *
           dxinv[1];
    dTdx = (-1.5_rt * prims(i, j, k, cls.QT) + 2.0_rt * prims(i + 1, j, k, cls.QT) -
            0.5_rt * prims(i + 2, j, k, cls.QT)) *
           dxinv[1];
  } else {
    dudx =
        (prims(i + 1, j, k, cls.QU) - prims(i - 1, j, k, cls.QU)) * 0.5_rt * dxinv[0];
    dvdx =
        (prims(i + 1, j, k, cls.QV) - prims(i - 1, j, k, cls.QV)) * 0.5_rt * dxinv[0];
    dwdx =
        (prims(i + 1, j, k, cls.QW) - prims(i - 1, j, k, cls.QW)) * 0.5_rt * dxinv[0];
    dTdx =
        (prims(i + 1, j, k, cls.QT) - prims(i - 1, j, k, cls.QT)) * 0.5_rt * dxinv[0];
  }

  // y direction //
  uy = prims(i, j, k, cls.QV);
  if (markers(i, j + 1, k, 0)) {  // the point to the top is solid then
                                  // one-sided derivative to bottom
    dudy = (1.5_rt * prims(i, j, k, cls.QU) - 2.0_rt * prims(i, j - 1, k, cls.QU) +
            0.5_rt * prims(i, j - 2, k, cls.QU)) *
           dxinv[1];
    dvdy = (1.5_rt * prims(i, j, k, cls.QV) - 2.0_rt * prims(i, j - 1, k, cls.QV) +
            0.5_rt * prims(i, j - 2, k, cls.QV)) *
           dxinv[1];
    dwdy = (1.5_rt * prims(i, j, k, cls.QW) - 2.0_rt * prims(i, j - 1, k, cls.QW) +
            0.5_rt * prims(i, j - 2, k, cls.QW)) *
           dxinv[1];
    dTdy = (1.5_rt * prims(i, j, k, cls.QT) - 2.0_rt * prims(i, j - 1, k, cls.QT) +
            0.5_rt * prims(i, j - 2, k, cls.QT)) *
           dxinv[1];
  } else if (markers(i, j - 1, k, 0)) {  // the point to the bottom is solid
                                         // then one-sided derivative to top
    dudy = (-1.5_rt * prims(i, j, k, cls.QU) + 2.0_rt * prims(i, j + 1, k, cls.QU) -
            0.5_rt * prims(i, j + 2, k, cls.QU)) *
           dxinv[1];
    dvdy = (-1.5_rt * prims(i, j, k, cls.QV) + 2.0_rt * prims(i, j + 1, k, cls.QV) -
            0.5_rt * prims(i, j + 2, k, cls.QV)) *
           dxinv[1];
    dwdy = (-1.5_rt * prims(i, j, k, cls.QW) + 2.0_rt * prims(i, j + 1, k, cls.QW) -
            0.5_rt * prims(i, j + 2, k, cls.QW)) *
           dxinv[1];
    dTdy = (-1.5_rt * prims(i, j, k, cls.QT) + 2.0_rt * prims(i, j + 1, k, cls.QT) -
            0.5_rt * prims(i, j + 2, k, cls.QT)) *
           dxinv[1];
  } else {
    dudy =
        (prims(i, j + 1, k, cls.QU) - prims(i, j - 1, k, cls.QU)) * 0.5_rt * dxinv[1];
    dvdy =
        (prims(i, j + 1, k, cls.QV) - prims(i, j - 1, k, cls.QV)) * 0.5_rt * dxinv[1];
    dwdy =
        (prims(i, j + 1, k, cls.QW) - prims(i, j - 1, k, cls.QW)) * 0.5_rt * dxinv[1];
    dTdy =
        (prims(i, j + 1, k, cls.QT) - prims(i, j - 1, k, cls.QT)) * 0.5_rt * dxinv[1];
  }

  // z direction
  uz = prims(i, j, k, cls.QW);
  if (markers(i, j, k - 1, 0)) {  // the point into the screen is solid then
                                  // one-sided derivative out of the screen
    dudz = (-1.5_rt * prims(i, j, k, cls.QU) + 2.0_rt * prims(i, j, k + 1, cls.QU) -
            0.5_rt * prims(i, j, k + 2, cls.QU)) *
           dxinv[2];
    dvdz = (-1.5_rt * prims(i, j, k, cls.QV) + 2.0_rt * prims(i, j, k + 1, cls.QV) -
            0.5_rt * prims(i, j, k + 2, cls.QV)) *
           dxinv[2];
    dwdz = (-1.5_rt * prims(i, j, k, cls.QW) + 2.0_rt * prims(i, j, k + 1, cls.QW) -
            0.5_rt * prims(i, j, k + 2, cls.QW)) *
           dxinv[2];
    dTdz = (-1.5_rt * prims(i, j, k, cls.QT) + 2.0_rt * prims(i, j, k + 1, cls.QT) -
            0.5_rt * prims(i, j, k + 2, cls.QT)) *
           dxinv[2];
  }

  if (markers(i, j, k + 1, 0)) {  // the point out of the screen is solid then
                                  // one-sided derivative into the screen
    dudz = (1.5_rt * prims(i, j, k, cls.QU) - 2.0_rt * prims(i, j, k - 1, cls.QU) +
            0.5_rt * prims(i, j, k - 2, cls.QU)) *
           dxinv[2];
    dvdz = (1.5_rt * prims(i, j, k, cls.QV) - 2.0_rt * prims(i, j, k - 1, cls.QV) +
            0.5_rt * prims(i, j, k - 2, cls.QV)) *
           dxinv[2];
    dwdz = (1.5_rt * prims(i, j, k, cls.QW) - 2.0_rt * prims(i, j, k - 1, cls.QW) +
            0.5_rt * prims(i, j, k - 2, cls.QW)) *
           dxinv[2];
    dTdz = (1.5_rt * prims(i, j, k, cls.QT) - 2.0_rt * prims(i, j, k - 1, cls.QT) +
            0.5_rt * prims(i, j, k - 2, cls.QT)) *
           dxinv[2];
  } else {
    uz = prims(i, j, k, cls.QW);
    dudz =
        (prims(i, j, k + 1, cls.QU) - prims(i, j, k - 1, cls.QU)) * 0.5_rt * dxinv[2];
    dvdz =
        (prims(i, j, k + 1, cls.QV) - prims(i, j, k - 1, cls.QV)) * 0.5_rt * dxinv[2];
    dwdz =
        (prims(i, j, k + 1, cls.QW) - prims(i, j, k - 1, cls.QW)) * 0.5_rt * dxinv[2];
    dTdz =
        (prims(i, j, k + 1, cls.QT) - prims(i, j, k - 1, cls.QT)) * 0.5_rt * dxinv[2];
  }

  // divergence
  Real div = dudx + dvdy + dwdz;

  // viscous fluxes
  Real tauxx = Real(2.0) * mu * (dudx - r1_3 * div);
  Real tauxy = mu * (dudy + dvdx);
  Real tauxz = mu * (dudz + dwdx);

  // tauxy = tauyx
  Real tauyy = Real(2.0) * mu * (dvdy - r1_3 * div);
  Real tauyz = mu * (dvdz + dwdy);

  // tauzx = tauxz;
  // tauzy = tauyz;
  Real tauzz = Real(2.0) * mu * (dwdz - r1_3 * div);

  // assemble fluxes on LHS
  fx(i, j, k, cls.URHO) = Real(0.0);
  fx(i, j, k, cls.UMX) = -tauxx;
  fx(i, j, k, cls.UMY) = -tauxy;
  fx(i, j, k, cls.UMZ) = -tauxz;
  fx(i, j, k, cls.UET) = -lambda * dTdx - tauxx * ux - tauxy * uy - tauxz * uz;

  fy(i, j, k, cls.URHO) = Real(0.0);
  fy(i, j, k, cls.UMX) = -tauxy;
  fy(i, j, k, cls.UMY) = -tauyy;
  fy(i, j, k, cls.UMZ) = -tauyz;
  fy(i, j, k, cls.UET) = -lambda * dTdy - tauxy * ux - tauyy * uy - tauyz * uz;

  fz(i, j, k, cls.URHO) = Real(0.0);
  fz(i, j, k, cls.UMX) = -tauxz;
  fz(i, j, k, cls.UMY) = -tauyz;
  fz(i, j, k, cls.UMZ) = -tauzz;
  fz(i, j, k, cls.UET) = -lambda * dTdz - tauxz * ux - tauyz * uy - tauzz * uz;
};

#endif