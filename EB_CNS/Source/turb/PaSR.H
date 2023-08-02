#ifndef PASR_H_
#define PASR_H_

#include "CNS.H"

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void unpack_pasr(
  int i, int j, int k, const Array4<Real>& snew_arr, Real& new_rho,
  const Array4<const Real>& sold_arr, const Array4<const Real>& rY,
  const Array4<const Real>& rYsrc, const Array4<const Real>& qarr,
  const Array4<const Real>& muarr, const GpuArray<Real, AMREX_SPACEDIM>& dx,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, amrex::Real dt)
{
  // Calculate laminar chemical source term
  Real omega[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    omega[n] =
      (rY(i, j, k, n) - sold_arr(i, j, k, UFS + n)) / dt - rYsrc(i, j, k, n);
  }

  // Calculate chemical timescale tau_chem = min(rY/|omega|)
  amrex::Real tau_chem = 1e10;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    tau_chem =
      amrex::min(rY(i, j, k, n) / std::max(std::abs(omega[n]),
                                           std::numeric_limits<Real>::denorm_min()),
                 tau_chem);
  }

  // Calculate mixing timescale tau_mix
  amrex::Real mu_T;
  amrex::Real delta =
    std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]), 1.0 / amrex::Real(AMREX_SPACEDIM));
  CNS::les_model->mu_T_cc(i, j, k, qarr, dxinv, delta, CNS::Cs, mu_T);
  amrex::Real tau_mix =
    CNS::Cm * delta * delta * qarr(i, j, k, QRHO) / (muarr(i, j, k) + mu_T);

  // Recompute rY
  amrex::Real rvfrac = 1.0 / (1.0 + tau_mix / tau_chem); // reactive volume frac
  for (int n = 0; n < NUM_SPECIES; ++n) {
    snew_arr(i, j, k, UFS + n) =
      sold_arr(i, j, k, UFS + n) + (rYsrc(i, j, k, n) + rvfrac * omega[n]) * dt;
    new_rho += snew_arr(i, j, k, UFS + n);
  }
}

#endif