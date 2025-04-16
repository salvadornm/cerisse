#ifndef PASR_H_
#define PASR_H_

#include "CNS.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void unpack_pasr(
  const int i, const int j, const int k, amrex::Real& new_rho,
  const amrex::Array4<const amrex::Real>& sold_arr,
  const amrex::Array4<amrex::Real>& rY,
  const amrex::Array4<const amrex::Real>& rYsrc,
  const amrex::Array4<const amrex::Real>& qarr,
  const Array4<const amrex::Real>& muarr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv, amrex::Real dt)
{
  using amrex::Real;

  // Calculate laminar chemical source term
  Real omega[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    omega[n] =
      (rY(i, j, k, n) - sold_arr(i, j, k, UFS + n)) / dt - rYsrc(i, j, k, n);
  }

  // Calculate chemical timescale tau_chem = min(rY/|omega|)
  Real tau_chem = 1e10;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    tau_chem =
      std::min(rY(i, j, k, n) / std::max(std::abs(omega[n]),
                                         std::numeric_limits<Real>::denorm_min()),
               tau_chem);
  }
  tau_chem = std::max(tau_chem, std::numeric_limits<Real>::epsilon());

  // Calculate mixing timescale tau_mix
  Real mu_T;
  Real delta = std::pow(AMREX_D_TERM(dxinv[0], *dxinv[1], *dxinv[2]),
                        -1.0 / Real(AMREX_SPACEDIM));
  CNS::les_model->mu_T_cc(i, j, k, qarr, dxinv, delta, CNS::Cs, mu_T);
  Real tau_mix =
    CNS::Cm * delta * delta * qarr(i, j, k, QRHO) / (muarr(i, j, k) + mu_T);

  // Recompute rY
  Real rvfrac = 1.0 / (1.0 + tau_mix / tau_chem); // reactive volume frac
  for (int n = 0; n < NUM_SPECIES; ++n) {
    rY(i, j, k, n) = std::max(0.0, sold_arr(i, j, k, UFS + n) +
                                     (rYsrc(i, j, k, n) + rvfrac * omega[n]) * dt);
    new_rho += rY(i, j, k, n);
  }
}

#endif