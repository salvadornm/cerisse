#ifndef CNS_RECON_H_
#define CNS_RECON_H_

#include <AMReX_FArrayBox.H>

#include "index_macros.H"

/**
 * \brief MUSCL/TVD/PLM.
 * @param s       stencils.
 * @param sm[out] the value of the interpolation on the left edge of the i cell.
 * @param sp[out] the value of the interpolation on the right edge of the i cell.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void muscl(amrex::Real s[3], amrex::Real& sm,
                                               amrex::Real& sp,
                                               amrex::Real plm_theta) noexcept
{
  using amrex::Real;

  // Calculate slope
  Real dlft, drgt;
  dlft = s[1] - s[0];
  drgt = s[2] - s[1];

  // Limiter
  Real dcen = Real(0.5) * (dlft + drgt);
  Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
  Real slop = plm_theta * amrex::min(amrex::Math::abs(dlft), amrex::Math::abs(drgt));
  Real dlim = (dlft * drgt >= Real(0.0)) ? slop : Real(0.0);
  Real dq = dsgn * amrex::min(dlim, amrex::Math::abs(dcen));

  sp = s[1] + 0.5 * dq;
  sm = s[1] - 0.5 * dq;
}

/// \brief WENO-Z 3rd order.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void wenoz3(const amrex::Real s[3],
                                                amrex::Real& sm,
                                                amrex::Real& sp) noexcept
{
  using amrex::Real;

  constexpr Real eps = 1e-6;
  Real vl[2];
  Real vr[2];
  Real beta[2];
  Real alpha[2];

  beta[1] = (s[0] - s[1]) * (s[0] - s[1]);
  beta[0] = (s[1] - s[2]) * (s[1] - s[2]);

  Real tau = amrex::Math::abs(beta[1] - beta[0]);

  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[1] = 2.0 * beta[1];
  alpha[0] = 1.0 * beta[0];
  Real alpha1 = 1.0 / (alpha[1] + alpha[0]);

  vl[1] = -s[0] + 3.0 * s[1];
  vl[0] = s[1] + s[2];

  sp = 0.5 * alpha1 * (alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[1] = (s[2] - s[1]) * (s[2] - s[1]);
  beta[0] = (s[1] - s[0]) * (s[1] - s[0]);

  tau = amrex::Math::abs(beta[1] - beta[0]);

  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[1] = 2.0 * beta[1];
  alpha[0] = 1.0 * beta[0];
  alpha1 = 1.0 / (alpha[1] + alpha[0]);

  vr[1] = 3.0 * s[1] - s[2];
  vr[0] = s[0] + s[1];

  sm = 0.5 * alpha1 * (alpha[1] * vr[1] + alpha[0] * vr[0]);
}

/// \brief WENO-JS 5th order. Ref https://doi.org/10.1006/jcph.1996.0130.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void wenojs5(amrex::Real s[5], amrex::Real& sm,
                                                 amrex::Real& sp) noexcept
{
  using amrex::Real;

  const Real eps = 1e-6;
  Real vl[3];
  Real vr[3];
  Real beta[3];
  Real alpha[3];
  Real alpha1;

  beta[2] = Real(13. / 12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
            0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13. / 12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
            0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13. / 12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
  beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
  beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
  vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
  vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

  sp = Real(1. / 6.) * alpha1 *
       (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[2] = Real(13. / 12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
            0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = Real(13. / 12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
            0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = Real(13. / 12.) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
  beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
  beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
  vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
  vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];

  sm = Real(1. / 6.) * alpha1 *
       (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

/// \brief WENO-Z 5th order. Ref https://doi.org/10.1016/j.jcp.2010.11.028.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void wenoz5(amrex::Real s[5], amrex::Real& sm,
                                                amrex::Real& sp) noexcept
{
  using amrex::Real;

  const Real eps = 1e-20;
  Real vl[3];
  Real vr[3];
  Real beta[3];
  Real alpha[3];
  Real alpha1, tau;

  beta[2] = Real(13. / 12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
            0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13. / 12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
            0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13. / 12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 10.0 * beta[1];
  alpha[0] = 5.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 3.0 * s[0] - 10.0 * s[1] + 15.0 * s[2];
  vl[1] = -s[1] + 6.0 * s[2] + 3.0 * s[3];
  vl[0] = 3.0 * s[2] + 6.0 * s[3] - s[4];

  sp = Real(1. / 8.) * alpha1 *
       (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[2] = Real(13. / 12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
            0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = Real(13. / 12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
            0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = Real(13. / 12.) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 10.0 * beta[1];
  alpha[0] = 5.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 3.0 * s[4] - 10.0 * s[3] + 15.0 * s[2];
  vr[1] = -s[3] + 6.0 * s[2] + 3.0 * s[1];
  vr[0] = 3.0 * s[2] + 6.0 * s[1] - s[0];

  sm = Real(1. / 8.) * alpha1 *
       (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void wenoz5(amrex::Real s[5], amrex::Real& sm,
//                                                 amrex::Real& sp) noexcept
// {
//   using amrex::Real;

//   const Real eps = 1e-20;
//   Real vl[3];
//   Real vr[3];
//   Real beta[3];
//   Real alpha[3];
//   Real alpha1, tau;

//   beta[2] = pow(s[0] - 2.0 * s[1] + s[2], 2) +
//             0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
//   beta[1] = pow(s[1] - 2.0 * s[3], 2) +
//             0.25 * (s[1] - s[3]) * (s[1] - s[3]);
//   beta[0] = pow(s[2] - 2.0 * s[3] + s[4], 2) +
//             0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

//   tau = amrex::Math::abs(beta[2] - beta[0]);
//   beta[2] = 1.0 + (tau / (eps + beta[2]));// * (tau / (eps + beta[2]));
//   beta[1] = 1.0 + (tau / (eps + beta[1]));// * (tau / (eps + beta[1]));
//   beta[0] = 1.0 + (tau / (eps + beta[0]));// * (tau / (eps + beta[0]));

//   alpha[2] = beta[2];
//   alpha[1] = 10.0 * beta[1];
//   alpha[0] = 5.0 * beta[0];
//   alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

//   vl[2] = 3.0 * s[0] - 10.0 * s[1] + 15.0 * s[2];
//   vl[1] = -s[1] + 6.0 * s[2] + 3.0 * s[3];
//   vl[0] = 3.0 * s[2] + 6.0 * s[3] - s[4];

//   sp = Real(1. / 8.) * alpha1 *
//        (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

//   beta[2] = pow(s[4] - 2.0 * s[3] + s[2], 2) +
//             0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
//   beta[1] = pow(s[3] - 2.0 * s[1], 2) +
//             0.25 * (s[3] - s[1]) * (s[3] - s[1]);
//   beta[0] = pow(s[2] - 2.0 * s[1] + s[0], 2) +
//             0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

//   tau = amrex::Math::abs(beta[2] - beta[0]);
//   beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
//   beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
//   beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

//   alpha[2] = beta[2];
//   alpha[1] = 10.0 * beta[1];
//   alpha[0] = 5.0 * beta[0];
//   alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

//   vr[2] = 3.0 * s[4] - 10.0 * s[3] + 15.0 * s[2];
//   vr[1] = -s[3] + 6.0 * s[2] + 3.0 * s[1];
//   vr[0] = 3.0 * s[2] + 6.0 * s[1] - s[0];

//   sm = Real(1. / 8.) * alpha1 *
//        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
// }

/// \brief TENO 5th order. Ref https://doi.org/10.1016/j.cpc.2018.10.009.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void teno5(amrex::Real s[5], amrex::Real& sm,
                                               amrex::Real& sp) noexcept
{
  using amrex::Real;

  // const Real eps = 1e-30; // the original 1e-40 is too close to machine zero
  constexpr Real eps = std::numeric_limits<Real>::epsilon();
  Real v[3];
  Real beta[3];
  Real alpha[3];
  Real xi[3];
  Real beta1, alpha1, tau;
  Real CT; // sharp cutoff smoothness indicator

  beta[2] = Real(13. / 12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
            0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13. / 12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
            0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13. / 12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  tau = amrex::Math::abs(tau - (1.0 / 6.0) * (beta[2] + 4 * beta[1] + beta[0]));
  xi[2] = tau / (eps + beta[2]);
  xi[1] = tau / (eps + beta[1]);
  xi[0] = tau / (eps + beta[0]);
  beta[2] = pow(1.0 + xi[2], int(6));
  beta[1] = pow(1.0 + xi[1], int(6));
  beta[0] = pow(1.0 + xi[0], int(6));
  beta1 = 1.0 / (beta[2] + beta[1] + beta[0]);  

  if (CNS::teno_cutoff <= 0) {
    // Adaptive TENO, Ref: https://doi.org/10.1016/j.jcp.2020.109902
    Real H = 0.1; // the stencil is considered discontinuous if max(xi) > 1/H
    amrex::Real theta = 1.0 / (1.0 + H * amrex::max(xi[0], xi[1], xi[2]));
    amrex::Real m = 4.0 + std::floor(theta * 6.0);
    CT = std::pow(10.0, -m);
  } else {
    CT = CNS::teno_cutoff;
  } 

  // TODO: consider optimised high spectral resolution weights 1.235341937, 5.065006634, 3.699651429?
  alpha[2] = beta[2] * beta1 < CT ? 0.0 : 1.0;
  alpha[1] = beta[1] * beta1 < CT ? 0.0 : 10.0;
  alpha[0] = beta[0] * beta1 < CT ? 0.0 : 5.0;
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  v[2] = 3.0 * s[0] - 10.0 * s[1] + 15.0 * s[2];
  v[1] = -s[1] + 6.0 * s[2] + 3.0 * s[3];
  v[0] = 3.0 * s[2] + 6.0 * s[3] - s[4];

  sp = Real(1. / 8.) * alpha1 *
       (alpha[2] * v[2] + alpha[1] * v[1] + alpha[0] * v[0]);

  beta[2] = Real(13. / 12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
            0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = Real(13. / 12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
            0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = Real(13. / 12.) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
            0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  tau = amrex::Math::abs(tau - (1.0 / 6.0) * (beta[2] + 4 * beta[1] + beta[0]));
  xi[2] = tau / (eps + beta[2]);
  xi[1] = tau / (eps + beta[1]);
  xi[0] = tau / (eps + beta[0]);
  beta[2] = pow(1.0 + xi[2], int(6));
  beta[1] = pow(1.0 + xi[1], int(6));
  beta[0] = pow(1.0 + xi[0], int(6));
  beta1 = 1.0 / (beta[2] + beta[1] + beta[0]);

  if (CNS::teno_cutoff <= 0) {
    // Adaptive TENO, Ref: https://doi.org/10.1016/j.jcp.2020.109902
    Real H = 0.1;
    amrex::Real theta = 1.0 / (1.0 + H * amrex::max(xi[0], xi[1], xi[2]));
    amrex::Real m = 4.0 + std::floor(theta * 6.0);
    CT = std::pow(10.0, -m);
  }

  alpha[2] = beta[2] * beta1 < CT ? 0.0 : 1.0;
  alpha[1] = beta[1] * beta1 < CT ? 0.0 : 10.0;
  alpha[0] = beta[0] * beta1 < CT ? 0.0 : 5.0;
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  v[2] = 3.0 * s[4] - 10.0 * s[3] + 15.0 * s[2];
  v[1] = -s[3] + 6.0 * s[2] + 3.0 * s[1];
  v[0] = 3.0 * s[2] + 6.0 * s[1] - s[0];

  sm = Real(1. / 8.) * alpha1 *
       (alpha[2] * v[2] + alpha[1] * v[1] + alpha[0] * v[0]);
}

/**
 * \brief Interpolate variables to cell faces.
 *
 * For numerical stability, we need to convert primitive variables into
 * characteristic variables, then convert back to primitive variables after
 * reconstruction.
 *
 * Note the indexing convention.
 * For cell averaged quantities: [  cell i-1  ]|[   cell i   ]|[  cell i+1  ]
 * For face averaged quantities:               ^face i        ^face i+1
 * For ql and qr:                         ql_i^ ^qr_i  ql_i+1^ ^qr_i+1
 *
 * @param i,j,k    x, y, z index.
 * @param n    the n-th variable in the array.
 * @param dir  direction, 0:x, 1:y, 2:z.
 * @param q    cell averaged variable.
 * @param[out] ql  left state at a face.
 * @param[out] qr  right state at a face.
 * @param recon_scheme reconstruction scheme.
 * @param plt_theta a paramter in MUSCL limiter, between 1.0 and 2.0.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_recon(
  int i, int j, int k, int n, int dir, amrex::Array4<amrex::Real const> const& q,
  amrex::Array4<amrex::Real> const& ql, amrex::Array4<amrex::Real> const& qr,
  int recon_scheme, amrex::Real plm_theta) noexcept
{
  using amrex::Real;

  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

  if (recon_scheme == 1) {
    // Simple Godunov
    qr(iv, n) = q(iv, n);
    ql(iv + iv_dir, n) = q(iv, n);

  } else if (recon_scheme == 2) {
    // MUSCL
    Real s[3];
    s[0] = q(iv - iv_dir, n);
    s[1] = q(iv, n);
    s[2] = q(iv + iv_dir, n);

    muscl(s, qr(iv, n), ql(iv + iv_dir, n), plm_theta);

  } else if (recon_scheme == 3) {
    // WENO-Z3
    Real s[3];
    s[0] = q(iv - iv_dir, n);
    s[1] = q(iv, n);
    s[2] = q(iv + iv_dir, n);

    wenoz3(s, qr(iv, n), ql(iv + iv_dir, n));

  } else if (recon_scheme == 4) {
    // WENO-JS5
    Real s[5];
    s[0] = q(iv - 2 * iv_dir, n);
    s[1] = q(iv - iv_dir, n);
    s[2] = q(iv, n);
    s[3] = q(iv + iv_dir, n);
    s[4] = q(iv + 2 * iv_dir, n);

    wenojs5(s, qr(iv, n), ql(iv + iv_dir, n));

  } else if (recon_scheme == 6) {
    // TENO-5
    Real s[5];
    s[0] = q(iv - 2 * iv_dir, n);
    s[1] = q(iv - iv_dir, n);
    s[2] = q(iv, n);
    s[3] = q(iv + iv_dir, n);
    s[4] = q(iv + 2 * iv_dir, n);

    teno5(s, qr(iv, n), ql(iv + iv_dir, n));

  } else {
    // Default: WENO-Z5
    Real s[5];
    s[0] = q(iv - 2 * iv_dir, n);
    s[1] = q(iv - iv_dir, n);
    s[2] = q(iv, n);
    s[3] = q(iv + iv_dir, n);
    s[4] = q(iv + 2 * iv_dir, n);

    wenoz5(s, qr(iv, n), ql(iv + iv_dir, n));
  }
}

#endif