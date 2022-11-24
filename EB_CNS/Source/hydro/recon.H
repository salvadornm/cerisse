#ifndef CNS_RECON_H_
#define CNS_RECON_H_

#include <AMReX_FArrayBox.H>

#include "index_macros.H"
// #include "parm.H"
// #include "PelePhysics.H"

// #include <cmath>

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// amrex::Real
// limiter (amrex::Real dlft, amrex::Real drgt, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

//     Real dcen = Real(0.5)*(dlft+drgt);
//     Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
//     Real slop = plm_theta * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
//     Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
//     return dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
muscl (amrex::Real s[3], amrex::Real& sm, amrex::Real& sp, amrex::Real plm_theta) noexcept
{
  using amrex::Real;
  
  // Calculate slope
  Real dlft, drgt;
  dlft = s[1] - s[0];
  drgt = s[2] - s[1];

  // Limiter
  Real dcen = Real(0.5)*(dlft+drgt);
  Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
  Real slop = plm_theta * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
  Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
  Real dq = dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
  
  sp = s[1] + 0.5*dq;
  sm = s[1] - 0.5*dq;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
wenojs5 (amrex::Real s[5], amrex::Real& sm, amrex::Real& sp) noexcept
{
  using amrex::Real;
    
  const Real eps = 1e-6;
  Real vl[3];
  Real vr[3];
  Real beta[3];
  Real alpha[3];
  Real alpha1;

  beta[2] = Real(13./12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
              0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
              0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
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

  sp = Real(1./6.) * alpha1 *
        (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[2] = Real(13./12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
              0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
              0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
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

  sm = Real(1./6.) * alpha1 *
        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
wenoz5 (amrex::Real s[5], amrex::Real& sm, amrex::Real& sp) noexcept
{
  using amrex::Real;

  const Real eps = 1e-6;
  Real vl[3];
  Real vr[3];
  Real beta[3];
  Real alpha[3];
  Real alpha1, tau;

  beta[2] = Real(13./12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
              0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
              0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
              0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
  vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
  vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

  sp = Real(1./6.) * alpha1 *
        (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

  beta[2] = Real(13./12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
              0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
              0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
              0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
  vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
  vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];

  sm = Real(1./6.) * alpha1 *
        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

/**
 * \brief Interpolate primitive variables to cell faces.
 * 
 * For numerical stability, we need to convert primitive variables into characteristic variables,
 * then convert back to primitive variables after reconstruction.
 * 
 * Notice the indexing convention. 
 * For cell averaged quantities: [  cell i-1  ]|[   cell i   ]|[  cell i+1  ]
 * For face averaged quantities:               ^face i        ^face i+1
 * For ql and qr:                         ql_i^ ^qr_i  ql_i+1^ ^qr_i+1
 * 
 * @param i    x index.
 * @param j    y index.
 * @param k    z index.
 * @param dir  direction, 0:x, 1:y, 2:z.
 * @param q    cell averaged primitive variable.
 * @param ql   left state at a face.
 * @param qr   right state at a face.
 * @param recon_scheme reconstruction scheme, 1:piecewise constant, 2:MUSCL, 3:WENO-JS5, 4:WENO-Z5.
 * @param plt_theta a paramter in MUSCL limiter, between 1.0 and 2.0.
 * \sa q2cha, cha2q
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_recon (int i, int j, int k, int n, int dir,
           amrex::Array4<amrex::Real const> const& q,
           amrex::Array4<amrex::Real> const& ql,
           amrex::Array4<amrex::Real> const& qr,
           int recon_scheme, amrex::Real plm_theta) noexcept
{
  using amrex::Real;

  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);
  
  if (recon_scheme == 1) {
    // Simple Godunov
    qr(iv,n) = q(iv,n);
    ql(iv+iv_dir,n) = q(iv,n);

  } else if (recon_scheme == 2) {
    // MUSCL    
    Real s[3];
    s[0] = q(iv-iv_dir,n);
    s[1] = q(iv,n);
    s[2] = q(iv+iv_dir,n);

    muscl(s, qr(iv,n), ql(iv+iv_dir,n), plm_theta);

  } else if (recon_scheme == 3) {
    // WENO-JS5
    Real s[5];    
    s[0] = q(iv-2*iv_dir,n);
    s[1] = q(iv-iv_dir,n);
    s[2] = q(iv,n);
    s[3] = q(iv+iv_dir,n);
    s[4] = q(iv+2*iv_dir,n);

    wenojs5(s, qr(iv,n), ql(iv+iv_dir,n));
  
  } else {
    // Default: WENO-Z5
    Real s[5];    
    s[0] = q(iv-2*iv_dir,n);
    s[1] = q(iv-iv_dir,n);
    s[2] = q(iv,n);
    s[3] = q(iv+iv_dir,n);
    s[4] = q(iv+2*iv_dir,n);

    wenoz5(s, qr(iv,n), ql(iv+iv_dir,n));
  }
}

#endif