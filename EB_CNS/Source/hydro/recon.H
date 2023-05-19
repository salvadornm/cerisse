#ifndef CNS_RECON_H_
#define CNS_RECON_H_

// #include <cmath>

#include <AMReX_FArrayBox.H>

#include "index_macros.H"
// #include "parm.H"
// #include "PelePhysics.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
wenoz3 (const amrex::Real s[3], amrex::Real& sm, amrex::Real& sp) noexcept
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
teno5 (amrex::Real s[5], amrex::Real& sm, amrex::Real& sp) noexcept
{
  using amrex::Real;

  const Real eps = 1e-40;
  Real vl[3];
  Real vr[3];
  Real beta[3];
  Real alpha[3];
  Real beta1, alpha1, tau;
  Real CT = 1e-4; //sharp cutoff smoothness indicator

  beta[2] = Real(13./12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
              0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
              0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
              0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  tau = amrex::Math::abs(tau - (1.0/6.0)*(beta[2] + 4*beta[1] + beta[0]));
  beta[2] = pow(1.0 + tau / (eps + beta[2]), 6);
  beta[1] = pow(1.0 + tau / (eps + beta[1]), 6);
  beta[0] = pow(1.0 + tau / (eps + beta[0]), 6);
  // beta[2] = pow(eps + beta[2], -6);
  // beta[1] = pow(eps + beta[1], -6);
  // beta[0] = pow(eps + beta[0], -6);
  beta1 = 1.0 / (beta[2] + beta[1] + beta[0]);

  alpha[2] = beta[2] * beta1 < CT ? 0.0 : 1.0;
  alpha[1] = beta[1] * beta1 < CT ? 0.0 : 6.0;
  alpha[0] = beta[0] * beta1 < CT ? 0.0 : 3.0;
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
  tau = amrex::Math::abs(tau - (1.0/6.0)*(beta[2] + 4*beta[1] + beta[0]));
  beta[2] = pow(1.0 + (tau / (eps + beta[2])), 6);
  beta[1] = pow(1.0 + (tau / (eps + beta[1])), 6);
  beta[0] = pow(1.0 + (tau / (eps + beta[0])), 6);
  // beta[2] = pow(eps + beta[2], -6);
  // beta[1] = pow(eps + beta[1], -6);
  // beta[0] = pow(eps + beta[0], -6);
  beta1 = 1.0 / (beta[2] + beta[1] + beta[0]);

  alpha[2] = beta[2] * beta1 < CT ? 0.0 : 1.0;
  alpha[1] = beta[1] * beta1 < CT ? 0.0 : 6.0;
  alpha[0] = beta[0] * beta1 < CT ? 0.0 : 3.0;
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
  vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
  vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];

  sm = Real(1./6.) * alpha1 *
        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pp_limiter_scalar (amrex::Real& ql, amrex::Real& qr, const amrex::Real& q, 
                   const amrex::Real m, const amrex::Real M = 1e10)
{
  using amrex::Real;
  Real min_q = amrex::min(ql, qr);
  Real max_q = amrex::max(ql, qr);

  if ((min_q < m) || (max_q > M)) {
    std::cout << "PP Limiter Scalar activated: ql=" << ql << " qr=" << qr << " q=" << q << "\n";
    
    Real t = amrex::Math::abs(q - m) / amrex::Math::abs(q - min_q);
    t = amrex::min(t, amrex::Math::abs(q - M) / amrex::Math::abs(q - max_q));
    AMREX_ASSERT(t <= 1);

    ql = t*ql + (1-t)*q;
    qr = t*qr + (1-t)*q;
    std::cout << "New ql=" << ql << " New qr=" << qr << "\n";
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pp_limiter_euler (const amrex::IntVect ivl, const amrex::IntVect iv, 
                  amrex::Array4<amrex::Real> const& ql, 
                  amrex::Array4<amrex::Real> const& qr, 
                  amrex::Array4<const amrex::Real> const& q, 
                  const amrex::Real mr, const amrex::Real mp)
{
  using amrex::Real;
  
  // Step 1: limit density as a scalar
  pp_limiter_scalar(ql(ivl, QRHO), qr(iv, QRHO), q(iv, QRHO), mr, 100);

  // Step 2: limit the whole state for positive pressure
  if ((ql(ivl, QPRES) < mp) || (qr(iv, QPRES) < mp)) {
    std::cout << iv << " PP Limiter Euler activated: pl=" << ql(ivl, QPRES) 
                                                << " pr=" << qr(iv, QPRES)
                                                << " p=" << q(iv, QPRES) << "\n";
    Real min_p = amrex::min(ql(ivl, QPRES), qr(iv, QPRES));
    Real t = amrex::Math::abs(q(iv, QPRES) - mp) / amrex::Math::abs(q(iv, QPRES) - min_p);
    // Real min_p = amrex::min(ql(ivl, QPRES)/ql(ivl, QRHO), qr(iv, QPRES)/qr(iv, QRHO));
    // Real t = amrex::Math::abs(q(iv, QPRES)/q(iv, QRHO) - mp/mr) 
    //        / amrex::Math::abs(q(iv, QPRES)/q(iv, QRHO) - min_p);
    AMREX_ASSERT(t <= 1);

    for (int n = QRHO; n <= QPRES; ++n) {
      std::cout << "Old ql=" << ql(ivl, n) << " Old qr=" << qr(iv, n);
      ql(ivl, n) = t*ql(ivl, n) + (1-t)*q(iv, n);
      qr(iv, n)  = t*qr(iv, n) + (1-t)*q(iv, n);
      std::cout << "  New ql=" << ql(ivl, n) << " New qr=" << qr(iv, n) << "\n";
    }

    Real T, Y[NUM_SPECIES];
    auto eos = pele::physics::PhysicsType::eos();
    for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = ql(ivl,QFS+n);
    eos.RYP2T(ql(ivl,QRHO), Y, ql(ivl, QPRES), T);
    std::cout << "New Tl=" << T;
    for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = qr(iv,QFS+n);
    eos.RYP2T(qr(iv,QRHO), Y, qr(iv, QPRES), T);
    std::cout << " New Tr=" << T << "\n";
  }
}

/**
 * \brief Interpolate variables to cell faces.
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
 * @param n    the n-th variable in the array.
 * @param dir  direction, 0:x, 1:y, 2:z.
 * @param q    cell averaged variable.
 * @param ql   left state at a face.
 * @param qr   right state at a face.
 * @param recon_scheme reconstruction scheme.
 * @param plt_theta a paramter in MUSCL limiter, between 1.0 and 2.0.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
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
    // WENO-Z3     
    Real s[3];
    s[0] = q(iv-iv_dir,n);
    s[1] = q(iv,n);
    s[2] = q(iv+iv_dir,n);

    wenoz3(s, qr(iv,n), ql(iv+iv_dir,n));

  } else if (recon_scheme == 4) {
    // WENO-JS5
    Real s[5];    
    s[0] = q(iv-2*iv_dir,n);
    s[1] = q(iv-iv_dir,n);
    s[2] = q(iv,n);
    s[3] = q(iv+iv_dir,n);
    s[4] = q(iv+2*iv_dir,n);

    wenojs5(s, qr(iv,n), ql(iv+iv_dir,n));
  
  } else if (recon_scheme == 6) {
    // TENO-5
    Real s[5];    
    s[0] = q(iv-2*iv_dir,n);
    s[1] = q(iv-iv_dir,n);
    s[2] = q(iv,n);
    s[3] = q(iv+iv_dir,n);
    s[4] = q(iv+2*iv_dir,n);

    teno5(s, qr(iv,n), ql(iv+iv_dir,n));
    
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