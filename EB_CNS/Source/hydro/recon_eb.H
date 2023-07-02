#ifndef CNS_RECON_EB_H_
#define CNS_RECON_EB_H_

#include "recon.H"

/**
 * \brief EB version of cns_recon.
 * 
 * @param i,j,k    x, y, z index.
 * @param n        the n-th variable in the array.
 * @param dir      direction, 0:x, 1:y, 2:z.
 * @param q        cell averaged variable.
 * @param[out] ql  left state at a face.
 * @param[out] qr  right state at a face.
 * @param recon_scheme reconstruction scheme.
 * @param plt_theta a paramter in MUSCL limiter, between 1.0 and 2.0.
 * @param flag     EB flags array.
 * @param eb_mode  0: fill covered stencil with values of the first cell from wall; 
 *                 1: automatically switch to lower order scheme.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_recon_eb (int i, int j, int k, int n, int dir,
              amrex::Array4<const amrex::Real> const& q,
              amrex::Array4<      amrex::Real> const& ql,
              amrex::Array4<      amrex::Real> const& qr,
              int recon_scheme, amrex::Real plm_theta,
              amrex::Array4<amrex::EBCellFlag const> const& flag,
              int eb_mode=1) noexcept
{
  using amrex::Real;

  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

  if (flag(iv).isCovered()) return; //EB covered cells do not need flux

  // Fill covered stencil points
  amrex::IntVect ivm2 = iv-2*iv_dir;
  amrex::IntVect ivm1 = iv-1*iv_dir;
  amrex::IntVect ivp1 = iv+1*iv_dir;
  amrex::IntVect ivp2 = iv+2*iv_dir;
  // Left
  bool cn1 = flag(iv).isConnected(-iv_dir);
  bool cn2 = flag(iv-iv_dir).isConnected(-iv_dir);      
  ivm1 = cn1 ? ivm1 : iv;
  ivm2 = cn1 && cn2 ? ivm2 : ivm1;
  // Right
  bool cp1 = flag(iv).isConnected(iv_dir);
  bool cp2 = flag(iv+iv_dir).isConnected(iv_dir);    
  ivp1 = cp1 ? ivp1 : iv;
  ivp2 = cp1 && cp2 ? ivp2 : ivp1;

  // Switch to second order scheme near EB
  if ((eb_mode == 1) && (recon_scheme > 2)) {
    if ((ivm2 == ivm1) || (ivp2 == ivp1)) {
      recon_scheme = 2;
    }
  }
  
  if (recon_scheme == 1) {
    // Simple Godunov
    qr(iv,n) = q(iv,n);
    ql(iv+iv_dir,n) = q(iv,n);

  } else if (recon_scheme == 2) {
    // MUSCL    
    Real s[3];
    s[0] = q(ivm1,n);
    s[1] = q(iv,n);
    s[2] = q(ivp1,n);

    muscl(s, qr(iv,n), ql(iv+iv_dir,n), plm_theta);

  } else if (recon_scheme == 3) {
    // WENO-Z3
    Real s[3];
    s[0] = q(ivm1,n);
    s[1] = q(iv,n);
    s[2] = q(ivp1,n);

    wenoz3(s, qr(iv,n), ql(iv+iv_dir,n));

  } else if (recon_scheme == 4) {
    // WENO-JS5
    Real s[5];
    s[0] = q(ivm2,n);
    s[1] = q(ivm1,n);
    s[2] = q(iv,n);
    s[3] = q(ivp1,n);
    s[4] = q(ivp2,n);

    wenojs5(s, qr(iv,n), ql(iv+iv_dir,n));
  
  } else if (recon_scheme == 6) {
    // TENO-5
    Real s[5];
    s[0] = q(ivm2,n);
    s[1] = q(ivm1,n);
    s[2] = q(iv,n);
    s[3] = q(ivp1,n);
    s[4] = q(ivp2,n);

    teno5(s, qr(iv,n), ql(iv+iv_dir,n));
    
  } else {
    // Default: WENO-Z5
    Real s[5];
    s[0] = q(ivm2,n);
    s[1] = q(ivm1,n);
    s[2] = q(iv,n);
    s[3] = q(ivp1,n);
    s[4] = q(ivp2,n);

    wenoz5(s, qr(iv,n), ql(iv+iv_dir,n));
  }
}

#endif