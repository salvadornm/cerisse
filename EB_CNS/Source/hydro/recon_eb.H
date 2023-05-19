#ifndef CNS_RECON_EB_H_
#define CNS_RECON_EB_H_

#include "recon.H"

/**
 * \brief Interpolate primitive variables to cell faces. EB version.
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
 * @param dir  direction.
 * @param q    cell averaged variable.
 * @param ql   left state at a face.
 * @param qr   right state at a face.
 * @param recon_scheme reconstruction scheme.
 * @param plt_theta a paramter in MUSCL limiter, between 1.0 and 2.0.
 * @param flag EB flags array.
 * @param eb_mode 0: fill covered stencil with values of the first cell from wall; 
 *                1: automatically switch to lower order scheme (not implemented).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_recon_eb (int i, int j, int k, int n, int dir,
              amrex::Array4<const amrex::Real> const& q,
              amrex::Array4<      amrex::Real> const& ql,
              amrex::Array4<      amrex::Real> const& qr,
              int recon_scheme, amrex::Real plm_theta,
              amrex::Array4<amrex::EBCellFlag const> const& flag,
              int eb_mode=0) noexcept
{
  using amrex::Real;

  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

  if (flag(iv).isCovered()) return;
  
  if (recon_scheme == 1) {
    // Simple Godunov
    qr(iv,n) = q(iv,n);
    ql(iv+iv_dir,n) = q(iv,n);

  } else if (recon_scheme == 2) {
    // MUSCL    
    // q(iv-iv_dir,n) = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // q(iv+iv_dir,n) = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n); 

    Real s[3];
    s[0] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    s[1] = q(iv,n);
    s[2] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);

    muscl(s, qr(iv,n), ql(iv+iv_dir,n), plm_theta);

  } else if (recon_scheme == 3) {
    // WENO-Z3
    Real s[3];
    s[0] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    s[1] = q(iv,n);
    s[2] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);

    wenoz3(s, qr(iv,n), ql(iv+iv_dir,n));

  } else if (recon_scheme == 4) {
    // WENO-JS5
    Real s[5];
    s[0] = flag(iv).isConnected(-iv_dir) ? 
            (flag(iv-iv_dir).isConnected(-iv_dir) ? 
              q(iv-2*iv_dir,n) : q(iv-iv_dir,n)) : q(iv,n);
    s[1] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    s[2] = q(iv,n);
    s[3] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    s[4] = flag(iv).isConnected(+iv_dir) ? 
            (flag(iv+iv_dir).isConnected(+iv_dir) ? 
              q(iv+2*iv_dir,n) : q(iv+iv_dir,n)) : q(iv,n);

    wenojs5(s, qr(iv,n), ql(iv+iv_dir,n));
  
  } else if (recon_scheme == 6) {
    // TENO-5
    Real s[5];    
    s[0] = flag(iv).isConnected(-iv_dir) ? 
            (flag(iv-iv_dir).isConnected(-iv_dir) ? 
              q(iv-2*iv_dir,n) : q(iv-iv_dir,n)) : q(iv,n);
    s[1] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    s[2] = q(iv,n);
    s[3] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    s[4] = flag(iv).isConnected(+iv_dir) ? 
            (flag(iv+iv_dir).isConnected(+iv_dir) ? 
              q(iv+2*iv_dir,n) : q(iv+iv_dir,n)) : q(iv,n);

    teno5(s, qr(iv,n), ql(iv+iv_dir,n));
    
  } else {
    // Default: WENO-Z5
    Real s[5];    
    s[0] = flag(iv).isConnected(-iv_dir) ? 
            (flag(iv-iv_dir).isConnected(-iv_dir) ? 
              q(iv-2*iv_dir,n) : q(iv-iv_dir,n)) : q(iv,n);
    s[1] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    s[2] = q(iv,n);
    s[3] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    s[4] = flag(iv).isConnected(+iv_dir) ? 
            (flag(iv+iv_dir).isConnected(+iv_dir) ? 
              q(iv+2*iv_dir,n) : q(iv+iv_dir,n)) : q(iv,n);

    wenoz5(s, qr(iv,n), ql(iv+iv_dir,n));
  }
}

#endif