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
 *                1: automatically switch to lower order scheme.
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
    // q(iv-iv_dir,n) = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // q(iv+iv_dir,n) = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n); 

    Real s[3];
    // s[0] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // s[1] = q(iv,n);
    // s[2] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);

    s[0] = q(ivm1,n);
    s[1] = q(iv,n);
    s[2] = q(ivp1,n);

    muscl(s, qr(iv,n), ql(iv+iv_dir,n), plm_theta);

  } else if (recon_scheme == 3) {
    // WENO-Z3
    Real s[3];
    // s[0] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // s[1] = q(iv,n);
    // s[2] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    s[0] = q(ivm1,n);
    s[1] = q(iv,n);
    s[2] = q(ivp1,n);

    wenoz3(s, qr(iv,n), ql(iv+iv_dir,n));

  } else if (recon_scheme == 4) {
    // WENO-JS5
    Real s[5];
    // s[0] = flag(iv).isConnected(-iv_dir) ? 
    //         (flag(iv-iv_dir).isConnected(-iv_dir) ? 
    //           q(iv-2*iv_dir,n) : q(iv-iv_dir,n)) : q(iv,n);
    // s[1] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // s[2] = q(iv,n);
    // s[3] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    // s[4] = flag(iv).isConnected(+iv_dir) ? 
    //         (flag(iv+iv_dir).isConnected(+iv_dir) ? 
    //           q(iv+2*iv_dir,n) : q(iv+iv_dir,n)) : q(iv,n);
    s[0] = q(ivm2,n);
    s[1] = q(ivm1,n);
    s[2] = q(iv,n);
    s[3] = q(ivp1,n);
    s[4] = q(ivp2,n);

    wenojs5(s, qr(iv,n), ql(iv+iv_dir,n));
  
  } else if (recon_scheme == 6) {
    // TENO-5
    Real s[5];    
    // s[0] = flag(iv).isConnected(-iv_dir) ? 
    //         (flag(iv-iv_dir).isConnected(-iv_dir) ? 
    //           q(iv-2*iv_dir,n) : q(iv-iv_dir,n)) : q(iv,n);
    // s[1] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // s[2] = q(iv,n);
    // s[3] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    // s[4] = flag(iv).isConnected(+iv_dir) ? 
    //         (flag(iv+iv_dir).isConnected(+iv_dir) ? 
    //           q(iv+2*iv_dir,n) : q(iv+iv_dir,n)) : q(iv,n);
    s[0] = q(ivm2,n);
    s[1] = q(ivm1,n);
    s[2] = q(iv,n);
    s[3] = q(ivp1,n);
    s[4] = q(ivp2,n);

    teno5(s, qr(iv,n), ql(iv+iv_dir,n));
    
  } else {
    // Default: WENO-Z5
    Real s[5];    
    // s[0] = flag(iv).isConnected(-iv_dir) ? 
    //         (flag(iv-iv_dir).isConnected(-iv_dir) ? 
    //           q(iv-2*iv_dir,n) : q(iv-iv_dir,n)) : q(iv,n);
    // s[1] = flag(iv).isConnected(-iv_dir) ? q(iv-iv_dir,n) : q(iv,n);
    // s[2] = q(iv,n);
    // s[3] = flag(iv).isConnected(+iv_dir) ? q(iv+iv_dir,n) : q(iv,n);
    // s[4] = flag(iv).isConnected(+iv_dir) ? 
    //         (flag(iv+iv_dir).isConnected(+iv_dir) ? 
    //           q(iv+2*iv_dir,n) : q(iv+iv_dir,n)) : q(iv,n);
    s[0] = q(ivm2,n);
    s[1] = q(ivm1,n);
    s[2] = q(iv,n);
    s[3] = q(ivp1,n);
    s[4] = q(ivp2,n);

    wenoz5(s, qr(iv,n), ql(iv+iv_dir,n));
  }
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
// is_valid_cell (amrex::IntVect const& iv, 
//                amrex::IntVect const& blo, 
//                amrex::IntVect const& bhi) 
// {
//   return AMREX_D_TERM( (blo[0] <= iv[0]) && (iv[0] <= bhi[0]),
//                     && (blo[1] <= iv[1]) && (iv[1] <= bhi[1]),
//                     && (blo[2] <= iv[2]) && (iv[2] <= bhi[2]) );
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// cns_fill_eb_q (int i, int j, int k, int dir,
//               amrex::Array4<amrex::Real> const& q,
//               amrex::Array4<amrex::EBCellFlag const> const& flag,
//               amrex::IntVect const& blo, amrex::IntVect const& bhi,
//               bool no_slip) //noexcept
// {
//   const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//   const amrex::IntVect idir = amrex::IntVect::TheDimensionVector(dir);

//   if ( !flag(iv).isCovered() ) return;

//   // Find the nearest non-covered cell
//   const bool cm2 = is_valid_cell(iv-2*idir, blo, bhi) ? !flag(iv-2*idir).isCovered() : false;
//   const bool cm1 = is_valid_cell(iv-  idir, blo, bhi) ? !flag(iv-  idir).isCovered() : false;
//   const bool cp2 = is_valid_cell(iv+2*idir, blo, bhi) ? !flag(iv+2*idir).isCovered() : false;
//   const bool cp1 = is_valid_cell(iv+  idir, blo, bhi) ? !flag(iv+  idir).isCovered() : false;

//   // const bool cm2 = !flag(iv-2*idir).isCovered();
//   // const bool cm1 = !flag(iv-  idir).isCovered();
//   // const bool cp2 = !flag(iv+2*idir).isCovered();
//   // const bool cp1 = !flag(iv+  idir).isCovered();

//   // Position to copy from (the nearest non-covered cell)
//   amrex::IntVect iv_cpy;
//   if (cp1) {
//     iv_cpy = iv + idir;
//   } else if (cm1) {
//     iv_cpy = iv - idir;
//   } else if (cp2) {
//     iv_cpy = iv + 2*idir;
//   } else if (cm2) {
//     iv_cpy = iv - 2*idir;
//   } else {
//     return;
//   }

//   q(iv,QRHO) = q(iv_cpy,QRHO);
//   AMREX_D_TERM(
//   q(iv,QU) = ((no_slip || dir == 0) ? -1 : 1) * q(iv_cpy,QU);,
//   q(iv,QV) = ((no_slip || dir == 1) ? -1 : 1) * q(iv_cpy,QV);,
//   q(iv,QW) = ((no_slip || dir == 2) ? -1 : 1) * q(iv_cpy,QW);)
//   // AMREX_D_TERM(
//   // q(iv,QU) = q(iv_cpy,QU);,
//   // q(iv,QV) = q(iv_cpy,QV);,
//   // q(iv,QW) = q(iv_cpy,QW);)
//   q(iv,QPRES) = q(iv_cpy,QPRES);
//   q(iv,QC) = q(iv_cpy,QC);
//   q(iv,QG) = q(iv_cpy,QG);
//   q(iv,QEINT) = q(iv_cpy,QEINT);
//   q(iv,QTEMP) = q(iv_cpy,QTEMP);
//   for (int n = 0; n < NUM_SPECIES; ++n)
//     q(iv,QFS+n) = q(iv_cpy,QFS+n);
// }



#endif