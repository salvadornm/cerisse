// #ifndef CNS_HYDRO_EB_K_H_
// #define CNS_HYDRO_EB_K_H_

// #include <AMReX_FArrayBox.H>

// #include "CNS.H"
// #include "CNS_hydro_K.H"

// #include <cmath>

// AMREX_GPU_DEVICE
// AMREX_FORCE_INLINE
// void
// cns_slope_eb (int i, int j, int k, int dir,
//               amrex::Array4<amrex::Real> const& dq,
//               amrex::Array4<amrex::Real const> const& q,
//               amrex::Array4<amrex::EBCellFlag const> const& flag,
//               int recon_scheme, amrex::Real plm_theta) noexcept
// {
//   using amrex::Real;

//   amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//   amrex::IntVect dir_iv = amrex::IntVect::TheDimensionVector(dir);

//   int QUN, QUT, QUTT;
//   if (dir == 0) {
//     QUN = QU;  QUT = QV;  QUTT = QW;
//   } else if (dir == 1) {
//     QUN = QV;  QUT = QU;  QUTT = QW;
//   } else {
//     QUN = QW;  QUT = QU;  QUTT = QV;
//   }

//   Real dlft[5+NUM_SPECIES] = {Real(0.)};
//   Real drgt[5+NUM_SPECIES] = {Real(0.)};

//   if (recon_scheme == 1) {
//     for (int n = 0; n < 5+NUM_SPECIES; ++n)
//       dq(iv,n) = Real(0.);
      
//   } else if (recon_scheme == 2) {
//     if (flag(iv).isConnected(-dir_iv))
//     {
//       dlft[0] = q(iv,QRHO) - q(iv-dir_iv,QRHO);
//       dlft[1] = q(iv,QUN) - q(iv-dir_iv,QUN);
//       dlft[2] = q(iv,QPRES) - q(iv-dir_iv,QPRES);
//       dlft[3] = q(iv,QUT) - q(iv-dir_iv,QUT);
//       dlft[4] = q(iv,QUTT) - q(iv-dir_iv,QUTT);
//       for (int n = 0; n < NUM_SPECIES; ++n)
//         dlft[5+n] = q(iv,QFS+n) - q(iv-dir_iv,QFS+n);
//     }

//     if (flag(iv).isConnected(dir_iv))
//     {
//       drgt[0] = q(iv+dir_iv,QRHO) - q(iv,QRHO);
//       drgt[1] = q(iv+dir_iv,QUN) - q(iv,QUN);
//       drgt[2] = q(iv+dir_iv,QPRES) - q(iv,QPRES);
//       drgt[3] = q(iv+dir_iv,QUT) - q(iv,QUT);
//       drgt[4] = q(iv+dir_iv,QUTT) - q(iv,QUTT);
//       for (int n = 0; n < NUM_SPECIES; ++n)
//         drgt[5+n] = q(iv+dir_iv,QFS+n) - q(iv,QFS+n);
//     }

//     for (int n = 0; n < 5+NUM_SPECIES; ++n)
//       dq(iv,n) = limiter(dlft[n], drgt[n], plm_theta);
//   }
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_eb_x (int i, int j, int k,
//                 amrex::Array4<amrex::Real> const& dq,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 int recon_scheme, amrex::Real plm_theta) noexcept
// {
//   cns_slope_eb(i, j, k, 0, dq, q, flag, recon_scheme, plm_theta);
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_eb_y (int i, int j, int k,
//                 amrex::Array4<amrex::Real> const& dq,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 int recon_scheme, amrex::Real plm_theta) noexcept
// {
//   cns_slope_eb(i, j, k, 1, dq, q, flag, recon_scheme, plm_theta);
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_eb_z (int i, int j, int k,
//                 amrex::Array4<amrex::Real> const& dq,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 int recon_scheme, amrex::Real plm_theta) noexcept
// {
//   cns_slope_eb(i, j, k, 2, dq, q, flag, recon_scheme, plm_theta);
// }

// #endif