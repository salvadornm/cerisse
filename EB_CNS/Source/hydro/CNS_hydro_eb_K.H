#ifndef CNS_HYDRO_EB_K_H_
#define CNS_HYDRO_EB_K_H_

#include <AMReX_FArrayBox.H>

#include "CNS.H"
#include "CNS_hydro_K.H"

#include <cmath>

// namespace {

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// amrex::Real limiter_eb (amrex::Real dlft, amrex::Real drgt, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

//     Real dcen = Real(0.5)*(dlft+drgt);
//     Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
//     Real slop = plm_theta * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
//     Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
//     return dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
// }

// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_eb_x (int i, int j, int k,
//                 amrex::Array4<amrex::Real> const& dq,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 int recon_scheme, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

//     Real dlft0 = Real(0.0); // rho
//     Real dlft1 = Real(0.0); // normal u
//     Real dlft2 = Real(0.0); // p
//     Real dlft3 = Real(0.0); // transverse u
//     Real dlft4 = Real(0.0); // transverse u
//     Real dlft5[NUM_SPECIES] = {Real(0.)}; // y

//     Real drgt0 = Real(0.0);
//     Real drgt1 = Real(0.0);
//     Real drgt2 = Real(0.0);
//     Real drgt3 = Real(0.0);
//     Real drgt4 = Real(0.0);
//     Real drgt5[NUM_SPECIES] = {Real(0.)}; // y

//     if (recon_scheme == 1) {
//         dq(i,j,k,0) = Real(0.);
//         dq(i,j,k,1) = Real(0.);
//         dq(i,j,k,2) = Real(0.);
//         dq(i,j,k,3) = Real(0.);
//         dq(i,j,k,4) = Real(0.);
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             dq(i,j,k,5+n) = Real(0.);
//         }
//     } else if (recon_scheme == 2) {
//         if (flag(i,j,k).isConnected(-1,0,0))
//         {
//             dlft0 = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
//             dlft1 = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO)) - (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
//             dlft2 = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
//             dlft3 = q(i,j,k,QV) - q(i-1,j,k,QV);
//             dlft4 = q(i,j,k,QW) - q(i-1,j,k,QW);
//             for (int n = 0; n < NUM_SPECIES; ++n)
//                 dlft5[n] = q(i,j,k,QFS+n) - q(i-1,j,k,QFS+n);
//         }        

//         if (flag(i,j,k).isConnected(1,0,0))
//         {
//             drgt0 = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
//             drgt1 = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO)) - (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
//             drgt2 = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
//             drgt3 = q(i+1,j,k,QV) - q(i,j,k,QV);
//             drgt4 = q(i+1,j,k,QW) - q(i,j,k,QW);
//             for (int n = 0; n < NUM_SPECIES; ++n)
//                 drgt5[n] = q(i+1,j,k,QFS+n) - q(i,j,k,QFS+n);
//         }

//         dq(i,j,k,0) = limiter(dlft0, drgt0, plm_theta);
//         dq(i,j,k,1) = limiter(dlft1, drgt1, plm_theta);
//         dq(i,j,k,2) = limiter(dlft2, drgt2, plm_theta);
//         dq(i,j,k,3) = limiter(dlft3, drgt3, plm_theta);
//         dq(i,j,k,4) = limiter(dlft4, drgt4, plm_theta);
//         for (int n = 0; n < NUM_SPECIES; ++n)
//             dq(i,j,k,5+n) = limiter(dlft5[n], drgt5[n], plm_theta);
//     }
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_eb_y (int i, int j, int k,
//                 amrex::Array4<amrex::Real> const& dq,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 int recon_scheme, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

//     Real dlft0 = Real(0.0);
//     Real dlft1 = Real(0.0);
//     Real dlft2 = Real(0.0);
//     Real dlft3 = Real(0.0);
//     Real dlft4 = Real(0.0);
//     Real dlft5[NUM_SPECIES] = {Real(0.)};

//     Real drgt0 = Real(0.0);
//     Real drgt1 = Real(0.0);
//     Real drgt2 = Real(0.0);
//     Real drgt3 = Real(0.0);
//     Real drgt4 = Real(0.0);
//     Real drgt5[NUM_SPECIES] = {Real(0.)};

//     if (recon_scheme == 1) {
//         dq(i,j,k,0) = Real(0.);
//         dq(i,j,k,1) = Real(0.);
//         dq(i,j,k,2) = Real(0.);
//         dq(i,j,k,3) = Real(0.);
//         dq(i,j,k,4) = Real(0.);
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             dq(i,j,k,5+n) = Real(0.);
//         }
//     } else if (recon_scheme == 2) {
//         if (flag(i,j,k).isConnected(0,-1,0))
//         {
//             dlft0 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
//             dlft1 = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO)) - (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
//             dlft2 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
//             dlft3 =q(i,j,k,QU) - q(i,j-1,k,QU);
//             dlft4 = q(i,j,k,QW) - q(i,j-1,k,QW);
//             for (int n = 0; n < NUM_SPECIES; ++n)
//                 dlft5[n] = q(i,j,k,QFS+n) - q(i,j-1,k,QFS+n);
//         }

//         if (flag(i,j,k).isConnected(0,1,0))
//         {
//             drgt0 = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
//             drgt1 = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO)) - (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
//             drgt2 = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
//             drgt3 = q(i,j+1,k,QU) - q(i,j,k,QU);
//             drgt4 = q(i,j+1,k,QW) - q(i,j,k,QW);
//             for (int n = 0; n < NUM_SPECIES; ++n)
//                 drgt5[n] = q(i,j+1,k,QFS+n) - q(i,j,k,QFS+n);
//         }

//         dq(i,j,k,0) = limiter(dlft0, drgt0, plm_theta);
//         dq(i,j,k,1) = limiter(dlft1, drgt1, plm_theta);
//         dq(i,j,k,2) = limiter(dlft2, drgt2, plm_theta);
//         dq(i,j,k,3) = limiter(dlft3, drgt3, plm_theta);
//         dq(i,j,k,4) = limiter(dlft4, drgt4, plm_theta);
//         for (int n = 0; n < NUM_SPECIES; ++n)
//             dq(i,j,k,5+n) = limiter(dlft5[n], drgt5[n], plm_theta);
//    }
// }

// #if (AMREX_SPACEDIM == 3)
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_eb_z (int i, int j, int k,
//                 amrex::Array4<amrex::Real> const& dq,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 int recon_scheme, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

//     Real dlft0 = Real(0.0);
//     Real dlft1 = Real(0.0);
//     Real dlft2 = Real(0.0);
//     Real dlft3 = Real(0.0);
//     Real dlft4 = Real(0.0);
//     Real dlft5[NUM_SPECIES] = {Real(0.)};

//     Real drgt0 = Real(0.0);
//     Real drgt1 = Real(0.0);
//     Real drgt2 = Real(0.0);
//     Real drgt3 = Real(0.0);
//     Real drgt4 = Real(0.0);
//     Real drgt5[NUM_SPECIES] = {Real(0.)};

//     if (recon_scheme == 1) {
//         dq(i,j,k,0) = Real(0.);
//         dq(i,j,k,1) = Real(0.);
//         dq(i,j,k,2) = Real(0.);
//         dq(i,j,k,3) = Real(0.);
//         dq(i,j,k,4) = Real(0.);
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             dq(i,j,k,5+n) = Real(0.);
//         }
//     } else if (recon_scheme == 2) {
//         if (flag(i,j,k).isConnected(0,0,-1))
//         {
//             dlft0 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
//             dlft1 = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO)) - (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
//             dlft2 = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
//             dlft3 = q(i,j,k,QU) - q(i,j,k-1,QU);
//             dlft4 = q(i,j,k,QV) - q(i,j,k-1,QV);
//             for (int n = 0; n < NUM_SPECIES; ++n)
//                 dlft5[n] = q(i,j,k,QFS+n) - q(i,j,k-1,QFS+n);
//         }

//         if (flag(i,j,k).isConnected(0,0,1))
//         {
//             drgt0 = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
//             drgt1 = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO)) - (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/(q(i,j,k,QCS)*q(i,j,k,QCS));
//             drgt2 = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
//             drgt3 = q(i,j,k+1,QU) - q(i,j,k,QU);
//             drgt4 = q(i,j,k+1,QV) - q(i,j,k,QV);
//             for (int n = 0; n < NUM_SPECIES; ++n)
//                 drgt5[n] = q(i,j,k+1,QFS+n) - q(i,j,k,QFS+n);
//         }

//         dq(i,j,k,0) = limiter(dlft0, drgt0, plm_theta);
//         dq(i,j,k,1) = limiter(dlft1, drgt1, plm_theta);
//         dq(i,j,k,2) = limiter(dlft2, drgt2, plm_theta);
//         dq(i,j,k,3) = limiter(dlft3, drgt3, plm_theta);
//         dq(i,j,k,4) = limiter(dlft4, drgt4, plm_theta);
//         for (int n = 0; n < NUM_SPECIES; ++n)
//             dq(i,j,k,5+n) = limiter(dlft5[n], drgt5[n], plm_theta);
//    }
// }
// #endif

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_slope_eb (int i, int j, int k, int dir,
              amrex::Array4<amrex::Real> const& dq,
              amrex::Array4<amrex::Real const> const& q,
              amrex::Array4<amrex::EBCellFlag const> const& flag,
              int recon_scheme, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    amrex::IntVect dir_iv = amrex::IntVect::TheDimensionVector(dir);

    int QUN, QUT, QUTT;
    if (dir == 0) {
        QUN = QU;  QUT = QV;  QUTT = QW;
    } else if (dir == 1) {
        QUN = QV;  QUT = QU;  QUTT = QW;
    } else {
        QUN = QW;  QUT = QU;  QUTT = QV;
    }

    Real dlft[5+NUM_SPECIES] = {Real(0.)};
    Real drgt[5+NUM_SPECIES] = {Real(0.)};

    if (recon_scheme == 1) {
        for (int n = 0; n < 5+NUM_SPECIES; ++n)
            dq(iv,n) = Real(0.);
        
    } else if (recon_scheme == 2) {
        if (flag(iv).isConnected(-dir_iv))
        {
            dlft[0] = q(iv,QRHO) - q(iv-dir_iv,QRHO);
            dlft[1] = q(iv,QUN) - q(iv-dir_iv,QUN);
            dlft[2] = q(iv,QPRES) - q(iv-dir_iv,QPRES);
            dlft[3] = q(iv,QUT) - q(iv-dir_iv,QUT);
            dlft[4] = q(iv,QUTT) - q(iv-dir_iv,QUTT);
            for (int n = 0; n < NUM_SPECIES; ++n)
                dlft[5+n] = q(iv,QFS+n) - q(iv-dir_iv,QFS+n);
        }

        if (flag(iv).isConnected(dir_iv))
        {
            drgt[0] = q(iv+dir_iv,QRHO) - q(iv,QRHO);
            drgt[1] = q(iv+dir_iv,QUN) - q(iv,QUN);
            drgt[2] = q(iv+dir_iv,QPRES) - q(iv,QPRES);
            drgt[3] = q(iv+dir_iv,QUT) - q(iv,QUT);
            drgt[4] = q(iv+dir_iv,QUTT) - q(iv,QUTT);
            for (int n = 0; n < NUM_SPECIES; ++n)
                drgt[5+n] = q(iv+dir_iv,QFS+n) - q(iv,QFS+n);
        }

        for (int n = 0; n < 5+NUM_SPECIES; ++n)
            dq(iv,n) = limiter(dlft[n], drgt[n], plm_theta);
   }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_eb_x (int i, int j, int k,
                amrex::Array4<amrex::Real> const& dq,
                amrex::Array4<amrex::Real const> const& q,
                amrex::Array4<amrex::EBCellFlag const> const& flag,
                int recon_scheme, amrex::Real plm_theta) noexcept
{
    cns_slope_eb(i, j, k, 0, dq, q, flag, recon_scheme, plm_theta);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_eb_y (int i, int j, int k,
                amrex::Array4<amrex::Real> const& dq,
                amrex::Array4<amrex::Real const> const& q,
                amrex::Array4<amrex::EBCellFlag const> const& flag,
                int recon_scheme, amrex::Real plm_theta) noexcept
{
    cns_slope_eb(i, j, k, 1, dq, q, flag, recon_scheme, plm_theta);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_eb_z (int i, int j, int k,
                amrex::Array4<amrex::Real> const& dq,
                amrex::Array4<amrex::Real const> const& q,
                amrex::Array4<amrex::EBCellFlag const> const& flag,
                int recon_scheme, amrex::Real plm_theta) noexcept
{
    cns_slope_eb(i, j, k, 2, dq, q, flag, recon_scheme, plm_theta);
}

#endif
