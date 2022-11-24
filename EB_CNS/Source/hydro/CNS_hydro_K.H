#ifndef CNS_HYDRO_K_H_
#define CNS_HYDRO_K_H_

#include <AMReX_FArrayBox.H>

#include "index_macros.H"
#include "parm.H"
#include "PelePhysics.H"

#include <cmath>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_ctoprim (int i, int j, int k,
             amrex::Array4<amrex::Real const> const& u,
             amrex::Array4<amrex::Real> const& q,
             Parm const& parm) noexcept
{
    amrex::Real rho = u(i,j,k,URHO);
    amrex::Real rhoinv = amrex::Real(1.0) / rho;
    AMREX_D_TERM(amrex::Real ux = u(i,j,k,UMX)*rhoinv; ,
                 amrex::Real uy = u(i,j,k,UMY)*rhoinv; ,
                 amrex::Real uz = u(i,j,k,UMZ)*rhoinv;);
    amrex::Real kineng = amrex::Real(0.5)*rho*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));
    amrex::Real ei = (u(i,j,k,UEDEN) - kineng) * rhoinv;
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
        q(i,j,k,QFS+n) = u(i,j,k,UFS+n) * rhoinv;
        massfrac[n] = q(i,j,k,QFS+n);
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T, p, cs, gamma;
    // eos.Y2WBAR(massfrac, wbar);
    eos.EY2T(ei, massfrac, T);
    T = amrex::max(T, 1e-10);
    eos.RTY2P(rho, T, massfrac, p);
    // eos.RTY2Cs(rho, T, massfrac, cs);
    // eos.RTY2G(rho, T, massfrac, gamma);
    // eos.RTY2dpde_dpdre(rho, T, massfrac, dpde, dpdr_e);

    q(i,j,k,QRHO) = rho;
    AMREX_D_TERM(q(i,j,k,QU) = ux; ,
                 q(i,j,k,QV) = uy; ,
                 q(i,j,k,QW) = uz;);
    q(i,j,k,QEINT) = ei;
    q(i,j,k,QPRES) = p;
    // q(i,j,k,QCS) = cs;
    q(i,j,k,QTEMP) = T; 
    // q(i,j,k,QG) = gamma;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_flux_to_dudt (int i, int j, int k, int n,
                  amrex::Array4<amrex::Real> const& dudt,
                  AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fx,
                               amrex::Array4<amrex::Real const> const& fy,
                               amrex::Array4<amrex::Real const> const& fz),
                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    dudt(i,j,k,n) = AMREX_D_TERM(dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n)),
                               + dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n)),
                               + dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n)));
}

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
limiter (amrex::Real dlft, amrex::Real drgt, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    Real dcen = Real(0.5)*(dlft+drgt);
    Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
    Real slop = plm_theta * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
    Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
}

}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_x (int i, int j, int k,
//              amrex::Array4<amrex::Real> const& dq,
//              amrex::Array4<amrex::Real const> const& q,
//              int recon_scheme, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

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
//         Real dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
//         Real drgt = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
//         dq(i,j,k,0) = limiter(dlft, drgt, plm_theta);

//         Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
//         dlft = (q(i,j,k,QRHO)-q(i-1,j,k,QRHO)) - (q(i,j,k,QPRES) - q(i-1,j,k,QPRES))/cs2;
//         drgt = (q(i+1,j,k,QRHO)-q(i,j,k,QRHO)) - (q(i+1,j,k,QPRES) - q(i,j,k,QPRES))/cs2;
//         dq(i,j,k,1) = limiter(dlft, drgt, plm_theta);

//         dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i-1,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QU) - q(i-1,j,k,QU));
//         drgt = Real(0.5)*(q(i+1,j,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i+1,j,k,QU) - q(i,j,k,QU));
//         dq(i,j,k,2) = limiter(dlft, drgt, plm_theta);

//         dlft = q(i,j,k,QV) - q(i-1,j,k,QV);
//         drgt = q(i+1,j,k,QV) - q(i,j,k,QV);
//         dq(i,j,k,3) = limiter(dlft, drgt, plm_theta);

//         dlft = q(i,j,k,QW) - q(i-1,j,k,QW);
//         drgt = q(i+1,j,k,QW) - q(i,j,k,QW);
//         dq(i,j,k,4) = limiter(dlft, drgt, plm_theta);

//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             dlft = q(i,j,k,QFS+n) - q(i-1,j,k,QFS+n);
//             drgt = q(i+1,j,k,QFS+n) - q(i,j,k,QFS+n);
//             dq(i,j,k,5+n) = limiter(dlft, drgt, plm_theta);
//         }
//     }
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_y (int i, int j, int k,
//              amrex::Array4<amrex::Real> const& dq,
//              amrex::Array4<amrex::Real const> const& q,
//              int recon_scheme, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

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
//         Real dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
//         Real drgt = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
//         dq(i,j,k,0) = limiter(dlft, drgt, plm_theta);

//         Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
//         dlft = (q(i,j,k,QRHO)-q(i,j-1,k,QRHO)) - (q(i,j,k,QPRES) - q(i,j-1,k,QPRES))/cs2;
//         drgt = (q(i,j+1,k,QRHO)-q(i,j,k,QRHO)) - (q(i,j+1,k,QPRES) - q(i,j,k,QPRES))/cs2;
//         dq(i,j,k,1) = limiter(dlft, drgt, plm_theta);

//         dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j-1,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QV) - q(i,j-1,k,QV));
//         drgt = Real(0.5)*(q(i,j+1,k,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j+1,k,QV) - q(i,j,k,QV));
//         dq(i,j,k,2) = limiter(dlft, drgt, plm_theta);

//         dlft = q(i,j,k,QU) - q(i,j-1,k,QU);
//         drgt = q(i,j+1,k,QU) - q(i,j,k,QU);
//         dq(i,j,k,3) = limiter(dlft, drgt, plm_theta);
        
//         dlft = q(i,j,k,QW) - q(i,j-1,k,QW);
//         drgt = q(i,j+1,k,QW) - q(i,j,k,QW);
//         dq(i,j,k,4) = limiter(dlft, drgt, plm_theta);

//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             dlft = q(i,j,k,QFS+n) - q(i,j-1,k,QFS+n);
//             drgt = q(i,j+1,k,QFS+n) - q(i,j,k,QFS+n);
//             dq(i,j,k,5+n) = limiter(dlft, drgt, plm_theta);
//         }
//     }
// }

// #if (AMREX_SPACEDIM == 3)
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_slope_z (int i, int j, int k,
//              amrex::Array4<amrex::Real> const& dq,
//              amrex::Array4<amrex::Real const> const& q,
//              int recon_scheme, amrex::Real plm_theta) noexcept
// {
//     using amrex::Real;

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
//         Real dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
//         Real drgt = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) - Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
//         dq(i,j,k,0) = limiter(dlft, drgt, plm_theta);

//         Real cs2 = q(i,j,k,QCS)*q(i,j,k,QCS);
//         dlft = (q(i,j,k,QRHO)-q(i,j,k-1,QRHO)) - (q(i,j,k,QPRES) - q(i,j,k-1,QPRES))/cs2;
//         drgt = (q(i,j,k+1,QRHO)-q(i,j,k,QRHO)) - (q(i,j,k+1,QPRES) - q(i,j,k,QPRES))/cs2;
//         dq(i,j,k,1) = limiter(dlft, drgt, plm_theta);

//         dlft = Real(0.5)*(q(i,j,k,QPRES)-q(i,j,k-1,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k,QW) - q(i,j,k-1,QW));
//         drgt = Real(0.5)*(q(i,j,k+1,QPRES)-q(i,j,k,QPRES))/q(i,j,k,QCS) + Real(0.5)*q(i,j,k,QRHO)*(q(i,j,k+1,QW) - q(i,j,k,QW));
//         dq(i,j,k,2) = limiter(dlft, drgt, plm_theta);

//         dlft = q(i,j,k,QU) - q(i,j,k-1,QU);
//         drgt = q(i,j,k+1,QU) - q(i,j,k,QU);
//         dq(i,j,k,3) = limiter(dlft, drgt, plm_theta);

//         dlft = q(i,j,k,QV) - q(i,j,k-1,QV);
//         drgt = q(i,j,k+1,QV) - q(i,j,k,QV);
//         dq(i,j,k,4) = limiter(dlft, drgt, plm_theta);

//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             dlft = q(i,j,k,QFS+n) - q(i,j,k-1,QFS+n);
//             drgt = q(i,j,k+1,QFS+n) - q(i,j,k,QFS+n);
//             dq(i,j,k,5+n) = limiter(dlft, drgt, plm_theta);
//         }
//     }
// }
// #endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope (int i, int j, int k, int dir,
           amrex::Array4<amrex::Real> const& dq,
           amrex::Array4<amrex::Real const> const& q,
           int recon_scheme, amrex::Real plm_theta) noexcept
{
    using amrex::Real;

    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect dir_iv = amrex::IntVect::TheDimensionVector(dir);

    int QUN, QUT, QUTT;
    if (dir == 0) {
        QUN = QU;  QUT = QV;  QUTT = QW;
    } else if (dir == 1) {
        QUN = QV;  QUT = QU;  QUTT = QW;
    } else {
        QUN = QW;  QUT = QU;  QUTT = QV;
    }

    if (recon_scheme == 1) {
        for (int n = 0; n < NVAR; ++n) {
            dq(iv,n) = Real(0.);
        }
    } else if (recon_scheme == 2) {
        // Real dlft = Real(0.5)*(q(iv,QPRES)-q(iv-dir_iv,QPRES))/q(iv,QCS) - Real(0.5)*q(iv,QRHO)*(q(iv,QUN) - q(iv-dir_iv,QUN));
        // Real drgt = Real(0.5)*(q(iv+dir_iv,QPRES)-q(iv,QPRES))/q(iv,QCS) - Real(0.5)*q(iv,QRHO)*(q(iv+dir_iv,QUN) - q(iv,QUN));
        // dq(iv,0) = limiter(dlft, drgt, plm_theta);

        // Real cs2 = q(iv,QCS)*q(iv,QCS);
        // dlft = (q(iv,QRHO)-q(iv-dir_iv,QRHO)) - (q(iv,QPRES) - q(iv-dir_iv,QPRES))/cs2;
        // drgt = (q(iv+dir_iv,QRHO)-q(iv,QRHO)) - (q(iv+dir_iv,QPRES) - q(iv,QPRES))/cs2;
        // dq(iv,1) = limiter(dlft, drgt, plm_theta);

        // dlft = Real(0.5)*(q(iv,QPRES)-q(iv-dir_iv,QPRES))/q(iv,QCS) + Real(0.5)*q(iv,QRHO)*(q(iv,QUN) - q(iv-dir_iv,QUN));
        // drgt = Real(0.5)*(q(iv+dir_iv,QPRES)-q(iv,QPRES))/q(iv,QCS) + Real(0.5)*q(iv,QRHO)*(q(iv+dir_iv,QUN) - q(iv,QUN));
        // dq(iv,2) = limiter(dlft, drgt, plm_theta);

        Real dlft = q(iv,QRHO) - q(iv-dir_iv,QRHO);
        Real drgt = q(iv+dir_iv,QRHO) - q(iv,QRHO);
        dq(iv,0) = limiter(dlft, drgt, plm_theta);

        dlft = q(iv,QUN) - q(iv-dir_iv,QUN);
        drgt = q(iv+dir_iv,QUN) - q(iv,QUN);
        dq(iv,1) = limiter(dlft, drgt, plm_theta);

        dlft = q(iv,QPRES) - q(iv-dir_iv,QPRES);
        drgt = q(iv+dir_iv,QPRES) - q(iv,QPRES);
        dq(iv,2) = limiter(dlft, drgt, plm_theta);

        dlft = q(iv,QUT) - q(iv-dir_iv,QUT);
        drgt = q(iv+dir_iv,QUT) - q(iv,QUT);
        dq(iv,3) = limiter(dlft, drgt, plm_theta);
            
        dlft = q(iv,QUTT) - q(iv-dir_iv,QUTT);
        drgt = q(iv+dir_iv,QUTT) - q(iv,QUTT);
        dq(iv,4) = limiter(dlft, drgt, plm_theta);

        for (int n = 0; n < NUM_SPECIES; ++n) {
            dlft = q(iv,QFS+n) - q(iv-dir_iv,QFS+n);
            drgt = q(iv+dir_iv,QFS+n) - q(iv,QFS+n);
            dq(iv,5+n) = limiter(dlft, drgt, plm_theta);
        }

        // for (int n = 0; n < NPRIM; ++n) {
        //     dlft = q(iv,n) - q(iv-dir_iv,n);
        //     drgt = q(iv+dir_iv,n) - q(iv,n);
        //     dq(iv,n) = limiter(dlft, drgt, plm_theta);
        // }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_x (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q,
             int recon_scheme, amrex::Real plm_theta) noexcept
{
    cns_slope(i, j, k, 0, dq, q, recon_scheme, plm_theta);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_y (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q,
             int recon_scheme, amrex::Real plm_theta) noexcept
{
    cns_slope(i, j, k, 1, dq, q, recon_scheme, plm_theta);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_slope_z (int i, int j, int k,
             amrex::Array4<amrex::Real> const& dq,
             amrex::Array4<amrex::Real const> const& q,
             int recon_scheme, amrex::Real plm_theta) noexcept
{
    cns_slope(i, j, k, 2, dq, q, recon_scheme, plm_theta);
}

namespace {

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// riemann (const amrex::Real gamma, const amrex::Real smallp, const amrex::Real /*smallr*/,
//          const amrex::Real rl, const amrex::Real ul, const amrex::Real pl, 
//          const amrex::Real ut1l, const amrex::Real ut2l, const amrex::Real yl[NUM_SPECIES],
//          const amrex::Real rr, const amrex::Real ur, const amrex::Real pr, 
//          const amrex::Real ut1r, const amrex::Real ut2r, const amrex::Real yr[NUM_SPECIES],
//          amrex::Real& flxrho, 
//          AMREX_D_DECL(amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt),
//          amrex::Real& flxe, 
//          amrex::Real flxrhoy[NUM_SPECIES]) noexcept
// {
//     using amrex::Real;
    
// #if (AMREX_SPACEDIM == 1)
//     amrex::ignore_unused(ut1l,ut1r);
// #endif
// #if (AMREX_SPACEDIM <= 2)
//     amrex::ignore_unused(ut2l,ut2r);
// #endif

//     constexpr Real weakwv = Real(1.e-3);
//     constexpr Real small = Real(1.e-6);

//     Real clsql = gamma*pl*rl;
//     Real clsqr = gamma*pr*rr;
//     Real wl = std::sqrt(clsql);
//     Real wr = std::sqrt(clsqr);
//     Real cleft = wl/rl;
//     Real cright = wr/rr;
//     Real ccsmall = small*(cleft+cright);

//     Real pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
//     pstar = amrex::max(pstar,smallp);
//     Real pstnm1 = pstar;

//     Real wlsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pl)+pstar)*rl;
//     Real wrsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pr)+pstar)*rr;

//     wl = std::sqrt(wlsq);
//     wr = std::sqrt(wrsq);
//     Real ustarp = ul - (pstar-pl)/wl;
//     Real ustarm = ur + (pstar-pr)/wr;

//     pstar = (wl*pr + wr*pl - wr*wl*(ur-ul))/(wl+wr);
//     pstar = amrex::max(pstar,smallp);

//     Real ustar;
//     for (int iter = 0; iter < 3; ++iter)
//     {
//         wlsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pl)+pstar)*rl;
//         wrsq = (Real(0.5)*(gamma-Real(1.))*(pstar+pr)+pstar)*rr;

//         wl = Real(1.)/std::sqrt(wlsq);
//         wr = Real(1.)/std::sqrt(wrsq);

//         Real ustnm1 = ustarm;
//         Real ustnp1 = ustarp;

//         ustarm = ur - (pr - pstar)*wr;
//         ustarp = ul + (pl - pstar)*wl;

//         Real dpditer = amrex::Math::abs(pstnm1-pstar);
//         Real zp = amrex::Math::abs(ustarp-ustnp1);
//         if (zp-weakwv*cleft < Real(0.0) ) {
//             zp = dpditer*wl;
//         }
//         Real zm = amrex::Math::abs(ustarm-ustnm1);
//         if (zm-weakwv*cright < Real(0.0) ) {
//             zm = dpditer*wr;
//         }

//         Real zz = zp+zm;
//         Real denom = dpditer/ amrex::max(zz,ccsmall);
//         pstnm1 = pstar;
//         pstar = pstar - denom*(ustarm-ustarp);
//         pstar = amrex::max(pstar,smallp);
//         ustar = Real(0.5)*(ustarm+ustarp);
//     }

//     Real ro, uo, po, sgnm AMREX_D_DECL(, utrans1, utrans2), yo[NUM_SPECIES];
//     if (ustar > Real(0.)) {
//         ro = rl;
//         uo = ul;
//         po = pl;
//         sgnm = Real(1.);
//         AMREX_D_TERM(, 
//             utrans1 = ut1l;,
//             utrans2 = ut2l;);
//         for (int n = 0; n < NUM_SPECIES; ++n) yo[n] = yl[n];
//     } else if (ustar < Real(0.)) {
//         ro = rr;
//         uo = ur;
//         po = pr;
//         sgnm = Real(-1.);
//         AMREX_D_TERM(, 
//             utrans1 = ut1r;,
//             utrans2 = ut2r;);
//         for (int n = 0; n < NUM_SPECIES; ++n) yo[n] = yr[n];
//     } else {
//         uo = Real(0.5)*(ur+ul);
//         po = Real(0.5)*(pr+pl);
//         ro = Real(2.)*(rl*rr)/(rl+rr);
//         sgnm = Real(1.);
//         AMREX_D_TERM(, 
//             utrans1 = Real(0.5)*(ut1l+ut1r);,
//             utrans2 = Real(0.5)*(ut2l+ut2r););
//         for (int n = 0; n < NUM_SPECIES; ++n) yo[n] = 0.5*(yl[n]+yr[n]);
//     }
//     Real wosq = (Real(0.5)*(gamma-Real(1.))*(pstar+po)+pstar)*ro;
//     Real co = std::sqrt(gamma * po / ro);
//     Real wo = std::sqrt(wosq);
//     Real dpjmp = pstar-po;
//     Real rstar = ro/(Real(1.)-ro*dpjmp/wosq);
//     Real cstar = std::sqrt(gamma * pstar / rstar);
//     Real spout = co-sgnm*uo;
//     Real spin = cstar - sgnm*uo;
//     if(pstar >= po) {
//         spin = wo/ro-sgnm*uo;
//         spout = spin;
//     }
//     Real ss = amrex::max(spout-spin, spout+spin);
//     Real frac = Real(0.5)*(Real(1.)+(spin+spout)/amrex::max(ss,ccsmall));
//     //const amrex::Real drho = (pstar - po) / (co * co);

//     Real rgdnv, ugdnv, pgdnv, ygdnv[NUM_SPECIES];
//     if (spout < Real(0.)) {
//         rgdnv = ro;
//         ugdnv = uo;
//         pgdnv = po;
//     } else if(spin >= Real(0.)) {
//         rgdnv = rstar;
//         ugdnv = ustar;
//         pgdnv = pstar;
//     } else {
//         rgdnv = frac*rstar + (Real(1.) - frac)* ro;
//         ugdnv = frac*ustar + (Real(1.) - frac)* uo;
//         pgdnv = frac*pstar + (Real(1.) - frac)* po;
//     }

//     flxrho = rgdnv*ugdnv;
//     AMREX_D_TERM(
//         flxu   = rgdnv*ugdnv*ugdnv + pgdnv;,
//         flxut  = rgdnv*ugdnv*utrans1;,
//         flxutt = rgdnv*ugdnv*utrans2;);
//     flxe = ugdnv*(Real(0.5)*rgdnv*(AMREX_D_TERM(ugdnv*ugdnv,+utrans1*utrans1,+utrans2*utrans2)) 
//           +pgdnv/(gamma -Real(1.)) + pgdnv);
//     for (int n = 0; n < NUM_SPECIES; ++n) 
//         flxrhoy[n] = rgdnv*ugdnv*yo[n];
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
riemann (const amrex::Real rl, 
         AMREX_D_DECL(const amrex::Real ul, const amrex::Real ut1l, const amrex::Real ut2l), 
         const amrex::Real pl, const amrex::Real cl, const amrex::Real el,
         const amrex::Real yl[NUM_SPECIES],
         const amrex::Real rr, 
         AMREX_D_DECL(const amrex::Real ur, const amrex::Real ut1r, const amrex::Real ut2r), 
         const amrex::Real pr, const amrex::Real cr, const amrex::Real er,
         const amrex::Real yr[NUM_SPECIES],
         amrex::Real& flxrho, 
         AMREX_D_DECL(amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt),
         amrex::Real& flxrhoe, 
         amrex::Real flxrhoy[NUM_SPECIES]) noexcept
{
    // HLLC Riemann solver   

    using amrex::Real;

    // Estimate wave speeds
    // Arithmetic average (Davis)
    Real sl = ul - cl;
    Real sr = ur + cr;
    // Roe-averaged speeds (Einfeldt)
    Real rp = std::sqrt(rr / rl);
    Real uroe = (ul + ur*rp)/(1. + rp);
    Real croe = (cl + cr*rp)/(1. + rp) + 0.5*rp/(1. + rp)/(1. + rp)
               *(AMREX_D_TERM((ul-ur)*(ul-ur), +(ut1l-ut1r)*(ut1l-ut1r), +(ut2l-ut2r)*(ut2l-ut2r)));
    sl = amrex::min(sl, uroe - croe);
    sr = amrex::max(sr, uroe + croe);

    if (sl >= 0) {
        //flx_l
        flxrho = 0.;
        AMREX_D_TERM(
            flxu = rl*ul*ul + pl; ,
            flxut = rl*ul*ut1l; ,
            flxutt = rl*ul*ut2l;);
        flxrhoe = ul*(rl*el + pl);
        for (int n = 0; n < NUM_SPECIES; ++n) {
            flxrhoy[n] = rl*ul*yl[n];
            flxrho += flxrhoy[n];
        }

    } else if (sr <= 0) {
        //flx_r
        flxrho = 0.;
        AMREX_D_TERM(
            flxu = rr*ur*ur + pr; ,
            flxut = rr*ur*ut1r; ,
            flxutt = rr*ur*ut2r;);
        flxrhoe = ur*(rr*er + pr);
        for (int n = 0; n < NUM_SPECIES; ++n) {
            flxrhoy[n] = rr*ur*yr[n];
            flxrho += flxrhoy[n];
        }

    } else {
        Real sstar = (pr - pl + rl*ul*(sl-ul) - rr*ur*(sr-ur)) / (rl*(sl-ul) - rr*(sr-ur)); // contact wave speed
        if (sstar >= 0) {
            //flx_l* = flx_l + sl * (q_l* - q_l)
            Real frac = (sl - ul) / (sl - sstar) - 1.;
            flxrho = 0.;
            AMREX_D_TERM(
                flxu   = rl*ul*ul + pl + sl*rl*((frac + 1.)*sstar - ul); ,
                flxut  = rl*ul*ut1l + sl*rl*frac*ut1l; ,
                flxutt = rl*ul*ut2l + sl*rl*frac*ut2l;);
            flxrhoe = ul*(rl*el + pl) 
                   + sl*rl*(frac*el + (sstar - ul)*(sstar + pl/rl/(sl-ul)));
            for (int n = 0; n < NUM_SPECIES; ++n) {
                flxrhoy[n] = rl*ul*yl[n] + sl*rl*frac*yl[n];
                flxrho += flxrhoy[n];
            }

        } else {
            //flx_r* = flx_r + sr * (q_r* - q_r)
            Real frac = (sr - ur) / (sr - sstar) - 1.;
            flxrho = 0.;
            AMREX_D_TERM(
                flxu   = rr*ur*ur + pr + sr*rr*((frac+1.)*sstar - ur); ,
                flxut  = rr*ur*ut1r + sr*rr*frac*ut1r; ,
                flxutt = rr*ur*ut2r + sr*rr*frac*ut2r;);
            flxrhoe = ur*(rr*er + pr)
                   + sr*rr*(frac*er + (sstar - ur)*(sstar + pr/rr/(sr-ur)));
            for (int n = 0; n < NUM_SPECIES; ++n) {
                flxrhoy[n] = rr*ur*yr[n] + sr*rr*frac*yr[n];
                flxrho += flxrhoy[n];
            }
        }
    }    
}

#if (AMREX_SPACEDIM > 1)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
compute_hyp_wallflux (const amrex::Real rho,
                      AMREX_D_DECL(const amrex::Real u, const amrex::Real v, const amrex::Real w),
                      const amrex::Real p,   const amrex::Real sp[NUM_SPECIES],
                      const amrex::Real axm, const amrex::Real axp,
                      const amrex::Real aym, const amrex::Real ayp,
#if (AMREX_SPACEDIM == 3)
                      const amrex::Real azm, const amrex::Real azp,
#endif
                      amrex::GpuArray<amrex::Real,NVAR>& fluxw, Parm const& parm) noexcept
{
    using amrex::Real;

    Real apnorm = std::sqrt( AMREX_D_TERM((axm-axp)*(axm-axp), + (aym-ayp)*(aym-ayp), + (azm-azp)*(azm-azp)) );
    Real apnorminv = 1./apnorm;
    Real un = AMREX_D_TERM(u*(axm-axp)*apnorminv, + v*(aym-ayp)*apnorminv, + w*(azm-azp)*apnorminv);
    Real AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2);
    Real y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) y[n] = sp[n];
    
    auto eos = pele::physics::PhysicsType::eos();
    Real e, cs;
    eos.RYP2E(rho, y, p, e);
    e += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
    eos.RPY2Cs(rho, p, y, cs);

    riemann(rho, AMREX_D_DECL( un, 0., 0.), p, cs, e, y,
            rho, AMREX_D_DECL(-un, 0., 0.), p, cs, e, y,
            fluxw[URHO], AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2), fluxw[UEDEN], &fluxw[UFS]);

    AMREX_D_TERM(fluxw[UMX] = (axm-axp) * momfluxn;,
                 fluxw[UMY] = (aym-ayp) * momfluxn;,
                 fluxw[UMZ] = (azm-azp) * momfluxn;);
}
#endif

} //namespace

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_riemann_x (int i, int j, int k,
//                amrex::Array4<amrex::Real> const& fx,
//                amrex::Array4<amrex::Real const> const& dq,
//                amrex::Array4<amrex::Real const> const& q,
//                Parm const& parm) noexcept
// {
//     using amrex::Real;

//     Real cspeed = q(i-1,j,k,QCS);
//     Real rl = q(i-1,j,k,QRHO) + Real(0.5) * ( (dq(i-1,j,k,0)+dq(i-1,j,k,2))/cspeed + dq(i-1,j,k,1));
//     rl = amrex::max(rl, parm.smallr);
//     Real ul = q(i-1,j,k,QU) + Real(0.5) * ( (dq(i-1,j,k,2)-dq(i-1,j,k,0))/q(i-1,j,k,QRHO));
//     Real pl = q(i-1,j,k,QPRES) + Real(0.5) *  (dq(i-1,j,k,0)+dq(i-1,j,k,2))*cspeed;
//     pl = amrex::max(pl, parm.smallp);
// #if (AMREX_SPACEDIM == 1)
//     Real ut1l = 0;
// #else
//     Real ut1l = q(i-1,j,k,QV) + Real(0.5) * dq(i-1,j,k,3);
// #endif
// #if (AMREX_SPACEDIM <= 2)
//     Real ut2l = 0.;
// #else
//     Real ut2l = q(i-1,j,k,QW) + Real(0.5) * dq(i-1,j,k,4);
// #endif
//     Real yl[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         yl[n] = q(i-1,j,k,UFS+n) + Real(0.5) * dq(i-1,j,k,5+n);
//     }

//     cspeed = q(i,j,k,QCS);
//     Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
//     rr = amrex::max(rr, parm.smallr);
//     Real ur = q(i,j,k,QU) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
//     Real pr = q(i,j,k,QPRES) - Real(0.5) * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
//     pr = amrex::max(pr, parm.smallp);
// #if (AMREX_SPACEDIM == 1)
//     Real ut1r = 0;
// #else
//     Real ut1r = q(i,j,k,QV) - Real(0.5) * dq(i,j,k,3);
// #endif
// #if (AMREX_SPACEDIM <= 2)
//     Real ut2r = 0.;
// #else
//     Real ut2r = q(i,j,k,QW) - Real(0.5) * dq(i,j,k,4);
// #endif
//     Real yr[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         yr[n] = q(i,j,k,UFS+n) - Real(0.5) * dq(i,j,k,5+n);
//     }

//     riemann(q(i,j,k,QG), parm.smallp, parm.smallr,
//             rl, ul, pl, ut1l, ut2l, yl,
//             rr, ur, pr, ut1r, ut2r, yr,
//             fx(i,j,k,URHO),
//             AMREX_D_DECL(fx(i,j,k,UMX), fx(i,j,k,UMY), fx(i,j,k,UMZ)), 
//             fx(i,j,k,UEDEN),
//             &fx(i,j,k,UFS));
// }

// #if (AMREX_SPACEDIM >= 2)
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_riemann_y (int i, int j, int k,
//                amrex::Array4<amrex::Real> const& fy,
//                amrex::Array4<amrex::Real const> const& dq,
//                amrex::Array4<amrex::Real const> const& q,
//                Parm const& parm) noexcept
// {
//     using amrex::Real;

//     Real cspeed = q(i,j-1,k,QCS);
//     Real rl = q(i,j-1,k,QRHO) + Real(0.5) * ( (dq(i,j-1,k,0)+dq(i,j-1,k,2))/cspeed + dq(i,j-1,k,1));
//     rl = amrex::max(rl, parm.smallr);
//     Real ul = q(i,j-1,k,QV) + Real(0.5) * ( (dq(i,j-1,k,2)-dq(i,j-1,k,0))/q(i,j-1,k,QRHO));
//     Real pl = q(i,j-1,k,QPRES) + Real(0.5) *  (dq(i,j-1,k,0)+dq(i,j-1,k,2))*cspeed;
//     pl = amrex::max(pl, parm.smallp);
//     Real ut1l = q(i,j-1,k,QU) + Real(0.5) * dq(i,j-1,k,3);
// #if (AMREX_SPACEDIM <= 2)
//     Real ut2l = 0.;
// #else
//     Real ut2l = q(i,j-1,k,QW) + Real(0.5) * dq(i,j-1,k,4);
// #endif
//     Real yl[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         yl[n] = q(i,j-1,k,UFS+n) + Real(0.5) * dq(i,j-1,k,5+n);
//     }

//     cspeed = q(i,j,k,QCS);
//     Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
//     rr = amrex::max(rr, parm.smallr);
//     Real ur = q(i,j,k,QV) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
//     Real pr = q(i,j,k,QPRES) - Real(0.5) * (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
//     pr = amrex::max(pr, parm.smallp);
//     Real ut1r = q(i,j,k,QU) - Real(0.5) * dq(i,j,k,3);
// #if (AMREX_SPACEDIM <= 2)
//     Real ut2r = 0.;
// #else
//     Real ut2r = q(i,j,k,QW) - Real(0.5) * dq(i,j,k,4);
// #endif
//     Real yr[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         yr[n] = q(i,j,k,UFS+n) - Real(0.5) * dq(i,j,k,5+n);
//     }

//     riemann(q(i,j,k,QG), parm.smallp, parm.smallr,
//             rl, ul, pl, ut1l, ut2l, yl,
//             rr, ur, pr, ut1r, ut2r, yr,
//             fy(i,j,k,URHO), 
//             AMREX_D_DECL(fy(i,j,k,UMY), fy(i,j,k,UMX), fy(i,j,k,UMZ)), 
//             fy(i,j,k,UEDEN),
//             &fy(i,j,k,UFS));
// }
// #endif

// #if (AMREX_SPACEDIM == 3)
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_riemann_z (int i, int j, int k,
//                amrex::Array4<amrex::Real> const& fz,
//                amrex::Array4<amrex::Real const> const& dq,
//                amrex::Array4<amrex::Real const> const& q,
//                Parm const& parm) noexcept
// {
//     using amrex::Real;

//     Real cspeed = q(i,j,k-1,QCS);
//     Real rl = q(i,j,k-1,QRHO) + Real(0.5) * ( (dq(i,j,k-1,0)+dq(i,j,k-1,2))/cspeed + dq(i,j,k-1,1));
//     rl = amrex::max(rl, parm.smallr);
//     Real ul = q(i,j,k-1,QW) + Real(0.5) * ( (dq(i,j,k-1,2)-dq(i,j,k-1,0))/q(i,j,k-1,QRHO));
//     Real pl = q(i,j,k-1,QPRES) + Real(0.5) *  (dq(i,j,k-1,0)+dq(i,j,k-1,2))*cspeed;
//     pl = amrex::max(pl, parm.smallp);
//     Real ut1l = q(i,j,k-1,QU) + Real(0.5) * dq(i,j,k-1,3);
//     Real ut2l = q(i,j,k-1,QV) + Real(0.5) * dq(i,j,k-1,4);
//     Real yl[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         yl[n] = q(i,j,k-1,UFS+n) + Real(0.5) * dq(i,j,k-1,5+n);
//     }

//     cspeed = q(i,j,k,QCS);
//     Real rr = q(i,j,k,QRHO) - Real(0.5) * ( (dq(i,j,k,0)+dq(i,j,k,2))/cspeed + dq(i,j,k,1));
//     rr = amrex::max(rr, parm.smallr);
//     Real ur = q(i,j,k,QW) - Real(0.5) * ( (dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO));
//     Real pr = q(i,j,k,QPRES) - Real(0.5) *  (dq(i,j,k,0)+dq(i,j,k,2))*cspeed;
//     pr = amrex::max(pr, parm.smallp);
// #if (AMREX_SPACEDIM == 1)
//     Real ut1r = 0.;
// #else
//     Real ut1r = q(i,j,k,QU) - Real(0.5) * dq(i,j,k,3);
// #endif
// #if (AMREX_SPACEDIM <= 2)
//     Real ut2r = 0.;
// #else
//     Real ut2r = q(i,j,k,QV) - Real(0.5) * dq(i,j,k,4);
// #endif

//     riemann(q(i,j,k,QG), parm.smallp, parm.smallr,
//             rl, ul, pl, ut1l, ut2l, yl,
//             rr, ur, pr, ut1r, ut2r, yr,
//             fz(i,j,k,URHO), 
//             fz(i,j,k,UMZ), fz(i,j,k,UMX), fz(i,j,k,UMY), 
//             fz(i,j,k,UEDEN),
//             &fz(i,j,k,UFS));
// }
// #endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_riemann (int i, int j, int k, int dir,
             amrex::Array4<amrex::Real> const& flx,
             amrex::Array4<amrex::Real const> const& dq,
             amrex::Array4<amrex::Real const> const& q,
             Parm const& parm) noexcept
{
    using amrex::Real;

    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

    int QUN, QUT, QUTT, UMN, UMT, UMTT;
    if (dir == 0) {
        QUN = QU;  QUT = QV;  QUTT = QW;
        UMN = UMX; UMT = UMY; UMTT = UMZ;  
    } else if (dir == 1) {
        QUN = QV;  QUT = QU;  QUTT = QW;
        UMN = UMY; UMT = UMX; UMTT = UMZ;  
    } else {
        QUN = QW;  QUT = QU;  QUTT = QV;
        UMN = UMZ; UMT = UMX; UMTT = UMY;  
    }
    
    Real cl;// = q(iv-iv_dir,QCS);
    Real rl = q(iv-iv_dir,QRHO) + Real(0.5) * dq(iv-iv_dir,0);//((dq(iv-iv_dir,0)+dq(iv-iv_dir,2))/cl + dq(iv-iv_dir,1));
    rl = amrex::max(rl, parm.smallr);
    AMREX_D_TERM(
        Real ul = q(iv-iv_dir,QUN) + Real(0.5) * dq(iv-iv_dir,1); //((dq(iv-iv_dir,2)-dq(iv-iv_dir,0))/q(iv-iv_dir,QRHO)); 
        Real ut1l = 0.; Real ut2l = 0.;
        ,
        ut1l = q(iv-iv_dir,QUT) + Real(0.5) * dq(iv-iv_dir,3); ,
        ut2l = q(iv-iv_dir,QUTT) + Real(0.5) * dq(iv-iv_dir,4););
    Real pl = q(iv-iv_dir,QPRES) + Real(0.5) * dq(iv-iv_dir,2); // (dq(iv-iv_dir,0)+dq(iv-iv_dir,2))*cl;
    pl = amrex::max(pl, parm.smallp);
    Real yl[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = q(iv-iv_dir,QFS+n) + Real(0.5) * dq(iv-iv_dir,5+n);
    }
    auto eos = pele::physics::PhysicsType::eos();
    Real el;
    eos.RPY2Cs(rl, pl, yl, cl);
    eos.RYP2E(rl, yl, pl, el);
    el += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));

    Real cr;// = q(i,j,k,QCS);
    Real rr = q(i,j,k,QRHO) - Real(0.5) * dq(i,j,k,0); //((dq(i,j,k,0)+dq(i,j,k,2))/cr + dq(i,j,k,1));
    rr = amrex::max(rr, parm.smallr);
    AMREX_D_TERM(
        Real ur = q(i,j,k,QUN) - Real(0.5) * dq(i,j,k,1); //((dq(i,j,k,2)-dq(i,j,k,0))/q(i,j,k,QRHO)); 
        Real ut1r = 0.; Real ut2r = 0.;
        ,
        ut1r = q(i,j,k,QUT) - Real(0.5) * dq(i,j,k,3); ,
        ut2r = q(i,j,k,QUTT) - Real(0.5) * dq(i,j,k,4););
    Real pr = q(i,j,k,QPRES) - Real(0.5) * dq(i,j,k,2); //(dq(i,j,k,0)+dq(i,j,k,2))*cr;
    pr = amrex::max(pr, parm.smallp);
    Real yr[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yr[n] = q(i,j,k,QFS+n) - Real(0.5) * dq(i,j,k,5+n);
    }
    Real er;
    eos.RPY2Cs(rr, pr, yr, cr);
    eos.RYP2E(rr, yr, pr, er);
    er += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));

    riemann(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
            rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
            flx(i,j,k,URHO), 
            AMREX_D_DECL(flx(i,j,k,UMN), flx(i,j,k,UMT), flx(i,j,k,UMTT)), 
            flx(i,j,k,UEDEN),
            &flx(i,j,k,UFS));
    // riemann(q(i,j,k,QG), parm.smallp, parm.smallr,
    //         rl, ul, pl, ut1l, ut2l, yl,
    //         rr, ur, pr, ut1r, ut2r, yr,
    //         flx(i,j,k,URHO), 
    //         AMREX_D_DECL(flx(i,j,k,UMN), flx(i,j,k,UMT), flx(i,j,k,UMTT)), 
    //         flx(i,j,k,UEDEN),
    //         &flx(i,j,k,UFS));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_riemann_x (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fx,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q,
               Parm const& parm) noexcept
{
    cns_riemann(i, j, k, 0, fx, dq, q, parm);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_riemann_y (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fy,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q,
               Parm const& parm) noexcept
{
    cns_riemann(i, j, k, 1, fy, dq, q, parm);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_riemann_z (int i, int j, int k,
               amrex::Array4<amrex::Real> const& fz,
               amrex::Array4<amrex::Real const> const& dq,
               amrex::Array4<amrex::Real const> const& q,
               Parm const& parm) noexcept
{
    cns_riemann(i, j, k, 2, fz, dq, q, parm);
}

#endif
