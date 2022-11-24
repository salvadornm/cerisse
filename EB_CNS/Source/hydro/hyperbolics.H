#ifndef CNS_HYPERBOLICS_H_
#define CNS_HYPERBOLICS_H_

#include <AMReX_FArrayBox.H>

#include "index_macros.H"
#include "parm.H"
#include "PelePhysics.H"
#include "recon.H"

#include <cmath>

/**
 * \brief Convert state variable to primitive variable.
 * 
 * @param i    x index.
 * @param j    y index.
 * @param k    z index.
 * @param u    state variables array.
 * @param q    primitive variables array.
 * @param parm parameters and constants (see parm.H).
 */
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
        massfrac[n] = u(i,j,k,UFS+n) * rhoinv;
        q(i,j,k,QFS+n) = massfrac[n];
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T, p, cs;
    // eos.Y2WBAR(massfrac, wbar);
    eos.EY2T(ei, massfrac, T); 
    eos.RTY2P(rho, T, massfrac, p);
    eos.RTY2Cs(rho, T, massfrac, cs);
    // eos.RTY2G(rho, T, massfrac, gamma);
    // eos.RTY2dpde_dpdre(rho, T, massfrac, dpde, dpdr_e);

    q(i,j,k,QRHO) = rho;
    AMREX_D_TERM(q(i,j,k,QU) = ux; ,
                 q(i,j,k,QV) = uy; ,
                 q(i,j,k,QW) = uz;);    
    q(i,j,k,QPRES) = p;
    q(i,j,k,QC)    = cs;
    q(i,j,k,QEINT) = ei;
    q(i,j,k,QTEMP) = T;
}

/**
 * \brief Convert primitive variables to characteristic variables.
 * 
 * @param i    x index.
 * @param j    y index.
 * @param k    z index.
 * @param dir  direction along characteristic decomposition.
 * @param q    primitive variables array.
 * @param w    characteristic variables array.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_ctochar (int i, int j, int k, int dir, 
             amrex::Array4<amrex::Real> const& q, 
             amrex::Array4<amrex::Real> const& w)
{
    using amrex::Real;

    int QUN, QUT, QUTT;
    if      (dir == 0) { QUN = QU;  QUT = QV;  QUTT = QW; } 
    else if (dir == 1) { QUN = QV;  QUT = QU;  QUTT = QW; } 
    else               { QUN = QW;  QUT = QU;  QUTT = QV; }
    
    Real rmpoc2 = q(i,j,k,QRHO) - q(i,j,k,QPRES)/q(i,j,k,QC)/q(i,j,k,QC);
    for (int n = 0; n < NUM_SPECIES; ++n) {
        w(i,j,k,WMAT+n) = rmpoc2 * q(i,j,k,QFS+n);
    }
    w(i,j,k,WACO)   = 0.5 * (q(i,j,k,QPRES)/q(i,j,k,QC) + q(i,j,k,QRHO)*q(i,j,k,QUN));
    w(i,j,k,WACO+1) = 0.5 * (q(i,j,k,QPRES)/q(i,j,k,QC) - q(i,j,k,QRHO)*q(i,j,k,QUN));
    AMREX_D_TERM(,
        w(i,j,k,WUT)   = q(i,j,k,QUT);,
        w(i,j,k,WUT+1) = q(i,j,k,QUTT););
}

/**
 * \brief Compute flux divergence. dSdt = -div(flux)*area/volume.
 * 
 * @param i     x index.
 * @param j     y index.
 * @param k     z index.
 * @param n     state variable array index.
 * @param dsdt  time derivatives of state variable.
 * @param f     fluxes in x, y, z direction.
 * @param dxinv inverse of cell size, 1/[dx,dy,dz].
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_div (int i, int j, int k, int n,
         amrex::Array4<amrex::Real> const& dsdt,
         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fx,
                      amrex::Array4<amrex::Real const> const& fy,
                      amrex::Array4<amrex::Real const> const& fz),
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    dsdt(i,j,k,n) = AMREX_D_TERM(dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n)),
                               + dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n)),
                               + dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n)));
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// wave_speed_est (const amrex::Real cl,)
//c, u, ut1, ut2, e

/**
 * \brief HLLC approximate Riemann solver.
 * 
 * @param rl    left density.
 * @param ul    left normal velocity.
 * @param ut1l  left tangential velocity.
 * @param ut2l  left tangential velocity.
 * @param pl    left pressure.
 * @param cl    left speed of sound.
 * @param el    left specific total energy.
 * @param yl    left mass fraction.
 * @param rr    right density.
 * @param ur    right normal velocity.
 * @param ut1r  right tangential velocity.
 * @param ut2r  right tangential velocity.
 * @param pr    right pressure.
 * @param cr    right speed of sound.
 * @param er    right specific total energy.
 * @param yr    right mass fraction.
 * @param flxrho  output density flux (= sum mass fraction fluxes).
 * @param flxu    output momentum fluxes.
 * @param flxrhoe output total energy fluxes.
 * @param flxrhoy output mass fraction fluxes.
 * @param parm  parameters and constants (see parm.H).
 */
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
         amrex::Real flxrhoy[NUM_SPECIES],
         Parm const& parm) noexcept
{
    // HLLC Riemann solver   

    using amrex::Real;

    // Estimate wave speeds
    // Arithmetic average (Davis), this alone won't work!
    Real sl = ul - cl;
    Real sr = ur + cr;
    // Roe-average (Einfeldt's HLLE)
    Real rp = sqrt(rr / rl);
    // AMREX_D_TERM(
    Real uroe = (ul + ur*rp)/(1. + rp); //,
    // Real vroe = (ut1l + ut1r*rp)/(1. + rp); ,
    // Real wroe = (ut2l + ut2r*rp)/(1. + rp););
    Real croe = (cl*cl + cr*cr*rp)/(1. + rp); + 0.5*rp/(1. + rp)/(1. + rp)
                *(AMREX_D_TERM((ul-ur)*(ul-ur), +(ut1l-ut1r)*(ut1l-ut1r), +(ut2l-ut2r)*(ut2l-ut2r)));
    croe = sqrt(croe); 
    // auto eos = pele::physics::PhysicsType::eos();
    // Real Tl, Tr, hl, hr;
    // Real ytmp[NUM_SPECIES];
    // for (int n = 0; n < NUM_SPECIES; ++n) { ytmp[n] = yl[n]; }
    // eos.EY2T(el, ytmp, Tl);
    // eos.TY2H(Tl, ytmp, hl);
    // hl += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));
    // for (int n = 0; n < NUM_SPECIES; ++n) { ytmp[n] = yr[n]; }
    // eos.EY2T(er, ytmp, Tr);
    // eos.TY2H(Tr, ytmp, hr);
    // hr += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
    // Real hroe = (hl + hr*rp)/(1. + rp);
    // hroe -= 0.5*(AMREX_D_TERM(uroe*uroe, +vroe*vroe, +wroe*wroe));
    // Real yroe[NUM_SPECIES];
    // Real rroe = 0.0;
    // Real Troe, croe;
    // for (int n = 0; n < NUM_SPECIES; ++n) {
    //     yroe[n] = sqrt(yl[n] * yr[n]);
    //     rroe += yroe[n];
    // }
    // eos.HY2T(hroe, yroe, Troe);
    // eos.RTY2Cs(rroe, Troe, yroe, croe);    
    // sl = uroe - croe;
    // sr = uroe + croe;
    sl = amrex::min(sl, uroe-croe);
    sr = amrex::max(sr, uroe+croe);

    if (sl > 0) {
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

    } else if (sr < 0) {
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
        Real sstar = (pr - pl + rl*ul*(sl-ul) - rr*ur*(sr-ur)) / (rl*(sl-ul) - rr*(sr-ur)); //contact wave speed
        
        if (sstar > parm.smallu*(amrex::Math::abs(ul) + amrex::Math::abs(ur))) {
            //flx_l* = flx_l + sl * (q_l* - q_l)
            Real frac = (sl - ul) / (sl - sstar) - 1.;
            flxrho = 0.;
            AMREX_D_TERM(
                flxu = rl*ul*ul + pl + sl*rl*((frac + 1.)*sstar - ul); ,
                flxut  = rl*ul*ut1l  + sl*rl*frac*ut1l; ,
                flxutt = rl*ul*ut2l  + sl*rl*frac*ut2l;);
            flxrhoe = ul*(rl*el + pl)
                   + sl*rl*(frac*el + (sstar - ul)*(sstar + pl/rl/(sl-ul)));
            for (int n = 0; n < NUM_SPECIES; ++n) {
                flxrhoy[n] = rl*ul*yl[n] + sl*rl*frac*yl[n];
                flxrho += flxrhoy[n];
            }

        } else if (sstar < -parm.smallu*(amrex::Math::abs(ul) + amrex::Math::abs(ur))) {
            //flx_r* = flx_r + sr * (q_r* - q_r)
            Real frac = (sr - ur) / (sr - sstar) - 1.;
            flxrho = 0.;
            AMREX_D_TERM(
                flxu = rr*ur*ur + pr + sr*rr*((frac + 1.)*sstar - ur); ,
                flxut  = rr*ur*ut1r  + sr*rr*frac*ut1r; ,
                flxutt = rr*ur*ut2r  + sr*rr*frac*ut2r;);
            flxrhoe = ur*(rr*er + pr)
                   + sr*rr*(frac*er + (sstar - ur)*(sstar + pr/rr/(sr-ur)));
            for (int n = 0; n < NUM_SPECIES; ++n) {
                flxrhoy[n] = rr*ur*yr[n] + sr*rr*frac*yr[n];
                flxrho += flxrhoy[n];
            }

        } else { 
            //sstar ~ 0, return to central scheme
            Real rc = 0.5 * (rl + rr);
            Real uc = 0.5 * (ul + ur);
            Real pc = 0.5 * (pl + pr);
            flxrho = 0.;
            AMREX_D_TERM(
                flxu = rc*uc*uc + pc; ,
                flxut  = rc*uc*0.5*(ut1l+ut1r); ,
                flxutt = rc*uc*0.5*(ut2l+ut2r););
            flxrhoe = uc*(rc*0.5*(el+er) + pc);
            for (int n = 0; n < NUM_SPECIES; ++n) {
                flxrhoy[n] = rc*ul*0.5*(yl[n]+yr[n]);
                flxrho += flxrhoy[n];
            }

        }
    }
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// riemann (const amrex::Real rl, 
//          AMREX_D_DECL(const amrex::Real ul, const amrex::Real ut1l, const amrex::Real ut2l), 
//          const amrex::Real pl, const amrex::Real cl, const amrex::Real el,
//          const amrex::Real yl[NUM_SPECIES],
//          const amrex::Real rr, 
//          AMREX_D_DECL(const amrex::Real ur, const amrex::Real ut1r, const amrex::Real ut2r), 
//          const amrex::Real pr, const amrex::Real cr, const amrex::Real er,
//          const amrex::Real yr[NUM_SPECIES],
//          amrex::Real& flxrho, 
//          AMREX_D_DECL(amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt),
//          amrex::Real& flxrhoe, 
//          amrex::Real flxrhoy[NUM_SPECIES],
//          Parm const& parm) noexcept
// {
//     // HLLC Riemann solver   

//     using amrex::Real;

//     // Estimate wave speeds
//     // Arithmetic average (Davis), this alone won't work!
//     Real sl = ul - cl;
//     Real sr = ur + cr;
//     // Roe-average (Einfeldt's HLLE)
//     Real rp = sqrt(rr / rl);
//     Real uroe = (ul + ur*rp)/(1. + rp);
//     Real croe = (cl*cl + cr*cr*rp)/(1. + rp); + 0.5*rp/(1. + rp)/(1. + rp)
//                 *(AMREX_D_TERM((ul-ur)*(ul-ur), +(ut1l-ut1r)*(ut1l-ut1r), +(ut2l-ut2r)*(ut2l-ut2r)));
//     croe = sqrt(croe); 
//     sl = amrex::min(sl, uroe-croe);
//     sr = amrex::max(sr, uroe+croe);

//     Real yout[NUM_SPECIES];
//     Real rhoout, AMREX_D_DECL(uout, ut1out, ut2out), pout;
//     if (sl > 0) {
//         //q_out = q_l
//         rhoout = rl;
//         AMREX_D_TERM(
//             uout   = ul; ,
//             ut1out = ut1l; ,
//             ut2out = ut2l;);
//         pout = pl;
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             yout[n] = yl[n];
//         }

//     } else if (sr < 0) {
//         //q_out = q_r
//         rhoout = rr;
//         AMREX_D_TERM(
//             uout   = ur; ,
//             ut1out = ut1r; ,
//             ut2out = ut2r;);
//         pout = pr;
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             yout[n] = yr[n];
//         }

//     } else {
//         Real sstar = (pr - pl + rl*ul*(sl-ul) - rr*ur*(sr-ur)) / (rl*(sl-ul) - rr*(sr-ur)); //contact wave speed

//         if (sstar > parm.smallu*(amrex::Math::abs(ul) + amrex::Math::abs(ur))) {
//             //q_out = q_l*
//             rhoout = rl*(sl-ul)/(sl-sstar);
//             AMREX_D_TERM(
//                 uout   = sstar; ,
//                 ut1out = ut1l; ,
//                 ut2out = ut2l;);
//             pout = pl + rl*(sl-ul)*(sstar-ul);
//             for (int n = 0; n < NUM_SPECIES; ++n) {
//                 yout[n] = yl[n];
//             }

//         } else if (sstar < -parm.smallu*(amrex::Math::abs(ul) + amrex::Math::abs(ur))) {
//             //q_out = q_r*
//             rhoout = rr*(sr-ur)/(sr-sstar);
//             AMREX_D_TERM(
//                 uout   = sstar; ,
//                 ut1out = ut1r; ,
//                 ut2out = ut2r;);
//             pout = pr + rr*(sr-ur)*(sstar-ur);
//             for (int n = 0; n < NUM_SPECIES; ++n) {
//                 yout[n] = yr[n];
//             }

//         } else { 
//             //sstar ~ 0, return to central scheme
//             rhoout = 0.5*(rl+rr);
//             AMREX_D_TERM(
//                 uout   = 0.5*(ul+ur); ,
//                 ut1out = 0.5*(ut1l+ut1r); ,
//                 ut2out = 0.5*(ut2l+ut2r););
//             pout = 0.5*(pl+pr);
//             for (int n = 0; n < NUM_SPECIES; ++n) {
//                 yout[n] = 0.5*(yl[n]+yr[n]);                
//             }

//         }
//     }

//     // flux(q_out)
//     auto eos = pele::physics::PhysicsType::eos();
//     Real eout;
//     eos.RYP2E(rhoout, yout, pout, eout);
//     eout += 0.5*(AMREX_D_TERM(uout*uout, +ut1out*ut1out, +ut2out*ut2out));

//     flxrho = 0.;
//     AMREX_D_TERM(
//         flxu = rhoout*uout*uout + pout; ,
//         flxut  = rhoout*uout*ut1out; ,
//         flxutt = rhoout*uout*ut2out;);
//     flxrhoe = uout*(rhoout*eout + pout);
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         flxrhoy[n] = rhoout*uout*yout[n];
//         flxrho += flxrhoy[n];
//     }
// }

// #if (AMREX_SPACEDIM > 1)
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// compute_hyp_wallflux (const amrex::Real rho,
//                       AMREX_D_DECL(const amrex::Real u, const amrex::Real v, const amrex::Real w),
//                       const amrex::Real p,   const amrex::Real y[NUM_SPECIES],
//                       const amrex::Real axm, const amrex::Real axp,
//                       const amrex::Real aym, const amrex::Real ayp,
// #if (AMREX_SPACEDIM == 3)
//                       const amrex::Real azm, const amrex::Real azp,
// #endif
//                       amrex::GpuArray<amrex::Real,NVAR>& fluxw, Parm const& parm) noexcept
// {
//     using amrex::Real;

//     Real apnorm = std::sqrt( AMREX_D_TERM((axm-axp)*(axm-axp), + (aym-ayp)*(aym-ayp), + (azm-azp)*(azm-azp)) );
//     Real apnorminv = 1./apnorm;
//     Real un = AMREX_D_TERM(u*(axm-axp)*apnorminv, + v*(aym-ayp)*apnorminv, + w*(azm-azp)*apnorminv);
//     Real AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2);

//     auto eos = pele::physics::PhysicsType::eos();
//     Real e, cs;
//     eos.RYP2E(rho, y, p, e);
//     e += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
//     eos.RPY2Cs(rho, p, y, cs);

//     riemann(rho, AMREX_D_DECL( un, 0., 0.), p, cs, e, y,
//             rho, AMREX_D_DECL(-un, 0., 0.), p, cs, e, y,
//             fluxw[URHO], AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2), fluxw[UEDEN], &fluxw[UFS]);

//     AMREX_D_TERM(fluxw[UMX] = (axm-axp) * momfluxn;,
//                  fluxw[UMY] = (aym-ayp) * momfluxn;,
//                  fluxw[UMZ] = (azm-azp) * momfluxn;);
// }
// #endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_riemann (int i, int j, int k, int dir,
             amrex::Array4<amrex::Real> const& flx,
             amrex::Array4<amrex::Real const> const& q,
             amrex::Array4<amrex::Real const> const& wl,
             amrex::Array4<amrex::Real const> const& wr,
             Parm const& parm) noexcept
{
    using amrex::Real;

    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

    int UMN, UMT, UMTT;
    if (dir == 0) {
        UMN = UMX; UMT = UMY; UMTT = UMZ;  
    } else if (dir == 1) {
        UMN = UMY; UMT = UMX; UMTT = UMZ;  
    } else {
        UMN = UMZ; UMT = UMX; UMTT = UMY;  
    }

    // Approximate face (rho, c, y)
    // Real rp = sqrt(q(iv,QRHO)/q(iv-iv_dir,QRHO));
    // Real cs = sqrt( (q(iv-iv_dir,QC)*q(iv-iv_dir,QC) + q(iv,QC)*q(iv,QC)*rp)/(1.+rp) );
    // Real rho = 1. + rp;
    // Real cs = 0.5*(q(iv-iv_dir,QC) + q(iv,QC));
    // Real invc = 1.0 / cs;
    // Real rho = 0.5*(q(iv-iv_dir,QRHO) + q(iv,QRHO));

    // Convert wl to ql (rl, ul, utl, pl, cl, el, yl)
    Real cs = q(iv-iv_dir,QC);
    Real invc = 1.0 / cs;
    Real rho = q(iv-iv_dir,QRHO);
    Real ryl[NUM_SPECIES];
    Real yl[NUM_SPECIES];
    Real rl = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        ryl[n] = wl(iv,WMAT+n) + (wl(iv,WACO) + wl(iv,WACO+1)) * q(iv-iv_dir,QFS+n) * invc;
        ryl[n] = amrex::max(ryl[n], parm.smallr); //enforce positive density
        rl += ryl[n];
    }
    Real invr = 1.0 / rl;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = ryl[n] * invr;
    }
    AMREX_D_TERM(
        Real ul   = (wl(iv,WACO) - wl(iv,WACO+1)) / rl; ,
        Real ut1l = wl(iv,WUT); ,
        Real ut2l = wl(iv,WUT+1););
    Real pl = (wl(iv,WACO) + wl(iv,WACO+1)) * cs;
    pl = amrex::max(pl, parm.smallp); //enforce positive pressure
    Real el, cl;
    auto eos = pele::physics::PhysicsType::eos();
    eos.RYP2E(rl, yl, pl, el);
    el += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));   
    eos.RPY2Cs(rl, pl, yl, cl);
    
    // Convert wr to qr
    cs = q(iv,QC);
    invc = 1.0 / cs;
    rho = q(iv,QRHO);
    Real ryr[NUM_SPECIES];
    Real yr[NUM_SPECIES];
    Real rr = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        ryr[n] = wr(iv,WMAT+n) + (wr(iv,WACO)+wr(iv,WACO+1)) * q(iv,QFS+n) * invc;
        ryr[n] = amrex::max(ryr[n], parm.smallr); //enforce positive mass
        rr += ryr[n];
    }
    invr = 1.0 / rr;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yr[n] = ryr[n] * invr;
    }
    AMREX_D_TERM(
        Real ur   = (wr(iv,WACO) - wr(iv,WACO+1)) / rr; ,
        Real ut1r = wr(iv,WUT); ,
        Real ut2r = wr(iv,WUT+1););
    Real pr = (wr(iv,WACO) + wr(iv,WACO+1)) * cs;
    pr = amrex::max(pr, parm.smallp); //enforce positive pressure
    Real er, cr;
    eos.RYP2E(rr, yr, pr, er);
    er += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
    eos.RPY2Cs(rr, pr, yr, cr);
       
    riemann(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
            rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
            flx(i,j,k,URHO), 
            AMREX_D_DECL(flx(i,j,k,UMN), flx(i,j,k,UMT), flx(i,j,k,UMTT)), 
            flx(i,j,k,UEDEN),
            &flx(i,j,k,UFS),
            parm);
}
#endif