#ifndef CNS_HYPERBOLICS_H_
#define CNS_HYPERBOLICS_H_

#include <AMReX_FArrayBox.H>

#include "index_macros.H"
#include "parm.H"
#include "PelePhysics.H"

#include <cmath>

/**
 * \brief Convert state variable to primitive variable
 * 
 * @param i    x index.
 * @param j    y index.
 * @param k    z index.
 * @param ns   starting component index for u.
 * @param u    state variables array.
 * @param q    primitive variables array.
 * @param parm parameters and constants (see parm.H).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_ctoprim (int i, int j, int k, int ns,
             amrex::Array4<amrex::Real const> const& u,
             amrex::Array4<amrex::Real> const& q,
             Parm const& parm) noexcept
{
    amrex::Real rho = u(i,j,k,ns+URHO);
    // if (rho <= 0.0) amrex::Print() << "rho corrected \n";
    // rho = amrex::max<amrex::Real>(rho, parm.smallr); //preserve positivity of mass
    amrex::Real rhoinv = amrex::Real(1.0) / rho;
    AMREX_D_TERM(amrex::Real ux = u(i,j,k,ns+UMX)*rhoinv; ,
                 amrex::Real uy = u(i,j,k,ns+UMY)*rhoinv; ,
                 amrex::Real uz = u(i,j,k,ns+UMZ)*rhoinv;);
    amrex::Real ei = u(i,j,k,ns+UEDEN) * rhoinv - amrex::Real(0.5)*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real sumRhoY = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = u(i,j,k,ns+UFS+n) * rhoinv;        
        // massfrac[n] = amrex::max(massfrac[n], 0.0); //preserve positivity of mass
        // check_rho += massfrac[n];
        sumRhoY += u(i,j,k,ns+UFS+n);
        q(i,j,k,QFS+n) = massfrac[n];
    }
    if (std::abs(sumRhoY - rho) > parm.smallr) {
        // if (std::abs(check_rho - 1.0) > 1e-4) {
        //     amrex::Print() << "i,j,k=" << i << "," << j << "," << k << ": sum(Y)-1 = " << check_rho - 1.0 << "!!!\n";
        // }
        if (std::abs(sumRhoY - rho) > 1e-5) {
            amrex::Print() << "i,j,k=" << i << "," << j << "," << k << ": sum(rY)-rho = " << sumRhoY - rho << "!!!\n";
        }
        // for (int n = 0; n < NUM_SPECIES; n++) {
        //     // massfrac[n] /= check_rho;
        //     // q(i,j,k,QFS+n) = massfrac[n];
        //     q(i,j,k,ns+URHO) = sumRhoY;
        // }
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T, p, cs, gamma;
    // eos.Y2WBAR(massfrac, wbar);
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2P(rho, T, massfrac, p);
    // if (T <= 0.0) amrex::Print() << "T corrected \n";
    // T = amrex::max<amrex::Real>(T, 1e-10); //clip temperature
    eos.RTY2Cs(rho, T, massfrac, cs);
    eos.RTY2G(rho, T, massfrac, gamma);
    // eos.RTY2dpde_dpdre(rho, T, massfrac, dpde, dpdr_e);
    

    q(i,j,k,QRHO) = rho;
    AMREX_D_TERM(q(i,j,k,QU) = ux; ,
                 q(i,j,k,QV) = uy; ,
                 q(i,j,k,QW) = uz;);    
    q(i,j,k,QPRES) = p; //amrex::max(p, parm.smallp); //preserve positivity of pressure
    q(i,j,k,QC)    = cs;
    q(i,j,k,QG)    = gamma;
    q(i,j,k,QEINT) = ei;
    q(i,j,k,QTEMP) = T;

    // if (T < 50.0) std::cout << "i=" << i << " j=" << j << " k=" << k <<
    //   " r=" << rho << " T=" << T << " p=" << p << " c=" << cs << std::endl;
    
    if (amrex::isnan(cs) || amrex::isnan(ux) || (rho < 0.0) || (p < 0.0)) {
        std::cout << "i=" << i << " j=" << j << " k=" << k << " ns=" << ns << ": rho=" << rho << " T=" << T << " p=" << p;// << "\n";
        std::cout << " rhoE=" << u(i,j,k,ns+UEDEN) << " ei=" << ei << "\n";
        for (int n = 0; n < NUM_SPECIES; n++) {
            std::cout << u(i,j,k,ns+UFS+n) << " ";
        }
        std::cout << "\n";
        // amrex::Abort("cs is nan");
    }
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
 * @param sys  gamma system or sos system.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_ctochar (int i, int j, int k, int dir, 
             amrex::Array4<const amrex::Real> const& q, 
             amrex::Array4<amrex::Real> const& w,
             const int sys)
{
    using amrex::Real;

    int QUN, QUT, QUTT;
    if      (dir == 0) { QUN = QU;  QUT = QV;  QUTT = QW; } 
    else if (dir == 1) { QUN = QV;  QUT = QU;  QUTT = QW; } 
    else               { QUN = QW;  QUT = QU;  QUTT = QV; }
    
    if (sys == 0) { // speed of sound system
        // Real rmpoc2 = q(i,j,k,QRHO) - q(i,j,k,QPRES)/q(i,j,k,QC)/q(i,j,k,QC); //rho minus p over c^2
        w(i,j,k,WRHO) = q(i,j,k,QRHO) - q(i,j,k,QPRES)/q(i,j,k,QC)/q(i,j,k,QC); //rho minus p over c^2
        w(i,j,k,WACO)   = 0.5 * (q(i,j,k,QPRES)/q(i,j,k,QC) + q(i,j,k,QRHO)*q(i,j,k,QUN));
        w(i,j,k,WACO+1) = 0.5 * (q(i,j,k,QPRES)/q(i,j,k,QC) - q(i,j,k,QRHO)*q(i,j,k,QUN));

    } else { // gamma system
        w(i,j,k,WRHO) = q(i,j,k,QRHO) * (1.0 - 1.0 / q(i,j,k,QG));
        w(i,j,k,WACO)   = 0.5 * (q(i,j,k,QPRES) + std::sqrt(q(i,j,k,QG)*q(i,j,k,QRHO)*q(i,j,k,QPRES))*q(i,j,k,QUN));
        w(i,j,k,WACO+1) = 0.5 * (q(i,j,k,QPRES) - std::sqrt(q(i,j,k,QG)*q(i,j,k,QRHO)*q(i,j,k,QPRES))*q(i,j,k,QUN));
    }

    // Passive scalars are the same for both systems
    for (int n = 0; n < NUM_SPECIES; ++n) {
        w(i,j,k,WY+n) = q(i,j,k,QFS+n);
    }
    w(i,j,k,WC) = q(i,j,k,QC);
    AMREX_D_TERM(,
        w(i,j,k,WUT)   = q(i,j,k,QUT);,
        w(i,j,k,WUT+1) = q(i,j,k,QUTT););
}

/**
 * \brief Convert characteristic variables to primitive variables.
 * 
 * This is done before the Riemann solver.
 * 
 * @param w1    characteristic variables WRHO.
 * @param w2    characteristic variables WACO.
 * @param w3    characteristic variables WACO+1.
 * @param c_est speed of sound estimate.
 * @param r_est density estimate.
 * @param g_est gamma estimate.
 * @param rho   unpacked density.
 * @param u     unpacked (normal) velocity.
 * @param p     unpacked pressure.
 * @param sys   gamma system (1) or sos system (0).
 * @param parm  parameters and constants (see parm.H).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_unpackchar (const amrex::Real w1, const amrex::Real w2, const amrex::Real w3, 
                const amrex::Real c, const amrex::Real /*r*/, const amrex::Real gamma, 
                amrex::Real& rho, amrex::Real& u, amrex::Real& p,
                const int sys, Parm const& parm)
{
    // using amrex::Real;
    
    if (sys == 0) { // speed of sound system
        rho = w1 + (w2 + w3) / c;        
        u = (w2 - w3) / rho;
        p = (w2 + w3) * c;

    } else { // gamma system
        p = w2 + w3;
        // rho = w1 + (w2 + w3) * r / gamma / p; 
        rho = w1 / (1.0 - 1.0 / gamma);
        u = (w2 - w3) / std::sqrt(gamma*rho*p);
    }

    // rho = amrex::max(rho, parm.smallr); //enforce positive density
    // p = amrex::max(p, parm.smallp); //enforce positive pressure
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
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
                               
    if (isnan(dsdt(i,j,k,n))) amrex::Print() << "div(" << i << "," << j << "," << k << ") is nan \n";
}

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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
hllc (const amrex::Real rl, 
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
    Real croe = (cl*cl + cr*cr*rp)/(1. + rp) + 0.5*rp/(1. + rp)/(1. + rp)
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
                flxrhoy[n] = rc*uc*0.5*(yl[n]+yr[n]);
                flxrho += flxrhoy[n];
            }

        }
    }
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// hllc (const amrex::Real rl, 
//       AMREX_D_DECL(const amrex::Real ul, const amrex::Real ut1l, const amrex::Real ut2l), 
//       const amrex::Real pl, const amrex::Real cl, const amrex::Real el,
//       const amrex::Real yl[NUM_SPECIES],
//       const amrex::Real rr, 
//       AMREX_D_DECL(const amrex::Real ur, const amrex::Real ut1r, const amrex::Real ut2r), 
//       const amrex::Real pr, const amrex::Real cr, const amrex::Real er,
//       const amrex::Real yr[NUM_SPECIES],
//       amrex::Real& flxrho, 
//       AMREX_D_DECL(amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt),
//       amrex::Real& pflxu, amrex::Real& pflxe,
//       amrex::Real& flxrhoe, 
//       amrex::Real flxrhoy[NUM_SPECIES],
//       Parm const& parm) noexcept
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
//     Real croe = (cl*cl + cr*cr*rp)/(1. + rp) + 0.5*rp/(1. + rp)/(1. + rp)
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
//             ut2out = ut2l;)
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
//             ut2out = ut2r;)
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
//                 ut2out = ut2l;)
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
//                 ut2out = ut2r;)
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
//                 ut2out = 0.5*(ut2l+ut2r);)
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
//         flxu   = rhoout*uout*uout + pout; ,
//         flxut  = rhoout*uout*ut1out; ,
//         flxutt = rhoout*uout*ut2out;)
//     flxrhoe = uout*(rhoout*eout + pout);
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//         flxrhoy[n] = rhoout*uout*yout[n];
//         flxrho += flxrhoy[n];
//     }
//     pflxu = pout;
//     pflxe = uout;
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

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_riemann (int i, int j, int k, int dir,
//              amrex::Array4<amrex::Real> const& flx,
//              amrex::Array4<amrex::Real const> const& q,
//              amrex::Array4<amrex::Real const> const& wl,
//              amrex::Array4<amrex::Real const> const& wr,
//              Parm const& parm) noexcept
// {
//     using amrex::Real;

//     const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//     const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

//     int UMN, UMT, UMTT;
//     if (dir == 0) {
//         UMN = UMX; UMT = UMY; UMTT = UMZ;  
//     } else if (dir == 1) {
//         UMN = UMY; UMT = UMX; UMTT = UMZ;  
//     } else {
//         UMN = UMZ; UMT = UMX; UMTT = UMY;  
//     }

//     // Convert wl to ql (rl, ul, utl, pl, cl, el, yl)
//     Real cl = q(iv-iv_dir,QC); //first guess
//     Real ryl[NUM_SPECIES];
//     Real yl[NUM_SPECIES];
//     Real rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, el;
//     Real invc, invr;
//     auto eos = pele::physics::PhysicsType::eos();

//     for (int iter = 0; iter < 1; ++iter) {
//         invc = 1.0 / cl;
//         rl = Real(0.0);
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             ryl[n] = wl(iv,WMAT+n) + (wl(iv,WACO) + wl(iv,WACO+1)) * q(iv-iv_dir,QFS+n) * invc;
//             ryl[n] = amrex::max(ryl[n], parm.smallr); //enforce positive density
//             rl += ryl[n];
//         }
//         invr = 1.0 / rl;
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             yl[n] = ryl[n] * invr;
//         }
//         AMREX_D_TERM(
//             ul   = (wl(iv,WACO) - wl(iv,WACO+1)) / rl; ,
//             ut1l = wl(iv,WUT); ,
//             ut2l = wl(iv,WUT+1););
//         pl = (wl(iv,WACO) + wl(iv,WACO+1)) * cl;
//         pl = amrex::max(pl, parm.smallp); //enforce positive pressure
//         eos.RYP2E(rl, yl, pl, el);
//         el += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));   
//         eos.RPY2Cs(rl, pl, yl, cl);
//     }
    
//     // Convert wr to qr
//     Real cr = q(iv,QC); //first guess
//     Real ryr[NUM_SPECIES];
//     Real yr[NUM_SPECIES];
//     Real rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, er;

//     for (int iter = 0; iter < 1; ++iter) {
//         invc = 1.0 / cr;    
//         rr = 0.0;
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             ryr[n] = wr(iv,WMAT+n) + (wr(iv,WACO)+wr(iv,WACO+1)) * q(iv,QFS+n) * invc;
//             ryr[n] = amrex::max(ryr[n], parm.smallr); //enforce positive mass
//             rr += ryr[n];
//         }
//         invr = 1.0 / rr;
//         for (int n = 0; n < NUM_SPECIES; ++n) {
//             yr[n] = ryr[n] * invr;
//         }
//         AMREX_D_TERM(
//             ur   = (wr(iv,WACO) - wr(iv,WACO+1)) / rr; ,
//             ut1r = wr(iv,WUT); ,
//             ut2r = wr(iv,WUT+1););
//         pr = (wr(iv,WACO) + wr(iv,WACO+1)) * cr;
//         pr = amrex::max(pr, parm.smallp); //enforce positive pressure
//         eos.RYP2E(rr, yr, pr, er);
//         er += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
//         eos.RPY2Cs(rr, pr, yr, cr);
//     }

//     hllc(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
//          rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
//          flx(i,j,k,URHO), 
//          AMREX_D_DECL(flx(i,j,k,UMN), flx(i,j,k,UMT), flx(i,j,k,UMTT)), 
//          flx(i,j,k,UEDEN),
//          &flx(i,j,k,UFS),
//          parm);
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_riemann (int i, int j, int k, int dir,
             amrex::Array4<amrex::Real> const& adv_flx,
             amrex::Array4<amrex::Real> const& p_flx,
             amrex::Array4<const amrex::Real> const& q,
             amrex::Array4<const amrex::Real> const& wl,
             amrex::Array4<const amrex::Real> const& wr,
             const int sys,
             Parm const& parm) noexcept
{
    using amrex::Real;

    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

    int UMN, UMT, UMTT;//, QUN;
    if (dir == 0) {
        UMN = UMX; UMT = UMY; UMTT = UMZ; //QUN = QU;
    } else if (dir == 1) {
        UMN = UMY; UMT = UMX; UMTT = UMZ; //QUN = QV;
    } else {
        UMN = UMZ; UMT = UMX; UMTT = UMY; //QUN = QW;
    }

    // Real r = 0.5*(q(iv,QRHO) + q(iv-iv_dir,QRHO));
    // Real r = sqrt(q(iv,QRHO) * q(iv-iv_dir,QRHO));
    Real g = 0.5*(q(iv,QG) + q(iv-iv_dir,QG));
    // Real rp = sqrt(q(iv,QRHO) / q(iv-iv_dir,QRHO));
    // Real c = (q(iv-iv_dir,QC) + q(iv,QC)*rp)/(1. + rp);
    // Real c = (q(iv-iv_dir,QC)*q(iv-iv_dir,QC) + q(iv,QC)*q(iv,QC)*rp)/(1. + rp)
    //     + 0.5*rp/(1. + rp)/(1. + rp)*((q(iv-iv_dir,QU)-q(iv,QU))*(q(iv-iv_dir,QU)-q(iv,QU)));
    // c = sqrt(c); 
    // Real c = sqrt(q(iv-iv_dir,QC)*q(iv,QC));
    auto eos = pele::physics::PhysicsType::eos();

    // Unpack wl to ql (rl, ul, utl, pl, cl, el, yl)
    Real r = q(iv-iv_dir,QRHO);
    // if (wl(iv,WC) < 0.0) std::cout << "cl < 0 = " << wl(iv,WC) << std::endl;
    // Real cl = amrex::max(1e-10, wl(iv,WC));
    Real cl = q(iv-iv_dir,QC);
    Real yl[NUM_SPECIES];
    Real rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, el;
    
    Real sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = amrex::max<amrex::Real>(0.0, wl(iv,WY+n));
        sumY += yl[n];
    }    
    if (std::abs(sumY - 1.0) > parm.smally) {
        Real invsumY = 1.0 / sumY;
        for (int n = 0; n < NUM_SPECIES; ++n) {
            yl[n] *= invsumY;
        }
    }
    AMREX_D_TERM(,
        ut1l = wl(iv,WUT); ,
        ut2l = wl(iv,WUT+1););

    cns_unpackchar(wl(iv,WRHO), wl(iv,WACO), wl(iv,WACO+1), cl, r, g, 
                   rl, ul, pl, sys, parm);
    
    // Convert wr to qr
    r = q(iv,QRHO);
    // if (wr(iv,WC) < 0.0) std::cout << "cr < 0 = " << wr(iv,WC) << std::endl;
    // Real cr = amrex::max(1e-10, wr(iv,WC));
    Real cr = q(iv,QC);
    Real yr[NUM_SPECIES];
    Real rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, er;

    sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yr[n] = amrex::max<amrex::Real>(0.0, wr(iv,WY+n));
        sumY += yr[n];
    }
    if (std::abs(sumY - 1.0) > parm.smally) {
        Real invsumY = 1.0 / sumY;
        for (int n = 0; n < NUM_SPECIES; ++n) {
            yr[n] *= invsumY;
        }
    }
    AMREX_D_TERM(,
        ut1r = wr(iv,WUT); ,
        ut2r = wr(iv,WUT+1););

    cns_unpackchar(wr(iv,WRHO), wr(iv,WACO), wr(iv,WACO+1), cr, r, g, 
                   rr, ur, pr, sys, parm);

    // Limit recon values: return to Godunov scheme if not admissable
    // if ((rl <= 0.0) || (pl <= 0.0) || (rr <= 0.0) || (pr <= 0.0)) {
    //     amrex::Print() << " Recon state should return to Godunov at " << iv << " "
    //                    << "rl=" << rl << " rr=" << rr << " pl=" << pl << " pr=" << pr << std::endl;
    //     rl = q(iv-iv_dir, QRHO);  rr = q(iv, QRHO);
    //     ul = q(iv-iv_dir, QUN);   ur = q(iv, QUN);
    //     pl = q(iv-iv_dir, QPRES); pr = q(iv, QPRES);
    //     for (int n = 0; n < NUM_SPECIES; ++n) {
    //         yl[n] = q(iv-iv_dir, QFS+n);  yr[n] = q(iv, QFS+n);
    //     }
    // }
    // if ((rl < 0.0) && (pl < 0.0)) { amrex::Print() << "Reverted r,pl at " << iv << "\n"; rl = -rl;  pl = -pl; }
    // if ((rr < 0.0) && (pr < 0.0)) { amrex::Print() << "Reverted r,pr at " << iv << "\n"; rr = -rr;  pr = -pr; }

    eos.RYP2E(rl, yl, pl, el);
    el += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));   
    eos.RPY2Cs(rl, pl, yl, cl);

    eos.RYP2E(rr, yr, pr, er);
    er += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
    eos.RPY2Cs(rr, pr, yr, cr);

    // if ((rl <= 0.0) || (pl <= 0.0) || (rr <= 0.0) || (pr <= 0.0) || (er <= 0.0) || (cr <= 0.0) || (el <= 0.0) || (cl <= 0.0)) {
    //     amrex::Print() << " rl=" << rl << " rr=" << rr << " pl=" << pl << " pr=" << pr << std::endl;
    //     amrex::Print() << " el=" << el << " er=" << er << " cl=" << cl << " cr=" << cr << std::endl;
    // }

    amrex::Real ry_flx[NUM_SPECIES] = {0.0};

    // This's wrong DONT USE!!!!!!!!
    // hllc(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
    //      rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
    //      adv_flx(i,j,k,URHO), 
    //      AMREX_D_DECL(adv_flx(i,j,k,UMN), adv_flx(i,j,k,UMT), adv_flx(i,j,k,UMTT)), 
    //      p_flx(i,j,k,UMN), p_flx(i,j,k,UEDEN),
    //      adv_flx(i,j,k,UEDEN),
    //      ry_flx, //&adv_flx(i,j,k,UFS),
    //      parm);

    hllc(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
         rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
         adv_flx(i,j,k,URHO), 
         AMREX_D_DECL(adv_flx(i,j,k,UMN), adv_flx(i,j,k,UMT), adv_flx(i,j,k,UMTT)), 
         adv_flx(i,j,k,UEDEN), ry_flx, //&adv_flx(i,j,k,UFS),
         parm);

    for (int n = 0; n < NUM_SPECIES; ++n) { adv_flx(i,j,k,UFS+n) = ry_flx[n]; }

    // if (std::abs(adv_flx(i,j,k,URHO)-adv_flx(i,j,k,UFS)-adv_flx(i,j,k,UFS+1))>1e-10) 
    //     amrex::Print() << i << j << k << dir << " " << adv_flx(i,j,k,URHO) << " " 
    //                                                 << adv_flx(i,j,k,UFS) << " " << adv_flx(i,j,k,UFS+1) << std::endl; //<< ", "
                                                    //<< rl << " " << yl[0] << " " << yl[1] << ", " << rr << " " << yr[0] << " " << yr[1] << std::endl;
}

// #if CNS_USE_EB
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// cns_riemann_eb (int i, int j, int k, int dir,
//                 amrex::Array4<amrex::Real> const& adv_flx,
//                 amrex::Array4<amrex::Real> const& p_flx,
//                 amrex::Array4<const amrex::Real> const& q,
//                 amrex::Array4<const amrex::Real> const& wl,
//                 amrex::Array4<const amrex::Real> const& wr,
//                 amrex::Array4<const amrex::EBCellFlag> const& flag,
//                 Parm const& parm) noexcept
// {
//     const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//     const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

//     if ( !flag(iv).isCovered() && !flag(iv-iv_dir).isCovered() ) {
//         cns_riemann(i, j, k, dir, adv_flx, p_flx, q, wl, wr, 0, parm);
//     }
// }
// #endif

// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//                                           EXPERIMENTAL
// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


// #include "recon.H"

// /**
//  * \param iv  (input) index of cell
//  * \param n   (input) which component of characteristic variable (0: material_wave, 1: +acoustic_wave, 2: -acoutic_wave)
//  * \param QUN (input) which component of q is u_normal
//  * \param q   (input) primitive variables array
//  * \param sys (input) indicating what thermodynamics system is used (0: speed of sound, 1: gamma)
//  * \param w   (output) chosen component of characteristic variable
//  */
// void q2w (const amrex::IntVect iv, const int n, const int QUN,
//           amrex::Array4<const amrex::Real> const& q, const int sys,
//           amrex::Real& w) 
// { 
//   if (sys == 0) { // speed of sound system
//     const Real cs = q(iv,QC);
//     if (n == 0) { // material_wave
//       w = q(iv,QRHO) - q(iv,QPRES)/cs/cs; 
//     } else if (n == 1) { // +acoustic_wave
//       // w = 0.5 * (q(iv,QPRES)/cs + q(iv,QRHO)*q(iv,QUN));
//       w = 0.5 * (q(iv,QPRES)/cs + q(iv,QRHO)*q(iv,QUN)) / cs;
//     } else if (n == 2) { // -acoustic_wave
//       // w = 0.5 * (q(iv,QPRES)/cs - q(iv,QRHO)*q(iv,QUN));
//       w = 0.5 * (q(iv,QPRES)/cs - q(iv,QRHO)*q(iv,QUN)) / cs;
//     }
//   } 
//   else { // specific heat system
//     const Real gamma = q(iv,QG);
//     if (n == 0) { // material_wave
//       w = q(iv,QRHO) - (1.0 - 1.0/gamma); 
//     } else if (n == 1) { // +acoustic_wave
//       w = 0.5*(q(iv,QPRES) + std::sqrt(gamma*q(iv,QRHO)*q(iv,QPRES))*q(iv,QUN));
//     } else if (n == 2) { // -acoustic_wave
//       w = 0.5*(q(iv,QPRES) - std::sqrt(gamma*q(iv,QRHO)*q(iv,QPRES))*q(iv,QUN));
//     }
//   }
// }

// // /**
// //  * \param q   (input) primitive variables [rho, u_normal, p]
// //  * \param thv (input) thermodynamics variable, can be c or gamma, consistent with system
// //  * \param sys (input) indicating what thermodynamics system is used (0: speed of sound, 1: gamma)
// //  * \param w   (output) characteristic variables [material_wave, +acoustic_wave, -acoutic_wave]
// //  */
// // void q2w (const amrex::Real q[3],
// //           const amrex::Real thv,
// //           const int sys,
// //           amrex::Real w[3])
// // {
// //   if (sys == 0) { // speed of sound system: thv = c
// //     w[0] = q[0] - q[2]/thv/thv;
// //     w[1] = 0.5 * (q[2]/thv + q[0]*q[1]);
// //     w[2] = 0.5 * (q[2]/thv - q[0]*q[1]);
// //   } 
// //   else { // specific heat system: thv = gamma
// //     w[0] = q[0] * (1.0 - 1.0/thv);
// //     w[1] = 0.5 * (q[2] + std::sqrt(thv*q[0]*q[2])*q[1]);
// //     w[2] = 0.5 * (q[2] - std::sqrt(thv*q[0]*q[2])*q[1]);
// //   }
// // }

// /**
//  * \param iv  (input) index of cell in output q
//  * \param QUN (input) which component of q is u_normal
//  * \param w   (input) characteristic variables [material_wave, +acoustic_wave, -acoutic_wave] 
//  * \param thv (input) thermodynamics variable, can be c or gamma, consistent with system
//  * \param sys (input) indicating what thermodynamics system is used (0: speed of sound, 1: gamma)
//  * \param q   (output) primitive variables array at face
//  */
// void w2q (const amrex::IntVect iv, const int QUN,
//           const amrex::Real w[3],
//           const amrex::Real thv,
//           const int sys,
//           amrex::Array4<amrex::Real> const& q)
// {
//   if (sys == 0) { // speed of sound system: thv = c
//     // q(iv,QRHO)  = w[0] + (w[1] + w[2]) / thv;        
//     // q(iv,QUN)   = (w[1] - w[2]) / q(iv,QRHO);
//     // q(iv,QPRES) = (w[1] + w[2]) * thv;
//     q(iv,QRHO)  = w[0] + (w[1] + w[2]) ;        
//     q(iv,QUN)   = (w[1] - w[2]) / q(iv,QRHO) * thv;
//     q(iv,QPRES) = (w[1] + w[2]) * thv * thv;
//   } 
//   else { // specific heat system: thv = gamma
//     q(iv,QPRES) = w[1] + w[2];
//     q(iv,QRHO)  = w[0] / (1.0 - (w[1] + w[2]) / thv / q(iv,QPRES));
//     q(iv,QUN)   = (w[1] - w[2]) / std::sqrt(thv*q(iv,QRHO)*q(iv,QPRES));
//   }
// }

// /** 
//  * \brief Reconstruct states at cell interface.
//  */
// void
// new_recons (int i, int j, int k, int dir,
//             amrex::Array4<const amrex::Real> const& q,
//             amrex::Array4<amrex::Real> const& ql,
//             amrex::Array4<amrex::Real> const& qr,
//             const int recon_scheme, 
//             const int sys, 
//             const bool do_pp_limiter,
//             const amrex::Real plm_theta,
//             Parm const& parm) noexcept
// {
//   using amrex::Real;

//   const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//   const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

//   int QUN, QUT, QUTT;
//   if (dir == 0) {
//     QUN = QU;  QUT = QV;  QUTT = QW;
//   } else if (dir == 1) {
//     QUN = QV;  QUT = QU;  QUTT = QW;
//   } else {
//     QUN = QW;  QUT = QU;  QUTT = QV;
//   }

//   int QTHV;
//   if (sys == 0) {
//     QTHV = QC;
//   } else {
//     QTHV = QG;
//   }

//   // Reconstruction
//   Real wm[3], wp[3];
//   Real thvm, thvp;

//   if (recon_scheme == 1) {                               ///////// Godunov scheme /////////
//     ql(iv+iv_dir,QRHO) = qr(iv,QRHO) = q(iv,QRHO);
//     AMREX_D_TERM(ql(iv+iv_dir,QUN)  = qr(iv,QUN)  = q(iv,QUN); ,
//                  ql(iv+iv_dir,QUT)  = qr(iv,QUT)  = q(iv,QUT); ,
//                  ql(iv+iv_dir,QUTT) = qr(iv,QUTT) = q(iv,QUTT););
//     ql(iv+iv_dir,QPRES) = qr(iv,QPRES) = q(iv,QPRES); 
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       ql(iv+iv_dir,QFS+n) = qr(iv,QFS+n) = q(iv,QFS+n);
//     }
//   } 
//   else if ((recon_scheme == 2) || (recon_scheme == 3)) { ///////// MUCSL / WENO-Z3 /////////
//     Real s[3]; //stencil

//     // Reconstruct euler equation variables
//     for (int n = 0; n < 3; ++n) {
//       q2w(iv-iv_dir, n, QUN, q, sys, s[0]);
//       q2w(iv,        n, QUN, q, sys, s[1]);
//       q2w(iv+iv_dir, n, QUN, q, sys, s[2]);
//       if (recon_scheme == 2) {
//         muscl(s, wm[n], wp[n], plm_theta);
//       } else {
//         wenoz3(s, wm[n], wp[n]);
//       }
//     }
//     // Reconstruct the thermodynamic variable for Jacobian
//     s[0] = q(iv-iv_dir,QTHV);
//     s[1] = q(iv,QTHV);
//     s[2] = q(iv+iv_dir,QTHV);
//     if (recon_scheme == 2) {
//       muscl(s, thvm, thvp, plm_theta);
//     } else {
//       wenoz3(s, thvm, thvp);
//     }
//     pp_limiter_scalar(thvm, thvp, s[1], 1e-4, 10*s[0]);
//     pp_limiter_scalar(thvm, thvp, s[1], 1e-4, 10*s[2]);

//     // Reconstruct other passive scalars
//     AMREX_D_TERM(, {   
//       s[0] = q(iv-iv_dir,QUT);
//       s[1] = q(iv,QUT);
//       s[2] = q(iv+iv_dir,QUT);
//       if (recon_scheme == 2) {
//         muscl(s, qr(iv,QUT), ql(iv+iv_dir,QUT), plm_theta);
//       } else {
//         wenoz3(s, qr(iv,QUT), ql(iv+iv_dir,QUT));
//       }
//     }, {
//       s[0] = q(iv-iv_dir,QUTT);
//       s[1] = q(iv,QUTT);
//       s[2] = q(iv+iv_dir,QUTT);
//       if (recon_scheme == 2) {
//         muscl(s, qr(iv,QUTT), ql(iv+iv_dir,QUTT), plm_theta);
//       } else {
//         wenoz3(s, qr(iv,QUTT), ql(iv+iv_dir,QUTT));
//       }
//     });
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       s[0] = q(iv-iv_dir,QFS+n);
//       s[1] = q(iv,QFS+n);
//       s[2] = q(iv+iv_dir,QFS+n);
      
//       if (recon_scheme == 2) {
//         muscl(s, qr(iv,QFS+n), ql(iv+iv_dir,QFS+n), plm_theta);
//       } else {
//         wenoz3(s, qr(iv,QFS+n), ql(iv+iv_dir,QFS+n));
//       }
//     }
//   } 
//   else if ((recon_scheme == 4) || (recon_scheme == 5)) { ///////// WENO-JS5 / WENO-Z5 /////////
//     Real s[5]; //stencil

//     // Reconstruct euler equation variables
//     for (int n = 0; n < 3; ++n) {
//       q2w(iv-2*iv_dir, n, QUN, q, sys, s[0]);
//       q2w(iv-iv_dir,   n, QUN, q, sys, s[1]);
//       q2w(iv,          n, QUN, q, sys, s[2]);
//       q2w(iv+iv_dir,   n, QUN, q, sys, s[3]);
//       q2w(iv+2*iv_dir, n, QUN, q, sys, s[4]);
//       if (recon_scheme == 4) {
//         wenojs5(s, wm[n], wp[n]);
//       } else {
//         wenoz5(s, wm[n], wp[n]);
//       }      
//     }
//     // Reconstruct the thermodynamic variable for Jacobian
//     s[0] = q(iv-2*iv_dir,QTHV);
//     s[1] = q(iv-iv_dir,QTHV);
//     s[2] = q(iv,QTHV);
//     s[3] = q(iv+iv_dir,QTHV);
//     s[4] = q(iv+2*iv_dir,QTHV);
//     if (recon_scheme == 4) {
//       wenojs5(s, thvm, thvp);
//     } else {
//       wenoz5(s, thvm, thvp);
//     }
//     pp_limiter_scalar(thvm, thvp, s[2], 1e-4, 10*s[1]);
//     pp_limiter_scalar(thvm, thvp, s[1], 1e-4, 10*s[3]);
//     // thvm = thvp = q(iv,QTHV);

//     // Reconstruct other passive scalars
//     AMREX_D_TERM(, {
//       s[0] = q(iv-2*iv_dir,QUT);
//       s[1] = q(iv-iv_dir,QUT);
//       s[2] = q(iv,QUT);
//       s[3] = q(iv+iv_dir,QUT);
//       s[4] = q(iv+2*iv_dir,QUT);
//       if (recon_scheme == 4) {
//         wenojs5(s, qr(iv,QUT), ql(iv+iv_dir,QUT));
//       } else {
//         wenoz5(s, qr(iv,QUT), ql(iv+iv_dir,QUT));
//       }
//     }, {
//       s[0] = q(iv-2*iv_dir,QUTT);
//       s[1] = q(iv-iv_dir,QUTT);
//       s[2] = q(iv,QUTT);
//       s[3] = q(iv+iv_dir,QUTT);
//       s[4] = q(iv+2*iv_dir,QUTT);
//       if (recon_scheme == 4) {
//         wenojs5(s, qr(iv,QUTT), ql(iv+iv_dir,QUTT));
//       } else {
//         wenoz5(s, qr(iv,QUTT), ql(iv+iv_dir,QUTT));
//       }
//     });
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       s[0] = q(iv-2*iv_dir,QFS+n);
//       s[1] = q(iv-iv_dir,QFS+n);
//       s[2] = q(iv,QFS+n);
//       s[3] = q(iv+iv_dir,QFS+n);
//       s[4] = q(iv+2*iv_dir,QFS+n);
//       if (recon_scheme == 4) {
//         wenojs5(s, qr(iv,QFS+n), ql(iv+iv_dir,QFS+n));
//       } else {
//         wenoz5(s, qr(iv,QFS+n), ql(iv+iv_dir,QFS+n));
//       }
//     }
//   } 
//   else {
//     amrex::Abort("Unknown reconstruction type");
//   } //end of recon

//   if (recon_scheme != 1) {
//     // Convert to primitive
//     w2q(iv, QUN, wm, thvm, sys, qr);
//     w2q(iv+iv_dir, QUN, wp, thvp, sys, ql);

//     // Enforce sum(Y) = 1
//     Real sumYl = 0.0, sumYr = 0.0;
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       ql(iv+iv_dir,QFS+n) = amrex::max<amrex::Real>(0.0, ql(iv+iv_dir,QFS+n));
//       sumYl += ql(iv+iv_dir,QFS+n);
      
//       qr(iv,QFS+n) = amrex::max<amrex::Real>(0.0, qr(iv,QFS+n));
//       sumYr += qr(iv,QFS+n);
//     }
//     if (std::abs(sumYl - 1.0) > parm.smally) {
//       Real invsumY = 1.0 / sumYl;
//       for (int n = 0; n < NUM_SPECIES; ++n) {
//         ql(iv+iv_dir,QFS+n) *= invsumY;
//       }
//     }
//     if (std::abs(sumYr - 1.0) > parm.smally) {
//       Real invsumY = 1.0 / sumYr;
//       for (int n = 0; n < NUM_SPECIES; ++n) {
//         qr(iv,QFS+n) *= invsumY;
//       }
//     }

//     // Prositivity preserving limiter
//     if (do_pp_limiter) {
//       pp_limiter_euler(iv+iv_dir, iv, ql, qr, q, parm.smallr, parm.smallp);
//     }
//   }
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pure_riemann (int i, int j, int k, int dir,
              amrex::Array4<amrex::Real> const& adv_flx,
              amrex::Array4<amrex::Real> const& /*p_flx*/,
              amrex::Array4<const amrex::Real> const& ql,
              amrex::Array4<const amrex::Real> const& qr,
              Parm const& parm) noexcept
{
  using amrex::Real;

  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

  int UMN, UMT, UMTT, QUN, QUT, QUTT;
  if (dir == 0) {
    UMN = UMX; UMT = UMY; UMTT = UMZ;  QUN = QU;  QUT = QV;  QUTT = QW;
  } else if (dir == 1) {
    UMN = UMY; UMT = UMX; UMTT = UMZ;  QUN = QV;  QUT = QU;  QUTT = QW;
  } else {
    UMN = UMZ; UMT = UMX; UMTT = UMY;  QUN = QW;  QUT = QU;  QUTT = QV;
  }

  auto eos = pele::physics::PhysicsType::eos();

  Real rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, yl[NUM_SPECIES], cl, el;
  rl = ql(iv,QRHO);
  AMREX_D_TERM(ul   = ql(iv,QUN); ,
               ut1l = ql(iv,QUT); ,
               ut2l = ql(iv,QUTT););
  pl = ql(iv,QPRES);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    yl[n] = ql(iv,QFS+n);
  }
  eos.RYP2E(rl, yl, pl, el);
  el += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));   
  eos.RPY2Cs(rl, pl, yl, cl);

  Real rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, yr[NUM_SPECIES], cr, er;
  rr = qr(iv,QRHO);
  AMREX_D_TERM(ur   = qr(iv,QUN); ,
               ut1r = qr(iv,QUT); ,
               ut2r = qr(iv,QUTT););
  pr = qr(iv,QPRES);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    yr[n] = qr(iv,QFS+n);
  }
  eos.RYP2E(rr, yr, pr, er);
  er += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
  eos.RPY2Cs(rr, pr, yr, cr);

  amrex::Real ry_flx[NUM_SPECIES] = {0.0};

  hllc(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
       rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
       adv_flx(i,j,k,URHO), 
       AMREX_D_DECL(adv_flx(i,j,k,UMN), adv_flx(i,j,k,UMT), adv_flx(i,j,k,UMTT)), 
       adv_flx(i,j,k,UEDEN), ry_flx,
       parm);

  for (int n = 0; n < NUM_SPECIES; ++n) { 
    adv_flx(i,j,k,UFS+n) = ry_flx[n]; 
  }

  if (isnan(adv_flx(i,j,k,URHO))) {
    amrex::Print() << dir << " flx(" << i << "," << j << "," << k << ") is nan \n";
    amrex::Print() << rl << " " << ul << " " << ut1l << " " << pl << " " << cl << " " << el << " " << yl[0] << " " <<
                      rr << " " << ur << " " << ut1r << " " << pr << " " << cr << " " << er << " " << yr[0] << "\n";
  }
}

// // /** 
// //  * \brief Reconstruction and Riemann solver all-in-one.
// //  */
// // void
// // recon_riemann_aio (int i, int j, int k, int dir,
// //                    amrex::Array4<amrex::Real const> const& q,
// //                    amrex::Array4<amrex::Real> const& adv_flx,
// //                    amrex::Array4<amrex::Real> const& p_flx,
// //                    const int recon_scheme, const amrex::Real plm_theta,
// //                    Parm const& parm) noexcept
// // {
// //     using amrex::Real;

// //     const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
// //     const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

// //     int UMN, UMT, UMTT, QUN, QUT, QUTT;
// //     if (dir == 0) {
// //         UMN = UMX; UMT = UMY; UMTT = UMZ;  
// //         QUN = QU;  QUT = QV;  QUTT = QW;
// //     } else if (dir == 1) {
// //         UMN = UMY; UMT = UMX; UMTT = UMZ;  
// //         QUN = QV;  QUT = QU;  QUTT = QW;
// //     } else {
// //         UMN = UMZ; UMT = UMX; UMTT = UMY;  
// //         QUN = QW;  QUT = QU;  QUTT = QV;
// //     }

// //     // Reconstruction
// //     Real rr, ur, utr, uttr, pr, er, yr[NUM_SPECIES], cr;
// //     Real rl, ul, utl, uttl, pl, el, yl[NUM_SPECIES], cl;

// //     Real rho, cs, pres; //jacobian variables
// //     rho = 0.5*(q(iv-iv_dir,QRHO) + q(iv,QRHO));
// //     pres = 0.5*(q(iv-iv_dir,QPRES) + q(iv,QPRES));
// //     cs = 0.5*(q(iv-iv_dir,QC) + q(iv,QC));
// //     // rho = sqrt(q(iv,QRHO) * q(iv-iv_dir,QRHO));
// //     // Real rp = sqrt(q(iv,QRHO) / q(iv-iv_dir,QRHO));
// //     // pres = (q(iv-iv_dir,QPRES) + q(iv,QPRES)*rp)/(1+rp);    
// //     // cs = (q(iv-iv_dir,QC) + q(iv,QC)*rp)/(1+rp);
    
// //     Real dummy;
// //     Real w[3];

// //     if (recon_scheme != 5) {
// //       amrex::Abort("Not implemented");
// //     } else {
// //         // Default: WENO-Z5
// //         Real s[5];   

// //         // Right
// //         for (int n = 0; n < 3; ++n) {
// //             q2w(iv-2*iv_dir, n, QUN, q, s[0], rho, cs, pres);
// //             q2w(iv-iv_dir, n, QUN, q, s[1], rho, cs, pres);
// //             q2w(iv, n, QUN, q, s[2], rho, cs, pres);
// //             q2w(iv+iv_dir, n, QUN, q, s[3], rho, cs, pres);       
// //             q2w(iv+2*iv_dir, n, QUN, q, s[4], rho, cs, pres);     
// //             wenoz5(s, w[n], dummy);
// //         }
// //         cns_unpackchar(w[0], w[1], w[2], cs, rho, dummy, rr, ur, pr, 0, parm);
// //         AMREX_D_TERM(,
// //         {   
// //             s[0] = q(iv-2*iv_dir,QUT);
// //             s[1] = q(iv-iv_dir,QUT);
// //             s[2] = q(iv,QUT);
// //             s[3] = q(iv+iv_dir,QUT);
// //             s[4] = q(iv+2*iv_dir,QUT);
// //             wenoz5(s, utr, dummy);
// //         }
// //         ,
// //         {
// //             s[0] = q(iv-2*iv_dir,QUTT);
// //             s[1] = q(iv-iv_dir,QUTT);
// //             s[2] = q(iv,QUTT);
// //             s[3] = q(iv+iv_dir,QUTT);
// //             s[4] = q(iv+2*iv_dir,QUTT);
// //             wenoz5(s, uttr, dummy);
// //         });
// //         for (int n = 0; n < NUM_SPECIES; ++n) {
// //             s[0] = q(iv-2*iv_dir,QFS+n);
// //             s[1] = q(iv-iv_dir,QFS+n);
// //             s[2] = q(iv,QFS+n);
// //             s[3] = q(iv+iv_dir,QFS+n);
// //             s[4] = q(iv+2*iv_dir,QFS+n);
// //             wenoz5(s, yr[n], dummy);
// //         }

// //         // Left
// //         for (int n = 0; n < 3; ++n) {
// //             q2w(iv-3*iv_dir, n, QUN, q, s[0], rho, cs, pres);
// //             q2w(iv-2*iv_dir, n, QUN, q, s[1], rho, cs, pres);
// //             q2w(iv-iv_dir, n, QUN, q, s[2], rho, cs, pres);
// //             q2w(iv, n, QUN, q, s[3], rho, cs, pres);       
// //             q2w(iv+iv_dir, n, QUN, q, s[4], rho, cs, pres);          
// //             wenoz5(s, dummy, w[n]);            
// //         }
// //         cns_unpackchar(w[0], w[1], w[2], cs, rho, dummy, rl, ul, pl, 0, parm);
// //         AMREX_D_TERM(,
// //         {   
// //             s[0] = q(iv-3*iv_dir,QUT);
// //             s[1] = q(iv-2*iv_dir,QUT);
// //             s[2] = q(iv-iv_dir,QUT);
// //             s[3] = q(iv,QUT);
// //             s[4] = q(iv+iv_dir,QUT);
// //             wenoz5(s, dummy, utl);
// //         }
// //         ,
// //         {
// //             s[0] = q(iv-3*iv_dir,QUTT);
// //             s[1] = q(iv-2*iv_dir,QUTT);
// //             s[2] = q(iv-iv_dir,QUTT);
// //             s[3] = q(iv,QUTT);
// //             s[4] = q(iv+iv_dir,QUTT);
// //             wenoz5(s, dummy, uttl);
// //         });
// //         for (int n = 0; n < NUM_SPECIES; ++n) {
// //             s[0] = q(iv-3*iv_dir,QFS+n);
// //             s[1] = q(iv-2*iv_dir,QFS+n);
// //             s[2] = q(iv-iv_dir,QFS+n);
// //             s[3] = q(iv,QFS+n);
// //             s[4] = q(iv+iv_dir,QFS+n);
// //             wenoz5(s, dummy, yl[n]);
// //         }

// //     } //end of recon

// //     auto eos = pele::physics::PhysicsType::eos();

// //     eos.RYP2E(rr, yr, pr, er);
// //     er += 0.5*(AMREX_D_TERM(ur*ur, +utr*utr, +uttr*uttr));
// //     eos.RPY2Cs(rr, pr, yr, cr);

// //     eos.RYP2E(rl, yl, pl, el);
// //     el += 0.5*(AMREX_D_TERM(ul*ul, +utl*utl, +uttl*uttl));
// //     eos.RPY2Cs(rl, pl, yl, cl);

// //     Real ry_flx[NUM_SPECIES] = {0.0};

// //     hllc(rl, AMREX_D_DECL(ul, utl, uttl), pl, cl, el, yl,
// //          rr, AMREX_D_DECL(ur, utr, uttr), pr, cr, er, yr,
// //          adv_flx(i,j,k,URHO), 
// //          AMREX_D_DECL(adv_flx(i,j,k,UMN), adv_flx(i,j,k,UMT), adv_flx(i,j,k,UMTT)), 
// //         //  p_flx(i,j,k,UMN), p_flx(i,j,k,UEDEN),
// //          adv_flx(i,j,k,UEDEN),
// //          ry_flx, //&adv_flx(i,j,k,UFS),
// //          parm);

// //     for (int n = 0; n < NUM_SPECIES; ++n) { adv_flx(i,j,k,UFS+n) = ry_flx[n]; }
// // }


#endif