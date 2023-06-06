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
 * \note We do not do clipping here because CNS::enforce_consistent_state will be called before this.
 * 
 * @param i,j,k   x, y, z index.
 * @param ns      starting component index for u.
 * @param u       state variables array.
 * @param[out] q  primitive variables array.
 * @param parm    parameters and constants (see parm.H).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_ctoprim (int i, int j, int k, int ns,
             amrex::Array4<amrex::Real const> const& u,
             amrex::Array4<amrex::Real> const& q,
             Parm const& parm) noexcept
{
    amrex::Real rho = u(i,j,k,ns+URHO);
    amrex::Real rhoinv = amrex::Real(1.0) / rho;
    AMREX_D_TERM(amrex::Real ux = u(i,j,k,ns+UMX)*rhoinv; ,
                 amrex::Real uy = u(i,j,k,ns+UMY)*rhoinv; ,
                 amrex::Real uz = u(i,j,k,ns+UMZ)*rhoinv;);
    amrex::Real ei = u(i,j,k,ns+UEDEN) * rhoinv - amrex::Real(0.5)*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));

    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = u(i,j,k,ns+UFS+n) * rhoinv;
        q(i,j,k,QFS+n) = massfrac[n];
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T, p, cs, gamma;
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2P(rho, T, massfrac, p);
    eos.RTY2Cs(rho, T, massfrac, cs);
    eos.RTY2G(rho, T, massfrac, gamma);
    
    q(i,j,k,QRHO) = rho;
    AMREX_D_TERM(q(i,j,k,QU) = ux; ,
                 q(i,j,k,QV) = uy; ,
                 q(i,j,k,QW) = uz;);    
    q(i,j,k,QPRES) = p; 
    q(i,j,k,QC)    = cs;
    q(i,j,k,QG)    = gamma;
    q(i,j,k,QEINT) = ei;
    q(i,j,k,QTEMP) = T;
    
    // Monitor problem cell
    if (amrex::isnan(cs) || amrex::isnan(ux) || (rho <= 0.0) || (p <= 0.0)) {
        std::cout << "@ " << i << "," << j << "," << k << ": ns=" << ns << ": rho=" << rho << " T=" << T << " p=" << p;// << "\n";
        std::cout << " rhoE=" << u(i,j,k,ns+UEDEN) << " ei=" << ei << "\n";
        for (int n = 0; n < NUM_SPECIES; n++) {
            std::cout << "rhoY" << n << "=" << u(i,j,k,ns+UFS+n) << " ";
        }
        std::cout << "\n";
    }
}

/**
 * \brief Convert primitive variables to characteristic variables.
 * 
 * @param i,j,k    x, y, z index.
 * @param dir      direction along characteristic decomposition. 0:x, 1:y, 2:z.
 * @param q        primitive variables array.
 * @param[out] w   characteristic variables array.
 * @param sys      sos system (0) or gamma system (1).
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
        w(i,j,k,WRHO)   = q(i,j,k,QRHO) - q(i,j,k,QPRES)/q(i,j,k,QC)/q(i,j,k,QC); //rho minus p over c^2
        w(i,j,k,WACO)   = 0.5 * (q(i,j,k,QPRES)/q(i,j,k,QC) + q(i,j,k,QRHO)*q(i,j,k,QUN));
        w(i,j,k,WACO+1) = 0.5 * (q(i,j,k,QPRES)/q(i,j,k,QC) - q(i,j,k,QRHO)*q(i,j,k,QUN));
    } 
    else if (sys == 1) { // gamma system
        w(i,j,k,WRHO)   = q(i,j,k,QRHO) * (1.0 - 1.0 / q(i,j,k,QG));
        w(i,j,k,WACO)   = 0.5 * (q(i,j,k,QPRES) + std::sqrt(q(i,j,k,QG)*q(i,j,k,QRHO)*q(i,j,k,QPRES))*q(i,j,k,QUN));
        w(i,j,k,WACO+1) = 0.5 * (q(i,j,k,QPRES) - std::sqrt(q(i,j,k,QG)*q(i,j,k,QRHO)*q(i,j,k,QPRES))*q(i,j,k,QUN));
    } else {
        amrex::Abort("Unrecognised thermodynamic system. Choose between sos (0) and gamma (1).");
    }

    // Passive scalars are the same for both systems
    for (int n = 0; n < NUM_SPECIES; ++n) {
        w(i,j,k,WY+n) = q(i,j,k,QFS+n);
    }
    AMREX_D_TERM(,
        w(i,j,k,WUT)   = q(i,j,k,QUT);,
        w(i,j,k,WUT+1) = q(i,j,k,QUTT););
}

/**
 * \brief Convert characteristic variables to primitive variables.
 * 
 * This is called before the Riemann solver.
 * 
 * @param w1    characteristic variables WRHO.
 * @param w2    characteristic variables WACO.
 * @param w3    characteristic variables WACO+1.
 * @param c     speed of sound estimate.
 * @param r     density estimate.
 * @param gamma gamma estimate.
 * @param[out] rho  unpacked density.
 * @param[out] u    unpacked (normal) velocity.
 * @param[out] p    unpacked pressure.
 * @param sys   sos system (0) or gamma system (1).
 * @param parm  parameters and constants (see parm.H).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_unpackchar (const amrex::Real w1, const amrex::Real w2, const amrex::Real w3, 
                const amrex::Real c, const amrex::Real /*r*/, const amrex::Real gamma, 
                amrex::Real& rho, amrex::Real& u, amrex::Real& p,
                const int sys, Parm const& parm)
{
    if (sys == 0) { // speed of sound system
        rho = w1 + (w2 + w3) / c;        
        u = (w2 - w3) / rho;
        p = (w2 + w3) * c;
    } 
    else if (sys == 1) { // gamma system
        p = w2 + w3;        
        rho = w1 / (1.0 - 1.0 / gamma); // rho = w1 + (w2 + w3) * r / gamma / p; 
        u = (w2 - w3) / std::sqrt(gamma*rho*p);
    } 
    else {
        amrex::Abort("Unrecognised thermodynamic system. Choose between sos (0) and gamma (1).");
    }
}

/**
 * \brief Compute flux divergence. dSdt = -div(flux)*area/volume.
 * 
 * @param i,j,k     x, y, z index.
 * @param n         state variable array index.
 * @param[out] dsdt output time derivatives of state variables.
 * @param fx,fy,fz  fluxes in x, y, z direction.
 * @param dxinv     array of inverse of cell size, 1/[dx,dy,dz].
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
    
    // Monitor problem cell
    if (isnan(dsdt(i,j,k,n))) 
        std::cout << "div(" << i << "," << j << "," << k << ") is nan \n";
}

/**
 * \brief HLLC approximate Riemann solver.
 * 
 * @param rl    left density.
 * @param ul    left normal velocity.
 * @param ut1l,ut2l left tangential velocities.
 * @param pl    left pressure.
 * @param cl    left speed of sound.
 * @param el    left specific total energy.
 * @param yl    left mass fraction array.
 * @param rr    right density.
 * @param ur    right normal velocity.
 * @param ut1r,ut2r  right tangential velocities.
 * @param pr    right pressure.
 * @param cr    right speed of sound.
 * @param er    right specific total energy.
 * @param yr    right mass fraction array.
 * @param[out] flxrho  output density flux (= sum mass fraction fluxes).
 * @param[out] flxu,flxut,flxutt  output momentum fluxes.
 * @param[out] flxrhoe output total energy fluxes.
 * @param[out] flxrhoy output mass fraction fluxes.
 * @param parm  parameters and constants (see parm.H).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
hllc (const amrex::Real rl, 
      AMREX_D_DECL(const amrex::Real ul, const amrex::Real ut1l, const amrex::Real ut2l), 
      const amrex::Real pl, const amrex::Real cl, const amrex::Real el,
      amrex::Real yl[NUM_SPECIES],
      const amrex::Real rr, 
      AMREX_D_DECL(const amrex::Real ur, const amrex::Real ut1r, const amrex::Real ut2r), 
      const amrex::Real pr, const amrex::Real cr, const amrex::Real er,
      amrex::Real yr[NUM_SPECIES],
      amrex::Real& flxrho, 
      AMREX_D_DECL(amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt),
      amrex::Real& flxrhoe, 
      amrex::Real flxrhoy[NUM_SPECIES],
      Parm const& parm) noexcept
{
    using amrex::Real;

    // Estimate wave speeds
    // Arithmetic average (Davis), this alone won't work!
    Real sl = ul - cl;
    Real sr = ur + cr;
    // Roe-average (Einfeldt's HLLE)
    Real rp = sqrt(rr / rl);
    // AMREX_D_TERM(
    Real uroe = (ul + ur*rp)/(1. + rp);//,
    // Real ut1roe = (ut1l + ut1r*rp)/(1. + rp);,
    // Real ut2roe = (ut2l + ut1r*rp)/(1. + rp);)

    // auto eos = pele::physics::PhysicsType::eos();
    // Real ei, T, gl, gr, hl, hr;
    // ei = el - 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));
    // eos.REY2T(rl, ei, yl, T);
    // eos.RTY2G(rl, T, yl, gl);
    // eos.TY2H(T, yl, hl);
    // hl += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));
    // ei = er - 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
    // eos.REY2T(rr, ei, yr, T);
    // eos.RTY2G(rr, T, yr, gr);
    // eos.TY2H(T, yr, hr);
    // hr += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
    // Real croe = 0.5*(gl+gr-2)*((hl + hr*rp)/(1. + rp) - 0.5*(AMREX_D_TERM(uroe*uroe, +ut1roe*ut1roe, +ut2roe*ut2roe)));
    Real croe = (cl*cl + cr*cr*rp)/(1. + rp) + 0.5*rp/(1. + rp)/(1. + rp)
                *(AMREX_D_TERM((ul-ur)*(ul-ur), +(ut1l-ut1r)*(ut1l-ut1r), +(ut2l-ut2r)*(ut2l-ut2r)));
    croe = sqrt(croe); 
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
        
        if (sstar >= 0) { //parm.smallu*(amrex::Math::abs(ul) + amrex::Math::abs(ur))) {
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

        } else { // if (sstar < -parm.smallu*(amrex::Math::abs(ul) + amrex::Math::abs(ur))) {
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

        // } else { 
        //     //sstar ~ 0, return to central scheme
        //     Real rc = 0.5 * (rl + rr);
        //     Real uc = 0.5 * (ul + ur);
        //     Real pc = 0.5 * (pl + pr);
        //     flxrho = 0.;
        //     AMREX_D_TERM(
        //         flxu = rc*uc*uc + pc; ,
        //         flxut  = rc*uc*0.5*(ut1l+ut1r); ,
        //         flxutt = rc*uc*0.5*(ut2l+ut2r););
        //     flxrhoe = uc*(rc*0.5*(el+er) + pc);
        //     for (int n = 0; n < NUM_SPECIES; ++n) {
        //         flxrhoy[n] = rc*uc*0.5*(yl[n]+yr[n]);
        //         flxrho += flxrhoy[n];
        //     }
        }
    }
    
    // // Dellacherie's low Mach correction
    // Real Ml = sqrt(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l)) / cl;
    // Real Mr = sqrt(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r)) / cr;
    // Real thetaM = amrex::min( amrex::max(Ml, Mr), 1.0 ) - 1.0; // = (Mmax - 1), Xie
    // // Real thetaM = amrex::min( amrex::Math::abs(uroe)/croe, 1.0 ) - 1.0; // Dellacherie
    // Real rroe = sqrt(rl*rr);
    // AMREX_D_TERM(
    //     flxu += thetaM*rroe*croe*0.5*(ul-ur); ,
    //     flxut  += thetaM*rroe*croe*0.5*(ut1l-ut1r); ,
    //     flxutt += thetaM*rroe*croe*0.5*(ut2l-ut2r););
}

/**
 * \brief Unpack reconstructed characteristic variables, then run the HLLC Riemann solver.
 * 
 * @param i,j,k         x, y, z index.
 * @param dir           direction, 0:x, 1:y, 2:z.
 * @param[out] adv_flx  output hyperbolic fluxes.
 * @param[out] p_flx    output separately pressure and energy fluxes.
 * @param parm          parameters and constants (@see Parm).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_riemann (int i, int j, int k, int dir,
             amrex::Array4<amrex::Real> const& adv_flx,
             amrex::Array4<amrex::Real> const& /*p_flx*/,
             amrex::Array4<const amrex::Real> const& q,
             amrex::Array4<const amrex::Real> const& wl,
             amrex::Array4<const amrex::Real> const& wr,
             const int sys,
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

    auto eos = pele::physics::PhysicsType::eos();

    // Unpack wl to ql (rl, ul, utl, pl, cl, el, yl)
    Real g = q(iv-iv_dir,QG);
    Real r = q(iv-iv_dir,QRHO);
    Real cl = q(iv-iv_dir,QC);
    Real yl[NUM_SPECIES];
    Real rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, el;
    
    Real sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = amrex::max<Real>(0.0, wl(iv,WY+n));
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
    g = q(iv,QG);
    r = q(iv,QRHO);
    Real cr = q(iv,QC);
    Real yr[NUM_SPECIES];
    Real rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, er;

    sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        yr[n] = amrex::max<Real>(0.0, wr(iv,WY+n));
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

    eos.RYP2E(rl, yl, pl, el);
    el += 0.5*(AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));   
    eos.RPY2Cs(rl, pl, yl, cl);

    eos.RYP2E(rr, yr, pr, er);
    er += 0.5*(AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));
    eos.RPY2Cs(rr, pr, yr, cr);

    Real ry_flx[NUM_SPECIES] = {0.0};

    hllc(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl,
         rr, AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr,
         adv_flx(i,j,k,URHO), 
         AMREX_D_DECL(adv_flx(i,j,k,UMN), adv_flx(i,j,k,UMT), adv_flx(i,j,k,UMTT)), 
         adv_flx(i,j,k,UEDEN), ry_flx,
         parm);

    for (int n = 0; n < NUM_SPECIES; ++n) { adv_flx(i,j,k,UFS+n) = ry_flx[n]; }
}

#endif