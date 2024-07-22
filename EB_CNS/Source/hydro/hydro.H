#ifndef CNS_HYDRO_H_
#define CNS_HYDRO_H_

#include <AMReX_FArrayBox.H>

#include <cmath>

#include "PelePhysics.H"
#include "index_macros.H"

// #define MONITOR_PROBLEM_CELL

/**
 * \brief Convert state variable to primitive variable
 *
 * \note We do not do clipping here because CNS::enforce_consistent_state will be
 * called before this.
 *
 * @param i,j,k   x, y, z index.
 * @param ns      starting component index for u.
 * @param u       state variables array.
 * @param[out] q  primitive variables array.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_ctoprim(int i, int j, int k, int ns, amrex::Array4<amrex::Real const> const& u,
            amrex::Array4<amrex::Real> const& q) noexcept
{
  amrex::Real rho = u(i, j, k, ns + URHO);
  amrex::Real rhoinv = amrex::Real(1.0) / rho;
  AMREX_D_TERM(amrex::Real ux = u(i, j, k, ns + UMX) * rhoinv;
               , amrex::Real uy = u(i, j, k, ns + UMY) * rhoinv;
               , amrex::Real uz = u(i, j, k, ns + UMZ) * rhoinv;);
  amrex::Real ei = u(i, j, k, ns + UEDEN) * rhoinv -
                   amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));

  amrex::Real Y[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y[n] = u(i, j, k, ns + UFS + n) * rhoinv;
    q(i, j, k, QFS + n) = Y[n];
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real T = u(i, j, k, ns + UTEMP); // Necessary for SRK EOS
  amrex::Real p, cs, gamma, wbar;
  eos.REY2T(rho, ei, Y, T);
  eos.RTY2P(rho, T, Y, p);
  eos.Y2WBAR(Y, wbar);
  eos.RTY2Cs(rho, T, Y, cs);
  eos.RTY2G(rho, T, Y, gamma);

  q(i, j, k, QRHO) = rho;
  AMREX_D_TERM(q(i, j, k, QU) = ux;, q(i, j, k, QV) = uy;, q(i, j, k, QW) = uz;);
  q(i, j, k, QPRES) = p;
  q(i, j, k, QC) = cs;
  q(i, j, k, QG) = gamma;
  q(i, j, k, QEINT) = ei;
  q(i, j, k, QTEMP) = T;
  q(i, j, k, QWBAR) = wbar;

#ifdef MONITOR_PROBLEM_CELL
  // Monitor problem cell
  // std::ostringstream ss;
  std::ostream& ss = std::cout;

  if ((T <= 0.0) || amrex::isnan(cs) || amrex::isnan(ux) || (rho <= 0.0) ||
      (p <= 0.0)) {
    ss << "cns_ctoprim @" << i << "," << j << "," << k << "(ns=" << ns
       << "): rho=" << rho << " T=" << u(i, j, k, ns + UTEMP) << "/" << T
       << " p=" << p << " cs=" << cs << " gamma=" << gamma << " u=" << ux
       << " rhoE=" << u(i, j, k, ns + UEDEN) << " ei=" << ei << '\n';
    for (int n = 0; n < NUM_SPECIES; ++n) {
      ss << " rhoY" << n << "=" << u(i, j, k, ns + UFS + n);
    }
    ss << '\n';

    amrex::Abort();
    // BL_BACKTRACE_PUSH(ss.str()); // PUSH takes std::string
  }
#endif
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
cns_ctochar(int i, int j, int k, int dir, amrex::Array4<const amrex::Real> const& q,
            amrex::Array4<amrex::Real> const& w, const int sys)
{
  using amrex::Real;

  int QUN, QUT, QUTT;
  if (dir == 0) {
    QUN = QU;
    QUT = QV;
    QUTT = QW;
  } else if (dir == 1) {
    QUN = QV;
    QUT = QU;
    QUTT = QW;
  } else {
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  if (sys == 0) { // speed of sound system
    w(i, j, k, WRHO) = q(i, j, k, QRHO) - q(i, j, k, QPRES) / q(i, j, k, QC) /
                                            q(i, j, k, QC); // rho minus p over c^2
    w(i, j, k, WACO) = 0.5 * (q(i, j, k, QPRES) / q(i, j, k, QC) +
                              q(i, j, k, QRHO) * q(i, j, k, QUN));
    w(i, j, k, WACO + 1) = 0.5 * (q(i, j, k, QPRES) / q(i, j, k, QC) -
                                  q(i, j, k, QRHO) * q(i, j, k, QUN));
    w(i, j, k, WTHER) = q(i, j, k, QC);
  } else if (sys == 1) { // gamma system
    w(i, j, k, WRHO) = q(i, j, k, QRHO) * (1.0 - 1.0 / q(i, j, k, QG));
    w(i, j, k, WACO) =
      0.5 * (q(i, j, k, QPRES) +
             std::sqrt(q(i, j, k, QG) * q(i, j, k, QRHO) * q(i, j, k, QPRES)) *
               q(i, j, k, QUN));
    w(i, j, k, WACO + 1) =
      0.5 * (q(i, j, k, QPRES) -
             std::sqrt(q(i, j, k, QG) * q(i, j, k, QRHO) * q(i, j, k, QPRES)) *
               q(i, j, k, QUN));
    w(i, j, k, WTHER) = q(i, j, k, QG);
  } else {
    amrex::Abort(
      "Unrecognised thermodynamic system. Choose between sos (0) and gamma (1).");
  }

  // Passive scalars are the same for both systems
  for (int n = 0; n < NUM_SPECIES; ++n) { w(i, j, k, WY + n) = q(i, j, k, QFS + n); }
  AMREX_D_TERM(, w(i, j, k, WUT) = q(i, j, k, QUT);
               , w(i, j, k, WUT + 1) = q(i, j, k, QUTT););
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
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_unpackchar(const amrex::Real w1, const amrex::Real w2, const amrex::Real w3,
               const amrex::Real c, const amrex::Real /*r*/, const amrex::Real gamma,
               amrex::Real& rho, amrex::Real& u, amrex::Real& p, const int sys)
{
  if (sys == 0) { // speed of sound system
    rho = w1 + (w2 + w3) / c;
    u = (w2 - w3) / rho;
    p = (w2 + w3) * c;
  } else if (sys == 1) { // gamma system
    p = w2 + w3;
    rho = w1 / (1.0 - 1.0 / gamma);
    u = (w2 - w3) / std::sqrt(gamma * rho * p);
  } else {
    amrex::Abort(
      "Unrecognised thermodynamic system. Choose between sos (0) and gamma (1).");
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
cns_div(int i, int j, int k, int n, amrex::Array4<amrex::Real> const& dsdt,
        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fx,
                     amrex::Array4<amrex::Real const> const& fy,
                     amrex::Array4<amrex::Real const> const& fz),
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  dsdt(i, j, k, n) = AMREX_D_TERM(dxinv[0] * (fx(i, j, k, n) - fx(i + 1, j, k, n)),
                                  +dxinv[1] * (fy(i, j, k, n) - fy(i, j + 1, k, n)),
                                  +dxinv[2] * (fz(i, j, k, n) - fz(i, j, k + 1, n)));

#ifdef MONITOR_PROBLEM_CELL
  // Monitor problem cell
  if (isnan(dsdt(i, j, k, n))) {
    std::cout << "div(" << i << "," << j << "," << k << "," << n << ") is nan\n";
    std::cout << "fx=" << fx(i, j, k, n) << " " << fx(i + 1, j, k, n) << '\n';
#if AMREX_SPACEDIM > 1
    std::cout << "fy=" << fy(i, j, k, n) << " " << fy(i, j + 1, k, n) << '\n';
#endif
#if AMREX_SPACEDIM == 3
    std::cout << "fz=" << fz(i, j, k, n) << " " << fz(i, j, k + 1, n) << '\n';
#endif
    amrex::Abort();
  }
#endif
}

/**
 * @brief Compute flux function for the Euler system.
 * @param[in] iv Cell indices.
 * @param[in] dir Direction of the interface.
 * @param[in] q Primitive variables.
 * @param[out] f Flux function.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_flux(amrex::IntVect iv, int dir, amrex::Array4<const amrex::Real> const& q,
         amrex::Real f[NVAR]) noexcept
{
  using amrex::Real;
  const int QUN = QU + dir;

  f[URHO] = q(iv, QRHO) * q(iv, QUN);
  f[UMX] = q(iv, QRHO) * q(iv, QU) * q(iv, QUN) + (dir == 0 ? q(iv, QPRES) : 0.0);
  if constexpr (amrex::SpaceDim > 1) {
    f[UMY] = q(iv, QRHO) * q(iv, QV) * q(iv, QUN) + (dir == 1 ? q(iv, QPRES) : 0.0);
  } else {
    f[UMY] = 0.0;
  }
  if constexpr (amrex::SpaceDim == 3) {
    f[UMZ] = q(iv, QRHO) * q(iv, QW) * q(iv, QUN) + (dir == 2 ? q(iv, QPRES) : 0.0);
  } else {
    f[UMZ] = 0.0;
  }
  Real ret =
    q(iv, QRHO) * (q(iv, QEINT) +
                   0.5 * (AMREX_D_TERM(q(iv, QU) * q(iv, QU), +q(iv, QV) * q(iv, QV),
                                       +q(iv, QW) * q(iv, QW))));
  f[UEDEN] = q(iv, QUN) * (q(iv, QPRES) + ret);
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    f[UFS + ns] = q(iv, QRHO) * q(iv, QUN) * q(iv, QFS + ns);
  }
  f[UTEMP] = 0.0;
}

/// \brief Original Riemann solver in CNS example (Colella & Glaz)
AMREX_GPU_DEVICE inline void colella_glaz(
  const amrex::Real rl, const amrex::Real ul, const amrex::Real ut1l,
  const amrex::Real ut2l, const amrex::Real pl, const amrex::Real cl,
  const amrex::Real /*el*/, const amrex::Real* /*yl[NUM_SPECIES]*/,
  const amrex::Real rr, const amrex::Real ur, const amrex::Real ut1r,
  const amrex::Real ut2r, const amrex::Real pr, const amrex::Real cr,
  const amrex::Real /*er*/, const amrex::Real* /*yr[NUM_SPECIES]*/,
  amrex::Real& flxrho, amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt,
  amrex::Real& flxe, amrex::Real flxy[NUM_SPECIES]) noexcept
{
  using amrex::Real;

  constexpr Real weakwv = Real(1.e-3);
  constexpr Real small = Real(1.e-6);
  constexpr Real gamma = Real(1.4);
  Real smallp = CNSConstants::smallp;

  Real wl = cl;
  Real wr = cr;
  Real cleft = wl / rl;
  Real cright = wr / rr;
  Real ccsmall = small * (cleft + cright);

  Real pstar = (wl * pr + wr * pl - wr * wl * (ur - ul)) / (wl + wr);
  pstar = amrex::max(pstar, smallp);
  Real pstnm1 = pstar;

  Real wlsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pl) + pstar) * rl;
  Real wrsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pr) + pstar) * rr;

  wl = std::sqrt(wlsq);
  wr = std::sqrt(wrsq);
  Real ustarp = ul - (pstar - pl) / wl;
  Real ustarm = ur + (pstar - pr) / wr;

  pstar = (wl * pr + wr * pl - wr * wl * (ur - ul)) / (wl + wr);
  pstar = amrex::max(pstar, smallp);

  Real ustar;
  for (int iter = 0; iter < 3; ++iter) {
    wlsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pl) + pstar) * rl;
    wrsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pr) + pstar) * rr;

    wl = Real(1.) / std::sqrt(wlsq);
    wr = Real(1.) / std::sqrt(wrsq);

    Real ustnm1 = ustarm;
    Real ustnp1 = ustarp;

    ustarm = ur - (pr - pstar) * wr;
    ustarp = ul + (pl - pstar) * wl;

    Real dpditer = std::abs(pstnm1 - pstar);
    Real zp = std::abs(ustarp - ustnp1);
    if (zp - weakwv * cleft < Real(0.0)) { zp = dpditer * wl; }
    Real zm = std::abs(ustarm - ustnm1);
    if (zm - weakwv * cright < Real(0.0)) { zm = dpditer * wr; }

    Real zz = zp + zm;
    Real denom = dpditer / amrex::max(zz, ccsmall);
    pstnm1 = pstar;
    pstar = pstar - denom * (ustarm - ustarp);
    pstar = amrex::max(pstar, smallp);
    ustar = Real(0.5) * (ustarm + ustarp);
  }

  Real ro, uo, po, sgnm, utrans1, utrans2;
  if (ustar > Real(0.)) {
    ro = rl;
    uo = ul;
    po = pl;
    sgnm = Real(1.);
    utrans1 = ut1l;
    utrans2 = ut2l;
  } else if (ustar < Real(0.)) {
    ro = rr;
    uo = ur;
    po = pr;
    sgnm = Real(-1.);
    utrans1 = ut1r;
    utrans2 = ut2r;
  } else {
    uo = Real(0.5) * (ur + ul);
    po = Real(0.5) * (pr + pl);
    ro = Real(2.) * (rl * rr) / (rl + rr);
    sgnm = Real(1.);
    utrans1 = Real(0.5) * (ut1l + ut1r);
    utrans2 = Real(0.5) * (ut2l + ut2r);
  }
  Real wosq = (Real(0.5) * (gamma - Real(1.)) * (pstar + po) + pstar) * ro;
  Real co = std::sqrt(gamma * po / ro);
  Real wo = std::sqrt(wosq);
  Real dpjmp = pstar - po;
  Real rstar = ro / (Real(1.) - ro * dpjmp / wosq);
  Real cstar = std::sqrt(gamma * pstar / rstar);
  Real spout = co - sgnm * uo;
  Real spin = cstar - sgnm * uo;
  if (pstar >= po) {
    spin = wo / ro - sgnm * uo;
    spout = spin;
  }
  Real ss = amrex::max(spout - spin, spout + spin);
  Real frac = Real(0.5) * (Real(1.) + (spin + spout) / amrex::max(ss, ccsmall));

  Real rgdnv, ugdnv, pgdnv;
  if (spout < Real(0.)) {
    rgdnv = ro;
    ugdnv = uo;
    pgdnv = po;
  } else if (spin >= Real(0.)) {
    rgdnv = rstar;
    ugdnv = ustar;
    pgdnv = pstar;
  } else {
    rgdnv = frac * rstar + (Real(1.) - frac) * ro;
    ugdnv = frac * ustar + (Real(1.) - frac) * uo;
    pgdnv = frac * pstar + (Real(1.) - frac) * po;
  }

  flxrho = rgdnv * ugdnv;
  flxu = rgdnv * ugdnv * ugdnv + pgdnv;
  flxut = rgdnv * ugdnv * utrans1;
  flxutt = rgdnv * ugdnv * utrans2;
  flxe = ugdnv * (Real(0.5) * rgdnv *
                    (ugdnv * ugdnv + utrans1 * utrans1 + utrans2 * utrans2) +
                  pgdnv / (gamma - Real(1.)) + pgdnv);
  for (int n = 0; n < NUM_SPECIES; ++n) flxy[n] = rgdnv * ugdnv;
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
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void hllc(
  const amrex::Real rl,
  AMREX_D_DECL(const amrex::Real ul, const amrex::Real ut1l, const amrex::Real ut2l),
  const amrex::Real pl, const amrex::Real cl, const amrex::Real el,
  amrex::Real yl[NUM_SPECIES], const amrex::Real rr,
  AMREX_D_DECL(const amrex::Real ur, const amrex::Real ut1r, const amrex::Real ut2r),
  const amrex::Real pr, const amrex::Real cr, const amrex::Real er,
  amrex::Real yr[NUM_SPECIES], amrex::Real& flxrho,
  AMREX_D_DECL(amrex::Real& flxu, amrex::Real& flxut, amrex::Real& flxutt),
  amrex::Real& flxrhoe, amrex::Real flxrhoy[NUM_SPECIES]) noexcept
{
  using amrex::Real;

  // Estimate wave speeds
  // Arithmetic average (Davis)
  Real sl = amrex::min(ul - cl, ur - cr);
  Real sr = amrex::max(ul + cl, ur + cr);
  // Roe-average (Einfeldt's HLLE)
  Real rp = sqrt(rr / rl);
  Real uroe = (ul + ur * rp) / (1. + rp);
  Real croe = (cl + cr * rp) / (1. + rp);
  sl = amrex::min(sl, uroe - croe);
  sr = amrex::max(sr, uroe + croe);

  if (sl > 0) {
    // flx_l
    flxrho = 0.;
    AMREX_D_TERM(flxu = rl * ul * ul + pl;, flxut = rl * ul * ut1l;
                 , flxutt = rl * ul * ut2l;);
    flxrhoe = ul * (rl * el + pl);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      flxrhoy[n] = rl * ul * yl[n];
      flxrho += flxrhoy[n];
    }

  } else if (sr < 0) {
    // flx_r
    flxrho = 0.;
    AMREX_D_TERM(flxu = rr * ur * ur + pr;, flxut = rr * ur * ut1r;
                 , flxutt = rr * ur * ut2r;);
    flxrhoe = ur * (rr * er + pr);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      flxrhoy[n] = rr * ur * yr[n];
      flxrho += flxrhoy[n];
    }

  } else {
    Real sstar = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) /
                 (rl * (sl - ul) - rr * (sr - ur)); // contact wave speed

    if (sstar >= 0) {
      // flx_l* = flx_l + sl * (q_l* - q_l)
      Real frac = (sl - ul) / (sl - sstar) - 1.;

      flxrho = 0.;
      AMREX_D_TERM(flxu = rl * ul * ul + pl + sl * rl * ((frac + 1.) * sstar - ul);
                   , flxut = rl * ul * ut1l + sl * rl * frac * ut1l;
                   , flxutt = rl * ul * ut2l + sl * rl * frac * ut2l;);
      flxrhoe = ul * (rl * el + pl) +
                sl * rl * (frac * el + (sstar - ul) * (sstar + pl / rl / (sl - ul)));
      for (int n = 0; n < NUM_SPECIES; ++n) {
        flxrhoy[n] = rl * ul * yl[n] + sl * rl * frac * yl[n];
        flxrho += flxrhoy[n];
      }

    } else {
      // flx_r* = flx_r + sr * (q_r* - q_r)
      Real frac = (sr - ur) / (sr - sstar) - 1.;

      flxrho = 0.;
      AMREX_D_TERM(flxu = rr * ur * ur + pr + sr * rr * ((frac + 1.) * sstar - ur);
                   , flxut = rr * ur * ut1r + sr * rr * frac * ut1r;
                   , flxutt = rr * ur * ut2r + sr * rr * frac * ut2r;);
      flxrhoe = ur * (rr * er + pr) +
                sr * rr * (frac * er + (sstar - ur) * (sstar + pr / rr / (sr - ur)));
      for (int n = 0; n < NUM_SPECIES; ++n) {
        flxrhoy[n] = rr * ur * yr[n] + sr * rr * frac * yr[n];
        flxrho += flxrhoy[n];
      }
    }
  }
}

/**
 * \brief Unpack reconstructed characteristic variables, then run the HLLC Riemann
 * solver.
 *
 * @param i,j,k         x, y, z index.
 * @param dir           direction, 0:x, 1:y, 2:z.
 * @param[out] adv_flx  output hyperbolic fluxes.
 * @param[out] p_flx    output separately pressure and energy fluxes.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann(
  int i, int j, int k, int dir, amrex::Array4<amrex::Real> const& adv_flx,
  //  amrex::Array4<amrex::Real> const& /*p_flx*/,
  amrex::Array4<const amrex::Real> const& q, amrex::Array4<amrex::Real> const& wl,
  amrex::Array4<amrex::Real> const& wr, const int sys,
  const bool recon_char_var) noexcept
{
  using amrex::Real;

  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

  int UMN, UMT, UMTT, QUN, QUT, QUTT;
  if (dir == 0) {
    UMN = UMX;
    UMT = UMY;
    UMTT = UMZ;
    QUN = QU;
    QUT = QV;
    QUTT = QW;
  } else if (dir == 1) {
    UMN = UMY;
    UMT = UMX;
    UMTT = UMZ;
    QUN = QV;
    QUT = QU;
    QUTT = QW;
  } else {
    UMN = UMZ;
    UMT = UMX;
    UMTT = UMY;
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  auto eos = pele::physics::PhysicsType::eos();

  if (recon_char_var &&
      (!std::isgreater(wl(iv, WTHER), 0.0) ||
       !std::isgreater(wr(iv, WTHER), 0.0))) { // isgreater() returns false with nan
    // Limit reconstructed thermodynamic variable
    wl(iv, WTHER) = (sys == 0) ? q(iv - iv_dir, QC) : q(iv - iv_dir, QG);
    wr(iv, WTHER) = (sys == 0) ? q(iv, QC) : q(iv, QG);
  }

  // Unpack wl to ql (rl, ul, utl, pl, cl, el, yl)
  Real g = (recon_char_var && sys == 1) ? wl(iv, WTHER) : q(iv - iv_dir, QG);
  Real cl = (recon_char_var && sys == 0) ? wl(iv, WTHER) : q(iv - iv_dir, QC);
  Real r = q(iv - iv_dir, QRHO);
  Real yl[NUM_SPECIES];
  Real rl, AMREX_D_DECL(ul, ut1l, ut2l), pl;

  Real sumY = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    yl[n] = amrex::max<Real>(0.0, wl(iv, WY + n));
    sumY += yl[n];
  }
  // if (std::abs(sumY - 1.0) > CNSConstants::smally) {
  Real invsumY = 1.0 / sumY;
  for (int n = 0; n < NUM_SPECIES; ++n) { yl[n] *= invsumY; }
  // }
  AMREX_D_TERM(, ut1l = wl(iv, WUT);, ut2l = wl(iv, WUT + 1););

  cns_unpackchar(wl(iv, WRHO), wl(iv, WACO), wl(iv, WACO + 1), cl, r, g, rl, ul, pl,
                 sys);

  // Convert wr to qr
  g = (recon_char_var && sys == 1) ? wr(iv, WTHER) : q(iv, QG);
  Real cr = (recon_char_var && sys == 0) ? wr(iv, WTHER) : q(iv, QC);
  r = q(iv, QRHO);
  Real yr[NUM_SPECIES];
  Real rr, AMREX_D_DECL(ur, ut1r, ut2r), pr;

  sumY = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    yr[n] = amrex::max<Real>(0.0, wr(iv, WY + n));
    sumY += yr[n];
  }
  // if (std::abs(sumY - 1.0) > CNSConstants::smally) {
  invsumY = 1.0 / sumY;
  for (int n = 0; n < NUM_SPECIES; ++n) { yr[n] *= invsumY; }
  // }
  AMREX_D_TERM(, ut1r = wr(iv, WUT);, ut2r = wr(iv, WUT + 1););

  cns_unpackchar(wr(iv, WRHO), wr(iv, WACO), wr(iv, WACO + 1), cr, r, g, rr, ur, pr,
                 sys);

  if (!std::isgreater(rl, 0.0) || !std::isgreater(pl, 0.0) ||
      !std::isgreater(rr, 0.0) || !std::isgreater(pr, 0.0) ||
      isnan(AMREX_D_TERM(ul + ur, +ut1l + ut1r, +ut2l + ut2r))) {
    // Revert to piecewise linear (Godunov) recon if any value is non-physical
    rl = q(iv - iv_dir, QRHO);
    AMREX_D_TERM(ul = q(iv - iv_dir, QUN);, ut1l = q(iv - iv_dir, QUT);
                 , ut2l = q(iv - iv_dir, QUTT););
    pl = q(iv - iv_dir, QPRES);
    for (int n = 0; n < NUM_SPECIES; ++n) { yl[n] = q(iv - iv_dir, QFS + n); }

    rr = q(iv, QRHO);
    AMREX_D_TERM(ur = q(iv, QUN);, ut1r = q(iv, QUT);, ut2r = q(iv, QUTT););
    pr = q(iv, QPRES);
    for (int n = 0; n < NUM_SPECIES; ++n) { yr[n] = q(iv, QFS + n); }
  }

  Real el, er;
  eos.RYP2E(rl, yl, pl, el);
  el += 0.5 * (AMREX_D_TERM(ul * ul, +ut1l * ut1l, +ut2l * ut2l));
  eos.RPY2Cs(rl, pl, yl, cl);

  eos.RYP2E(rr, yr, pr, er);
  er += 0.5 * (AMREX_D_TERM(ur * ur, +ut1r * ut1r, +ut2r * ut2r));
  eos.RPY2Cs(rr, pr, yr, cr);

  if (isnan(cl) || isnan(cr)) {
    cl = q(iv - iv_dir, QC);
    cr = q(iv, QC);
  }

  Real ry_flx[NUM_SPECIES] = {0.0};
  hllc(rl, AMREX_D_DECL(ul, ut1l, ut2l), pl, cl, el, yl, rr,
       AMREX_D_DECL(ur, ut1r, ut2r), pr, cr, er, yr, adv_flx(i, j, k, URHO),
       AMREX_D_DECL(adv_flx(i, j, k, UMN), adv_flx(i, j, k, UMT),
                    adv_flx(i, j, k, UMTT)),
       adv_flx(i, j, k, UEDEN), ry_flx);
  // colella_glaz(rl, ul, ut1l, 0.0, pl, cl, el, yl, rr, ur, ut1r, 0.0, pr, cr, er,
  // yr,
  //              adv_flx(i, j, k, URHO), adv_flx(i, j, k, UMN), adv_flx(i, j, k,
  //              UMT), adv_flx(i, j, k, UMTT), adv_flx(i, j, k, UEDEN), ry_flx);

  for (int n = 0; n < NUM_SPECIES; ++n) { adv_flx(i, j, k, UFS + n) = ry_flx[n]; }

#ifdef MONITOR_PROBLEM_CELL
  // Monitor problem cell
  for (int n = 0; n < NVAR; ++n) {
    if (isnan(adv_flx(i, j, k, n))) {
      std::cout << "flx_" << dir << "(" << i << "," << j << "," << k << "," << n
                << ") is nan\n";
    }
  }
#endif
}

// Add higher order derivatives correction term to achieve real high order
// Ref: Alternative finite difference WENO https://doi.org/10.1137/17M115757X
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_afd_correction(
  int i, int j, int k, int dir, amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& adv_flx) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

  constexpr int ng = 3;
  Real f[2 * ng][NVAR];
  for (int m = 0; m < 2 * ng; ++m) {
    cns_flux(iv + (m - ng) * iv_dir, dir, q, f[m]);
  }

  for (int n = 0; n < NVAR; ++n) {
    adv_flx(i, j, k, n) +=
      -1.0 / 1152.0 *
        (-5.0 * f[0][n] + 39.0 * f[1][n] - 34.0 * f[2][n] - 34.0 * f[3][n] +
         39.0 * f[4][n] - 5.0 * f[5][n]) // dx^2*fxx/24
      + 7.0 / 11520.0 *
          (f[0][n] - 3.0 * f[1][n] + 2.0 * f[2][n] + 2.0 * f[3][n] - 3.0 * f[4][n] +
           f[5][n]); // dx^4*fxxxx*7/5760
  }
}

#endif