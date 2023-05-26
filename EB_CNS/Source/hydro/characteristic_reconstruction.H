#ifndef CHARACTERISTIC_RECON_H
#define CHARACTERISTIC_RECON_H

using namespace amrex;

/**
 * \brief Convert primitive variables to characteristic.
 * \param rho (input) density
 * \param un  (input) noraml velocity
 * \param p   (input) pressure
 * \param sys (input) indicating what thermodynamic system to use (0: speed of sound, 1: gamma, 2: both)
 * \param cs  (input) speed of sound
 * \param ga  (input) gamma
 * \param u   (input) averaged u
 * \param n   (input) which component of characteristic variable (0: material_wave, 1: +acoustic_wave, 2: -acoutic_wave)
 * \param w   (output) chosen component of characteristic variables 
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
q2w (const Real rho, const Real  un,  const Real p,
     const int  sys, const Real  cs,  const Real ga,  const Real u, 
     const int  n,         Real& w) noexcept
{
  if (sys == 0) { // speed of sound system (thv == cs)
    if (n == 0) {
      w = rho - p/cs/cs; // material_wave
    } else if (n == 1) {
      w = 0.5 * (p/cs + rho*un) / cs; // +acoustic_wave
    } else if (n == 2) {
      w = 0.5 * (p/cs - rho*un) / cs; // -acoustic_wave
    }
  }
  else if (sys == 1) { // specific heat system (thv == gamma)
    if (n == 0) {
      w = rho - (1.0 - 1.0/ga); // material_wave
    } else if (n == 1) {
      w = 0.5*(p + std::sqrt(ga*rho*p)*un); // +acoustic_wave
    } else if (n == 2) {
      w = 0.5*(p - std::sqrt(ga*rho*p)*un); // -acoustic_wave
    }
  }
  else if (sys == 2) {
    const Real ooc2 = 1.0/cs/cs;
    const Real ga12 = 0.5*(ga-1);
    const Real uoc = u/cs;
    const Real rhoE = p/(ga-1) + 0.5*rho*un*un;

    if (n == 0) {
      w = 0.5*uoc*(ga12*uoc+1)*rho - (ga12*uoc+0.5)/cs*rho*un + ga12*ooc2*rhoE;
    } else if (n == 1) {
      w = (1-ga12*uoc*uoc)*rho + (2*ga12*u*ooc2)*rho*un - (ga-1)*ooc2*rhoE;
    } else if (n == 2) {
      w = 0.5*uoc*(ga12*uoc-1)*rho - (ga12*uoc-0.5)/cs*rho*un + ga12*ooc2*rhoE;
    }
  }
  else {
    Abort("Unrecognised thermodynamic system");
  }
}

/**
 * \brief Convert characteristic variables to primitive.
 * \param w   (input) characteristic variable array
 * \param sys (input) indicating what thermodynamic system to use (0: speed of sound, 1: gamma, 2: both)
 * \param cs  (input) speed of sound
 * \param ga  (input) gamma
 * \param u   (input) averaged u
 * \param iv  (input) InVect object (= {i,j,k})
 * \param QUN (input) index of un in q
 * \param q   (output) chosen component of characteristic variables 
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
w2q (const Real    w[3], const int sys, const Real cs, const Real ga, const Real u, 
     const IntVect iv,   const int QUN, Array4<Real> const& q) noexcept
{
  if (sys == 0) { // speed of sound system (thv == cs)
    q(iv,QRHO)  = w[0] + (w[1] + w[2]);
    q(iv,QPRES) = (w[1] + w[2]) * cs * cs;       
    q(iv,QUN)   = (w[1] - w[2]) / q(iv,QRHO) * cs;
  } 
  else if (sys == 1) { // specific heat system (thv == gamma)
    q(iv,QRHO)  = w[0] / (1.0 - 1.0 / ga);
    q(iv,QPRES) = w[1] + w[2];
    q(iv,QUN)   = (w[1] - w[2]) / std::sqrt(ga*q(iv,QRHO)*q(iv,QPRES));
  }
  else if (sys == 2) {
    const Real u2o2 = 0.5*u*u;
    const Real c2og = cs*cs/(ga-1);

    q(iv,QRHO) = w[0] + w[1] + w[2];    
    const Real u2 = (u-cs)*w[0] + u*w[1] + (u+cs)*w[2];
    const Real u3 = (u2o2+c2og-u*cs)*w[0] + u2o2*w[1] + (u2o2+c2og+u*cs)*w[2];
    
    q(iv,QUN) = u2 / q(iv,QRHO);
    q(iv,QPRES) = (u3 - 0.5*q(iv,QRHO)*q(iv,QUN)*q(iv,QUN))*(ga-1);
  }
  else {
    Abort("Unrecognised thermodynamic system");
  }
}

/**
 * \brief Reconstruct values on both sides of a face using TVD MUSCL. 
 * \param s  (input) stencil array [iv-3, iv-2, iv-1, iv, iv+1, iv+2]
 * \param sr (output) reconstructed state on the right of face iv (left of cell iv)
 * \param sl (output) reconstructed state on the left of face iv (right of cell iv-1)
 * \param plm_theta (input) parameter for MUSCL limiter
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
face_muscl (const Real s[4], Real& sr, Real& sl, const Real plm_theta) noexcept
{
  // Calculate slope
  Real dlft, drgt;
  dlft = s[2] - s[1];
  drgt = s[3] - s[2];

  Real dcen = Real(0.5)*(dlft+drgt);
  Real dsgn = amrex::Math::copysign(Real(1.0), dcen);
  Real slop = plm_theta * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
  Real dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
  Real dq = dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
  
  sr = s[2] - 0.5*dq;

  //
  dlft = s[1] - s[0];
  drgt = s[2] - s[1];

  dcen = Real(0.5)*(dlft+drgt);
  dsgn = amrex::Math::copysign(Real(1.0), dcen);
  slop = plm_theta * amrex::min(amrex::Math::abs(dlft),amrex::Math::abs(drgt));
  dlim = (dlft*drgt >= Real(0.0)) ? slop : Real(0.0);
  dq = dsgn * amrex::min(dlim,amrex::Math::abs(dcen));
  
  sl = s[1] + 0.5*dq;
}

/**
 * \brief Reconstruct values on both sides of a face using WENO-Z5. 
 * \param s  (input) stencil array [iv-3, iv-2, iv-1, iv, iv+1, iv+2]
 * \param sr (output) reconstructed state on the right of face iv (left of cell iv)
 * \param sl (output) reconstructed state on the left of face iv (right of cell iv-1)
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
face_wenoz5 (const Real s[6], Real& sr, Real& sl) noexcept
{
  const Real eps = 1e-10;

  Real vl[3];
  Real vr[3];
  Real beta[3];
  Real alpha[3];
  Real alpha1, tau;

  beta[2] = Real(13./12.) * pow(s[5] - 2.0 * s[4] + s[3], 2) +
              0.25 * pow(s[5] - 4.0 * s[4] + 3.0 * s[3], 2);
  beta[1] = Real(13./12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
              0.25 * (s[4] - s[2]) * (s[4] - s[2]);
  beta[0] = Real(13./12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
              0.25 * pow(3.0 * s[3] - 4.0 * s[2] + s[1], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[3] - 7.0 * s[4] + 2.0 * s[5];
  vr[1] = -s[4] + 5.0 * s[3] + 2.0 * s[2];
  vr[0] = 2.0 * s[3] + 5.0 * s[2] - s[1];

  sr = Real(1./6.) * alpha1 *
        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);

  beta[2] = Real(13./12.) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
              0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
              0.25 * (s[1] - s[3]) * (s[1] - s[3]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
              0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
  vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
  vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

  sl = Real(1./6.) * alpha1 *
        (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
wenoz5_m (Real s[5], Real& sm) noexcept
{
  const Real eps = 1e-6;
  Real vl[3], vr[3], beta[3], alpha[3], alpha1, tau;

  beta[2] = Real(13./12.) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
              0.25 * pow(s[4] - 4.0 * s[3] + 3.0 * s[2], 2);
  beta[1] = Real(13./12.) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
              0.25 * (s[3] - s[1]) * (s[3] - s[1]);
  beta[0] = Real(13./12.) * pow(s[2] - 2.0 * s[1] + s[0], 2) +
              0.25 * pow(3.0 * s[2] - 4.0 * s[1] + s[0], 2);

  tau = amrex::Math::abs(beta[2] - beta[0]);
  beta[2] = 1.0 + (tau / (eps + beta[2])) * (tau / (eps + beta[2]));
  beta[1] = 1.0 + (tau / (eps + beta[1])) * (tau / (eps + beta[1]));
  beta[0] = 1.0 + (tau / (eps + beta[0])) * (tau / (eps + beta[0]));

  alpha[2] = beta[2];
  alpha[1] = 6.0 * beta[1];
  alpha[0] = 3.0 * beta[0];
  alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

  vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
  vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
  vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];

  sm = Real(1./6.) * alpha1 *
        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}

/**
 * \brief Positivity-preserving limiter from Zhang and Shu 2011
*/
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pp_limiter_scalar (Real& qm, Real& qp, const Real& q, 
                   const Real m, const Real /*M*/)
{
  Real w = 1/6; // first and last Gauss-Lobatto weight of order (k+3)/2
  Real mid_p = (q - w*qm - w*qp) / (1 - 2*w); // some point of the polynomial in the cell
  Real min_q = amrex::min(amrex::min(qm, qp), mid_p);
  // Real max_q = amrex::max(amrex::max(qm, qp), mid_p);

  if ((min_q < m) /*|| (max_q > M)*/) {
    std::cout << "PP Limiter Scalar activated: q^+(i-1/2)=" << qm << " q^-(i+1/2)=" << qp << " q=" << q << "\n";
    
    Real t = amrex::Math::abs((q - m) / (q - min_q));
    // t = amrex::min(t, amrex::Math::abs((q - M) / (q - max_q)));
    AMREX_ALWAYS_ASSERT(t <= 1);

    qm = t*qm + (1-t)*q;
    qp = t*qp + (1-t)*q;
    std::cout << "New q^+(i-1/2)=" << qm << " q^-(i+1/2)=" << qp << "\n";
  }
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// pp_limiter_euler (const IntVect ivl, const IntVect iv, 
//                   Array4<Real> const& ql, 
//                   Array4<Real> const& qr, 
//                   Array4<const Real> const& q, 
//                   const Real mr, const Real mp)
// {
//   // Step 1: limit density as a scalar
//   pp_limiter_scalar(ql(ivl, QRHO), qr(iv, QRHO), q(iv, QRHO), mr, 100);

//   // Step 2: limit the whole state for positive pressure
//   if ((ql(ivl, QPRES) < mp) || (qr(iv, QPRES) < mp)) {
//     std::cout << iv << " PP Limiter Euler activated: pl=" << ql(ivl, QPRES) 
//                     << " pr=" << qr(iv, QPRES) << " p=" << q(iv, QPRES) << "\n";
//     Real min_p = amrex::min(ql(ivl, QPRES), qr(iv, QPRES));
//     Real t = amrex::Math::abs(q(iv, QPRES) - mp) / amrex::Math::abs(q(iv, QPRES) - min_p);
//     // Real min_p = amrex::min(ql(ivl, QPRES)/ql(ivl, QRHO), qr(iv, QPRES)/qr(iv, QRHO));
//     // Real t = amrex::Math::abs(q(iv, QPRES)/q(iv, QRHO) - mp/mr) 
//     //        / amrex::Math::abs(q(iv, QPRES)/q(iv, QRHO) - min_p);
//     AMREX_ASSERT(t <= 1);

//     for (int n = QRHO; n <= QPRES; ++n) {
//       std::cout << "Old ql=" << ql(ivl, n) << " Old qr=" << qr(iv, n);
//       ql(ivl, n) = t*ql(ivl, n) + (1-t)*q(iv, n);
//       qr(iv, n)  = t*qr(iv, n) + (1-t)*q(iv, n);
//       std::cout << "  New ql=" << ql(ivl, n) << " New qr=" << qr(iv, n) << "\n";
//     }

//     Real T, Y[NUM_SPECIES];
//     auto eos = pele::physics::PhysicsType::eos();
//     for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = ql(ivl,QFS+n);
//     eos.RYP2T(ql(ivl,QRHO), Y, ql(ivl, QPRES), T);
//     std::cout << "New Tl=" << T;
//     for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = qr(iv,QFS+n);
//     eos.RYP2T(qr(iv,QRHO), Y, qr(iv, QPRES), T);
//     std::cout << " New Tr=" << T << "\n";
//   }
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
pp_limiter_pres (const IntVect iv, 
                 Real& rhol, Real& ul, Real& pl, 
                 Real& rhor, Real& ur, Real& pr, 
                 Real  rho,  Real  u,  Real  p, 
                 const Real min_r, const Real min_p)
{
  // Step 1: limit density as a scalar
  pp_limiter_scalar(rhol, rhor, rho, min_r, 1000);

  // Step 2: limit the whole state for positive pressure
  if ((pl < min_p) || (pr < min_p)) {
    std::cout << iv << " PP Limiter Pres activated: pl=" << pl << " pr=" << pr << " p=" << p << "\n";
    Real min_p_lr = amrex::min(pl, pr);
    Real t = (min_p - p) / (min_p_lr - p);
    AMREX_ALWAYS_ASSERT(t <= 1);

    rhol = rho + t * (rhol - rho);
    ul = u + t * (ul - u);
    pl = p + t * (pl - p);
    std::cout << "  New rhol=" << rhol << " ul=" << ul << " pl=" << pl << "\n";

    rhor = rho + t * (rhor - rho);
    ur = u + t * (ur - u);
    pr = p + t * (pr - p);
    std::cout << "  New rhor=" << rhor << " ur=" << ur << " pr=" << pr << "\n";   
  }
}

/**
 * \brief Positivity-preserving limiter by returning invalid reconstruction to first order
 * \param ql    (in/output) reconstructed left state
 * \param ql    (in/output) reconstructed right state
 * \param qL    (input) left cell state
 * \param qR    (input) right cell state
 * \param min_q (input) minimum value for limiter to activate
*/
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
reduce_order_limiter (Real& ql, Real& qr, 
                      Real  qL, Real  qR,
                      const Real min_q) 
{
  if ((ql < min_q) || (qr < min_q)) {
    // Return to Godunov reconstruction
    ql = qL;
    qr = qR;
  }
}

/** 
 * \brief Reconstruct states at cell interface. To be used with pure_riemann().
 * \param i   (input) x-index
 * \param j   (input) y-index
 * \param k   (input) z-index
 * \param dir (input) direction along which reconstruction takes place (0: x, 1: y, 2: z)
 * \param q   (input) primitive variables array
 * \param ql  (output) reconstructive primitive variables on the left of the interface
 * \param qr  (output) reconstructive primitive variables on the right of the interface
 * \param recon_scheme (input) chosen reconstruction scheme
 * \param plm_theta (input) parameter in MUSCL limiter
 * \param sys (input) indicating what thermodynamics system to use (0: speed of sound, 1: gamma)
 * \param parm (input) Parm object
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
char_recon (int i, int j, int k, int dir,
            Array4<const Real> const& q,
            Array4<Real> const& ql, Array4<Real> const& qr,
            const int recon_scheme, 
            const Real plm_theta,
            const int sys,
            // const bool do_pp_limiter,
            Parm const& parm,
            Array4<const EBCellFlag> const* flag = nullptr) noexcept
{
  // Spatial indices
  const IntVect iv(AMREX_D_DECL(i, j, k));
  const IntVect iv_dir = IntVect::TheDimensionVector(dir);

  IntVect ivm3 = iv-3*iv_dir;
  IntVect ivm2 = iv-2*iv_dir;
  IntVect ivm1 = iv-1*iv_dir;
  IntVect ivp1 = iv+1*iv_dir;
  IntVect ivp2 = iv+2*iv_dir;
  if (flag != nullptr) {
    // Do not work on EB interfaces 
    if ((*flag)(iv).isCovered() || (*flag)(ivm1).isCovered()) return;

    // Check EB connectivity
    // Left
    bool cn1 = (*flag)(iv).isConnected(-iv_dir);
    bool cn2 = (*flag)(iv-iv_dir).isConnected(-iv_dir);
    bool cn3 = (*flag)(iv-2*iv_dir).isConnected(-iv_dir);        
    ivm1 = cn1 ? ivm1 : iv;
    ivm2 = cn1 && cn2 ? ivm2 : ivm1;
    ivm3 = cn1 && cn2 && cn3 ? ivm3 : ivm2;
    // Right
    bool cp1 = (*flag)(iv).isConnected(iv_dir);
    bool cp2 = (*flag)(iv+iv_dir).isConnected(iv_dir);    
    ivp1 = cp1 ? ivp1 : iv;
    ivp2 = cp1 && cp2 ? ivp2 : ivp1;    
  }

  // Indices of normal and transverse velocity
  int QUN, QUT, QUTT; 
  if      (dir == 0) { QUN = QU;  QUT = QV;  QUTT = QW; } 
  else if (dir == 1) { QUN = QV;  QUT = QU;  QUTT = QW; } 
  else if (dir == 2) { QUN = QW;  QUT = QU;  QUTT = QV; }
  
  // Averaged thermodynamic variable at face
  Real cs, ga, u; 
  // Simple avg for gamma (it is almost constant)
  ga = 0.5*(q(iv,QG) + q(ivm1,QG)); 
  // Roe avg for cs and u
  Real rp = sqrt(q(iv,QRHO) / q(ivm1,QRHO));
  // if (isnan(rp)) { amrex::Print() << "rp is nan rhop=" << q(iv,QRHO) << " rhom=" << q(ivm1,QRHO) << "\n"; rp = 1.0; }
  u  = (q(ivm1,QUN) + q(iv,QUN)*rp)/(1. + rp);
  cs = sqrt((q(ivm1,QC)*q(ivm1,QC) + q(iv,QC)*q(iv,QC)*rp)/(1. + rp));
  // ga = (q(ivm1,QG) + q(iv,QG)*rp)/(1. + rp);
        // + 0.5*rp/(1. + rp)/(1. + rp)*(AMREX_D_TERM((q(iv,QU)-q(ivm1,QU))*(q(iv,QU)-q(ivm1,QU)), 
        //                                           +(q(iv,QV)-q(ivm1,QV))*(q(iv,QV)-q(ivm1,QV)), 
        //                                           +(q(iv,QW)-q(ivm1,QW))*(q(iv,QW)-q(ivm1,QW)))));
  // Real u2 = (q(ivm1,QUN)*q(ivm1,QUN) + q(iv,QUN)*q(iv,QUN)*rp)/(1. + rp);
  // cs = sqrt((q(ivm1,QC)*q(ivm1,QC) + q(iv,QC)*q(iv,QC)*rp)/(1. + rp) + 0.5*(ga-1)*(u2 - u*u));
  if (isnan(rp)) amrex::Print() << "rp is nan. r"<<iv<<"=" << q(iv,QRHO) << " r(iv-1)=" << q(ivm1,QRHO) << " \n";
  if (isnan(ga)) amrex::Print() << "ga is nan. g"<<iv<<"=" << q(iv,QG) << " g(iv-1)=" << q(ivm1,QG) << " \n";
  if (isnan(u )) amrex::Print() << "u  is nan. u"<<iv<<"=" << q(iv,QUN) << " u(iv-1)=" << q(ivm1,QUN) << " \n";
  if (isnan(cs)) amrex::Print() << "cs is nan. c"<<iv<<"=" << q(iv,QC) << " c(iv-1)=" << q(ivm1,QC) << " \n";

  // Reconstruction q -> wl, wr
  Real wl[3], wr[3]; // arrays to hold w+ and w-

  if (recon_scheme == 1) { // Godunov (no need characteristics)
    ql(iv,QRHO) = q(ivm1,QRHO);  qr(iv,QRHO) = q(iv,QRHO);
    AMREX_D_TERM(ql(iv,QUN) = q(ivm1,QUN);   qr(iv,QUN) = q(iv,QUN); ,
                 ql(iv,QUT) = q(ivm1,QUT);   qr(iv,QUT) = q(iv,QUT); ,
                 ql(iv,QUTT) = q(ivm1,QUTT); qr(iv,QUTT) = q(iv,QUTT););
    ql(iv,QPRES) = q(ivm1,QPRES);  qr(iv,QPRES) = q(iv,QPRES); 
    for (int n = 0; n < NUM_SPECIES; ++n) {
      ql(iv,QFS+n) = q(ivm1,QFS+n);  qr(iv,QFS+n) = q(iv,QFS+n);
    }
  } 
  else if (recon_scheme == 2) { // MUSCL
    Real s[4]; //stencil

    for (int n = 0; n < 3; ++n) {
      // Convert q -> w
      Real un = (ivm2 == ivm1) ? -q(ivm1,QUN) : q(ivm2,QUN);
      q2w(q(ivm2,QRHO), un, q(ivm2,QPRES), sys, cs, ga, u, n, s[0]);
      un = (ivm1 == iv) ? -q(iv,QUN) : q(ivm1,QUN);
      q2w(q(ivm1,QRHO), un, q(ivm1,QPRES), sys, cs, ga, u, n, s[1]);      
      q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, cs, ga, u, n, s[2]);
      un = (ivp1 == iv) ? -q(iv,QUN) : q(ivp1,QUN);
      q2w(q(ivp1,QRHO), un, q(ivp1,QPRES), sys, cs, ga, u, n, s[3]);

      // q2w(q(ivm2,QRHO), q(ivm2,QUN), q(ivm2,QPRES), sys, cs, ga, u, n, s[0]);
      // q2w(q(ivm1,QRHO), q(ivm1,QUN), q(ivm1,QPRES), sys, cs, ga, u, n, s[1]);
      // q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, cs, ga, u, n, s[2]);
      // q2w(q(ivp1,QRHO), q(ivp1,QUN), q(ivp1,QPRES), sys, cs, ga, u, n, s[3]);
    
      // Reconstruct to i-1/2
      face_muscl(s, wr[n], wl[n], plm_theta);
    }

    // Reconstruct other passive scalars
    s[0] = q(ivm2,QUT);
    s[1] = q(ivm1,QUT);
    s[2] = q(iv,QUT);
    s[3] = q(ivp1,QUT);
    face_muscl(s, qr(iv,QUT), ql(iv,QUT), plm_theta);

    s[0] = q(ivm2,QUTT);
    s[1] = q(ivm1,QUTT);
    s[2] = q(iv,QUTT);
    s[3] = q(ivp1,QUTT);
    face_muscl(s, qr(iv,QUTT), ql(iv,QUTT), plm_theta);

    for (int n = 0; n < NUM_SPECIES; ++n) {
      s[0] = q(ivm2,QFS+n);
      s[1] = q(ivm1,QFS+n);
      s[2] = q(iv,QFS+n);
      s[3] = q(ivp1,QFS+n);
      face_muscl(s, qr(iv,QFS+n), ql(iv,QFS+n), plm_theta);
    }
  }
  else if (recon_scheme == 5) { // WENO-Z5
    Real s[6]; //stencil
  
    for (int n = 0; n < 3; ++n) {
      // Convert q -> w
      Real un = (ivm3 == ivm2) ? -q(ivm2,QUN) : q(ivm3,QUN);
      q2w(q(ivm3,QRHO), un, q(ivm3,QPRES), sys, cs, ga, u, n, s[0]);
      un = (ivm2 == ivm1) ? -q(ivm1,QUN) : q(ivm2,QUN);
      q2w(q(ivm2,QRHO), un, q(ivm2,QPRES), sys, cs, ga, u, n, s[1]);
      un = (ivm1 == iv) ? -q(iv,QUN) : q(ivm1,QUN);
      q2w(q(ivm1,QRHO), un, q(ivm1,QPRES), sys, cs, ga, u, n, s[2]);      
      q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, cs, ga, u, n, s[3]);
      un = (ivp1 == iv) ? -q(iv,QUN) : q(ivp1,QUN);
      q2w(q(ivp1,QRHO), un, q(ivp1,QPRES), sys, cs, ga, u, n, s[4]);
      un = (ivp2 == ivp1) ? -q(ivp1,QUN) : q(ivp2,QUN);
      q2w(q(ivp2,QRHO), un, q(ivp2,QPRES), sys, cs, ga, u, n, s[5]);

      // Real un = (ivm1 == iv) ? -q(iv,QUN) : q(iv-3*iv_dir,QUN);
      // q2w(q(ivm3,QRHO), un, q(ivm3,QPRES), sys, cs, ga, u, n, s[0]);
      // un = (ivm1 == iv) ? -q(iv,QUN) : q(iv-2*iv_dir,QUN);
      // q2w(q(ivm2,QRHO), un, q(ivm2,QPRES), sys, cs, ga, u, n, s[1]);
      // un = (ivm1 == iv) ? -q(iv,QUN) : q(iv-iv_dir,QUN);
      // q2w(q(ivm1,QRHO), un, q(ivm1,QPRES), sys, cs, ga, u, n, s[2]);      
      // q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, cs, ga, u, n, s[3]);
      // un = (ivp1 == iv) ? -q(iv,QUN) : q(iv+iv_dir,QUN);
      // q2w(q(ivp1,QRHO), un, q(ivp1,QPRES), sys, cs, ga, u, n, s[4]);
      // un = (ivp1 == iv) ? -q(iv,QUN) : q(iv+2*iv_dir,QUN);
      // q2w(q(ivp2,QRHO), un, q(ivp2,QPRES), sys, cs, ga, u, n, s[5]);

      // q2w(q(ivm3,QRHO), q(ivm3,QUN), q(ivm3,QPRES), sys, cs, ga, u, n, s[0]);
      // q2w(q(ivm2,QRHO), q(ivm2,QUN), q(ivm2,QPRES), sys, cs, ga, u, n, s[1]);
      // q2w(q(ivm1,QRHO), q(ivm1,QUN), q(ivm1,QPRES), sys, cs, ga, u, n, s[2]);
      // q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, cs, ga, u, n, s[3]);
      // q2w(q(ivp1,QRHO), q(ivp1,QUN), q(ivp1,QPRES), sys, cs, ga, u, n, s[4]);
      // q2w(q(ivp2,QRHO), q(ivp2,QUN), q(ivp2,QPRES), sys, cs, ga, u, n, s[5]);
    
      // Reconstruct to i-1/2
      face_wenoz5(s, wr[n], wl[n]);
    }

    // Reconstruct other passive scalars
    s[0] = q(ivm3,QUT);
    s[1] = q(ivm2,QUT);
    s[2] = q(ivm1,QUT);
    s[3] = q(iv,QUT);
    s[4] = q(ivp1,QUT);
    s[5] = q(ivp2,QUT);
    face_wenoz5(s, qr(iv,QUT), ql(iv,QUT));

    s[0] = q(ivm3,QUTT);
    s[1] = q(ivm2,QUTT);
    s[2] = q(ivm1,QUTT);
    s[3] = q(iv,QUTT);
    s[4] = q(ivp1,QUTT);
    s[5] = q(ivp2,QUTT);
    face_wenoz5(s, qr(iv,QUTT), ql(iv,QUTT));

    for (int n = 0; n < NUM_SPECIES; ++n) {      
      s[0] = q(ivm3,QFS+n);
      s[1] = q(ivm2,QFS+n);
      s[2] = q(ivm1,QFS+n);
      s[3] = q(iv,QFS+n);
      s[4] = q(ivp1,QFS+n);
      s[5] = q(ivp2,QFS+n);
      face_wenoz5(s, qr(iv,QFS+n), ql(iv,QFS+n));
    }
  } 
  else {
    Abort("Unknown reconstruction type");
  } 

  // Convert back wl, wr -> ql, qr
  if (recon_scheme != 1) {
    w2q(wl, sys, cs, ga, u, iv, QUN, ql);
    w2q(wr, sys, cs, ga, u, iv, QUN, qr);

    // Enforce sum(Y) = 1
    Real sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      ql(iv,QFS+n) = max<Real>(0.0, ql(iv,QFS+n));
      sumY += ql(iv,QFS+n);
    }
    if (std::abs(sumY - 1.0) > parm.smally) {
      Real invsumY = 1.0 / sumY;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        ql(iv+iv_dir,QFS+n) *= invsumY;
      }
    }
    sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      qr(iv,QFS+n) = max<Real>(0.0, qr(iv,QFS+n));
      sumY += qr(iv,QFS+n);
    }
    if (std::abs(sumY - 1.0) > parm.smally) {
      Real invsumY = 1.0 / sumY;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        qr(iv,QFS+n) *= invsumY;
      }
    }

    // Prositivity preserving limiter
    // if (q(iv,QRHO) > 0 && q(iv-iv_dir,QRHO) > 0) {
    if (iv != ivm1) { 
      // pp_limiter_scalar(qr(iv,QRHO), ql(iv,QRHO), 0.5*(q(iv,QRHO)+q(ivm1,QRHO)), parm.smallr, 1e10);
      // pp_limiter_scalar(qr(iv,QPRES), ql(iv,QPRES), 0.5*(q(iv,QPRES)+q(ivm1,QRHO)), parm.smallp, 1e20);
      // pp_limiter_pres(iv, 
      //                 ql(iv,QRHO), ql(iv,QUN), ql(iv,QPRES), 
      //                 qr(iv,QRHO), qr(iv,QUN), qr(iv,QPRES), 
      //                 0.5*(q(iv,QRHO)+q(ivm1,QRHO)), 0.5*(q(iv,QUN)+q(ivm1,QUN)), 0.5*(q(iv,QPRES)+q(ivm1,QRHO)), 
      //                 parm.smallr, parm.smallp);
      reduce_order_limiter(qr(iv,QRHO), ql(iv,QRHO), q(iv,QRHO), q(ivm1,QRHO), parm.smallr);
      reduce_order_limiter(qr(iv,QPRES), ql(iv,QPRES), q(iv,QPRES), q(ivm1,QPRES), parm.smallr);
    
      if (isnan(qr(iv,QRHO))) amrex::Print() << dir << " qr(" << iv << ",QRHO) is nan. \n";
      else if (qr(iv,QRHO) <= -parm.smallr) amrex::Print() << dir << " qr(" << iv << ",QRHO) <= 0. \n";
      if (isnan(ql(iv,QRHO))) amrex::Print() << dir << " ql(" << iv << ",QRHO) is nan. \n";
      else if (ql(iv,QRHO) <= -parm.smallr) amrex::Print() << dir << " ql(" << iv << ",QRHO) <= 0. \n";
      if (isnan(qr(iv,QPRES))) amrex::Print() << dir << " qr(" << iv << ",QPRES) is nan. \n";
      else if (qr(iv,QPRES) <= -parm.smallp) amrex::Print() << dir << " qr(" << iv << ",QPRES) <= 0. \n";
      if (isnan(ql(iv,QPRES))) amrex::Print() << dir << " ql(" << iv << ",QPRES) is nan. \n";
      else if (ql(iv,QPRES) <= -parm.smallp) amrex::Print() << dir << " ql(" << iv << ",QPRES) <= 0. \n";
    }
  }
}

// /** 
//  * \brief Reconstruct states at cell interface. To be used with pure_riemann().
//  * \param i   (input) x-index
//  * \param j   (input) y-index
//  * \param k   (input) z-index
//  * \param dir (input) direction along which reconstruction takes place (0: x, 1: y, 2: z)
//  * \param q   (input) primitive variables array
//  * \param ql  (output) reconstructive primitive variables on the left of the interface
//  * \param qr  (output) reconstructive primitive variables on the right of the interface
//  * \param recon_scheme (input) chosen reconstruction scheme
//  * \param plm_theta (input) parameter in MUSCL limiter
//  * \param sys (input) indicating what thermodynamics system to use (0: speed of sound, 1: gamma)
//  * \param do_pp_limiter (input) whether to use positivity limiter
//  * \param parm (input) Parm object
//  */
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// char_pp_recon (int i, int j, int k, int dir,
//               Array4<const Real> const& q,
//               Array4<Real> const& ql, Array4<Real> const& qr,
//               const int recon_scheme, 
//               const Real plm_theta,
//               const int sys,
//               Parm const& parm) noexcept
// {
//   const IntVect iv(AMREX_D_DECL(i, j, k));
//   const IntVect iv_dir = IntVect::TheDimensionVector(dir);

//   int QUN, QUT, QUTT; // indices of normal and transverse velocity
//   if (dir == 0) {
//     QUN = QU;  QUT = QV;  QUTT = QW;
//   } 
//   else if (dir == 1) {
//     QUN = QV;  QUT = QU;  QUTT = QW;
//   } 
//   else {
//     QUN = QW;  QUT = QU;  QUTT = QV;
//   }

//   Real thvl, thvr; // averaged thermodynamic variable
//   if (sys == 0) {
//     thvl = 0.5*(q(iv,QC) + q(iv-iv_dir,QC));
//     thvr = 0.5*(q(iv,QC) + q(iv+iv_dir,QC));
//   }
//   else if (sys == 1) {
//     thvl = 0.5*(q(iv,QG) + q(iv-iv_dir,QG));
//     thvr = 0.5*(q(iv,QC) + q(iv+iv_dir,QC));
//   }
//   else{
//     Abort("Unrecognised system: choose between 0: speed of sound and 1: gamma");
//   }

//   // Reconstruction q -> wl, wr
//   Real wl[3], wr[3]; // arrays to hold w+ and w-

//   if (recon_scheme == 1) { // Godunov (no need characteristics)
//     ql(iv+iv_dir,QRHO) = qr(iv,QRHO) = q(iv,QRHO);
//     AMREX_D_TERM(ql(iv+iv_dir,QUN)  = qr(iv,QUN)  = q(iv,QUN); ,
//                  ql(iv+iv_dir,QUT)  = qr(iv,QUT)  = q(iv,QUT); ,
//                  ql(iv+iv_dir,QUTT) = qr(iv,QUTT) = q(iv,QUTT););
//     ql(iv+iv_dir,QPRES) = qr(iv,QPRES) = q(iv,QPRES); 
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       ql(iv+iv_dir,QFS+n) = qr(iv,QFS+n) = q(iv,QFS+n);
//     }
//   } 
//   else if (recon_scheme == 5) { // WENO-Z5
//     Real s[5]; //stencil    

//     // q^+_{i-1/2}
//     for (int n = 0; n < 3; ++n) {
//       // Convert q -> w
//       q2w(q(iv-2*iv_dir,QRHO), q(iv-2*iv_dir,QUN), q(iv-2*iv_dir,QPRES), sys, thvl, n, s[0]);
//       q2w(q(iv-iv_dir,QRHO), q(iv-iv_dir,QUN), q(iv-iv_dir,QPRES), sys, thvl, n, s[1]);
//       q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, thvl, n, s[2]);
//       q2w(q(iv+iv_dir,QRHO), q(iv+iv_dir,QUN), q(iv+iv_dir,QPRES), sys, thvl, n, s[3]);
//       q2w(q(iv+2*iv_dir,QRHO), q(iv+2*iv_dir,QUN), q(iv+2*iv_dir,QPRES), sys, thvl, n, s[4]);

//       // Reconstruct
//       wenoz5_m(s, wl[n]);
//     }

//     // q^-_{i+1/2}
//     for (int n = 0; n < 3; ++n) {
//       // Convert q -> w
//       q2w(q(iv-2*iv_dir,QRHO), q(iv-2*iv_dir,QUN), q(iv-2*iv_dir,QPRES), sys, thvr, n, s[4]);
//       q2w(q(iv-iv_dir,QRHO), q(iv-iv_dir,QUN), q(iv-iv_dir,QPRES), sys, thvr, n, s[3]);
//       q2w(q(iv,QRHO), q(iv,QUN), q(iv,QPRES), sys, thvr, n, s[2]);
//       q2w(q(iv+iv_dir,QRHO), q(iv+iv_dir,QUN), q(iv+iv_dir,QPRES), sys, thvr, n, s[1]);
//       q2w(q(iv+2*iv_dir,QRHO), q(iv+2*iv_dir,QUN), q(iv+2*iv_dir,QPRES), sys, thvr, n, s[0]);

//       // Reconstruct
//       wenoz5_m(s, wr[n]);
//     }

//     // Reconstruct other passive scalars
//     s[0] = q(iv-2*iv_dir,QUT);
//     s[1] = q(iv-iv_dir,QUT);
//     s[2] = q(iv,QUT);
//     s[3] = q(iv+iv_dir,QUT);
//     s[4] = q(iv+2*iv_dir,QUT);
//     wenoz5(s, qr(iv,QUT), ql(iv+iv_dir,QUT));

//     s[0] = q(iv-2*iv_dir,QUTT);
//     s[1] = q(iv-iv_dir,QUTT);
//     s[2] = q(iv,QUTT);
//     s[3] = q(iv+iv_dir,QUTT);
//     s[4] = q(iv+2*iv_dir,QUTT);
//     wenoz5(s, qr(iv,QUTT), ql(iv+iv_dir,QUTT));

//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       s[0] = q(iv-2*iv_dir,QFS+n);
//       s[1] = q(iv-iv_dir,QFS+n);
//       s[2] = q(iv,QFS+n);
//       s[3] = q(iv+iv_dir,QFS+n);
//       s[4] = q(iv+2*iv_dir,QFS+n);
//       wenoz5(s, qr(iv,QFS+n), ql(iv+iv_dir,QFS+n));
//     }
//   } 
//   else {
//     Abort("Unknown reconstruction type");
//   } 

//   // Convert back wl, wr -> ql, qr
//   if (recon_scheme != 1) {
//     w2q(wl, sys, thvl, iv, QUN, qr);
//     w2q(wr, sys, thvr, iv+iv_dir, QUN, ql);

//     // Enforce sum(Y) = 1
//     Real sumY = 0.0;
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       ql(iv+iv_dir,QFS+n) = max<Real>(0.0, ql(iv+iv_dir,QFS+n));
//       sumY += ql(iv+iv_dir,QFS+n);
//     }
//     if (std::abs(sumY - 1.0) > parm.smally) {
//       Real invsumY = 1.0 / sumY;
//       for (int n = 0; n < NUM_SPECIES; ++n) {
//         ql(iv+iv_dir,QFS+n) *= invsumY;
//       }
//     }
//     sumY = 0.0;
//     for (int n = 0; n < NUM_SPECIES; ++n) {
//       qr(iv,QFS+n) = max<Real>(0.0, qr(iv,QFS+n));
//       sumY += qr(iv,QFS+n);
//     }
//     if (std::abs(sumY - 1.0) > parm.smally) {
//       Real invsumY = 1.0 / sumY;
//       for (int n = 0; n < NUM_SPECIES; ++n) {
//         qr(iv,QFS+n) *= invsumY;
//       }
//     }

//     // Prositivity preserving limiter
//     if (q(iv,QRHO) > 0) {
//       // pp_limiter_euler(iv+iv_dir, iv, ql, qr, q, parm.smallr, parm.smallp);
//       pp_limiter_scalar(qr(iv,QRHO), ql(iv+iv_dir,QRHO), q(iv,QRHO), parm.smallr, 1e10);
//       pp_limiter_scalar(qr(iv,QPRES), ql(iv+iv_dir,QPRES), q(iv,QPRES), parm.smallp, 1e10);
//     }
//   }
// }

#endif