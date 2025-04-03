#ifndef CNS_DIFFUSION_EB_K_H_
#define CNS_DIFFUSION_EB_K_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

#include <cmath>

#include "index_macros.H"
#include "wall_model.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_species_eb(amrex::IntVect const& iv, const int idir,
                    amrex::Array4<amrex::Real const> const& q,
                    amrex::Array4<amrex::Real const> const& coeffs,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
                    amrex::Real flx[NVAR])
{
#if (NUM_SPECIES > 1) // no point doing diffusion for single species
  using amrex::Real;

  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(idir));

  // Get massfrac, molefrac, enthalpy
  Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
  Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    mass1[n] = q(iv, QFS + n);
    mass2[n] = q(ivm, QFS + n);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);

  // Compute species and enthalpy fluxes for ideal EOS
  // Get species/enthalpy diffusion, compute correction vel
  Real T = q(iv, QTEMP);
  eos.RTY2Hi(q(iv, QRHO), T, mass1, hi1);
  T = q(ivm, QTEMP);
  eos.RTY2Hi(q(ivm, QRHO), T, mass2, hi2);

  Real Vc = 0.0;
  const amrex::Real dpdx = (q(iv, QPRES) - q(ivm, QPRES)) * dxinv[0];
  const amrex::Real dlnp = dpdx / (0.5 * (q(iv, QPRES) + q(ivm, QPRES)));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real Xface = 0.5 * (mole1[n] + mole2[n]);
    const amrex::Real Yface = 0.5 * (mass1[n] + mass2[n]);
    const amrex::Real hface = 0.5 * (hi1[n] + hi2[n]);
    const amrex::Real dXdx = (mole1[n] - mole2[n]) * dxinv[0];
    const amrex::Real Vd = -0.5 * (coeffs(iv, CRHOD + n) + coeffs(ivm, CRHOD + n)) *
                           (dXdx + (Xface - Yface) * dlnp);
    Vc += Vd;
    flx[UFS + n] += Vd;
    flx[UEDEN] += Vd * hface;
  }
  // Add correction velocity to fluxes so sum(Vd) = 0
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    const amrex::Real Yface = 0.5 * (mass1[ns] + mass2[ns]);
    const amrex::Real hface = 0.5 * (hi1[ns] + hi2[ns]);
    flx[UFS + ns] -= Yface * Vc;
    flx[UEDEN] -= Yface * hface * Vc;
  }
#endif
}

/**
 * @brief EB version of cns_diff.
 *
 * @param i,j,k  x, y, z index.
 * @param dir    direction, 0:x, 1:y, 2:z.
 * @param q      primitive variables.
 * @param coeffs transport coefficients.
 * @param flag   EB status flags.
 * @param dxinv  1/dx.
 * @param[out] flx_arr  output diffusion fluxes.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_eb(const amrex::IntVect& iv, const int d1,
            amrex::Array4<const amrex::Real> const& q,
            amrex::Array4<const amrex::Real> const& coeffs,
            amrex::Array4<const amrex::EBCellFlag> const& flag,
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
            amrex::Real flx[NVAR]) noexcept
{
  using amrex::Real;
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(d1));
  if (flag(iv).isCovered() || flag(ivm).isCovered()) { return; }

  constexpr Real weights[3] = {0.0, 1.0, 0.5}; // EB weights

  const int d2 = d1 == 0 ? 1 : 0;
  const int d3 = d1 == 2 ? 1 : 2;
  AMREX_D_TERM(const int QU1 = QU + d1;, const int QU2 = QU + d2;
               , const int QU3 = QU + d3;)
  AMREX_D_TERM(const int UM1 = UMX + d1;, const int UM2 = UMX + d2;
               , const int UM3 = UMX + d3;)

  // Aij = dA_i/dx_j
  const Real dTdn = (q(iv, QTEMP) - q(ivm, QTEMP)) * dxinv[d1];
  const Real u11 = (q(iv, QU1) - q(ivm, QU1)) * dxinv[d1];
  // EB only allows for 2D and 3D
  const Real u21 = (q(iv, QU2) - q(ivm, QU2)) * dxinv[d1];
  const amrex::IntVect ivd2(amrex::IntVect::TheDimensionVector(d2));
  int hip = +(flag(iv).isConnected(ivd2) ? 1 : 0);
  int him = -(flag(iv).isConnected(-ivd2) ? 1 : 0);
  int lop = +(flag(ivm).isConnected(ivd2) ? 1 : 0);
  int lom = -(flag(ivm).isConnected(-ivd2) ? 1 : 0);
  Real whi = weights[hip - him];
  Real wlo = weights[lop - lom];
  const Real u12 = (0.5 * dxinv[d2]) *
                   ((q(iv + hip * ivd2, QU1) - q(iv + him * ivd2, QU1)) * whi +
                    (q(ivm + lop * ivd2, QU1) - q(ivm + lom * ivd2, QU1)) * wlo);
  const Real u22 = (0.5 * dxinv[d2]) *
                   ((q(iv + hip * ivd2, QU2) - q(iv + him * ivd2, QU2)) * whi +
                    (q(ivm + lop * ivd2, QU2) - q(ivm + lom * ivd2, QU2)) * wlo);
  // Real whip = hip == 1 ? 1 : -1;
  // Real whim = him == -1 ? 1 : -1;
  // Real wlop = lop == 1 ? 1 : -1;
  // Real wlom = lom == -1 ? 1 : -1;
  // const Real u12 = (0.25 * dxinv[d2]) *
  //                  (q(iv + hip * ivd2, QU1) * whip - q(iv + him * ivd2, QU1) * whim +
  //                   q(ivm + lop * ivd2, QU1) * wlop - q(ivm + lom * ivd2, QU1) * wlom);
  // const Real u22 = (0.25 * dxinv[d2]) *
  //                  (q(iv + hip * ivd2, QU2) * whip - q(iv + him * ivd2, QU2) * whim +
  //                   q(ivm + lop * ivd2, QU2) * wlop - q(ivm + lom * ivd2, QU2) * wlom);
#if (AMREX_SPACEDIM == 3)
  const Real u31 = (q(iv, QU3) - q(ivm, QU3)) * dxinv[d1];
  const amrex::IntVect ivd3(amrex::IntVect::TheDimensionVector(d3));
  hip = +(flag(iv).isConnected(ivd3) ? 1 : 0);
  him = -(flag(iv).isConnected(-ivd3) ? 1 : 0);
  lop = +(flag(ivm).isConnected(ivd3) ? 1 : 0);
  lom = -(flag(ivm).isConnected(-ivd3) ? 1 : 0);
  whi = weights[hip - him];
  wlo = weights[lop - lom];
  const Real u13 = (0.5 * dxinv[d3]) *
                   ((q(iv + hip * ivd3, QU1) - q(iv + him * ivd3, QU1)) * whi +
                    (q(ivm + lop * ivd3, QU1) - q(ivm + lom * ivd3, QU1)) * wlo);
  const Real u33 = (0.5 * dxinv[d3]) *
                   ((q(iv + hip * ivd3, QU3) - q(iv + him * ivd3, QU3)) * whi +
                    (q(ivm + lop * ivd3, QU3) - q(ivm + lom * ivd3, QU3)) * wlo);
  // whip = hip == 1 ? 1 : -1;
  // whim = him == -1 ? 1 : -1;
  // wlop = lop == 1 ? 1 : -1;
  // wlom = lom == -1 ? 1 : -1;
  // const Real u13 = (0.25 * dxinv[d3]) *
  //                  (q(iv + hip * ivd3, QU1) * whip - q(iv + him * ivd3, QU1) * whim +
  //                   q(ivm + lop * ivd3, QU1) * wlop - q(ivm + lom * ivd3, QU1) * wlom);
  // const Real u33 = (0.25 * dxinv[d3]) *
  //                  (q(iv + hip * ivd3, QU3) * whip - q(iv + him * ivd3, QU3) * whim +
  //                   q(ivm + lop * ivd3, QU3) * wlop - q(ivm + lom * ivd3, QU3) * wlom);

#endif
  const Real divu = AMREX_D_TERM(u11, +u22, +u33);
  const Real muf = 0.5 * (coeffs(iv, CMU) + coeffs(ivm, CMU));
  const Real xif = 0.5 * (coeffs(iv, CXI) + coeffs(ivm, CXI));
  AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - (2.0 / 3.0) * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

  AMREX_D_TERM(flx[UM1] -= tau11;, flx[UM2] -= tau12;, flx[UM3] -= tau13;)
  flx[UEDEN] -= 0.5 * (AMREX_D_TERM((q(iv, QU1) + q(ivm, QU1)) * tau11,
                                        +(q(iv, QU2) + q(ivm, QU2)) * tau12,
                                        +(q(iv, QU3) + q(ivm, QU3)) * tau13) +
                           (coeffs(iv, CLAM) + coeffs(ivm, CLAM)) * dTdn);

  // Species transport
  cns_diff_species_eb(iv, d1, q, coeffs, dxinv, flx);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_interp1d(amrex::Real cym, amrex::Real cy0, amrex::Real cyp,
                 amrex::Array1D<amrex::Real, 1, 3>& v)
{
  Real x = cym * v(1) + cy0 * v(2) + cyp * v(3);
  return x;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real compute_interp2d(
  amrex::Real cym, amrex::Real cy0, amrex::Real cyp, amrex::Real czm,
  amrex::Real cz0, amrex::Real czp, amrex::Array2D<amrex::Real, 1, 3, 1, 3>& v)
{
  Real x = czm * (cym * v(1, 1) + cy0 * v(2, 1) + cyp * v(3, 1)) +
           cz0 * (cym * v(1, 2) + cy0 * v(2, 2) + cyp * v(3, 2)) +
           czp * (cym * v(1, 3) + cy0 * v(2, 3) + cyp * v(3, 3));
  return x;
}

#if (AMREX_SPACEDIM == 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_diff_wallflux(
  int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
  amrex::Array4<amrex::Real const> const& coefs,
  amrex::Array4<amrex::Real const> const& bcent, const amrex::Real axm,
  const amrex::Real axp, const amrex::Real aym, const amrex::Real ayp,
  const amrex::GpuArray<amrex::Real, 2> dxinv,
  amrex::GpuArray<amrex::Real, NVAR>& viscw, bool eb_no_slip, bool eb_isothermal,
  amrex::Real eb_wall_temp, const bool do_wall_model) noexcept
{
  //  This implementation assumes adiabatic or isothermal walls

  using amrex::Real;

  for (int n = 0; n < NVAR; n++) viscw[n] = 0.;

  Real dapx = axp - axm;
  Real dapy = ayp - aym;

  Real apnorm = std::sqrt(dapx * dapx + dapy * dapy); // this is cut-face area / dx
  if (apnorm <= 0.0) amrex::Abort("compute_diff_wallflux: we are in trouble.");

  Real apnorminv = 1.0 / apnorm;
  Real anrmx = -dapx * apnorminv; // unit vector pointing toward the wall
  Real anrmy = -dapy * apnorminv;

  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> bct;

  // The center of the wall
  bct[0] = bcent(i, j, k, 0);
  bct[1] = bcent(i, j, k, 1);

  Real u1, v1, T1, u2, v2, T2, d1, d2;

  if (std::abs(anrmx) >= std::abs(anrmy)) {
    // y line: x = const
    // the equation for the line:  x = bct(0) - d*anrmx
    //                             y = bct(1) - d*anrmy
    Real s = -anrmx > 0 ? 1. : -1.;
    int is = static_cast<int>(round(s));

    //
    // the line intersects the y line (x = s) at ...
    //
    d1 = (bct[0] - s) *
         (1.0 / anrmx); // this is also the distance from wall to intersection
    Real yit = bct[1] - d1 * anrmy;
    int iyit = j + static_cast<int>(round(yit));
    yit = yit - static_cast<int>(round(
                  yit)); // shift so that the center of the nine cells are (0.,0.)

    // coefficients for quadratic interpolation

    Real cym = 0.50 * yit * (yit - 1.0);
    Real cy0 = 1.0 - yit * yit;
    Real cyp = 0.50 * yit * (yit + 1.0);

    // interpolation

    Array1D<Real, 1, 3> v;

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(i + is, iyit - 2 + ii, k, QU); }
    u1 = compute_interp1d(cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(i + is, iyit - 2 + ii, k, QV); }
    v1 = compute_interp1d(cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(i + is, iyit - 2 + ii, k, QTEMP); }
    T1 = compute_interp1d(cym, cy0, cyp, v);

    //
    // the line intersects the y-z plane (x = 2*s) at ...
    //

    d2 = (bct[0] - 2.0 * s) * (1.0 / anrmx);
    yit = bct[1] - d2 * anrmy;
    iyit = j + static_cast<int>(round(yit));
    yit = yit - static_cast<int>(round(
                  yit)); // shift so that the center of the nine cells are (0.,0.)

    // coefficients for quadratic interpolation

    cym = 0.5 * yit * (yit - 1.);
    cy0 = 1. - yit * yit;
    cyp = 0.5 * yit * (yit + 1.);

    // interpolation

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(i + 2 * is, iyit - 2 + ii, k, QU); }
    u2 = compute_interp1d(cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(i + 2 * is, iyit - 2 + ii, k, QV); }
    v2 = compute_interp1d(cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      v(ii) = q(i + 2 * is, iyit - 2 + ii, k, QTEMP);
    }
    T2 = compute_interp1d(cym, cy0, cyp, v);

    // Tangent pressure gradient for wall model

    // Real p1 = q(i + 2 * is, iyit - 1, QPRES);
    // Real p2 = q(i + 2 * is, iyit + 1, QPRES);
    // dpdy = (p2 + p1) * yit + 0.5 * (p2 - p1);

  } else {
    // x line
    Real s = -anrmy > 0. ? 1. : -1.;
    int is = static_cast<int>(round(s));

    d1 = (bct[1] - s) / anrmy;
    Real xit = bct[0] - d1 * anrmx;

    int ixit = i + static_cast<int>(round(xit));
    xit = xit - static_cast<int>(round(xit));

    Real cxm = 0.5 * xit * (xit - 1.);
    Real cx0 = 1. - xit * xit;
    Real cxp = 0.5 * xit * (xit + 1.);

    Array1D<Real, 1, 3> v;

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(ixit - 2 + ii, j + is, k, QU); }
    u1 = compute_interp1d(cxm, cx0, cxp, v);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(ixit - 2 + ii, j + is, k, QV); }
    v1 = compute_interp1d(cxm, cx0, cxp, v);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(ixit - 2 + ii, j + is, k, QTEMP); }
    T1 = compute_interp1d(cxm, cx0, cxp, v);

    d2 = (bct[1] - 2.0 * s) * (1.0 / anrmy);
    xit = bct[0] - d2 * anrmx;
    ixit = i + static_cast<int>(round(xit));
    xit = xit - static_cast<int>(round(xit));

    cxm = 0.5 * xit * (xit - 1.);
    cx0 = 1. - xit * xit;
    cxp = 0.5 * xit * (xit + 1.);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(ixit - 2 + ii, j + 2 * is, k, QU); }
    u2 = compute_interp1d(cxm, cx0, cxp, v);

    for (int ii = 1; ii < 4; ii++) { v(ii) = q(ixit - 2 + ii, j + 2 * is, k, QV); }
    v2 = compute_interp1d(cxm, cx0, cxp, v);

    for (int ii = 1; ii < 4; ii++) {
      v(ii) = q(ixit - 2 + ii, j + 2 * is, k, QTEMP);
    }
    T2 = compute_interp1d(cxm, cx0, cxp, v);

    // Tangent pressure gradient for wall model

    // Real p1 = q(ixit - 1, j + 2 * is, k, QPRES);
    // Real p2 = q(ixit + 1, j + 2 * is, k, QPRES);
    // dpdx = ((p2 + p1) * xit + 0.5 * (p2 - p1));
  }
  // Real dpdt = - dpdx * anrmy - dpdy * anrmx;

  // get density and transport coefficients 
  Real rho, mu, xi, lam;
  Real p = q(i, j, k, QPRES), Y[NUM_SPECIES]; // assume dp/dn, dY/dn ~ 0
  for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = q(i, j, k, QFS + n); }
  if (eb_isothermal) {
    auto trans = pele::physics::PhysicsType::transport();
    auto const* ltransparm = CNS::trans_parms.device_trans_parm();
    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(p, Y, eb_wall_temp, rho);

    const bool get_xi = true, get_mu = true, get_lam = true, get_Ddiag = false,
               get_chi = false;
    trans.transport(get_xi, get_mu, get_lam, get_Ddiag, get_chi, eb_wall_temp, rho,
                    Y, nullptr, nullptr, mu, xi, lam, ltransparm);
  } else {
    // Adiabatic wall dT/dn ~ 0
    rho = q(i, j, k, QRHO);
    mu = coefs(i, j, k, CMU);
    xi = coefs(i, j, k, CXI);
    lam = coefs(i, j, k, CLAM);
  }

  Real ddinv = 1. / (d1 * d2 * (d2 - d1));

  if (eb_no_slip) {
    // compute derivatives on the wall given that velocity is zero on the wall
    // note that the normal vector points toward the wall
    Real dudn = -ddinv * (d2 * d2 * u1 - d1 * d1 * u2);
    Real dvdn = -ddinv * (d2 * d2 * v1 - d1 * d1 * v2);

    // transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
    Real dudx = dudn * anrmx;
    Real dudy = dudn * anrmy;
    Real dvdx = dvdn * anrmx;
    Real dvdy = dvdn * anrmy;
    Real divu = dudx + dvdy;

    Real tautmp = (xi - (2.0 / 3.0) * mu) * divu;
    Real tauxx = mu * 2. * dudx + tautmp;
    Real tauyy = mu * 2. * dvdy + tautmp;
    Real tauxy = mu * (dudy + dvdx);

    // assumes dx == dy == dz
    viscw[UMX] = (dapx * tauxx + dapy * tauxy);
    viscw[UMY] = (dapx * tauxy + dapy * tauyy);
  }

  if (eb_isothermal) {
    // normal (into cell) temp gradient
    Real dTdn =
      ddinv * (d2 * d2 * T1 - d1 * d1 * T2 + (d1 * d1 - d2 * d2) * eb_wall_temp);

    viscw[UEDEN] = apnorm * lam * dTdn;
  }

  // Modify if using a wall model
  if (do_wall_model) {
    // Check y+
    const Real dx = 1.0 / dxinv[0]; // dx == dy
    const Real cfarea = apnorm * dx; // cut-face area / dx
    Real tau = std::sqrt(viscw[UMX] * viscw[UMX] + viscw[UMY] * viscw[UMY]) / cfarea;
    Real utau = std::sqrt(tau / rho);
    Real yplus = rho * utau / mu * dx;

    if (yplus > 11.0) {
      // Wall model
      // LawOfTheWall wm;
      EquilibriumODE wm;
      Real u_parallel = u2 * anrmy - v2 * anrmx; // u dot tangent
      Real ts = (u_parallel > 0.0) ? 1.0 : -1.0;
      u_parallel *= ts; // abs(u)
      if (u_parallel > 1e-3 && T2 > 90.0 && T2 < 5000.0 && rho > 0.0) {
        Real T_wall = eb_isothermal ? eb_wall_temp : -1.0;
        Real h = d2 * dx;
        Real tauw, qw;
        wm.parallel_wall_stress(u_parallel, T2, rho, Y, h, mu, lam, T_wall, tauw, qw);
        // Real dudy = -ddinv * (d2 * d2 * u1 - d1 * d1 * u2);
        // std::cout << u_parallel << " " << h << " " << mu*dudy/dx << " " << tau << "  ";
        // std::cout << "y+: " << yplus << " -> " << rho * std::sqrt(tauw / rho) / mu * dx;
        // std::cout << "  viscous wall stress: " << tau * cfarea << " (" << viscw[UMX] << ", " << viscw[UMY]
        //           << ") -> ";

        viscw[UMX] = cfarea * ts * anrmy * tauw;
        viscw[UMY] = -cfarea * ts * anrmx * tauw;
        // std::cout << tauw * cfarea << " (" << viscw[UMX] << ", " << viscw[UMY] << ")";
        // viscw[UEDEN] = cfarea * qw; // TODO: no wall model heat flux for now
        // std::cout << "  " << apnorm << " (" << anrmx << ", " << anrmy << ")\n";
      }
    } 
    // else {
    //   std::cout << "y+: " << yplus << std::endl;
    // }
  }
}

#else // 3d version below

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_diff_wallflux(
  int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
  amrex::Array4<amrex::Real const> const& coefs,
  amrex::Array4<amrex::Real const> const& bcent, const amrex::Real axm,
  const amrex::Real axp, const amrex::Real aym, const amrex::Real ayp,
  const amrex::Real azm, const amrex::Real azp,
  const amrex::GpuArray<amrex::Real, 3> dxinv,
  amrex::GpuArray<amrex::Real, NVAR>& viscw, bool eb_no_slip, bool eb_isothermal,
  amrex::Real eb_wall_temp, const bool do_wall_model) noexcept
{
  //  This implementation assumes adiabatic or isothermal walls

  using amrex::Real;

  AMREX_D_TERM(Real dapx = axp - axm;, Real dapy = ayp - aym;
               , Real dapz = azp - azm;);

  Real apnorm = std::sqrt(dapx * dapx + dapy * dapy + dapz * dapz);
  if (apnorm == 0.0) amrex::Abort("compute_diff_wallflux: we are in trouble.");

  Real apnorminv = 1.0 / apnorm;
  AMREX_D_TERM(Real anrmx = -dapx * apnorminv;
               , // unit vector pointing toward the wall
               Real anrmy = -dapy * apnorminv;
               , Real anrmz = -dapz * apnorminv;);

  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> bct;

  // The center of the wall
  AMREX_D_TERM(bct[0] = bcent(i, j, k, 0);, bct[1] = bcent(i, j, k, 1);
               , bct[2] = bcent(i, j, k, 2););

  Real u1, v1, w1, T1, u2, v2, w2, T2, d1, d2;

  if ((std::abs(anrmx) >= std::abs(anrmy)) && (std::abs(anrmx) >= std::abs(anrmz))) {
    // y-z plane: x = const
    // the equation for the line:  x = bct(0) - d*anrmx
    //                             y = bct(1) - d*anrmy
    //                             z = bct(2) - d*anrmz
    Real s = -anrmx > 0 ? 1. : -1.;
    int is = static_cast<int>(round(s));

    //
    // the line intersects the y-z plane (x = s) at ...
    //
    d1 = (bct[0] - s) *
         (1.0 / anrmx); // this is also the distance from wall to intersection
    Real yit = bct[1] - d1 * anrmy;
    int iyit = j + static_cast<int>(round(yit));
    yit = yit - static_cast<int>(round(
                  yit)); // shift so that the center of the nine cells are (0.,0.)

    Real zit = bct[2] - d1 * anrmz;
    int izit = k + static_cast<int>(round(zit));
    zit = zit - static_cast<int>(round(zit));

    // coefficients for quadratic interpolation

    Real cym = 0.50 * yit * (yit - 1.0);
    Real cy0 = 1.0 - yit * yit;
    Real cyp = 0.50 * yit * (yit + 1.0);

    Real czm = 0.50 * zit * (zit - 1.0);
    Real cz0 = 1.0 - zit * zit;
    Real czp = 0.50 * zit * (zit + 1.0);

    // interpolation

    Array2D<Real, 1, 3, 1, 3> v;

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + is, iyit - 2 + ii, izit - 2 + jj, QU);
    }
    u1 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + is, iyit - 2 + ii, izit - 2 + jj, QV);
    }
    v1 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + is, iyit - 2 + ii, izit - 2 + jj, QW);
    }
    w1 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + is, iyit - 2 + ii, izit - 2 + jj, QTEMP);
    }
    T1 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    //
    // the line intersects the y-z plane (x = 2*s) at ...
    //

    d2 = (bct[0] - 2.0 * s) * (1.0 / anrmx);
    yit = bct[1] - d2 * anrmy;
    zit = bct[2] - d2 * anrmz;
    iyit = j + static_cast<int>(round(yit));
    izit = k + static_cast<int>(round(zit));
    yit = yit - static_cast<int>(round(
                  yit)); // shift so that the center of the nine cells are (0.,0.)
    zit = zit - static_cast<int>(round(zit));

    // coefficients for quadratic interpolation

    cym = 0.5 * yit * (yit - 1.);
    cy0 = 1. - yit * yit;
    cyp = 0.5 * yit * (yit + 1.);
    czm = 0.5 * zit * (zit - 1.);
    cz0 = 1. - zit * zit;
    czp = 0.5 * zit * (zit + 1.);

    // interpolation

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + 2 * is, iyit - 2 + ii, izit - 2 + jj, QU);
    }
    u2 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + 2 * is, iyit - 2 + ii, izit - 2 + jj, QV);
    }
    v2 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + 2 * is, iyit - 2 + ii, izit - 2 + jj, QW);
    }
    w2 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(i + 2 * is, iyit - 2 + ii, izit - 2 + jj, QTEMP);
    }
    T2 = compute_interp2d(cym, cy0, cyp, czm, cz0, czp, v);

  } else if (std::abs(anrmy) > std::abs(anrmx) &&
             std::abs(anrmy) > std::abs(anrmz)) {
    // z-x plane
    Real s = -anrmy > 0. ? 1. : -1.;
    int is = static_cast<int>(round(s));

    AMREX_D_TERM(d1 = (bct[1] - s) / anrmy;, Real xit = bct[0] - d1 * anrmx;
                 , Real zit = bct[2] - d1 * anrmz;);
    int ixit = i + static_cast<int>(round(xit));
    xit = xit - static_cast<int>(round(xit));

    Real cxm = 0.5 * xit * (xit - 1.);
    Real cx0 = 1. - xit * xit;
    Real cxp = 0.5 * xit * (xit + 1.);

    int izit = k + static_cast<int>(round(zit));
    zit = zit - static_cast<int>(round(zit));
    Real czm = 0.5 * zit * (zit - 1.);
    Real cz0 = 1. - zit * zit;
    Real czp = 0.5 * zit * (zit + 1.);

    Array2D<Real, 1, 3, 1, 3> v;

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + is, izit - 2 + jj, QU);
    }
    u1 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + is, izit - 2 + jj, QV);
    }
    v1 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + is, izit - 2 + jj, QW);
    }
    w1 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + is, izit - 2 + jj, QTEMP);
    }
    T1 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    d2 = (bct[1] - 2.0 * s) * (1.0 / anrmy);
    xit = bct[0] - d2 * anrmx;
    ixit = i + static_cast<int>(round(xit));
    xit = xit - static_cast<int>(round(xit));

    cxm = 0.5 * xit * (xit - 1.);
    cx0 = 1. - xit * xit;
    cxp = 0.5 * xit * (xit + 1.);

    zit = bct[2] - d2 * anrmz;
    izit = k + static_cast<int>(round(zit));
    zit = zit - static_cast<int>(round(zit));
    czm = 0.5 * zit * (zit - 1.);
    cz0 = 1. - zit * zit;
    czp = 0.5 * zit * (zit + 1.);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + 2 * is, izit - 2 + jj, QU);
    }
    u2 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + 2 * is, izit - 2 + jj, QV);
    }
    v2 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + 2 * is, izit - 2 + jj, QW);
    }
    w2 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, j + 2 * is, izit - 2 + jj, QTEMP);
    }
    T2 = compute_interp2d(cxm, cx0, cxp, czm, cz0, czp, v);

  } else {
    // x-y plane

    Real s = -anrmz > 0 ? 1. : -1.;
    int is = static_cast<int>(round(s));

    d1 = (bct[2] - s) * (1.0 / anrmz);
    Real xit = bct[0] - d1 * anrmx;
    Real yit = bct[1] - d1 * anrmy;
    int ixit = i + static_cast<int>(round(xit));
    int iyit = j + static_cast<int>(round(yit));
    xit = xit - static_cast<int>(round(xit));
    yit = yit - static_cast<int>(round(yit));

    Real cxm = 0.5 * xit * (xit - 1.);
    Real cx0 = 1. - xit * xit;
    Real cxp = 0.5 * xit * (xit + 1.);
    Real cym = 0.5 * yit * (yit - 1.);
    Real cy0 = 1. - yit * yit;
    Real cyp = 0.5 * yit * (yit + 1.);

    Array2D<Real, 1, 3, 1, 3> v;

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + is, QU);
    }
    u1 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + is, QV);
    }
    v1 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + is, QW);
    }
    w1 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + is, QTEMP);
    }
    T1 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    d2 = (bct[2] - 2. * s) * (1.0 / anrmz);
    xit = bct[0] - d2 * anrmx;
    yit = bct[1] - d2 * anrmy;
    ixit = i + static_cast<int>(round(xit));
    iyit = j + static_cast<int>(round(yit));
    xit = xit - static_cast<int>(round(xit));
    yit = yit - static_cast<int>(round(yit));

    cxm = 0.5 * xit * (xit - 1.);
    cx0 = 1. - xit * xit;
    cxp = 0.5 * xit * (xit + 1.);
    cym = 0.5 * yit * (yit - 1.);
    cy0 = 1. - yit * yit;
    cyp = 0.5 * yit * (yit + 1.);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + 2 * is, QU);
    }
    u2 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + 2 * is, QV);
    }
    v2 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + 2 * is, QW);
    }
    w2 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);

    for (int ii = 1; ii < 4; ii++) {
      for (int jj = 1; jj < 4; jj++)
        v(ii, jj) = q(ixit - 2 + ii, iyit - 2 + jj, k + 2 * is, QTEMP);
    }
    T2 = compute_interp2d(cxm, cx0, cxp, cym, cy0, cyp, v);
  }

  // compute derivatives on the wall given that velocity is zero on the wall.

  Real ddinv = 1. / (d1 * d2 * (d2 - d1));

  // get density and transport coefficients
  Real rho, mu, xi, lam;
  Real p = q(i, j, k, QPRES), Y[NUM_SPECIES]; // assume dp/dn, dY/dn ~ 0
  for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = q(i, j, k, QFS + n); }
  // if (eb_isothermal) {
  //   auto trans = pele::physics::PhysicsType::transport();
  //   auto const* ltransparm = CNS::trans_parms.device_trans_parm();
  //   Real dummy2;

  //   auto eos = pele::physics::PhysicsType::eos();
  //   eos.PYT2R(p, Y, eb_wall_temp, rho); // assume dp/dn ~ 0

  //   trans.transport(true, true, true, false, false, eb_wall_temp, rho, Y, nullptr,
  //                   nullptr, mu, xi, lam, ltransparm);
  // } else 
  {
    // Adiabatic wall dT/dn ~ 0
    rho = q(i, j, k, QRHO);
    mu = coefs(i, j, k, CMU);
    xi = coefs(i, j, k, CXI);
    lam = coefs(i, j, k, CLAM);
  }

// #ifdef MONITOR_PROBLEM_CELL
  // if (!std::isgreater(rho, 0.0) || !std::isgreater(mu, 0.0) ||
  //     !std::isgreater(lam, 0.0)) {
  //   std::cout << "p: " << p << ", Tw: " << eb_wall_temp << ", Y:";
  //   for (int n = 0; n < NUM_SPECIES; ++n) { std::cout << Y[n] << " "; }
  //   std::cout << " rho: " << rho << ", mu: " << mu << ", lam: " << lam << std::endl;
  // }
  // if (std::isnan(v1) || std::isnan(v2) || std::isnan(u1) || std::isnan(u2) ||
  //     std::isnan(w1) || std::isnan(w2) || 
  //     !std::isgreater(T1, 0.0) || !std::isgreater(T2, 0.0)) {
  //   std::cout << "d1: " << d1 << ", d2: " << d2 << " v1: " << v1 << ", v2: " << v2
  //             << ", u1: " << u1 << ", u2: " << u2 << ", w1: " << w1 << ", w2: " << w2
  //             << ", T1: " << T1 << ", T2: " << T2 << std::endl;
  // }
// #endif

  {
    Real tauxx, tauyy, tauxy, tauzz, tauxz, tauyz; 
    Real dudn, dvdn, dwdn;
  if (eb_no_slip) {
    dudn =
      -ddinv * (d2 * d2 * u1 -
                d1 * d1 * u2); // note that the normal vector points toward the wall
    dvdn = -ddinv * (d2 * d2 * v1 - d1 * d1 * v2);
    dwdn = -ddinv * (d2 * d2 * w1 - d1 * d1 * w2);

    // transform them to d/dx, d/dy and d/dz given transverse derivatives are zero

    Real dudx = dudn * anrmx;
    Real dvdx = dvdn * anrmx;
    Real dwdx = dwdn * anrmx;

    Real dudy = dudn * anrmy;
    Real dvdy = dvdn * anrmy;
    Real dwdy = dwdn * anrmy;

    Real dudz = dudn * anrmz;
    Real dvdz = dvdn * anrmz;
    Real dwdz = dwdn * anrmz;

    Real divu = dudx + dvdy + dwdz;

    Real tautmp = (xi - (2.0 / 3.0) * mu) * divu;
    tauxx = mu * 2. * dudx + tautmp;
    tauyy = mu * 2. * dvdy + tautmp;
    tauxy = mu * (dudy + dvdx);

    tauzz = mu * 2. * dwdz + tautmp;
    tauxz = mu * (dudz + dwdx);
    tauyz = mu * (dwdy + dvdz);

    // assumes dx == dy == dz
    viscw[UMX] = (dapx * tauxx + dapy * tauxy + dapz * tauxz);
    viscw[UMY] = (dapx * tauxy + dapy * tauyy + dapz * tauyz);
    viscw[UMZ] = (dapx * tauxz + dapy * tauyz + dapz * tauzz);
  }

  if (eb_isothermal) {
    Real dTdn = ddinv * (d2 * d2 * T1 - d1 * d1 * T2 +
                         (d1 * d1 - d2 * d2) *
                           eb_wall_temp); // normal (into cell) temp gradient

    viscw[UEDEN] = apnorm * lam * dTdn;
  }

#ifdef MONITOR_PROBLEM_CELL
 if (std::isnan(viscw[UMX]) || std::isnan(viscw[UMY]) || std::isnan(viscw[UMZ]) ||
     std::isnan(viscw[UEDEN])) {
     std::cout << "viscous wall flux: " << viscw[UMX] << ", " << viscw[UMY] << ", "
             << viscw[UMZ] << ", " << viscw[UEDEN] << std::endl;
     std::cout << "dapx = " << dapx << "tauxx = " << tauxx << "dapy =" << dapy << "tauxy =" << tauxy
                  << "dapz = " << dapz << "tauxy = " << tauxy << "tauyy = " << tauyy << "tauyz ="<< tauyz << "tauxz ="
                  << tauxz << "tauyz " << tauyz << "tauzz " << tauzz << std::endl;

    std::cout << "dudn = " << dudn << " dvdn = " << dvdn << " dwdn = " << dwdn << std::endl; 
     amrex::Abort("Problem with wall flux");
  } 
#endif
  }

  // Modify if using a wall model
  if (do_wall_model) {
    // Check y+
    const Real dx = 1.0 / dxinv[0]; // dx == dy == dz
    const Real cfarea = apnorm * dx; // cut-face area / dx
    Real tau = std::sqrt(viscw[UMX] * viscw[UMX] + viscw[UMY] * viscw[UMY] +
                         viscw[UMZ] * viscw[UMZ]) / cfarea;
    Real utau = std::sqrt(tau / rho);
    Real yplus = rho * utau / mu * dx;

    if (yplus > 11.0) {
      BL_PROFILE("EB_WallModel");
      // Wall model
      // LawOfTheWall wm;
      EquilibriumODE wm;
      // Tangent vectors (t1.n = 0, t2.t1 = 0, u.t2 = 0)
      Real t1[3];
      Real udotn = u1 * anrmx + v1 * anrmy + w1 * anrmz;
      t1[0] = u1 - udotn * anrmx;
      t1[1] = v1 - udotn * anrmy;
      t1[2] = w1 - udotn * anrmz;
      Real tnorm = std::sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
      t1[0] /= tnorm;
      t1[1] /= tnorm;
      t1[2] /= tnorm;
      // if(std::abs(t1[0] * anrmx + t1[1] * anrmy + t1[2] * anrmz) > 1.e-10) {
      //   std::cout << t1[0] << " " << t1[1] << " " << t1[2] << "; " 
      //             << anrmx << " " << anrmy << " " << anrmz << "; "
      //             << u2 << " " << v2 << " " << w2 << "; "
      //             << std::abs(t1[0] * anrmx + t1[1] * anrmy + t1[2] * anrmz)
      //             << std::endl;
      //   amrex::Abort("Tangent vector is not orthogonal to the normal vector");
      // }
      Real u_parallel = u1 * t1[0] + v1 * t1[1] + w1 * t1[2]; // u dot tangent
      Real ts = (u_parallel > 0.0) ? 1.0 : -1.0;
      u_parallel *= ts; // abs(u)
      if (u_parallel > 100.0 && T1 > 90.0 && T1 < 3000.0 && rho > 0.0) {
        Real T_wall = eb_isothermal ? eb_wall_temp : -1.0;
        Real h = d1 * dx;
        Real tauw, qw;
        wm.parallel_wall_stress(u_parallel, T1, rho, Y, h, mu, lam, T_wall, tauw, qw);
        // std::cout << "y+: " << yplus << " -> "
        //           << rho * std::sqrt(tauw / rho) / mu * dx;
        // std::cout << "  viscous wall stress: " << tau * cfarea << " (" << viscw[UMX]
        //           << ", " << viscw[UMY] << ", " << viscw[UMY] << ") -> ";
        if (!isnan(tauw) && !isnan(qw)) {
        viscw[UMX] = cfarea * ts * t1[0] * tauw;
        viscw[UMY] = cfarea * ts * t1[1] * tauw;
        viscw[UMZ] = cfarea * ts * t1[2] * tauw;
        // std::cout << tauw * cfarea << " (" << viscw[UMX] << ", " << viscw[UMY]
        //           << ", " << viscw[UMZ] << ")";
        viscw[UEDEN] = cfarea * qw;
        // std::cout << "  " << apnorm << " (" << anrmx << ", " << anrmy << ", "
        //           << anrmz << ")\n";
        
        // if ((i == 100 && j == 16 && k == 23))
        //   std::cout << "(" << i <<","<<j<<","<<k<<")" << ": tauw=" << tauw << ", qw=" << qw 
        //             << ", t=" << ts*t1[0] <<"," << ts*t1[1] <<","<< ts*t1[2] << ", cfarea=" << cfarea 
        //             << ", guess y+=" << yplus << std::endl;
        }
      }
    }
  }
}
#endif
#endif
