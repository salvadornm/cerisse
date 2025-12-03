#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
  auto eos = pele::physics::PhysicsType::eos();

  Real M, p0 = -1.0_rt, T0 = -1.0_rt, p = -1.0_rt, T = -1.0_rt, mom_ratio, M_j = 1.0_rt;
  {
    ParmParse pp("prob");
    pp.get("M", M);                 // inflow Mach number
    pp.query("p0", p0);             // total pressure [Ba]
    pp.query("T0", T0);             // inflow total temperature [K]
    pp.query("p", p);               // static pressure [Ba]
    pp.query("T", T);               // inflow static temperature [K]
    pp.get("mom_ratio", mom_ratio); // jet to inflow momentum ratio
    pp.get("T_j", CNS::h_prob_parm->T_j);
    pp.query("A", CNS::h_prob_parm->A);
    pp.query("theta0", CNS::h_prob_parm->theta0);
    pp.query("M_j", M_j);
    pp.query("r_j", CNS::h_prob_parm->r_j);
    pp.query("do_spark", CNS::h_prob_parm->do_spark);
    pp.query("record_statistics", CNS::h_prob_parm->record_statistics);
    pp.query("clean_aux_on_restart", CNS::h_prob_parm->clean_aux_on_restart);
  }
  if constexpr (NUM_AUX != 13) {
    if (CNS::h_prob_parm->record_statistics)
      amrex::Abort("Please compile with NUM_AUX=13 to record statistics");
  }
  if (!((p0 > 0.0_rt && T0 > 0.0_rt && p < 0.0_rt && T < 0.0_rt) ||
        (p0 < 0.0_rt && T0 < 0.0_rt && p > 0.0_rt && T > 0.0_rt))) {
    amrex::Abort("Please specify either (p0, T0) or (p, T)");
  }

  // Calculate inflow conditions
  CNS::h_prob_parm->Y[O2_ID] = 0.26_rt;
  CNS::h_prob_parm->Y[N2_ID] = 0.74_rt;
  Real rho, gam = 1.4_rt;
  int iter = 0;
  if (p0 > 0.0_rt && T0 > 0.0_rt && p < 0.0_rt && T < 0.0_rt) {
    // Total T, p: Iterate to find T, p, rho, gam
    Real gam_old = 1.0e10_rt;
    while (std::abs(gam - gam_old) > 1.0e-4_rt && iter < 20) {
      iter += 1;
      gam_old = gam;

      // Isentropic relations
      T = T0 / (1.0_rt + 0.5_rt * (gam - 1.0_rt) * M * M);
      p = p0 * std::pow(1.0_rt + 0.5_rt * (gam - 1.0_rt) * M * M, -gam / (gam - 1.0_rt));
      eos.PYT2R(p, CNS::h_prob_parm->Y.begin(), T, rho);
      eos.RTY2G(rho, T, CNS::h_prob_parm->Y.begin(), gam);
    }}
  else if (p0 < 0.0_rt && T0 < 0.0_rt && p > 0.0_rt && T > 0.0_rt) {
    // Static T, p: Find T0, p0, rho, gam directly
    eos.PYT2R(p, CNS::h_prob_parm->Y.begin(), T, rho);
    eos.RTY2G(rho, T, CNS::h_prob_parm->Y.begin(), gam);
    T0 = T * (1.0_rt + 0.5_rt * (gam - 1.0_rt) * M * M);
    p0 = p * std::pow(1.0_rt + 0.5_rt * (gam - 1.0_rt) * M * M, gam / (gam - 1.0_rt));
  }
  Real ei;
  eos.RTY2E(rho, T, CNS::h_prob_parm->Y.begin(), ei);
  Real cs = std::sqrt(gam * p / rho);
  Real u = M * cs;

  CNS::h_prob_parm->ei_inf = ei;
  CNS::h_prob_parm->rho_inf = rho;
  CNS::h_prob_parm->T_inf = T;
  CNS::h_prob_parm->u_inf = u;
  amrex::Print() << "Inflow (gamma, T0, p0, rho, T, p, u, ei) = " << gam << ", "
                 << T0 << ", " << p0 << ", " << rho << ", " << T << ", " << p << ", "
                 << u << ", " << ei << " converged in " << iter << " iterations\n";

  // Calculate jet conditions
  CNS::h_prob_parm->Y_j[H2_ID] = 1.0_rt;
  Real cs_j;
  eos.RTY2Cs(1.0_rt, CNS::h_prob_parm->T_j, CNS::h_prob_parm->Y_j.begin(), cs_j); // cs should be independent of rho
  Real u_j = M_j * cs_j;
  Real rho_j = mom_ratio * rho * u / u_j;
  Real ei_j, gam_j, p_j;
  eos.RTY2E(rho_j, CNS::h_prob_parm->T_j, CNS::h_prob_parm->Y_j.begin(), ei_j);
  eos.RTY2G(rho_j, CNS::h_prob_parm->T_j, CNS::h_prob_parm->Y_j.begin(), gam_j);
  eos.RTY2P(rho_j, CNS::h_prob_parm->T_j, CNS::h_prob_parm->Y_j.begin(), p_j);
  Real T0_j = CNS::h_prob_parm->T_j * (1.0_rt + 0.5_rt * (gam_j - 1.0_rt) * M_j * M_j);
  Real p0_j = p_j * std::pow(1.0_rt + 0.5_rt * (gam_j - 1.0_rt) * M_j * M_j, gam_j / (gam_j - 1.0_rt));

  CNS::h_prob_parm->ei_j = ei_j;
  CNS::h_prob_parm->rho_j = rho_j;
  CNS::h_prob_parm->u_j = u_j;
  amrex::Print() << "Jet    (gamma, T0, p0, rho, T, p, u, ei) = " << gam_j << ", " 
                 << T0_j << ", " << p0_j << ", " << rho_j << ", " << CNS::h_prob_parm->T_j 
                 << ", " << p_j << ", " << u_j << ", " << ei_j << "\n";

  // Report some global numbers
  auto tp = pele::physics::PhysicsType::transport();
  auto* tparm = &CNS::trans_parms.host_trans_parm();
  const bool wtr_get_xi = false;
  const bool wtr_get_mu = true;
  const bool wtr_get_lam = false;
  const bool wtr_get_Ddiag = false;
  const bool wtr_get_chi = false;
  Real muloc, xiloc, lamloc;
  Real *Ddiag = nullptr, *chi_mix = nullptr;
  tp.transport(wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T,
               rho, CNS::h_prob_parm->Y.begin(), Ddiag, chi_mix, muloc, xiloc,
               lamloc, tparm);
  Real Re_inf = rho * u / muloc;
  tp.transport(wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi,
               CNS::h_prob_parm->T_j, rho_j, CNS::h_prob_parm->Y_j.begin(), Ddiag,
               chi_mix, muloc, xiloc, lamloc, tparm);
  Real Re_j = rho_j * u_j / muloc;
  Real A = 100.0_rt;
  Real A_j = M_PI * CNS::h_prob_parm->r_j * CNS::h_prob_parm->r_j;
  Real ER = rho_j * u_j * A_j * CNS::h_prob_parm->Y_j[H2_ID] /
            (rho * u * A * CNS::h_prob_parm->Y[O2_ID]) * 8.0_rt;
  amrex::Print() << "Re_inf = " << Re_inf << ", " << "Re_jet = " << Re_j << ", "
                 << "ER = " << ER << "\n";

  Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
            CNS::d_prob_parm);
}
}
