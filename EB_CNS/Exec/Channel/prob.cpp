#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  amrex::Real Pr = 0.7;
  {
    amrex::ParmParse pp("prob");
    pp.query("H", CNS::h_prob_parm->H);
    pp.query("Re_b", CNS::h_prob_parm->Re_b);
    pp.query("M_b", CNS::h_prob_parm->M_b);
    pp.query("Re_tau", CNS::h_prob_parm->Re_tau);
    pp.query("u_tau", CNS::h_prob_parm->u_tau);
    pp.query("Tw", CNS::h_prob_parm->T_w);
    pp.query("muw", CNS::h_prob_parm->mu_w);
    pp.query("Pr", Pr);
  }

  // Wall quantities
  CNS::h_prob_parm->rho_w = CNS::h_prob_parm->Re_tau * CNS::h_prob_parm->mu_w /
                            CNS::h_prob_parm->u_tau / CNS::h_prob_parm->H;

  // Bulk quantities
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real csw;
  eos.RTY2Cs(CNS::h_prob_parm->rho_w, CNS::h_prob_parm->T_w,
             CNS::h_prob_parm->massfrac.begin(), csw);
  CNS::h_prob_parm->u_b = CNS::h_prob_parm->M_b * csw;
  CNS::h_prob_parm->rho_b = CNS::h_prob_parm->Re_b * CNS::h_prob_parm->mu_w /
                            CNS::h_prob_parm->u_b / CNS::h_prob_parm->H;

  // Forcing per unit density
  CNS::h_prob_parm->f_x = CNS::h_prob_parm->rho_w * CNS::h_prob_parm->u_tau *
                          CNS::h_prob_parm->u_tau / CNS::h_prob_parm->H /
                          CNS::h_prob_parm->rho_b;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();

  // Use Sutherland law to approximate mu/mu_w = (T/T_w)^0.7
  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  trans_parm.viscosity_mu_ref = CNS::h_prob_parm->mu_w;
  trans_parm.viscosity_T_ref = CNS::h_prob_parm->T_w;
  trans_parm.viscosity_S =
    CNS::h_prob_parm->T_w *
    0.463; // use this to approximate mu ~ T^0.7 in T/Tw between 1 and 1.25
  trans_parm.Prandtl_number = Pr;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;

  CNS::trans_parms.sync_to_device();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src, ProbParm const& pp)
{
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    ext_src(i, j, k, nf * NVAR + UMX) += state(i, j, k, nf * NVAR + URHO) * pp.f_x;
    ext_src(i, j, k, nf * NVAR + UEDEN) += state(i, j, k, nf * NVAR + UMX) * pp.f_x;
  }
}