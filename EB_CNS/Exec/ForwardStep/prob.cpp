#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real e, c;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p_l, massfrac, CNS::h_prob_parm->T_l,
             CNS::h_prob_parm->rho_l, e);
  CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l * e;

  eos.PYT2RE(CNS::h_prob_parm->p_r, massfrac, CNS::h_prob_parm->T_r,
             CNS::h_prob_parm->rho_r, e);
  CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r * e;

  eos.RTY2Cs(CNS::h_prob_parm->rho_l, CNS::h_prob_parm->T_l, massfrac, c);
  CNS::h_prob_parm->u_l = 3.0 * c;

  amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                   CNS::d_prob_parm);

  // auto& trans_parm = CNS::trans_parms.host_trans_parm();
  // trans_parm.const_bulk_viscosity = 0.0;
  // trans_parm.const_diffusivity = 0.0;
  // trans_parm.const_viscosity = CNS::h_prob_parm->rho_l * CNS::h_prob_parm->u_l /
  // 1e4; trans_parm.const_conductivity = trans_parm.const_viscosity * 1.005e7 / 0.7;
  // CNS::trans_parms.sync_to_device();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& /*state*/,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& pp)
{
}