#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real e, T_r;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(CNS::h_prob_parm->rho_l, massfrac, CNS::h_prob_parm->p_l, e);
  // eos.EY2T(e, massfrac, CNS::h_prob_parm->T_l);
  CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l * e;

  eos.RYP2E(CNS::h_prob_parm->rho_r, massfrac, CNS::h_prob_parm->p_r, e);
  eos.EY2T(e, massfrac, T_r);
  CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r * e;
  eos.RTY2Cs(CNS::h_prob_parm->rho_r, T_r, massfrac, CNS::h_prob_parm->c_r);

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& /*state*/,
                       amrex::Array4<amrex::Real> const& ext_src, ProbParm const& pp)
{
}