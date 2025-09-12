#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real e, T;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(CNS::h_prob_parm->rho_l, massfrac, CNS::h_prob_parm->p_l, e);
  // eos.EY2T(e, massfrac, CNS::h_prob_parm->T_l);
  CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l * e;

  eos.RYP2E(CNS::h_prob_parm->rho_r, massfrac, CNS::h_prob_parm->p_r, e);
  eos.EY2T(e, massfrac, T);
  CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r * e;
  eos.RTY2Cs(CNS::h_prob_parm->rho_r, T, massfrac, CNS::h_prob_parm->c_r);

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}