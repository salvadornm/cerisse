#include "prob.H"

#include "CNS.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  molefrac[H2_ID] = 2.0;
  molefrac[O2_ID] = 1.0;
  molefrac[AR_ID] = 7.0;

  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, CNS::h_prob_parm->massfrac_l.begin());
  eos.X2Y(molefrac, CNS::h_prob_parm->massfrac_r.begin());

  eos.RYP2E(CNS::h_prob_parm->rho_l, CNS::h_prob_parm->massfrac_l.begin(),
            CNS::h_prob_parm->p_l, CNS::h_prob_parm->e_l);
  eos.RYP2E(CNS::h_prob_parm->rho_r, CNS::h_prob_parm->massfrac_r.begin(),
            CNS::h_prob_parm->p_r, CNS::h_prob_parm->e_r);

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}