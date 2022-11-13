#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit (const int* /*init*/,
                     const int* /*name*/,
                     const int* /*namelen*/,
                     const amrex_real* /*problo*/,
                     const amrex_real* /*probhi*/)
{
    amrex::Real massfrac[NUM_SPECIES] = {1.0};
    amrex::Real e;

    auto eos = pele::physics::PhysicsType::eos();
    eos.RYP2E(CNS::h_prob_parm->rho_1, massfrac,
              CNS::h_prob_parm->p_1, e);
    eos.EY2T(e, massfrac, CNS::h_prob_parm->T_1);
    CNS::h_prob_parm->rhoe_1 = CNS::h_prob_parm->rho_1 * e;

    eos.RYP2E(CNS::h_prob_parm->rho_2, massfrac,
              CNS::h_prob_parm->p_2, e);
    eos.EY2T(e, massfrac, CNS::h_prob_parm->T_2);
    CNS::h_prob_parm->rhoe_2 = CNS::h_prob_parm->rho_2 * e;

    eos.RYP2E(CNS::h_prob_parm->rho_3, massfrac,
              CNS::h_prob_parm->p_3, e);
    eos.EY2T(e, massfrac, CNS::h_prob_parm->T_3);
    CNS::h_prob_parm->rhoe_3 = CNS::h_prob_parm->rho_3 * e;

    eos.RYP2E(CNS::h_prob_parm->rho_4, massfrac,
              CNS::h_prob_parm->p_4, e);
    eos.EY2T(e, massfrac, CNS::h_prob_parm->T_4);
    CNS::h_prob_parm->rhoe_4 = CNS::h_prob_parm->rho_4 * e;

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm+1,
                     CNS::d_prob_parm);
}
}
