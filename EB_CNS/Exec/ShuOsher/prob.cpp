#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  // amrex::ParmParse pp("prob");
  // pp.query("p_l",   CNS::h_prob_parm->p_l);
  // pp.query("u_l",   CNS::h_prob_parm->u_l);
  // pp.query("rho_l", CNS::h_prob_parm->rho_l);
  // pp.query("sd_u_l", CNS::h_prob_parm->sd_u_l);
  // pp.query("p_r",   CNS::h_prob_parm->p_r);
  // pp.query("u_r",   CNS::h_prob_parm->u_r);
  // pp.query("rho_r", CNS::h_prob_parm->rho_r);

  // amrex::Real e_l, e_r;
  // auto eos = pele::physics::PhysicsType::eos();
  // eos.RYP2E(CNS::h_prob_parm->rho_l, CNS::h_prob_parm->massfrac_l,
  //           CNS::h_prob_parm->p_l, e_l);
  // CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l * e_l;

  // eos.RYP2E(CNS::h_prob_parm->rho_r, CNS::h_prob_parm->massfrac_r,
  //           CNS::h_prob_parm->p_r, e_r);
  // CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r * e_r;

  // amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm+1,
  // CNS::d_prob_parm);
}
}