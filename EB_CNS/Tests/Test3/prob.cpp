#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real e;

  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p_1", CNS::h_prob_parm->p_1);
    pp.query("u_1", CNS::h_prob_parm->u_1);
    pp.query("v_1", CNS::h_prob_parm->v_1);
    pp.query("rho_1", CNS::h_prob_parm->rho_1);

    pp.query("p_2", CNS::h_prob_parm->p_2);
    pp.query("u_2", CNS::h_prob_parm->u_2);
    pp.query("v_2", CNS::h_prob_parm->v_2);
    pp.query("rho_2", CNS::h_prob_parm->rho_2);

    pp.query("p_3", CNS::h_prob_parm->p_3);
    pp.query("u_3", CNS::h_prob_parm->u_3);
    pp.query("v_3", CNS::h_prob_parm->v_3);
    pp.query("rho_3", CNS::h_prob_parm->rho_3);

    pp.query("p_4", CNS::h_prob_parm->p_4);
    pp.query("u_4", CNS::h_prob_parm->u_4);
    pp.query("v_4", CNS::h_prob_parm->v_4);
    pp.query("rho_4", CNS::h_prob_parm->rho_4);
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(CNS::h_prob_parm->rho_1, massfrac, CNS::h_prob_parm->p_1, e);
  eos.EY2T(e, massfrac, CNS::h_prob_parm->T_1);
  CNS::h_prob_parm->rhoe_1 = CNS::h_prob_parm->rho_1 * e;

  eos.RYP2E(CNS::h_prob_parm->rho_2, massfrac, CNS::h_prob_parm->p_2, e);
  eos.EY2T(e, massfrac, CNS::h_prob_parm->T_2);
  CNS::h_prob_parm->rhoe_2 = CNS::h_prob_parm->rho_2 * e;

  eos.RYP2E(CNS::h_prob_parm->rho_3, massfrac, CNS::h_prob_parm->p_3, e);
  eos.EY2T(e, massfrac, CNS::h_prob_parm->T_3);
  CNS::h_prob_parm->rhoe_3 = CNS::h_prob_parm->rho_3 * e;

  eos.RYP2E(CNS::h_prob_parm->rho_4, massfrac, CNS::h_prob_parm->p_4, e);
  eos.EY2T(e, massfrac, CNS::h_prob_parm->T_4);
  CNS::h_prob_parm->rhoe_4 = CNS::h_prob_parm->rho_4 * e;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& /*state*/,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& pp)
{
}