#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::Real T_1, T_2;
  {
    amrex::ParmParse pp("prob");
    pp.query("p_1", CNS::h_prob_parm->p_1);
    pp.query("u_1", CNS::h_prob_parm->u_1);
    pp.query("v_1", CNS::h_prob_parm->v_1);
    pp.query("T_1", T_1);

    pp.query("p_2", CNS::h_prob_parm->p_2);
    pp.query("u_2", CNS::h_prob_parm->u_2);
    pp.query("v_2", CNS::h_prob_parm->v_2);
    pp.query("T_2", T_2);
  }

  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real X2[NUM_SPECIES]; // burnt
  X2[H2_ID] = 2.68915464e-01;
  X2[O2_ID] = 6.82036376e-06;
  X2[OH_ID] = 1.40332795e-03;
  X2[H2O_ID] = 1.27873575e-01;
  X2[CO_ID] = 2.82651366e-01;
  X2[CO2_ID] = 1.95262241e-02;
  X2[AR_ID] = 2.88751438e-01;
  eos.X2Y(X2, CNS::h_prob_parm->massfrac.begin());

  amrex::Real X[NUM_SPECIES]; // unburnt
  X[C3H8_ID] = 0.18;
  X[O2_ID] = 0.403;
  X[AR_ID] = 0.516;
  eos.X2Y(X, CNS::h_prob_parm->massfrac_2.begin());

  amrex::Real sumY = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) { sumY += CNS::h_prob_parm->massfrac[n]; }
  amrex::Print() << "Sum Y = " << sumY << std::endl;

  amrex::Real e;
  eos.PYT2RE(CNS::h_prob_parm->p_1, CNS::h_prob_parm->massfrac.begin(), T_1,
             CNS::h_prob_parm->rho_1, e);
  CNS::h_prob_parm->rhoe_1 = CNS::h_prob_parm->rho_1 * e;

  eos.PYT2RE(CNS::h_prob_parm->p_2, CNS::h_prob_parm->massfrac_2.begin(), T_2,
             CNS::h_prob_parm->rho_2, e);
  CNS::h_prob_parm->rhoe_2 = CNS::h_prob_parm->rho_2 * e;

  amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                   CNS::d_prob_parm);
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& /*state*/,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& pp)
{
}