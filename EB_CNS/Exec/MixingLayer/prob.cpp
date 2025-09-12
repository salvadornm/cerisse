#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("T1", CNS::h_prob_parm->T1);
    pp.query("u1", CNS::h_prob_parm->u1);
    pp.query("T2", CNS::h_prob_parm->T2);
    pp.query("u2", CNS::h_prob_parm->u2);
    pp.query("p", CNS::h_prob_parm->p);
    pp.query("vorticity_thickness", CNS::h_prob_parm->theta_w);
    pp.query("record_statistics", CNS::h_prob_parm->record_statistics);
    pp.query("clean_aux_on_restart", CNS::h_prob_parm->clean_aux_on_restart);
  }

  // Fuel stream
  CNS::h_prob_parm->massfrac1[H2_ID] = 0.05;
  CNS::h_prob_parm->massfrac1[N2_ID] = 0.95;
  // Oxidiser stream
  CNS::h_prob_parm->massfrac2[O2_ID] = 0.278;
  CNS::h_prob_parm->massfrac2[H2O_ID] = 0.17;
  CNS::h_prob_parm->massfrac2[H_ID] = 5.6e-7;
  CNS::h_prob_parm->massfrac2[O_ID] = 1.55e-4;
  CNS::h_prob_parm->massfrac2[OH_ID] = 1.83e-3;
  CNS::h_prob_parm->massfrac2[HO2_ID] = 5.1e-6;
  CNS::h_prob_parm->massfrac2[H2O2_ID] = 2.5e-7;
  CNS::h_prob_parm->massfrac2[N2_ID] = 0.55;

  // Initial density and ei
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p, CNS::h_prob_parm->massfrac1.begin(),
             CNS::h_prob_parm->T1, CNS::h_prob_parm->rho1, CNS::h_prob_parm->ei1);
  eos.PYT2RE(CNS::h_prob_parm->p, CNS::h_prob_parm->massfrac2.begin(),
             CNS::h_prob_parm->T2, CNS::h_prob_parm->rho2, CNS::h_prob_parm->ei2);

  amrex::Real c1, c2;
  eos.RTY2Cs(CNS::h_prob_parm->rho1, CNS::h_prob_parm->T1,
             CNS::h_prob_parm->massfrac1.begin(), c1);
  eos.RTY2Cs(CNS::h_prob_parm->rho2, CNS::h_prob_parm->T2,
             CNS::h_prob_parm->massfrac2.begin(), c2);
  CNS::h_prob_parm->uc =
    (c1 * CNS::h_prob_parm->u1 + c2 * CNS::h_prob_parm->u2) / (c1 + c2);

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}