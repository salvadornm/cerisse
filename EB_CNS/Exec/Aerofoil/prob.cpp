#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
  Real M = 0.729, Re = 6.5e6, aoa = 2.31, Pr = 0.72;
  {
    ParmParse pp("prob");
    pp.query("M", M);
    pp.query("Re", Re);
    pp.query("AoA", aoa);
  }

  Real e, c;
  Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p, massfrac, CNS::h_prob_parm->T,
             CNS::h_prob_parm->rho, e);
  CNS::h_prob_parm->rhoe = CNS::h_prob_parm->rho * e;
  eos.RTY2Cs(CNS::h_prob_parm->rho, CNS::h_prob_parm->T, massfrac, c);
  CNS::h_prob_parm->u = M * c * cos(aoa / 180.0 * M_PI);
  CNS::h_prob_parm->v = M * c * sin(aoa / 180.0 * M_PI);

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();

  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  trans_parm.viscosity_mu_ref =
    CNS::h_prob_parm->rho * CNS::h_prob_parm->u / Re; // length of foil = 1.0cm
  trans_parm.viscosity_T_ref = CNS::h_prob_parm->T;
  trans_parm.viscosity_S = 110.4;
  trans_parm.Prandtl_number = Pr;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;

  CNS::trans_parms.sync_to_device();
}
}

void CNS::fill_ext_src(int i, int j, int k, Real time, GeometryData const& geomdata,
                       Array4<const Real> const& /*state*/,
                       Array4<Real> const& ext_src, ProbParm const& pp)
{
}