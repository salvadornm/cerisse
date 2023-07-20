#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("mach", CNS::h_prob_parm->mach);
    pp.query("beta", CNS::h_prob_parm->beta);
  }

  // Initial density and velocity
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p0, massfrac, CNS::h_prob_parm->T0,
             CNS::h_prob_parm->rho0, eint);
  eos.RTY2Cs(CNS::h_prob_parm->rho0, CNS::h_prob_parm->T0, massfrac, cs);
  eos.TY2Cp(CNS::h_prob_parm->T0, massfrac, cp);
  CNS::h_prob_parm->v0 = CNS::h_prob_parm->mach * cs;

  // gamma = 1.4

  // Output IC
  amrex::Print() << "v0, p0, rho0, T0, beta = " << CNS::h_prob_parm->v0 << ", "
                 << CNS::h_prob_parm->p0 << ", " << CNS::h_prob_parm->rho0 << ", "
                 << CNS::h_prob_parm->T0 << ", " << CNS::h_prob_parm->beta
                 << std::endl;

  amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                   CNS::d_prob_parm);
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real t,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& prob_parm)
{
}