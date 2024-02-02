#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real M_s = 1.5;
  {
    amrex::ParmParse pp("prob");
    pp.query("M_s", M_s);
  }

  // Compute downstream conditions
  amrex::Real g = pele::physics::Constants::gamma;
  CNS::h_prob_parm->rho_r = CNS::h_prob_parm->rho_l * (g + 1.0) * M_s * M_s / (2.0 + (g - 1.0) * M_s * M_s);
  CNS::h_prob_parm->u_r = CNS::h_prob_parm->u_l * (2.0 + (g - 1.0) * M_s * M_s) / ((g + 1.0) * M_s * M_s);
  CNS::h_prob_parm->p_r = CNS::h_prob_parm->p_l * (1.0 + (2.0 * g / (g + 1.0)) * (M_s * M_s - 1.0));

  Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
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