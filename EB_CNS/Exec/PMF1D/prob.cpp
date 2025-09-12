#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  CNS::pmf_data.initialize(); // read pmf dat file

  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("x_offset", CNS::h_prob_parm->x_offset);
  pp.query("u_offset", CNS::h_prob_parm->u_offset);
  Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
            CNS::d_prob_parm);
}
}