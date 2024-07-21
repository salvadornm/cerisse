#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  amrex::Real e_l, e_r;

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2E(CNS::h_prob_parm->rho_l, massfrac, CNS::h_prob_parm->p_l, e_l);
  eos.EY2T(e_l, massfrac, CNS::h_prob_parm->T_l);
  CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l * e_l;

  eos.RYP2E(CNS::h_prob_parm->rho_r, massfrac, CNS::h_prob_parm->p_r, e_r);
  eos.EY2T(e_r, massfrac, CNS::h_prob_parm->T_r);
  CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r * e_r;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}

/**
 * \brief Fill external source term.
 *
 * \warning This function will overwrite dsdt fab. Always use += instead of =.
 *
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param time      time.
 * @param geomdata  domain geometry data.
 * @param state     state data.
 * @param ext_src   external source term.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 */
void CNS::fill_ext_src(int i, int j, int k, amrex::Real /*time*/,
                       amrex::GeometryData const& /*geomdata*/,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src,
                       ProbParm const& /*prob_parm*/)
{
  // Example
  // const Real g = 9.81 * 100;
  // ext_src(i, j, k, UMZ)   += g * state(i, j, k, URHO);
  // ext_src(i, j, k, UEDEN) += g * state(i, j, k, UMZ);
}