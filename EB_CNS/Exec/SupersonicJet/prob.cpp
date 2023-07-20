#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("record_statistics", CNS::h_prob_parm->record_statistics);
    pp.query("clean_aux_on_restart", CNS::h_prob_parm->clean_aux_on_restart);
  }

  auto eos = pele::physics::PhysicsType::eos();

  // Vitiated air inlet
  amrex::Real X[NUM_SPECIES];
  X[O2_ID] = 0.201;
  X[N2_ID] = 0.544;
  X[H2O_ID] = 0.255;
  eos.X2Y(X, CNS::h_prob_parm->Y1.begin());
  // Fuel inlet
  CNS::h_prob_parm->Y2[H2_ID] = 1.0;
  // Atmospheric conditions
  CNS::h_prob_parm->Y_a[N2_ID] = 0.77;
  CNS::h_prob_parm->Y_a[O2_ID] = 0.23;

  // Initial density and ei
  eos.PYT2RE(CNS::h_prob_parm->p1, CNS::h_prob_parm->Y1.begin(),
             CNS::h_prob_parm->T1, CNS::h_prob_parm->rho1, CNS::h_prob_parm->ei1);
  eos.PYT2RE(CNS::h_prob_parm->p2, CNS::h_prob_parm->Y2.begin(),
             CNS::h_prob_parm->T2, CNS::h_prob_parm->rho2, CNS::h_prob_parm->ei2);
  eos.PYT2RE(CNS::h_prob_parm->p_a, CNS::h_prob_parm->Y_a.begin(),
             CNS::h_prob_parm->T_a, CNS::h_prob_parm->rho_a, CNS::h_prob_parm->ei_a);

  amrex::Real c1, c2;
  eos.RTY2Cs(CNS::h_prob_parm->rho1, CNS::h_prob_parm->T1,
             CNS::h_prob_parm->Y1.begin(), c1);
  eos.RTY2Cs(CNS::h_prob_parm->rho2, CNS::h_prob_parm->T2,
             CNS::h_prob_parm->Y2.begin(), c2);

  amrex::Print() << "Confirm inlet Mach number: M_air = "
                 << CNS::h_prob_parm->u1 / c1
                 << " M_fuel = " << CNS::h_prob_parm->u2 / c2 << std::endl;

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