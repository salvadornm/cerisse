#include "prob.H"

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("reynolds", CNS::h_prob_parm->reynolds);
    pp.query("mach", CNS::h_prob_parm->mach);
    pp.query("prandtl", CNS::h_prob_parm->prandtl);
    pp.query("convecting", CNS::h_prob_parm->convecting);
    pp.query("omega_x", CNS::h_prob_parm->omega_x);
    pp.query("omega_y", CNS::h_prob_parm->omega_y);
    pp.query("omega_z", CNS::h_prob_parm->omega_z);
  }

  // Define the length scale
  CNS::h_prob_parm->L = 1.0 / M_PI;
  // CNS::h_prob_parm->L_x = probhi[0] - problo[0];
  // CNS::h_prob_parm->L_y = probhi[1] - problo[1];
  // CNS::h_prob_parm->L_z = probhi[2] - problo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    CNS::h_prob_parm->p0, massfrac, CNS::h_prob_parm->T0,
    CNS::h_prob_parm->rho0, eint);
  eos.RTY2Cs(
    CNS::h_prob_parm->rho0, CNS::h_prob_parm->T0, massfrac,
    cs);
  eos.TY2Cp(CNS::h_prob_parm->T0, massfrac, cp);
  CNS::h_prob_parm->v0 = CNS::h_prob_parm->mach * cs;

  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity =
    CNS::h_prob_parm->rho0 * CNS::h_prob_parm->v0 * CNS::h_prob_parm->L / CNS::h_prob_parm->reynolds;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / CNS::h_prob_parm->prandtl;
  CNS::trans_parms.sync_to_device();

  // // Output IC
  amrex::Print() << "V0 = " << CNS::h_prob_parm->v0 << " p0 = " << CNS::h_prob_parm->p0 
              << " rho0 = " << CNS::h_prob_parm->rho0 << " T0 = " << CNS::h_prob_parm->T0 << std::endl;

  // std::ofstream ofs("ic.txt", std::ofstream::out);
  // amrex::Print(ofs) << "L, rho0, v0, p0, T0, gamma, mu, k, c_s0, Reynolds, "
  //                      "Mach, Prandtl, omega_x, omega_y, omega_z"
  //                   << std::endl;
  // amrex::Print(ofs).SetPrecision(17)
  //   << CNS::h_prob_parm->L << "," << CNS::h_prob_parm->rho0
  //   << "," << CNS::h_prob_parm->v0 << ","
  //   << CNS::h_prob_parm->p0 << "," << CNS::h_prob_parm->T0
  //   << "," << eos.gamma << "," << trans_parm.const_viscosity << ","
  //   << trans_parm.const_conductivity << "," << cs << ","
  //   << CNS::h_prob_parm->reynolds << ","
  //   << CNS::h_prob_parm->mach << ","
  //   << CNS::h_prob_parm->prandtl << ","
  //   << CNS::h_prob_parm->omega_x << ","
  //   << CNS::h_prob_parm->omega_y << ","
  //   << CNS::h_prob_parm->omega_z << std::endl;
  // ofs.close();

  amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm+1, CNS::d_prob_parm);
}
}

void
CNS::fill_ext_src (int i, int j, int k, 
                   amrex::Real /*time*/,
                   amrex::GeometryData const& /*geomdata*/, 
                   amrex::Array4<const amrex::Real> const& state, 
                   amrex::Array4<amrex::Real> const& ext_src, 
                   Parm const& /*parm*/,
                   ProbParm const& /*prob_parm*/)
{
}