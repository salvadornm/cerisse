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
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  const amrex::Real rsq =
    AMREX_D_TERM((x - 4 * pp.theta_w) * (x - 4 * pp.theta_w), +y * y, +z * z);
  const amrex::Real theta_sq = pp.theta_w * pp.theta_w;

  if (rsq >= theta_sq) {
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      ext_src(i, j, k, nf * NVAR + UMY) +=
        0.01 * pp.uc * exp(-rsq / theta_sq) *
        sin(1e-6 * time); // Yuri says amplitude is 0.1*Uc ??
    }
  }
}