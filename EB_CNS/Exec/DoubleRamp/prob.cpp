#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
  {
    ParmParse pp("prob");
    pp.query("M", CNS::h_prob_parm->M);
  }

  CNS::h_prob_parm->massfrac[O2_ID] = 0.23;
  CNS::h_prob_parm->massfrac[N2_ID] = 0.77;

  Real e, c;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p, CNS::h_prob_parm->massfrac.begin(),
             CNS::h_prob_parm->T, CNS::h_prob_parm->rho, e);
  CNS::h_prob_parm->rhoe = CNS::h_prob_parm->rho * e;

  eos.RTY2Cs(CNS::h_prob_parm->rho, CNS::h_prob_parm->T,
             CNS::h_prob_parm->massfrac.begin(), c);
  CNS::h_prob_parm->u = CNS::h_prob_parm->M * c;

  Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
            CNS::d_prob_parm);
}
}

void CNS::fill_ext_src(int i, int j, int k, Real time, GeometryData const& geomdata,
                       Array4<const Real> const& /*state*/,
                       Array4<Real> const& ext_src, Parm const& /*parm*/,
                       ProbParm const& pp)
{
}

void DoubleRamp::build(const Geometry& geom, const int max_coarsening_level)
{
  Real ang1 = 15.0, ang2 = 45.0; //, h = 20.0;
  {
    ParmParse pp("double_ramp");
    pp.query("angle1", ang1);
    pp.query("angle2", ang2);
    // pp.query("height", h);
  }

  assert(ang1 < 90.0 && ang1 >= 0.0);
  assert(ang2 < 90.0 && ang2 >= 0.0);

  // Convert angles to radian
  ang1 *= M_PI / 180.0;
  ang2 *= M_PI / 180.0;

  Array<Real, AMREX_SPACEDIM> norm = {0.0};
  Array<Real, AMREX_SPACEDIM> point = {0.0};

  // Ramp plane 1
  point[0] = 0.0;
  point[1] = 0.0;

  norm[0] = sin(ang1);
  norm[1] = -cos(ang1);

  EB2::PlaneIF ramp1(point, norm);

  // Ramp plane 2
  // so that ramp 1 and 2 have the same length
  point[0] = 23.0; // h * cos(ang1) / (sin(ang1) + sin(ang2));
  point[1] = point[0] * tan(ang1);

  norm[0] = sin(ang2);
  norm[1] = -cos(ang2);

  EB2::PlaneIF ramp2(point, norm);

  // Bottom plane
  EB2::PlaneIF plane0({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(0.0, 1.0, 0.0)});

  // Rear plane
  EB2::PlaneIF plane1({AMREX_D_DECL(100.0, 0.0, 0.0)},
                      {AMREX_D_DECL(-1.0, 0.0, 0.0)});

  // auto double_ramp_IF = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF,
  // EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(
  //                                           ramp1, ramp2, plane0, plane1, plane2,
  //                                           plane3, plane4);

  auto r1 = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(
    ramp1, plane0, plane1);
  auto r2 = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(
    ramp2, plane0, plane1);
  auto r12 = EB2::makeUnion(r1, r2);
  auto extrude_r12 = EB2::extrude(r12, 2); // along z

  auto bx = EB2::BoxIF({AMREX_D_DECL(0.0, 0.0, -15.0)},
                       {AMREX_D_DECL(100.0, 20.0, 15.0)}, 0);

  auto double_ramp_IF = EB2::makeIntersection(extrude_r12, bx);

  auto gshop = EB2::makeShop(double_ramp_IF);
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}