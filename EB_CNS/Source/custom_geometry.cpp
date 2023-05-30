#include "custom_geometry.H"

//
// Need to define build outside the class for Register to find the class. Why??

#if AMREX_SPACEDIM == 3
void 
Combustor::build (const Geometry& geom, const int max_coarsening_level) 
{
  ParmParse pp("combustor");

  Real fwl;
  pp.get("far_wall_loc",fwl);

  EB2::PlaneIF farwall({AMREX_D_DECL(fwl,0.,0.)},
                       {AMREX_D_DECL(1. ,0.,0.)});

  Vector<Real> pl1pt, pl2pt, pl2nm, pl3pt;
  pp.getarr("ramp_plane1_point", pl1pt);
  pp.getarr("ramp_plane2_point", pl2pt);
  pp.getarr("ramp_plane2_normal", pl2nm);
  pp.getarr("ramp_plane3_point", pl3pt);

  auto ramp = EB2::makeIntersection(EB2::PlaneIF({pl1pt[0], pl1pt[1], 0.},
                                                 {      0.,      -1., 0.}),
                                    EB2::PlaneIF({pl2pt[0], pl2pt[1], 0.},
                                                 {pl2nm[0], pl2nm[1], 0.}),
                                    EB2::PlaneIF({pl3pt[0], pl3pt[1], 0.},
                                                 {      1.,       0., 0.}));

  Vector<Real> pipelo, pipehi;
  pp.getarr("pipe_lo", pipelo);
  pp.getarr("pipe_hi", pipehi);

  EB2::BoxIF pipe({pipelo[0], pipelo[1], -1.}, {pipehi[0], pipehi[1], 1.}, false);

  // where does plane 1 and plane 2 intersect?
  Real k2 = std::abs(pl2nm[0]/pl2nm[1]);
  Real secty = pl2pt[1] + k2*(pl3pt[0]-pl2pt[0]);
  // How much do we cut?
  Real dx = geom.CellSize(0);
  Real dycut = 4.*(1.+max_coarsening_level)*std::min(dx, k2*dx);
  EB2::BoxIF flat_corner({pl3pt[0], 0., -1.}, {1.e10, secty+dycut, 1.}, false);

  auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);

  Real lenx = geom.ProbLength(0);
  Real leny = geom.ProbLength(1);
  auto pr = EB2::translate(EB2::lathe(polys), {lenx*0.5, leny*0.5, 0.});

  auto gshop = EB2::makeShop(pr);
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}

void
ConvergingNozzle::build (const Geometry& geom, const int max_coarsening_level) 
{
  ParmParse pp("converging-nozzle");

  Real d_inlet, l_inlet, d_exit, l_nozzle;
  pp.get("d_inlet", d_inlet);
  pp.get("l_inlet", l_inlet);
  pp.get("d_exit", d_exit);
  pp.get("l_nozzle", l_nozzle);

  EB2::CylinderIF main(0.5 * d_inlet, 0,
    {AMREX_D_DECL(static_cast<amrex::Real>(0.5 * l_inlet), 0, 0)}, true);

  amrex::Real slope_nozzle = (0.5 * d_inlet - 0.5 * d_exit) / l_nozzle;
  amrex::Real norm = -1.0 / slope_nozzle;
  amrex::Real nmag = std::sqrt(1 + 1 / (norm * norm));
  amrex::EB2::PlaneIF nozzle_plane(
    {AMREX_D_DECL(0, 0, 0)},
    {AMREX_D_DECL(
      static_cast<amrex::Real>(1.0 / nmag), slope_nozzle / nmag, 0.0)},
    true);
  auto nozzle = amrex::EB2::translate(
    amrex::EB2::rotate(amrex::EB2::lathe(nozzle_plane), 90 * M_PI / 180, 1),
    {AMREX_D_DECL(
      l_inlet + static_cast<amrex::Real>(0.5 * d_inlet / slope_nozzle),
      0, 0)});

  amrex::EB2::CylinderIF exit(
    0.5 * d_exit, 0, {AMREX_D_DECL(l_inlet + l_nozzle, 0, 0)},
    true);
  auto nozzle_exit = amrex::EB2::makeIntersection(nozzle, exit);

  auto polys = amrex::EB2::makeUnion(main, nozzle_exit);
  auto gshop = amrex::EB2::makeShop(polys);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 6, false);
}
#endif

void 
Triangles::build (const Geometry& geom, const int max_coarsening_level) 
{
  // setting some constants
  // the polygon is triangle
  // we can only do a maximum of 5 triangles (change if needed)
  const int npts_in_tri = 3;
  const int max_tri = 5;

  // number of user defined triangles
  int num_tri;

  amrex::ParmParse pp("triangles");
  amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> alltri(
    npts_in_tri * max_tri);

  // initalize all triangles with some dummy values
  // that fall outside of the domain
  const amrex::Real* problo;
  const amrex::Real* probhi;
  amrex::Real maxlen;

  problo = geom.ProbLo();
  probhi = geom.ProbHi();

  maxlen = amrex::max<amrex::Real>(
    amrex::max<amrex::Real>(geom.ProbLength(0), geom.ProbLength(1)),
    geom.ProbLength(2));

  // setting all triangles to be waaay outside the domain initially
  for (int itri = 0; itri < max_tri; itri++) {
    alltri[npts_in_tri * itri + 0][0] = problo[0] + 100.0 * maxlen;
    alltri[npts_in_tri * itri + 0][1] = problo[1] - 100.0 * maxlen;
    alltri[npts_in_tri * itri + 0][2] = 0.0;

    alltri[npts_in_tri * itri + 1][0] = probhi[0] + 100.0 * maxlen;
    alltri[npts_in_tri * itri + 1][1] = problo[1] - 100.0 * maxlen;
    alltri[npts_in_tri * itri + 1][2] = 0.0;

    alltri[npts_in_tri * itri + 2][0] = probhi[0] + 100.0 * maxlen;
    alltri[npts_in_tri * itri + 2][1] = problo[1] + 100.0 * maxlen;
    alltri[npts_in_tri * itri + 2][2] = 0.0;
  }

  // get user defined number of triangles
  pp.get("num_tri", num_tri);

  for (int itri = 0; itri < num_tri; itri++) {
    amrex::Array<amrex::Real, AMREX_SPACEDIM> point{
      AMREX_D_DECL(0.0, 0.0, 0.0)};

    for (int ipt = 0; ipt < npts_in_tri; ipt++) {
      std::string pointstr = "tri_" + std::to_string(itri) + "_point_" + std::to_string(ipt);
      amrex::Vector<amrex::Real> vecpt;
      pp.getarr(pointstr.c_str(), vecpt, 0, AMREX_SPACEDIM);
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        point[dir] = vecpt[dir];
      }
      alltri[npts_in_tri * itri + ipt] = point;
    }
  }

  // intersection of the 3 planes in a triangle for all triangles
  amrex::Vector<std::unique_ptr<amrex::EB2::IntersectionIF<
    amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>>
    impfunc_triangles(max_tri);

  for (int itri = 0; itri < max_tri; itri++) {
    // make sure points are in anti clockwise direction to set the inside of
    // the triangle as solid phase correctly
    amrex::Array<amrex::Real, AMREX_SPACEDIM> norm0;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> norm1;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> norm2;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> point0;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> point1;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> point2;

    point0 = alltri[npts_in_tri * itri + 0];
    point1 = alltri[npts_in_tri * itri + 1];
    point2 = alltri[npts_in_tri * itri + 2];

    norm0[0] = -(point1[1] - point0[1]);
    norm0[1] = (point1[0] - point0[0]);
    norm0[2] = 0.0;

    norm1[0] = -(point2[1] - point1[1]);
    norm1[1] = (point2[0] - point1[0]);
    norm1[2] = 0.0;

    norm2[0] = -(point0[1] - point2[1]);
    norm2[1] = (point0[0] - point2[0]);
    norm2[2] = 0.0;

    // normalize so that magnitude is 1
    amrex::Real norm = sqrt(norm0[0] * norm0[0] + norm0[1] * norm0[1]);
    norm0[0] = norm0[0] / norm;
    norm0[1] = norm0[1] / norm;

    // normalize so that magnitude is 1
    norm = sqrt(norm1[0] * norm1[0] + norm1[1] * norm1[1]);
    norm1[0] = norm1[0] / norm;
    norm1[1] = norm1[1] / norm;

    // normalize so that magnitude is 1
    norm = sqrt(norm2[0] * norm2[0] + norm2[1] * norm2[1]);
    norm2[0] = norm2[0] / norm;
    norm2[1] = norm2[1] / norm;

    amrex::EB2::PlaneIF plane0(point0, norm0);
    amrex::EB2::PlaneIF plane1(point1, norm1);
    amrex::EB2::PlaneIF plane2(point2, norm2);

    impfunc_triangles[itri] = std::make_unique<amrex::EB2::IntersectionIF<
      amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>(

      plane0, plane1, plane2);
  }

  auto alltri_IF = amrex::EB2::makeUnion(
    *impfunc_triangles[0], *impfunc_triangles[1], *impfunc_triangles[2],
    *impfunc_triangles[3], *impfunc_triangles[4]);

  auto alltri_extrude_IF = amrex::EB2::extrude(alltri_IF, 2); // along z

  auto gshop = amrex::EB2::makeShop(alltri_extrude_IF);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}