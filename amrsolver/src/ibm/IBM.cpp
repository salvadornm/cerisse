#include <AMReX_ParmParse.H>
#include <IBM.h>

using namespace amrex;
using namespace IBM;

IBFab::IBFab(const Box& b, int ncomp, bool alloc, bool shared, Arena* ar)
    : BaseFab<bool>(b, ncomp, alloc, shared, ar) {
  // gpData = (gpData_t*)The_Managed_Arena()->alloc(sizeof(gpData_t));
}

IBFab::IBFab(const IBFab& rhs, MakeType make_type, int scomp, int ncomp)
    : BaseFab<bool>(rhs, make_type, scomp, ncomp) {}

IBFab::~IBFab() {
  // The_Managed_Arena()->free(gpData);
}

IBMultiFab::IBMultiFab(const BoxArray& bxs, const DistributionMapping& dm,
                       const int nvar, const int ngrow, const MFInfo& info,
                       const FabFactory<IBFab>& factory)
    : FabArray<IBFab>(bxs, dm, nvar, ngrow, info, factory) {}
IBMultiFab::~IBMultiFab() {}

IBMultiFab::IBMultiFab(IBMultiFab&& rhs) noexcept
    : FabArray<IBFab>(std::move(rhs)) {
#ifdef AMREX_MEM_PROFILING
  ++num_multifabs;
  num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

// for a single level
void IBMultiFab::copytoRealMF(MultiFab& mf, int ibcomp, int mfcomp) {
  for (MFIter mfi(*this, false); mfi.isValid(); ++mfi) {
    // const Box& ibbox = mfi.fabbox(); // box with ghost points
    const Box& ibbox = mfi.validbox();  // box without ghost points

    IBM::IBFab& ibfab = this->get(mfi);
    FArrayBox& realfab = mf.get(mfi);
    Array4<bool> ibMarkers = ibfab.array();    // boolean array
    Array4<Real> realfield = realfab.array();  // real array

    amrex::ParallelFor(
        ibbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          realfield(i, j, k, mfcomp) = ibMarkers(i, j, k, ibcomp);
          realfield(i, j, k, mfcomp + 1) = ibMarkers(i, j, k, ibcomp + 1);
        });
  }
}

// constructor and destructor
IB::IB() {}
IB::~IB() {
  // clear memory
  for (int ii = 0; ii < IB::NGEOM; ii++) {
    delete VtreePtr.at(ii);
    delete VInOutFunc.at(ii);
  };
}

// initialise IB
void IB::init(Amr* pointer_amr, int nghost) {
  IB::NGHOST_IB = nghost;
  IB::pamr = pointer_amr;  // store pointer to main Amr class object's instance
  IB::ref_ratio = pamr->refRatio();
  IB::MAX_LEVEL = pamr->maxLevel();
  IB::ibMFa.resize(MAX_LEVEL + 1);
  IB::lsMFa.resize(MAX_LEVEL + 1);

  // TODO make the next 15 lines simpler
  IB::cellSizes.resize(IB::MAX_LEVEL + 1);
  IB::disIM.resize(IB::MAX_LEVEL + 1);
  cellSizes[0] = pamr->Geom(0).CellSizeArray();
  for (int i = 1; i <= IB::MAX_LEVEL; i++) {
    for (int j = 0; j < AMREX_SPACEDIM; j++) {
      cellSizes[i][j] = cellSizes[i - 1][j] / ref_ratio[i - 1][j];
    }
  }

  // TODO make distance ip a parameter
  for (int i = 0; i <= IB::MAX_LEVEL; i++) {
    IB::disIM[i] = 0.6_rt * sqrt(cellSizes[i][0] * cellSizes[i][0] +
                                 cellSizes[i][1] * cellSizes[i][1] +
                                 cellSizes[i][2] * cellSizes[i][2]);
  }

  // read geometry
  readGeom();
}
// create IBMultiFabs at a level and store pointers to it
void IB::buildMFs(const BoxArray& bxa, const DistributionMapping& dm, int lev) {
  ibMFa[lev] = new IBMultiFab(bxa, dm, NVAR_IB, NGHOST_IB);
  lsMFa[lev].define(bxa, dm, 1, NGHOST_IB);
}

void IB::destroyMFs(int lev) {
  if (!ibMFa.empty()) {
    delete ibMFa.at(lev);
  }
  if (!lsMFa.empty()) {
    lsMFa[lev].clear();
  }
}

void IB::computeMarkers(int lev) {
  IBMultiFab& mfab = *ibMFa[lev];
  int nhalo =
      mfab.nGrow(0);  // assuming same number of ghost points in all directions
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = pamr->Geom(lev).ProbLoArray();

  for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
    IBM::IBFab& ibFab = mfab.get(mfi);
    const Box& bx = mfi.tilebox();
    const IntVect& lo = bx.smallEnd();
    const IntVect& hi = bx.bigEnd();
    auto const& ibMarkers = mfab.array(mfi);  // boolean array

    // compute sld markers (including at ghost points) - cannot use ParallelFor
    // - CGAL call causes problems
    for (int k = lo[2] - nhalo; k <= hi[2] + nhalo; ++k) {
      for (int j = lo[1] - nhalo; j <= hi[1] + nhalo; ++j) {
        for (int i = lo[0] - nhalo; i <= hi[0] + nhalo; ++i) {
          ibMarkers(i, j, k, 0) = false;  // initialise to false
          ibMarkers(i, j, k, 1) = false;  // initialise to false

          Real x = prob_lo[0] + (0.5_rt + Real(i)) * cellSizes[lev][0];
          Real y = prob_lo[1] + (0.5_rt + Real(j)) * cellSizes[lev][1];
          Real z = prob_lo[2] + (0.5_rt + Real(k)) * cellSizes[lev][2];
          Point gridpoint(x, y, z);

          for (int ii = 0; ii < IB::NGEOM; ii++) {
            // Print() << i << " " << j << " " << k << " " << ii << std::endl;
            inside_t& inside = *VInOutFunc[ii];
            CGAL::Bounded_side result = inside(gridpoint);
            AMREX_ASSERT_WITH_MESSAGE((result != CGAL::ON_BOUNDARY),
                                      "Grid point on IB surface");
            // if point inside any IB geometry, mark as solid, move on to
            // another point. This minimises the number of inout testing
            // (expensive) calls.
            if (int(result) == int(CGAL::ON_BOUNDED_SIDE)) {
              ibMarkers(i, j, k, 0) = true;
              break;
            }
          }
        }
      }
    };

    // compute ghost markers
    ibFab.gpData.ngps = 0;
    int nextra = 1;
    for (int k = lo[2] - nextra; k <= hi[2] + nextra; ++k) {
      for (int j = lo[1] - nextra; j <= hi[1] + nextra; ++j) {
        for (int i = lo[0] - nextra; i <= hi[0] + nextra; ++i) {
          bool ghost = false;
          if (ibMarkers(i, j, k, 0)) {
            for (int l = -1; l <= 1; l = l + 2) {
              ghost = ghost || (!ibMarkers(i + l, j, k, 0));
              ghost = ghost || (!ibMarkers(i, j + l, k, 0));
              ghost = ghost || (!ibMarkers(i, j, k + l, 0));
            }
            ibMarkers(i, j, k, 1) = ghost;
            ibFab.gpData.ngps += ghost;

            if (ghost) {
              // store GP index
              ibFab.gpData.gp_ijk.push_back(
                  Array1D<int, 0, AMREX_SPACEDIM - 1>{i, j, k});
            } else {
              ibMarkers(i, j, k, 1) = false;
            }
          }
        }
      };
    }
  }
}

void IB::initialiseGPs(int lev) {
  IBMultiFab& mfab = *ibMFa[lev];
  int nhalo =
      mfab.nGrow(0);  // assuming same number of ghost points in all directions
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = pamr->Geom(lev).ProbLoArray();

  for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
    IBM::IBFab& ibFab = mfab.get(mfi);
    gpData_t& gpData = ibFab.gpData;
    const Box& bxg = mfi.growntilebox(nhalo);
    // const Box& bx = mfi.tilebox();
    const IntVect& lo = bxg.smallEnd();
    const IntVect& hi = bxg.bigEnd();
    auto const& ibMarkers = mfab.array(mfi);  // boolean array
    auto const idxCube = ibFab.gpData.indexCube.data();

    // we need a CPU loop here (cannot be GPU loop) as CGAL tree seach for
    // closest element to a point needs to be called. instead of looping through
    // previously indexed gps, we loop through the whole ghost point field as it
    // is available on GPU and CPU at all times. Unlike the gp indexes, which
    // are only stored on GPU memory. Array1D<int,0,AMREX_SPACEDIM-1>& idx =
    // ibFab.gpData.gp_ijk[ii];

    for (int k = lo[2]; k <= hi[2]; ++k) {
      for (int j = lo[1]; j <= hi[1]; ++j) {
        for (int i = lo[0]; i <= hi[0]; ++i) {
          // for each ghost point
          if (ibMarkers(i, j, k, 1)) {
            Real x = prob_lo[0] + (0.5_rt + i) * cellSizes[lev][0];
            Real y = prob_lo[1] + (0.5_rt + j) * cellSizes[lev][1];
            Real z = prob_lo[2] + (0.5_rt + k) * cellSizes[lev][2];
            Point gp(x, y, z);

            // find and store geometery index for this GP. This index is used
            // for searching appropriate tree in initialiseGPs to find the
            // closest element and searching all trees.

            // in out test for each geometry
            int igeom;
            for (int ii = 0; ii < IB::NGEOM; ii++) {
              inside_t& inside = *VInOutFunc[ii];
              CGAL::Bounded_side result = inside(gp);

              if (int(result) == int(CGAL::ON_BOUNDED_SIDE)) {
                igeom = ii;
                ibFab.gpData.geomIdx.push_back(igeom);
                // if (ii == 1) {Print() << ii << std::endl;}
                break;
              }
              // TODO: assert geometries do not overlap?
            }

            // closest surface point and face --------------------------
            Point_and_primitive_id closest_elem =
                VtreePtr[igeom]->closest_point_and_primitive(gp);

            // store
            gpData.closest_cgal.push_back(closest_elem);

            // This closest point (cp) is between the face plane and the gp
            Point cp = closest_elem.first;
            Polyhedron::Face_handle face = closest_elem.second;

            // amrex::Print() << "------------------- " << std::endl;
            // Print() << "closest surface point: " << cp << std::endl;
            // Print() << "closest triangle: ( "
            //           << face->halfedge()->vertex()->point() << " , "
            //           << face->halfedge()->next()->vertex()->point() << " , "
            //           << face->halfedge()->next()->next()->vertex()->point()
            //           << " )"
            //           << std::endl;
            // Print() << "Normal " << Vfnormals[igeom][face] <<std::endl;
            // Print() << "cp-gp " << cp - gp << std::endl; // should be in the
            // direction of normal Print() << "Plane " << face->plane().a() << "
            // " << face->plane().b() << " "  << face->plane().c() << " " <<
            // face->plane().d() << std::endl;

            // IB point -------------------------------------------
            Vector_CGAL imp_gp(gp, cp);
            Real disGP = sqrt(CGAL::squared_distance(gp, cp));
            AMREX_ASSERT_WITH_MESSAGE(
                disGP < 1.0 * sqrt(cellSizes[lev][0] * cellSizes[lev][0] +
                                   cellSizes[lev][1] * cellSizes[lev][1] +
                                   cellSizes[lev][2] * cellSizes[lev][2]),
                "Ghost point and IB point distance larger than mesh diagonal");

            //*store*
            gpData.disGP.push_back(disGP);
            Array1D<Real, 0, AMREX_SPACEDIM - 1> norm = {
                Vfnormals[igeom][face][0], Vfnormals[igeom][face][1],
                Vfnormals[igeom][face][2]};

            Point p1 = face->halfedge()->vertex()->point();
            Point p2 = face->halfedge()->next()->vertex()->point();
            Vector_CGAL tan1_not_unit(p1, p2);
            Real len = sqrt(CGAL::squared_distance(p1, p2));
            Array1D<Real, 0, AMREX_SPACEDIM - 1> tan1;
            tan1(0) = tan1_not_unit[0] / len;
            tan1(1) = tan1_not_unit[1] / len;
            tan1(2) = tan1_not_unit[2] / len;

            // norm x tan1
            Array1D<Real, 0, AMREX_SPACEDIM - 1> tan2 = {
                norm(1) * tan1(2) - norm(2) * tan1(1),
                norm(2) * tan1(0) - norm(0) * tan1(2),
                norm(0) * tan1(1) - norm(1) * tan1(0)};

            // norm.tan1 and norm.tan2 == 0
            AMREX_ASSERT_WITH_MESSAGE(
                (norm(0) * tan1(0) + norm(1) * tan1(1) + norm(2) * tan1(2) +
                 norm(0) * tan2(0) + norm(1) * tan2(1) + norm(2) * tan2(2)) <
                    1.0e-9,
                "norm.tan1 or norm.tan2 not orthogonal");

            gpData.normal.push_back(norm);
            gpData.tangent1.push_back(tan1);
            gpData.tangent2.push_back(tan2);

            // IM points -------------------------------------------
            Array2D<Real, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
            Array2D<int, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;

            // find image point and the bottom left point closest to the image
            // point
            //  In 2D, same idea in 3D.
            //     i,j+1 (2) ---------------------     i+1,j+1 (3)
            //     |                                  |
            //     |         P                        |
            //     |                                  |
            //     |                                  |
            //     i,j  (1) ----------------------      i+1,j  (4)
            for (int jj = 0; jj < NIMPS; jj++) {
              for (int kk = 0; kk < AMREX_SPACEDIM; kk++) {
                imp_xyz(jj, kk) = cp[kk] + Real(jj + 1) * disIM[lev] *
                                               Vfnormals[igeom][face][kk];
                imp_ijk(jj, kk) =
                    floor((imp_xyz(jj, kk) - prob_lo[kk]) / cellSizes[lev][kk] -
                          0.5_rt);
              }

              AMREX_ASSERT_WITH_MESSAGE(
                  bxg.contains(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2)),
                  "Interpolation point outside fab");
            }
            // DEBUGGING //////////////
            // Print() << "disGP, disIM (from IB) " << fab.gpArray[ii].disGP <<
            // " " << IB::disIM[lev] << std::endl; for (int jj=0; jj<NIMPS;
            // jj++) {
            //   Print() << "imp" << jj+1 << " "<< fab.gpArray[ii].imps[jj] <<
            //   std::endl; Print() << "impCell" << jj+1 << " "<<
            //   fab.gpArray[ii].impscell[jj] << std::endl;
            // }
            ///////////////////////////
            // *store*
            gpData.disIM.push_back(disIM[lev]);
            gpData.imp_xyz.push_back(imp_xyz);
            gpData.imp_ijk.push_back(imp_ijk);

            // Interpolation points' (ips) weights for each image point
            Array2D<Real, 0, NIMPS - 1, 0, 7> ipweights;
            computeIPweights(ipweights, imp_xyz, imp_ijk, prob_lo,
                             cellSizes[lev], ibMarkers, idxCube);
            // *store*
            gpData.imp_ipweights.push_back(ipweights);

            //  if (k==35) {
            //   int jj = 0;
            //   Print() << "--- " << std::endl;
            //   Print() << "gp array idx " << gpData.normal.size() - 1 <<
            //   std::endl; Print() << "bxg " << bxg << std::endl; Print() <<
            //   "gp_ijk " << i << " " << j << " " << k << std::endl; Print() <<
            //   "im_ijk " << imp_ijk(jj,0) << " " << imp_ijk(jj,1) << " " <<
            //   imp_ijk(jj,2) << std::endl; Print() << "gp_xyz " << x << " " <<
            //   y << " " << z << std::endl; Print() << "ib_xyz " << cp[0] << "
            //   " << cp[1] << " " << cp[2] << std::endl; Print() << "ip_xyz "
            //   << imp_xyz(jj,0) << " " << imp_xyz(jj,1) << " " <<
            //   imp_xyz(jj,2) << std::endl; Print() << "norm " << norm(0) << "
            //   " << norm(1) << " " << norm(2) << std::endl; Print() << "tan1 "
            //   << tan1(0) << " " << tan1(1) << " " << tan1(2) << std::endl;
            //   Print() << "tan2 " << tan2(0) << " " << tan2(1) << " " <<
            //   tan2(2) << std::endl; Print() << "dx  " << cellSizes[lev][0] <<
            //   " " << cellSizes[lev][1] << " " << cellSizes[lev][2] <<
            //   std::endl; Real temp = sqrt(cellSizes[lev][0]*cellSizes[lev][0]
            //   + cellSizes[lev][1]*cellSizes[lev][1] +
            //   cellSizes[lev][2]*cellSizes[lev][2]); Print() << "dx diag " <<
            //   temp << std::endl; Print() << "disIM "  << disIM[lev] <<
            //   std::endl; Print() << "disGP " << disGP << std::endl; Print()
            //   << "weights " << ipweights(jj,0) << " " << ipweights(jj,1) << "
            //   " << ipweights(jj,2) << " " << ipweights(jj,3) << " " <<
            //   ipweights(jj,4) << " " << ipweights(jj,5) << " " <<
            //   ipweights(jj,6) << " " << ipweights(jj,7) << std::endl;
            //   }
          }
        }
      }
    }
  }
}

void IB::computeIPweights(
    Array2D<Real, 0, NIMPS - 1, 0, 7>& weights,
    Array2D<Real, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
    Array2D<int, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    GpuArray<Real, AMREX_SPACEDIM>& dxyz, const Array4<bool> ibFab,
    auto const idxCube) {
  for (int iim = 0; iim < NIMPS; iim++) {
    int i = imp_ijk(iim, 0);
    int j = imp_ijk(iim, 1);
    int k = imp_ijk(iim, 2);
    // tri-linear interpolation
    // note xl,xr, ...etc do not have prob_lo added to them. This does not
    // matter as we only need the relative distances between the points.
    Real xl = prob_lo[0] +
              (Real(i) + 0.5_rt) * dxyz[0];  // bottom left corner of cell
    Real xr = xl + dxyz[0];
    Real yl = prob_lo[1] + (Real(j) + 0.5_rt) * dxyz[1];
    Real yr = yl + dxyz[1];
    Real zl = prob_lo[2] + (Real(k) + 0.5_rt) * dxyz[2];
    Real zr = zl + dxyz[2];

    Real xd = (imp_xyz(iim, 0) - xl) / (xr - xl);
    Real yd = (imp_xyz(iim, 1) - yl) / (yr - yl);
    Real zd = (imp_xyz(iim, 2) - zl) / (zr - zl);

    int sumfluid = 0;
    Real sumweights = 0.0_rt;
    // zd = 0
    int ii = i + idxCube[0](0);
    int jj = j + idxCube[0](1);
    int kk = k + idxCube[0](2);
    int fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 0) = (1.0_rt - xd) * (1.0_rt - yd) * (1.0_rt - zd) * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 0);

    ii = i + idxCube[1](0);
    jj = j + idxCube[1](1);
    kk = k + idxCube[1](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 1) = (1.0_rt - xd) * yd * (1.0_rt - zd) * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 1);

    ii = i + idxCube[2](0);
    jj = j + idxCube[2](1);
    kk = k + idxCube[2](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 2) = xd * yd * (1.0_rt - zd) * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 2);

    ii = i + idxCube[3](0);
    jj = j + idxCube[3](1);
    kk = k + idxCube[3](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 3) = xd * (1.0_rt - yd) * (1.0_rt - zd) * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 3);

    // zd = 2
    ii = i + idxCube[4](0);
    jj = j + idxCube[4](1);
    kk = k + idxCube[4](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 4) = (1.0_rt - xd) * (1.0_rt - yd) * zd * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 4);

    ii = i + idxCube[5](0);
    jj = j + idxCube[5](1);
    kk = k + idxCube[5](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 5) = (1.0_rt - xd) * yd * zd * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 5);

    ii = i + idxCube[6](0);
    jj = j + idxCube[6](1);
    kk = k + idxCube[6](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 6) = xd * yd * zd * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 6);

    ii = i + idxCube[7](0);
    jj = j + idxCube[7](1);
    kk = k + idxCube[7](2);
    fluid = !ibFab(ii, jj, kk, 0);
    weights(iim, 7) = xd * (1.0_rt - yd) * zd * fluid;
    sumfluid += fluid;
    sumweights += weights(iim, 7);

    AMREX_ASSERT_WITH_MESSAGE(
        sumfluid >= 2, "Less than 2 interpolation points are fluid points");

    // re-normalise
    for (int ll = 0; ll < 8; ll++) {
      weights(iim, ll) = weights(iim, ll) / sumweights;
    }

    AMREX_ASSERT_WITH_MESSAGE(
        std::abs(weights(iim, 0) + weights(iim, 1) + weights(iim, 2) +
                 weights(iim, 3) + weights(iim, 4) + weights(iim, 5) +
                 weights(iim, 6) + weights(iim, 7) - Real(1.0)) < Real(1.e-9),
        "Interpolation point weights do not sum to 1.0");
  }
}

// Taylor expansion around IB point
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void extrapolateGP(
    Array2D<Real, 0, NIMPS + 1, 0, NPRIM - 1>& state, Real dgp, Real dim) {
  Real sgn_dgp = -dgp;  // negative sign as taylor expansion is around IB point,
                        // IM and GP are in opposite directions
  for (int kk = 0; kk < NPRIM; kk++) {
    // Linear
    // Real c1 = state(1,kk);
    // Real c2 = (state(2,kk) - state(1,kk))/dim;
    // state(0,kk) = c1 + c2*sgn_dgp;

    // Linear + Van Leer limiter
    Real c1 = state(1, kk);
    Real ratio =
        (state(3, kk) - state(2, kk)) / (state(2, kk) - state(1, kk) + 1.0e-12);
    Real phi = (ratio + std::abs(ratio)) / (1.0_rt + std::abs(ratio));
    Real c2 = phi * (state(2, kk) - state(1, kk)) / dim;
    state(0, kk) = c1 + c2 * sgn_dgp;

    // Real c1    =  state(1,kk);
    // Real c2    =  (-1.5_rt*state(1,kk) + 2.0_rt*state(2,kk) -
    // 0.5_rt*state(3,kk))/sgn_dgp; Real c3    =  (0.5_rt*state(1,kk) -
    // state(2,kk) + 0.5_rt*state(3,kk))/sgn_dgp/sgn_d; state(0,kk) = c1 +
    // c2*dim + c3*dim*dim;
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void global2local(
    int jj, Array2D<Real, 0, NIMPS + 1, 0, NPRIM - 1>& primStateNormal,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& norm,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tan1,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tan2) {
  Array1D<Real, 0, AMREX_SPACEDIM - 1> vel;
  vel(0) = primStateNormal(jj, QU);
  vel(1) = primStateNormal(jj, QV);
  vel(2) = primStateNormal(jj, QW);

  primStateNormal(jj, QU) =
      vel(0) * norm(0) + vel(1) * norm(1) + vel(2) * norm(2);
  primStateNormal(jj, QV) =
      vel(0) * tan1(0) + vel(1) * tan1(1) + vel(2) * tan1(2);
  primStateNormal(jj, QW) =
      vel(0) * tan2(0) + vel(1) * tan2(1) + vel(2) * tan2(2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void local2global(
    int jj, Array2D<Real, 0, NIMPS + 1, 0, NPRIM - 1>& primStateNormal,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& norm,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tan1,
    const Array1D<Real, 0, AMREX_SPACEDIM - 1>& tan2) {
  Array1D<Real, 0, AMREX_SPACEDIM - 1> vel;
  vel(0) = primStateNormal(jj, QU);
  vel(1) = primStateNormal(jj, QV);
  vel(2) = primStateNormal(jj, QW);

  primStateNormal(jj, QU) =
      vel(0) * norm(0) + vel(1) * tan1(0) + vel(2) * tan2(0);
  primStateNormal(jj, QV) =
      vel(0) * norm(1) + vel(1) * tan1(1) + vel(2) * tan2(1);
  primStateNormal(jj, QW) =
      vel(0) * norm(2) + vel(1) * tan1(2) + vel(2) * tan2(2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void ComputeGPState(
    int ii, auto const gp_ijk, auto const imp_ijk, auto const weights,
    auto const norm, auto const tan1, auto const tan2, const Real disGP,
    const Real disIM, const Array4<bool> ibFab, Array4<Real> primFab,
    Array4<Real> conFab, auto const idxCube,
    const PROB::ProbClosures& closures) {
  // if (norm(2)> 0.9999_rt) {
  // Print() << "-------" << std::endl;
  // Print() << "computeGP ii,i,j,k " << ii << " " << gp_ijk(0) << " " <<
  // gp_ijk(1) << " " << gp_ijk(2)  << std::endl; Print() << "imp_ijk " <<
  // imp_ijk(0,0) << " " << imp_ijk(0,1) << " " << imp_ijk(0,2) << std::endl;
  // Print() << "norm " << norm(0) << " " << norm(1) << " " << norm(2) <<
  // std::endl; Print() << "weights " << weights(0,0) << " " << weights(0,1) <<
  // " " << weights(0,2) << " " << weights(0,3) << " " << weights(0,4) << " " <<
  // weights(0,5) << " " << weights(0,6) << " " << weights(0,7) << std::endl;
  // }

  // interpolate IM ////////////////
  Array2D<Real, 0, NIMPS + 1, 0, NPRIM - 1> primStateNormal = {0.0};
  // for each image point
  for (int jj = 0; jj < NIMPS; jj++) {
    // GP (jj=0),IB (jj=1), IM1 (jj=2),IM2 (jj=3)...
    // for each IP point
    for (int ll = 0; ll < 8; ll++) {
      int itemp = imp_ijk(jj, 0) + idxCube[ll](0);
      int jtemp = imp_ijk(jj, 1) + idxCube[ll](1);
      int ktemp = imp_ijk(jj, 2) + idxCube[ll](2);
      // for each primitive variable
      for (int kk = 0; kk < NPRIM; kk++) {
        primStateNormal(jj + 2, kk) +=
            primFab(itemp, jtemp, ktemp, kk) * weights(jj, ll);
      }
    }
  }

  // Convert ux,uy,uz to un,ut1,ut2 at IM points only
  for (int jj = 2; jj < 2 + NIMPS; jj++) {
    global2local(jj, primStateNormal, norm, tan1, tan2);
  }

  // interpolate IB (enforcing BCs) -- TODO:generalise to IBM BCs and make user
  // defined no-slip velocity (in local coordinates)
  primStateNormal(1, QU) = 0.0_rt;                  // un
  primStateNormal(1, QV) = primStateNormal(2, QV);  // ut1
  primStateNormal(1, QW) = primStateNormal(2, QW);  // ut2
  // zerograd temperature and pressure
  primStateNormal(1, QPRES) = primStateNormal(2, QPRES);
  primStateNormal(1, QT) = primStateNormal(2, QT);
  // ensure thermodynamic consistency
  primStateNormal(1, QRHO) =
      primStateNormal(1, QPRES) / (primStateNormal(1, QT) * closures.Rspec);

  // extrapolate GP
  extrapolateGP(primStateNormal, disGP, disIM);

  // limiting p and T
  // primStateNormal(0,QPRES) = max(primStateNormal(0,QPRES),1.0);
  // primStateNormal(0,QT)    = max(primStateNormal(0,QT),50.0);

  // thermodynamic consistency
  primStateNormal(0, QRHO) =
      primStateNormal(0, QPRES) / (primStateNormal(0, QT) * closures.Rspec);

  // Convert un,ut1,ut2 to ux,uy,uz at GP point only, we don't need to transform
  // at other points
  local2global(0, primStateNormal, norm, tan1, tan2);

  // insert primitive variables into primsFab
  int i = gp_ijk(0);
  int j = gp_ijk(1);
  int k = gp_ijk(2);
  for (int kk = 0; kk < NPRIM; kk++) {
    primFab(i, j, k, kk) = primStateNormal(0, kk);
  }

  // primFab(i,j,k,jj)
  // AMREX_ASSERT_WITH_MESSAGE( primFab(i,j,k,QPRES)>50,"P<50 at GP");

  // insert conservative ghost state into consFab
  conFab(i, j, k, URHO) = primStateNormal(0, QRHO);
  conFab(i, j, k, UMX) = primStateNormal(0, QRHO) * primStateNormal(0, QU);
  conFab(i, j, k, UMY) = primStateNormal(0, QRHO) * primStateNormal(0, QV);
  conFab(i, j, k, UMZ) = primStateNormal(0, QRHO) * primStateNormal(0, QW);
  Real ek = 0.5_rt * (primStateNormal(0, QU) * primStateNormal(0, QU) +
                      primStateNormal(0, QV) * primStateNormal(0, QV) +
                      primStateNormal(0, QW) * primStateNormal(0, QW));
  conFab(i, j, k, UET) = primStateNormal(0, QPRES) / (closures.gamma - 1.0_rt) +
                         primStateNormal(0, QRHO) * ek;
}

void IB::computeGPs(const Array4<const Real>& cons, const Array4<Real>& prims,
                    const IBFab& ibFab, const PROB::ProbClosures& cls) {
  // for GP data
  auto const gp_ijk = ibFab.gpData.gp_ijk.data();
  auto const imp_ijk = ibFab.gpData.imp_ijk.data();
  auto const ipweights = ibFab.gpData.imp_ipweights.data();
  auto const idxCube = ibFab.gpData.indexCube.data();
  auto const disGP = ibFab.gpData.disGP.data();
  auto const disIM = ibFab.gpData.disIM.data();
  auto const norm = ibFab.gpData.normal.data();
  auto const tan1 = ibFab.gpData.tangent1.data();
  auto const tan2 = ibFab.gpData.tangent2.data();

  // if (WM) or can use templates
  //  ComputeGPStateWM()
  // else
  amrex::ParallelFor(ibFab.gpData.ngps, [=] AMREX_GPU_DEVICE(int ii) {
    ComputeGPState(ii, gp_ijk[ii], imp_ijk[ii], ipweights[ii], norm[ii],
                   tan1[ii], tan2[ii], disGP[ii], disIM[ii], ibFabArr, cons,
                   prims, idxCube, cls);
  });
}

void IB::compute_plane_equations(Polyhedron::Facet& f) {
  Polyhedron::Halfedge_handle h = f.halfedge();
  f.plane() =
      Polyhedron::Plane_3(h->opposite()->vertex()->point(),
                          h->vertex()->point(), h->next()->vertex()->point());
};

void IB::readGeom() {
  ParmParse pp;
  Vector<std::string> Vfilename;
  pp.getarr("ib.filename", Vfilename);

  IB::NGEOM = Vfilename.size();
  VGeom.resize(NGEOM);
  VtreePtr.resize(NGEOM);
  Vfnormals.resize(NGEOM);
  VInOutFunc.resize(NGEOM);

  namespace PMP = CGAL::Polygon_mesh_processing;
  Print() << "----------------------------------" << std::endl;
  for (int i = 0; i < NGEOM; i++) {
    Print() << "----------------------------------" << std::endl;
    if (!PMP::IO::read_polygon_mesh(Vfilename[i], IB::VGeom[i])) {
      std::cerr << "Invalid geometry filename" << std::endl;
      exit(1);
    }
    Print() << "Geometry (i=" << i << ") " << Vfilename[i] << " read"
            << std::endl;
    Print() << "Is geometry only made of triangles? "
            << VGeom[i].is_pure_triangle() << std::endl;
    Print() << "Number of facets " << VGeom[i].size_of_facets() << std::endl;

    // constructs AABB tree and computes internal KD-tree
    // data structure to accelerate distance queries
    VtreePtr[i] =
        new Tree(faces(VGeom[i]).first, faces(VGeom[i]).second, VGeom[i]);
    Print() << "AABB tree constructed" << std::endl;

    PMP::compute_face_normals(VGeom[i],
                              boost::make_assoc_property_map(Vfnormals[i]));
    Print() << "Face normals computed" << std::endl;

    // plane class also computes orthogonal direction to the face. However, the
    // orthogonal vector is not normalised.
    std::for_each(VGeom[i].facets_begin(), VGeom[i].facets_end(),
                  compute_plane_equations);
    Print() << "Plane equations per face computed" << std::endl;

    // make inside/outside function for each geometry
    for (int ii = 0; ii < IB::NGEOM; ii++) {
      VInOutFunc[ii] = new inside_t(IB::VGeom[ii]);
    }
    Print() << "In out testing function constructed" << std::endl;

    // create face to displacement map //
    // auto temp = boost::make_assoc_property_map(fdisplace);
    // for(face_descriptor f : faces(geom))
    // {
    //   Vector_CGAL vec;
    //   put(temp, f, vec);
    //   // std::cout << "face plane " << f->plane() << "\n";
    // }

    // create face to surfdata map //
    // auto map = boost::make_assoc_property_map(face2state);
    // for(face_descriptor f : faces(geom))
    // {
    //   surfdata data;
    //   put(map, f, data);
    //   // std::cout << "face plane" << f->plane() << "\n";
    // }
  }
  Print() << "----------------------------------" << std::endl;
  Print() << "----------------------------------" << std::endl;

  //  || CGAL::is_empty(IB::VGeom[i]) || !CGAL::is_triangle_mesh(IB::geom)
}

// void IB::moveGeom() {
//   // Displace verticies //
//   // For each vertex its position p is translated by the displacement vector
//   (di) for each ith face. Each vertex has nfaces, for a triangular closed
//   mesh this equals the number of edges at a vertex. This is called degree of
//   the vertex by CGGAL. for (Polyhedron::Facet_handle fh :
//   geom.facet_handles())
//   {
//     // Print() << "New face" << " \n";
//     Polyhedron::Halfedge_handle start = fh->halfedge(), h = start;
//     do {
//       int nfaces = h->vertex()->degree();
//       CGAL::Point_3 p = h->vertex()->point();
//       // std::cout << "Vertex degree = " << nfaces  << "\n";
//       // std::cout << "Vertex before = " << p << "\n";
//       face_descriptor f = fh->halfedge()->face();
//       Array<Real,AMREX_SPACEDIM> dis = face2state[f].displace;

//       CGAL::Vector_3<K2> di(dis[0]/nfaces,dis[1]/nfaces,dis[2]/nfaces);
//       CGAL::Aff_transformation_3<K2> translate(CGAL::TRANSLATION,di);
//       p = p.transform(translate);

//       // std::cout << "Vertex after = " << p << "\n";

//       h = h->next();
//     } while(h!=start);
//   }

//   // apply boundary conditions //
//   // using a map of boundary nodes?

//   // Misc geometry things //
//   // rebuild tree
//   treePtr->rebuild(faces(geom).first,faces(geom).second,geom);
//   Print() << "Tree rebuilt" << std::endl;

//   CGAL::Polygon_mesh_processing::compute_face_normals(geom,
//   boost::make_assoc_property_map(fnormals)); Print() << "Face normals
//   recomputed" << std::endl;

//   // plane class also computes orthogonal direction to the face. However, the
//   orthogonal vector is not normalised. std::for_each( geom.facets_begin(),
//   geom.facets_end(),compute_plane_equations); Print() << "Plane equations
//   recomputed" << std::endl;

//   // Geometry fair?
// }

void IB::computeSurf(int lev) {
  // for each level

  exit(0);
}