#ifndef eib_H_
#define eib_H_

#include <ratio>

#include <IBMultiFab.h>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_Derive.H>
// basic CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polyhedron_3.h>
// CGAL headers for AABB tree for closest point
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
// CGAL headers for AABB tree for surface data
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
// #include <CGAL/Unique_hash_map.h>
// #include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>
// CGAL header for inout testing
#include <CGAL/Side_of_triangle_mesh.h>


// CGAL types ------------------------------------------------------------------
typedef CGAL::Simple_cartesian<Real> K2;
typedef K2::Point_3 Point;
typedef CGAL::Polyhedron_3<K2> Polyhedron;

typedef K2::FT FT;
typedef K2::Segment_3 Segment;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K2, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

typedef K2::Vector_3 Vector_CGAL;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, K2> inside_t;

// typedef Polyhedron::Vertex_iterator      vertexIter;
// typedef Polyhedron::Face_iterator        faceIter;
// CGAL types ------------------------------------------------------------------

///
/// \brief Class to store ghost point arrays
/// \param eorder_tparm Number of image points (integer) 
///
template <int eorder_tparm>
struct gpData_t {
  gpData_t() {}

  // CPU only attributes
  static constexpr int IDIM=AMREX_SPACEDIM - 1;
  int ngps;
  // closest surface point (ib point) and face
  Vector<Point_and_primitive_id> closest_cgal;  

  // GPU/CPU attributes
  Gpu::ManagedVector<Array1D<int,  0, IDIM>> gp_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, IDIM>> normal, tangent1,
      tangent2, ib_xyz;
  Gpu::ManagedVector<Real> disGP, disIM;
  Gpu::ManagedVector<int> geomIdx;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>> imp_xyz;
  Gpu::ManagedVector<Array2D<int, 0, eorder_tparm - 1, 0, IDIM>> imp_ijk;
  // for imp1 [(i,j,k), (i+1,j,k), (i,j+1,k), (i,j,k+1),
  //  ... ]
  Gpu::ManagedVector<Array3D<int,0,eorder_tparm-1,0,7,0,AMREX_SPACEDIM-1>> imp_ip_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, 7>> imp_ipweights;
};

// main class ------------------------------------------------------------------

///
/// \brief eib_t is explicit geometry (triangulation based) immersed boundary method
/// class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
/// the geometry
///
/// \param iorder integer interpolation order
/// \param eorder integer extrapolation order
/// \param cip real interpolation distance factor (relative to mesh diagonal)
/// \param cls_t Number of image points (type)
///
template <int iorder_tparm, int eorder_tparm, typename cim_tparm, typename cls_t>
class eib_t
{
public:
  // constant factor for image point
  Real cim = Real(cim_tparm::num) / cim_tparm::den;

  // image point distance per level
  Vector<Real> di_a;

  // number of geometries
  int ngeom;

  // pointer to Amr class instance
  Amr* amr_p;

  // vector of refinement ratio per level in each direction
  Vector<IntVect> rratio_a;

  // vector of cell sizes per level in each direction
  Vector<GpuArray<Real, AMREX_SPACEDIM>> dx_a; // dx array

  // Immersed boundary MultiFab array
  Vector<IBMultiFab<bool, gpData_t<iorder_tparm>>*> bmf_a; //(bool multifab array)

  // Level set MultiFab array -- not needed for now
  // Vector<MultiFab> rmf_a; //(real multifab array)

  // IB explicit geometry
  // Polyhedron geom;
  Vector<Polyhedron> geom_a;

  // AABB tree
  // Tree* treePtr;
  Vector<Tree*> tree_pa; //

  // Instead of std::map you may use std::unordered_map, boost::unordered_map
  // or CGAL::Unique_hash_map
  // CGAL::Unique_hash_map<face_descriptor,Vector> fnormals;
  // boost::unordered_map<vertex_descriptor,Vector> vnormals;
  // face element normals
  // std::map<face_descriptor,Vector_CGAL> fnormals;
  Vector<std::map<face_descriptor, Vector_CGAL>> fnorm_a;

  // face state data map
  // std::map<face_descriptor,surfdata> face2state;
  // std::map stores information in a binary tree, it has log(N) complexity for
  // key-value pair insertion and value retrieval for a key. could also fit
  // normals into this -- however might need to modify compute_normals routine

  // face displacement map
  // std::map<face_descriptor,Vector_CGAL> fdisplace;
  // Vector or Real[3] don't work

  // in out testing function
  Vector<inside_t*> inout_fa; // function array

  // methods which are to be called from the solver class
  ~eib_t()
  {
    // clear memory
    for (int ii = 0; ii < ngeom; ii++) {
      delete tree_pa.at(ii);
      delete inout_fa.at(ii);
    };
  }

  // initialise IB
  void init(Amr* pointer_amr)
  {
    amr_p = pointer_amr;
    rratio_a = amr_p->refRatio();
    int lmax = amr_p->maxLevel();
    bmf_a.resize(lmax + 1);
    // rmf_a.resize(lmax + 1);

    dx_a.resize(lmax + 1);
    dx_a[0] = amr_p->Geom(0).CellSizeArray();
    for (int i = 1; i <= lmax; i++) {
      for (int j = 0; j < AMREX_SPACEDIM; j++) {
        dx_a[i][j] = dx_a[i - 1][j] / rratio_a[i - 1][j];
      }
    }

    di_a.resize(lmax + 1);
    for (int i = 0; i <= lmax; i++) {
      di_a[i] = cim*sqrt(pow(dx_a[i][0], 2) + pow(dx_a[i][1], 2) + pow(dx_a[i][2], 2));
    }

    read_geom();
  };

  // create IBMultiFabs at a level and store pointers to it
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    bmf_a[lev] =
      new IBMultiFab<bool, gpData_t<eorder_tparm>>(bxa, dm, 2, cls_t::NGHOST);
    // lsMFa[lev].define(bxa, dm, 1, NGHOST_IB);
  }

  void destroy_mf(int lev)
  {
    if (!bmf_a.empty()) { delete bmf_a.at(lev); }
    // if (!lsMFa.empty()) { lsMFa[lev].clear(); }
  }

  void computeMarkers(int lev)
  {
    auto& mfab = *bmf_a[lev];
    // assuming same number of ghost points in all directions
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      auto& ibFab = mfab.get(mfi);
      const Box& bx = mfi.tilebox();
      const IntVect& lo = bx.smallEnd();
      const IntVect& hi = bx.bigEnd();
      const auto& ibMarkers = mfab.array(mfi); // boolean array

      // compute sld markers (including at ghost points) - cannot use ParallelFor
      // - CGAL call causes problems
      for (int k = lo[2] - cls_t::NGHOST; k <= hi[2] + cls_t::NGHOST; ++k) {
        for (int j = lo[1] - cls_t::NGHOST; j <= hi[1] + cls_t::NGHOST; ++j) {
          for (int i = lo[0] - cls_t::NGHOST; i <= hi[0] + cls_t::NGHOST; ++i) {
            ibMarkers(i, j, k, 0) = false; // initialise to false
            ibMarkers(i, j, k, 1) = false; // initialise to false

            Real x = prob_lo[0] + (0.5_rt + Real(i)) * dx_a[lev][0];
            Real y = prob_lo[1] + (0.5_rt + Real(j)) * dx_a[lev][1];
            Real z = prob_lo[2] + (0.5_rt + Real(k)) * dx_a[lev][2];
            Point gridpoint(x, y, z);

            for (int ii = 0; ii < eib_t::ngeom; ii++) {
              inside_t& inside = *inout_fa[ii];
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

void initialiseGPs(int lev) {
  auto& mfab = *bmf_a[lev];
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

  for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
    auto& ibFab = mfab.get(mfi);
    auto& gpData = ibFab.gpData;
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    // const Box& bx = mfi.tilebox();
    const IntVect& lo = bxg.smallEnd();
    const IntVect& hi = bxg.bigEnd();
    auto const ibMarkers = mfab.array(mfi);  // boolean array

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
            Real x = prob_lo[0] + (0.5_rt + i) * dx_a[lev][0];
            Real y = prob_lo[1] + (0.5_rt + j) * dx_a[lev][1];
            Real z = prob_lo[2] + (0.5_rt + k) * dx_a[lev][2];
            Point gp(x, y, z);

            // find and store geometery index for this GP. This index is used
            // for searching appropriate tree in initialiseGPs to find the
            // closest element and searching all trees.

            // in out test for each geometry
            int igeom;
            for (int ii = 0; ii < eib_t::ngeom; ii++) {
              inside_t& inside = *inout_fa[ii];
              CGAL::Bounded_side result = inside(gp);

              if (int(result) == int(CGAL::ON_BOUNDED_SIDE)) {
                igeom = ii;
                ibFab.gpData.geomIdx.push_back(igeom);
                break;
              }
              // TODO: assert geometries do not overlap?
            }

            // closest surface point and face --------------------------
            Point_and_primitive_id closest_elem =
                tree_pa[igeom]->closest_point_and_primitive(gp);

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
            // Print() << "Normal " << fnorm_a[igeom][face] <<std::endl;
            // Print() << "cp-gp " << cp - gp << std::endl; // should be in the
            // direction of normal Print() << "Plane " << face->plane().a() << "
            // " << face->plane().b() << " "  << face->plane().c() << " " <<
            // face->plane().d() << std::endl;

            // IB point -------------------------------------------
            Vector_CGAL imp_gp(gp, cp);
            Real disGP = sqrt(CGAL::squared_distance(gp, cp));
            AMREX_ASSERT_WITH_MESSAGE(
                disGP < 1.0 * sqrt(dx_a[lev][0] * dx_a[lev][0] +
                                   dx_a[lev][1] * dx_a[lev][1] +
                                   dx_a[lev][2] * dx_a[lev][2]),
                "Ghost point and IB point distance larger than mesh diagonal");

            //*store*
            gpData.disGP.push_back(disGP);
            Array1D<Real, 0, AMREX_SPACEDIM - 1> norm = {
                fnorm_a[igeom][face][0], fnorm_a[igeom][face][1],
                fnorm_a[igeom][face][2]};

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
            Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
            Array2D<int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;

            // find image point and the bottom left point closest to the image
            // point
            //  In 2D, same idea in 3D.
            //     i,j+1 (2) ---------------------     i+1,j+1 (3)
            //     |                                  |
            //     |         P                        |
            //     |                                  |
            //     |                                  |
            //     i,j  (1) ----------------------      i+1,j  (4)
            for (int jj = 0; jj <= eorder_tparm-1; jj++) {
              for (int kk = 0; kk < AMREX_SPACEDIM; kk++) {
                imp_xyz(jj, kk) = cp[kk] + Real(jj + 1) * di_a[lev] *
                                               fnorm_a[igeom][face][kk];
                imp_ijk(jj, kk) =
                    floor((imp_xyz(jj, kk) - prob_lo[kk]) / dx_a[lev][kk] -
                          0.5_rt);
              }

              AMREX_ASSERT_WITH_MESSAGE(
                  bxg.contains(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2)),
                  "Interpolation point outside fab");
            }

            // store
            gpData.disIM.push_back(di_a[lev]);
            gpData.imp_xyz.push_back(imp_xyz);
            gpData.imp_ijk.push_back(imp_ijk);

            // Interpolation points' (ips) weights for each image point
            Array2D<Real,0,eorder_tparm-1,0,7> ipweights;
            Array3D<int,0,eorder_tparm-1,0,7,0,AMREX_SPACEDIM-1> ip_ijk;
            computeIPweights(ipweights,ip_ijk,imp_xyz, imp_ijk, prob_lo, dx_a[lev], ibMarkers);
            // *store*
            gpData.imp_ipweights.push_back(ipweights);
            gpData.imp_ip_ijk.push_back(ip_ijk);

        }
      }
    }
  }
}
}

// called per fab
void computeGPs(const MFIter& mfi, const Array4<Real>& cons, const Array4<Real>& prims, const cls_t* cls, int& lev) {
  
  auto& mfab = *bmf_a[lev];
  const auto& ibFab = mfab.get(mfi);

  // for GP data
  auto const gp_ijk = ibFab.gpData.gp_ijk.data();
  auto const imp_ijk= ibFab.gpData.imp_ijk.data();
  auto const imp_ipweights = ibFab.gpData.imp_ipweights.data();
  auto const imp_ip_ijk = ibFab.gpData.imp_ip_ijk.data();
  // TODO: combine disGP and disIM into 2D array
  auto const disGP = ibFab.gpData.disGP.data();
  auto const disIM = ibFab.gpData.disIM.data();
  // TODO: combine norm,tan1,tan2 into matrix (2d array)
  auto const norm = ibFab.gpData.normal.data();
  auto const tan1 = ibFab.gpData.tangent1.data();
  auto const tan2 = ibFab.gpData.tangent2.data();

    ParallelFor(ibFab.gpData.ngps, [=,copy=this] AMREX_GPU_DEVICE (int ii)
    {

      Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1> primsNormal={0.0};
      copy->interpolateIMs(imp_ip_ijk[ii],imp_ipweights[ii],prims,primsNormal);

      // transform velocity to local coordinates for image points only
      for (int iip=2; iip<2+eorder_tparm; iip++) {
        copy->global2local(iip, primsNormal, norm[ii], tan1[ii], tan2[ii]);
      }
      copy->computeIB(primsNormal,cls);
      copy->extrapolate(primsNormal, disGP[ii], disIM[ii]);
      // thermodynamic consistency
      primsNormal(0,cls_t::QRHO)  = primsNormal(0,cls_t::QPRES)/(primsNormal(0, cls_t::QT)*cls->Rspec);

      // only transform velocity back to global coordinates for gp only
      int idx=0;
      copy->local2global(idx,primsNormal,norm[ii],tan1[ii],tan2[ii]);

      // limiting p and T
      // primsNormal(0,QPRES) = max(primsNormal(0,QPRES),1.0);
      // primsNormal(0,QT)    = max(primsNormal(0,QT),50.0);

      // insert primitive variables into primsFab
      int i=gp_ijk[ii](0); int j=gp_ijk[ii](1); int k = gp_ijk[ii](2);
      for (int nn=0; nn<cls_t::NPRIM; nn++) {
        prims(i,j,k,nn) = primsNormal(0,nn);
      }

      // AMREX_ASSERT_WITH_MESSAGE( prims(i,j,k,QPRES)>50,"P<50 at GP");

      // insert conservative ghost state into consFab
      cons(i,j,k,cls_t::URHO) = primsNormal(0,cls_t::QRHO);
      cons(i,j,k,cls_t::UMX)  = primsNormal(0,cls_t::QRHO)*primsNormal(0,cls_t::QU);
      cons(i,j,k,cls_t::UMY)  = primsNormal(0,cls_t::QRHO)*primsNormal(0,cls_t::QV);
      cons(i,j,k,cls_t::UMZ)  = primsNormal(0,cls_t::QRHO)*primsNormal(0,cls_t::QW);
      Real ek   = 0.5_rt*(primsNormal(0,cls_t::QU)*primsNormal(0,cls_t::QU) + primsNormal(0,cls_t::QV)* primsNormal(0,cls_t::QV) + primsNormal(0,cls_t::QW)*primsNormal(0,cls_t::QW));
      cons(i,j,k,cls_t::UET) = primsNormal(0,cls_t::QPRES)/(cls->gamma-1.0_rt) + primsNormal(0,cls_t::QRHO)*ek;
      });
};


private:
    // Taylor expansion around IB point
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void extrapolate(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& stateNormal, Real dgp, Real dim) {
      Real sgn_dgp = -dgp ; // negative sign as taylor expansion is around IB point, IM and GP are in opposite directions
      for (int kk=0; kk<cls_t::NPRIM; kk++) {
        // Linear
        Real c1 = stateNormal(1,kk);
        Real c2 = (stateNormal(2,kk) - stateNormal(1,kk))/dim;
        stateNormal(0,kk) = c1 + c2*sgn_dgp;
      }
    }

    // note runs on CPU
    void computeIPweights(Array2D<Real,0,eorder_tparm-1,0,7>&weights,Array3D<int,0,eorder_tparm-1,0,7,0,AMREX_SPACEDIM-1>&ip_ijk, Array2D<Real,0,eorder_tparm-1,0,AMREX_SPACEDIM-1>&imp_xyz, Array2D<int,0,eorder_tparm-1,0,AMREX_SPACEDIM-1>& imp_ijk, const GpuArray<Real,AMREX_SPACEDIM>& prob_lo, GpuArray<Real,AMREX_SPACEDIM>& dxyz, const Array4<bool>& ibFab) {
    // Array2D<int,0,7,0,AMREX_SPACEDIM-1> indexCube={{0,0,0},{0,1,0},{1,1,0},{1,0,0},{0,0,1},{0,1,1},{1,1,1},{1,0,1}}; // index cube for image point interpolation
    // Anti-clockwise order k=0 plane first, then k=1 plane.
    for (int iim=0; iim<=eorder_tparm-1; iim++) {
      int i = imp_ijk(iim,0); int j = imp_ijk(iim,1); int k = imp_ijk(iim,2); 
      // note xl,xr, ...etc do not have prob_lo added to them. This does not matter as we only need the relative distances between the points.
      Real xl = prob_lo[0] + (Real(i)+0.5_rt) * dxyz[0];  // bottom left corner of cell
      Real xr = xl + dxyz[0];
      Real yl = prob_lo[1] + (Real(j)+0.5_rt) * dxyz[1];
      Real yr = yl + dxyz[1];
      Real zl = prob_lo[2] + (Real(k)+0.5_rt) * dxyz[2];
      Real zr = zl + dxyz[2];

      Real xd =  (imp_xyz(iim,0) - xl )/(xr-xl);
      Real yd =  (imp_xyz(iim,1) - yl )/(yr-yl);
      Real zd =  (imp_xyz(iim,2) - zl )/(zr-zl);

      int sumfluid = 0;
      Real sumweights = 0.0_rt;
      // zd = 0
      int ii = i;
      int jj = j;
      int kk = k;
      int iip= 0;
      int fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = (1.0_rt - xd) *(1.0_rt - yd)*(1.0_rt-zd)*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,0);

      ii = i + 0;
      jj = j + 1;
      kk = k + 0;
      iip= 1;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = (1.0_rt - xd) *yd*(1.0_rt-zd)*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,iip);

      ii = i + 1;
      jj = j + 1;
      kk = k + 0;
      iip= 2;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = xd*yd*(1.0_rt-zd)*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,iip);

      ii = i + 1;
      jj = j + 0;
      kk = k + 0;
      iip= 3;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = xd*(1.0_rt - yd)*(1.0_rt-zd)*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,iip);

      // zd = 2
      ii = i + 0;
      jj = j + 0;
      kk = k + 1;
      iip= 4;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = (1.0_rt - xd) *(1.0_rt - yd)*zd*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,iip);

      ii = i + 0;
      jj = j + 1;
      kk = k + 1;
      iip= 5;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = (1.0_rt - xd) *yd*zd*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,5);

      ii = i + 1;
      jj = j + 1;
      kk = k + 1;
      iip= 6;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = xd*yd*zd*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,iip);

      ii = i + 1;
      jj = j + 0;
      kk = k + 1;
      iip= 7;
      fluid = !ibFab(ii,jj,kk, 0);
      weights(iim,iip) = xd*(1.0_rt - yd)*zd*fluid;
      ip_ijk(iim,iip,0) = ii; ip_ijk(iim,iip,1) = jj; ip_ijk(iim,iip,2) = kk;
      sumfluid += fluid; sumweights += weights(iim,iip);

      AMREX_ASSERT_WITH_MESSAGE( sumfluid >= 2,"Less than 2 interpolation points are fluid points");

      // re-normalise
      for (int ll=0; ll<8; ll++) {
        weights(iim,ll) = weights(iim,ll)/sumweights;
      }

      AMREX_ASSERT_WITH_MESSAGE(std::abs(weights(iim,0) + weights(iim,1) + weights(iim,2) + weights(iim,3) + weights(iim,4) + weights(iim,5) + weights(iim,6) + weights(iim,7) - Real(1.0)) < Real(1.e-9),"Interpolation point weights do not sum to 1.0");
    }
  }

  // general interpolation routine -- over a given stencil points and weights
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
  void interpolateIMs(const Array3D<int,0,eorder_tparm-1,0,7,0,AMREX_SPACEDIM-1>& imp_ip_ijk, const Array2D<Real,0,eorder_tparm-1,0,7>& ipweights, const Array4<Real>& prims, Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& primsNormal){
    // for each image point
    for (int iim=0; iim<eorder_tparm; iim++) {
      // for each IP point
      for (int iip=0; iip<8; iip++) {
        int ii = imp_ip_ijk(iim,iip,0);
        int jj = imp_ip_ijk(iim,iip,1);
        int kk = imp_ip_ijk(iim,iip,2);
        // for each primitive variable
        for (int nn=0; nn<cls_t::NPRIM; nn++) {
          // GP (jj=0),IB (jj=1), IM1 (jj=2),IM2 (jj=3)...
          primsNormal(iim+2,nn) += prims(ii,jj,kk,nn)*ipweights(iim,iip);
        }
      }
    }
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
  void global2local( int iip, Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& primsNormal, const Array1D<Real,0,AMREX_SPACEDIM-1>& norm, const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1, const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2) {

    Array1D<Real,0,AMREX_SPACEDIM-1> vel;
    vel(0) = primsNormal(iip,cls_t::QU); vel(1) = primsNormal(iip,cls_t::QV); vel(2) = primsNormal(iip,cls_t::QW);

    primsNormal(iip,cls_t::QU) = vel(0)*norm(0) + vel(1)*norm(1) + vel(2)*norm(2);
    primsNormal(iip,cls_t::QV) = vel(0)*tan1(0) + vel(1)*tan1(1) + vel(2)*tan1(2);
    primsNormal(iip,cls_t::QW) = vel(0)*tan2(0) + vel(1)*tan2(1) + vel(2)*tan2(2);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
  void local2global (int jj, Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& primsNormal, const Array1D<Real,0,AMREX_SPACEDIM-1>& norm, const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1, const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2) {
    
    Array1D<Real,0,AMREX_SPACEDIM-1> vel;
    vel(0) = primsNormal(jj,cls_t::QU); vel(1) = primsNormal(jj,cls_t::QV); vel(2) = primsNormal(jj,cls_t::QW);
    primsNormal(jj,cls_t::QU) = vel(0)*norm(0) + vel(1)*tan1(0) + vel(2)*tan2(0);
    primsNormal(jj,cls_t::QV) = vel(0)*norm(1) + vel(1)*tan1(1) + vel(2)*tan2(1);
    primsNormal(jj,cls_t::QW) = vel(0)*norm(2) + vel(1)*tan1(2) + vel(2)*tan2(2);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void computeIB(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& primsNormal, const cls_t* cls) {
    // no-slip/slip velocity (in local coordinates)
    primsNormal(1,cls_t::QU) = 0.0_rt; // un
    primsNormal(1,cls_t::QV) = primsNormal(2,cls_t::QV); // ut1
    primsNormal(1,cls_t::QW) = primsNormal(2,cls_t::QW); // ut2
    // zerograd temperature and pressure
    primsNormal(1,cls_t::QPRES) = 2.0_rt*(2.0_rt*primsNormal(2,cls_t::QPRES) - 0.5_rt*primsNormal(3,cls_t::QPRES))/3.0_rt;
    primsNormal(1,cls_t::QT)    = primsNormal(2,cls_t::QT);
    // ensure thermodynamic consistency
    primsNormal(1,cls_t::QRHO)  = primsNormal(1,cls_t::QPRES)/(primsNormal(1,cls_t::QT)*cls->Rspec); 
    }

  void read_geom()
  {
    ParmParse pp;
    Vector<std::string> files_a;
    pp.getarr("ib.filename", files_a);

    eib_t::ngeom = files_a.size();
    eib_t::geom_a.resize(ngeom);
    eib_t::tree_pa.resize(ngeom);
    eib_t::fnorm_a.resize(ngeom);
    eib_t::inout_fa.resize(ngeom);

    namespace PMP = CGAL::Polygon_mesh_processing;
    Print() << "----------------------------------" << std::endl;
    for (int i = 0; i < ngeom; i++) {
      Print() << "----------------------------------" << std::endl;
      if (!PMP::IO::read_polygon_mesh(files_a[i], geom_a[i])) {
        std::cerr << "Invalid geometry filename" << std::endl;
        exit(0);
      }
      Print() << "Geometry (i=" << i << ") " << files_a[i] << " read" << std::endl;
      Print() << "Is geometry only made of triangles? "
              << geom_a[i].is_pure_triangle() << std::endl;
      Print() << "Number of facets " << geom_a[i].size_of_facets() << std::endl;

      // constructs AABB tree and computes internal KD-tree
      // data structure to accelerate distance queries
      tree_pa[i] =
        new Tree(faces(geom_a[i]).first, faces(geom_a[i]).second, geom_a[i]);
      Print() << "AABB tree constructed" << std::endl;

      PMP::compute_face_normals(geom_a[i],
                                boost::make_assoc_property_map(fnorm_a[i]));
      Print() << "Face normals computed" << std::endl;

      // plane class also computes orthogonal direction to the face. However, the
      // orthogonal vector is not normalised.
      std::for_each(geom_a[i].facets_begin(), geom_a[i].facets_end(),
                    eib_t::compute_plane_equations);
      Print() << "Plane equations per face computed" << std::endl;

      // make inside/outside function for each geometry
      for (int ii = 0; ii < ngeom; ii++) { inout_fa[ii] = new inside_t(geom_a[ii]); }
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

    //  || CGAL::is_empty(IB::geom_a[i]) || !CGAL::is_triangle_mesh(IB::geom)
  }

  void static compute_plane_equations(Polyhedron::Facet& f)
  {
    Polyhedron::Halfedge_handle h = f.halfedge();
    f.plane() =
      Polyhedron::Plane_3(h->opposite()->vertex()->point(), h->vertex()->point(),
                          h->next()->vertex()->point());
  };

};
#endif