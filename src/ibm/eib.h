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
#include <CGAL/AABB_traits.h> // deprecated
//#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
// CGAL headers for AABB tree for surface data
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
// #include <CGAL/Unique_hash_map.h>
// #include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>
// CGAL header for inout testing
#include <CGAL/Side_of_triangle_mesh.h>

// #include <ib_walltypes.h>


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

///
/// \brief Class to store surface data
/// \param eorder_tparm Number of image points (integer) 
///
template <int eorder_tparm>
struct surfData_t{
  surfData_t() {}

  Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
  Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;
  Array2D<Real,0,eorder_tparm-1,0,7> ipweights;
  Array3D<int,0,eorder_tparm-1,0,7,0,AMREX_SPACEDIM-1> ip_ijk;
  Array1D<Real, 0, AMREX_SPACEDIM - 1> ib_xyz;
  Array1D<Real, 0, AMREX_SPACEDIM - 1> norm;

  int  ifab,core,lev,iface,igeom;
  amrex::Real area;
  amrex::Real o_dis,pressure, temperature, dTdn;
  bool pointfound;

};


/// 
template <typename FacetHandle>
double compute_face_area(FacetHandle fd) {
    //using Point = typename FacetHandle::value_type::Point_3;
    using Point = CGAL::Simple_cartesian<double>::Point_3;
    using FT = typename Point::FT;

    auto h = fd->facet_begin();
    std::vector<Point> vertices;

    auto end = h;
    do {
        vertices.push_back(h->vertex()->point());
        ++h;
    } while (h != end);

    // Use first vertex as fan origin
    const Point& origin = vertices[0];
    FT total_area = FT(0);

    for (size_t i = 1; i + 1 < vertices.size(); ++i) {
        const Point& p1 = vertices[i];
        const Point& p2 = vertices[i + 1];
        total_area += CGAL::sqrt(CGAL::squared_area(origin, p1, p2));
    }

    return CGAL::to_double(total_area);
}
//


// main class ------------------------------------------------------------------

///
/// \brief eib_t is explicit geometry (triangulation based) immersed boundary method
/// class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
/// the geometry
///
template <typename wallmodel, typename param, typename cls_t>
class eib_t
{
public:
  // constant factor for image point
  //Real cim = Real(cim_tparm::num) / cim_tparm::den;
  //Real cim = 0.5;

  static const int iorder_tparm = param::interp_order;
  static const int eorder_tparm = param::extrap_order;  
  static constexpr Real cim = param::alpha;

  // image point distance per level
  Vector<Real> di_a;

  // number of geometries
  int ngeom=1;

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

  // surface data (independent of number of geometries)
  //Vector<std::map<face_descriptor,  surfData_t<iorder_tparm> >> surfdata_a;
  Vector<surfData_t<iorder_tparm> > surfdata_a;
  // Faces stored in per fab
  Vector<Vector<Polyhedron::Facet_const_handle>> faces_in_fab;
  // faces integers per fab and per level
  Vector<Vector<int>> intfaces_in_fab;
  Vector<Vector<int>> intfaces_in_lev;
  int ntotalfaces=0; //across all geometries


  //----


  // face state data map 
  //std::map<face_descriptor,surfdata> face2state;
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

      //delete fnorm_a.at(ii);
      //delete surfdata_a.at(ii); 

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
      const auto& ibMarkers = mfab.array(mfi); // boolean array

      // compute sld markers (including at ghost points) - cannot use ParallelFor - CGAL call causes problems
      // const IntVect& lo = bx.smallEnd();
      // const IntVect& hi = bx.bigEnd();
      // for (int k = lo[2] - cls_t::NGHOST; k <= hi[2] + cls_t::NGHOST; ++k) {
      //   for (int j = lo[1] - cls_t::NGHOST; j <= hi[1] + cls_t::NGHOST; ++j) {
      //     for (int i = lo[0] - cls_t::NGHOST; i <= hi[0] + cls_t::NGHOST; ++i) {
      amrex::LoopOnCpu(amrex::grow(bx, cls_t::NGHOST), [&](int i, int j, int k) {
            // NOTE: ibMarkers are accessed on CPU here always. This relies on amrex.the_arena_is_managed=1 option in the inputs. 
            // TODO: remove this and transfer bool for all i for given j,k to GPU in async arrays.

            // initialise to false
            ibMarkers(i, j, k, 0) = false;
            ibMarkers(i, j, k, 1) = false;

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
                ibFab.gpData.geomIdx.push_back(ii);
                break;
              }
            }
        //   }
        // }
      });

      // compute ghost markers -- move to GPU.
      ibFab.gpData.ngps = 0;
      int nextra = 1;
      // for (int k = lo[2] - nextra; k <= hi[2] + nextra; ++k) {
      //   for (int j = lo[1] - nextra; j <= hi[1] + nextra; ++j) {
      //     for (int i = lo[0] - nextra; i <= hi[0] + nextra; ++i) {
      amrex::LoopOnCpu(amrex::grow(bx, nextra), [&](int i, int j, int k) {
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
        //   }
        // }
      });
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
    auto const ibMarkers = mfab.array(mfi);  // boolean array

    // we need a CPU loop here (cannot be GPU loop) as CGAL tree seach for
    // closest element to a point needs to be called. instead of looping through
    // previously indexed gps, we loop through the whole ghost point field as it
    // is available on GPU and CPU at all times. Unlike the gp indexes, which
    // are only stored on GPU memory. Array1D<int,0,AMREX_SPACEDIM-1>& idx =
    // ibFab.gpData.gp_ijk[ii];

    // const IntVect& lo = bxg.smallEnd();
    // const IntVect& hi = bxg.bigEnd();
    // for (int k = lo[2]; k <= hi[2]; ++k) {
    //   for (int j = lo[1]; j <= hi[1]; ++j) {
    //     for (int i = lo[0]; i <= hi[0]; ++i) {
    amrex::LoopOnCpu(bxg, [&](int i, int j, int k) {
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
            int igeom = 0 ;
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

            //ib_xyz
            Array1D<Real, 0, AMREX_SPACEDIM - 1> ib_xyz = {cp[0],cp[1],cp[2]};
            gpData.ib_xyz.push_back(ib_xyz); //SNM

            // IM points -------------------------------------------
            Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
            Array2D<int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;


             //    imp_xyz(0, kk) = cp[kk] ??? 

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

            // check locations physical coordinates
            // printf("GP point: %f %f %f\n", gp[0], gp[1], gp[2]);
            // printf("IB point: %f %f %f\n", cp[0], cp[1], cp[2]);
            // printf("Normal: %f %f %f\n", norm(0), norm(1), norm(2));
            // printf("disGP: %f\n", disGP);
            // printf("disIM: %f\n", di_a[lev]);
            // printf("IM point: %f %f %f\n", imp_xyz(0,0), imp_xyz(0,1), imp_xyz(0,2));
        }
    //   }
    // }
  });
}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief compute surface indexes and store them (core, fab, lev)
// is independet of geomenries, it will store faces 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_surface_index(int lev) {

  amrex::Print() << " Compute Surface Index at LEVEL " << lev << std::endl;

  printf(" dx_a[lev] = %f \n ",dx_a[lev][0]);

  auto& mfab = *bmf_a[lev];
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

  int faces_notfound=0;
  
  ntotalfaces=0;

  int iface = -1;    // face counter

  // clear surfdata??

  for (int ii = 0; ii < ngeom; ii++) {    
    const Polyhedron& mesh   = geom_a[ii];

    // loop over the faces of geometry
      for (auto fd : faces(mesh)) {
        iface++;ntotalfaces++;

        // create a surface point if first level (otherwise get the surfdata)

        surfData_t<iorder_tparm> surf_dat; 


        if (lev==0)
        {
          surf_dat.pointfound = false;
        }
        else 
        { 
          auto& surf_dat = surfdata_a[iface];
        }
  
        // extract face normal norm
        Array1D<Real, 0, AMREX_SPACEDIM - 1> norm = {
                fnorm_a[ii][fd][0], fnorm_a[ii][fd][1],
                fnorm_a[ii][fd][2]};

        // extract face center point
        Point face_center = CGAL::centroid(fd->halfedge()->vertex()->point(), 
                                         fd->halfedge()->next()->vertex()->point(),
                                         fd->halfedge()->next()->next()->vertex()->point());
                                         
        Array1D<Real, 0, AMREX_SPACEDIM - 1> ib_xyz = {face_center[0],face_center[1],face_center[2]};

        // is the face inside the domain?
        // find closest  point i,j,k  to the face_center
        //const int i1 = int((face_center.x() - amr_p->Geom(lev).ProbLo()[0]) / dx_a[lev][0] - 0.5_rt);
        //const int j1 = int((face_center.y() - amr_p->Geom(lev).ProbLo()[1]) / dx_a[lev][1] - 0.5_rt);
        //const int k1 = int((face_center.z() - amr_p->Geom(lev).ProbLo()[2]) / dx_a[lev][2] - 0.5_rt);

        const int i1 = int((face_center.x() - prob_lo[0]) / dx_a[lev][0] - 0.5_rt);
        const int j1 = int((face_center.y() - prob_lo[1]) / dx_a[lev][1] - 0.5_rt);
        const int k1 = int((face_center.z() - prob_lo[2]) / dx_a[lev][2] - 0.5_rt);
         
        const bool is_inside = amr_p->Geom(lev).Domain().contains(i1, j1, k1);

       /// printf(" ifac=%d closest point i1 j1 k1 = %d %d %d \n",iface,i1,j1,k1);
        
        if (is_inside)
        {

          // compute interpolation weights imp_ipweights
          Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
          Array2D<int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;

          // find mirror point and the bottom left point closest to the image
          Real o_dis = 0.0;
          const int jj = 0;
          for (int kk = 0; kk < AMREX_SPACEDIM; kk++) {
            imp_xyz(jj, kk) = face_center[kk] + Real(jj + 1) * di_a[lev] * norm(kk);
            imp_ijk(jj, kk) = floor((imp_xyz(jj, kk) - prob_lo[kk]) / dx_a[lev][kk] - 0.5_rt);        
            o_dis += (face_center[kk] - imp_xyz(jj, kk))*(face_center[kk] - imp_xyz(jj, kk));          
          }                  
          o_dis = 1.0/sqrt(o_dis);

        
          double area = compute_face_area(fd);
         
          
          // mirror point 
          int i =  imp_ijk(jj, 0);int j =  imp_ijk(jj, 1);int k = imp_ijk(jj, 2);
                            
          int nfab = mfab.local_size(); 
          faces_in_fab.resize(nfab);  // resize fab Vector to store faces

          intfaces_in_fab.resize(nfab); 

          // locate the mirror points, loop over fabs
          bool face_present_core = false; 
          for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {             
    
            const int ifab= mfi.index();
            auto& ibFab = mfab.get(mfi);                
            const Box& bx = mfi.tilebox();
            const Box& bxg = mfi.growntilebox(cls_t::NGHOST);

            // point located
            if (bx.contains(i, j, k)) {

              auto const ibMarkers = mfab.array(mfi);
              // Interpolation points' (ips) weights for each image point
              Array2D<Real,0,eorder_tparm-1,0,7> ipweights;
              Array3D<int,0,eorder_tparm-1,0,7,0,AMREX_SPACEDIM-1> ip_ijk;

              computeIPweights(ipweights,ip_ijk,imp_xyz, imp_ijk, prob_lo, dx_a[lev], ibMarkers);

              face_present_core = true; 
              // store ip_ijk, ipweights, ifab, core  <------
                                
              surf_dat.pointfound = true; 
              // 
              if (surf_dat.pointfound)
              {                
                surf_dat.ip_ijk= ip_ijk;
                surf_dat.ipweights= ipweights;              
                surf_dat.ifab = ifab;
                surf_dat.lev= lev;    
                surf_dat.iface = iface;
                surf_dat.igeom = ii;      
                surf_dat.o_dis   = o_dis;
                surf_dat.imp_xyz = imp_xyz;
                surf_dat.imp_ijk = imp_ijk;
                surf_dat.ib_xyz  = ib_xyz;
                surf_dat.norm    = norm;   
                surf_dat.area  = area;        
                // push into arrays
                faces_in_fab[ifab].push_back(fd);             // Store face for fab index
                intfaces_in_fab[ifab].push_back(iface);       // Store face number for fab index              
                surfdata_a.emplace_back(surf_dat);             // store surface data 
              }  
            
              
            }
          } // end loop over fabs

          if (!face_present_core) faces_notfound++;

        }
        else // face outside domain
        {
          // stoes 0
          faces_notfound++;
        }
      }  // end loop faces
  } // end loop geometries

  amrex::Print() << " Faces not found  (outside domain/lost) " << faces_notfound << " out of " << ntotalfaces << std::endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief compute surface properties for each surface face
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_surface_props(MultiFab& stateprops,const cls_t* cls,int lev) {
  

  amrex::Print() << " Compute Surface Properties at LEVEL" << lev << std::endl;

  // loop over mfi
  auto& mfab = *bmf_a[lev];
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();
  for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) { 
    
    const int ifab= mfi.index();
    auto& ibFab = mfab.get(mfi);                
    const Box& bx = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    auto const ibMarkers = mfab.array(mfi);  // boolean array

    // pointer to  U
    Array4<Real> const& cons = stateprops.array(mfi);
    // primitives array (store a local copy)
    FArrayBox primf(bxg, cls_t::NPRIM, The_Async_Arena());
    Array4<Real> const& prims= primf.array();
    // convert to Q (local copy)
    cls->cons2prims(mfi, cons, prims); 
    // ........................................
    for (int ii = 0; ii < ngeom; ii++) {    

      //
      const Polyhedron& mesh   = geom_a[ii];
            
      //for (auto fd : faces_in_fab[ifab]){
      for (int j = 0; j < intfaces_in_fab[ifab].size(); ++j){
       // const auto& fd = faces_in_fab[ifab][j];
        
        int iface =  intfaces_in_fab[ifab][j];

        auto& surf_dat = surfdata_a[iface];

        // extract arrays from surface data
        auto const norm      = surf_dat.norm;
        auto const ip_ijk    = surf_dat.ip_ijk;
        auto const ipweights = surf_dat.ipweights;
        auto const ib_xyz    = surf_dat.ib_xyz;     
        auto const o_dis     = surf_dat.o_dis;
                      
        Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1> primsNormal={0.0};                

        // calculate surface properties    
        interpolateIMs(ip_ijk,ipweights,prims,primsNormal);
    
        wallmodel::compute_surfIB(ib_xyz,norm,primsNormal,cls);   

        // compute one-sided gradients dT/dn du/dn
        Real dTdn = (primsNormal(2,cls_t::QT) - primsNormal(1,cls_t::QT))*o_dis;

        // store values
        surf_dat.pressure       = primsNormal(1,cls_t::QPRES); 
        surf_dat.temperature    = primsNormal(1,cls_t::QT); 
        surf_dat.dTdn           = dTdn;


        // if (isnan(surf_dat.pressure)  ) {
        // //if (surf_dat.temperature > 5000) {

        //   printf(" iface= %d ifab=%d\n",iface,ifab);
        //   printf(" TEMP prims1= %f prims2= %f \n",primsNormal(1,cls_t::QT),primsNormal(2,cls_t::QT));
        //   printf(" PRES prims1= %f prims2= %f \n",primsNormal(1,cls_t::QPRES),primsNormal(2,cls_t::QPRES));
        //   printf(" QT=%d QPRES=%d \n",cls_t::QT,cls_t::QPRES);
        //   printf(" norm= %f %f %f \n",norm(0),norm(1),norm(2));


        //   for (int iip=0;iip<8;iip++) {
        //     int i1 = ip_ijk(0,iip,0);
        //     int j1 = ip_ijk(0,iip,1);
        //     int k1 = ip_ijk(0,iip,2);
           
        //     printf(" iip=%d i=%d j=%d k=%d  \n", iip,i1,j1,k1);
        //     printf(" weights= %f\n ibMarker0 =%d ",ipweights(0,iip),ibMarkers(i1,j1,k1,0));
            
        //     printf(" --------------------- PRIMS \n");
        //     for (int nv=0;nv<cls_t::NPRIM;nv++){
        //       printf(" nv= %d Q =%f \n",nv,prims(i1,j1,k1,nv));
        //     }
        //     printf(" --------------------- CONS \n");
        //     for (int nv=0;nv<cls_t::NCONS;nv++){
        //       printf(" nv= %d U =%f \n",nv,cons(i1,j1,k1,nv));
        //     }
        //     printf(" ------------------------- \n");

        //   }

        //   amrex::Abort();
        // } 

    
      } //end loop faces          
    } //end loop geometry 
    //...........................................
  } // end looping mfi
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief computeGPs computes ghost points (GPs) for each fab
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  // surface coordinates
  auto const ib_xyz =  ibFab.gpData.ib_xyz.data();

  // snm surface coordinates (per GP) (no only relevant IM)
  //auto const imp_xyz = ibFab.gpData.imp_xyz.data();

  // create a copy of prims to use local prims0 (needed?)
  const Box& bxg = mfi.growntilebox(cls->NGHOST);
  FArrayBox primf(bxg, cls_t::NPRIM, The_Async_Arena());
  Array4<Real> const& prims0= primf.array();
  ParallelFor(bxg, cls_t::NPRIM, [=] 
      AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      prims0(i,j,k,n) = prims(i,j,k,n);    
    });
  //
  

  ParallelFor(ibFab.gpData.ngps, [=,copy=this] AMREX_GPU_DEVICE (int ii)
  {

// #if TEST_IB_INTERPOLATION
    // printf("GP %d -------\n",ii);
    // printf("GP ijk = (%d, %d, %d)\n", gp_ijk[ii](0), gp_ijk[ii](1), gp_ijk[ii](2) );
    // printf("norm: (%f, %f, %f)\n", norm[ii](0), norm[ii](1), norm[ii](2) );

    // // for each image point
    // for (int iim=0; iim<eorder_tparm; iim++) {
    //   // for each IP point
    //   printf("IM %d -------\n",iim);
    //   printf("IM ijk = (%d, %d, %d)\n", imp_ijk[ii](iim,0), imp_ijk[ii](iim,1), imp_ijk[ii](iim,2) );
    //   for (int iip=0; iip<8; iip++) {
    //     int iii = imp_ip_ijk[ii](iim,iip,0);
    //     int jjj = imp_ip_ijk[ii](iim,iip,1);
    //     int kkk = imp_ip_ijk[ii](iim,iip,2);
    //     printf("IP %d, ijk=(%d, %d, %d), weight= %f \n",iip,iii,jjj,kkk,imp_ipweights[ii](iim,iip));
    //   }
    //   printf("------- \n");
    // }

// #endif

      Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1> primsNormal={0.0};

      // check put prims0 instead of prims
      copy->interpolateIMs(imp_ip_ijk[ii],imp_ipweights[ii],prims0,primsNormal);

      // transform velocity to local coordinates for image points only >2
      for (int iip=2; iip<2+eorder_tparm; iip++) {
        copy->global2local(iip, primsNormal, norm[ii], tan1[ii], tan2[ii]);
      }
     
      // copy->computeIB(primsNormal,cls);                      // only P,T,Y set 

      // compute surface values based on wallmodel  (set u,P,T,Y)      
      wallmodel::compute_surfIB(ib_xyz[ii],norm[ii],primsNormal,cls);   
          
      // only extrapolate to GP the values of  u,P,T,Y
      copy->extrapolate(primsNormal, disGP[ii], disIM[ii]);  
      
      // transform velocity back to global coordinates (GP only)
      int idx=0;
      copy->local2global(idx,primsNormal,norm[ii],tan1[ii],tan2[ii]);      

      ///.  copy primsNormal  -> Q
      Real P,T,Y[NUM_SPECIES]={0.0};
      P = primsNormal(0,cls_t::QPRES);
      T = primsNormal(0,cls_t::QT);
#if NUM_SPECIES > 1    
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n]   =  primsNormal(0,cls_t::QFS+n);
      }
#endif 
      Real ux =  primsNormal(0,cls_t::QU);
      Real uy =  primsNormal(0,cls_t::QV);
      Real uz =  primsNormal(0,cls_t::QW);
            
      // ensure Thermodynamic consistency and that the prims array is filled      
      Real Q[cls_t::NPRIM];
      cls->ensurePTYfillq(P, T, Y, ux,uy,uz,Q); 
            
      // insert primitive variables into primsFab
      int i=gp_ijk[ii](0); int j=gp_ijk[ii](1); int k = gp_ijk[ii](2);
      for (int nn=0; nn<cls_t::NPRIM; nn++) {
        prims(i,j,k,nn) = Q[nn];       
      }      
    });
};


  ////////////////////////////////////////////////////////////////
  //  \brief plot surface mesh to file
  //  \param igeom geometry index
  //  \param filename output file name
  //  \note uses CGAL Polygon_mesh_processing IO functions
  void plot_surface(const int igeom, const std::string& filename) {
    const Polyhedron& mesh   = geom_a[igeom];
    const auto& face_normals = fnorm_a[igeom];

    std::ofstream out(filename);
    out << "# vtk DataFile Version 3.0\n";
    out << "CGAL Polyhedron\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    // Step 1: Write vertices
    std::map<Polyhedron::Vertex_const_handle, int> vmap;
    int vidx = 0;
    out << "POINTS " << mesh.size_of_vertices() << " float\n";
    for (auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
        const auto& p = vit->point();
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
        vmap[vit] = vidx++;
    }

    // Step 2: Write faces
    int num_faces = mesh.size_of_facets();
    out << "POLYGONS " << num_faces << " " << num_faces * 4 << "\n";
    for (auto fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit) {
        auto h = fit->halfedge();
        out << "3 "
            << vmap[h->vertex()]
            << " " << vmap[h->next()->vertex()]
            << " " << vmap[h->next()->next()->vertex()] << "\n";
    }

    out << "CELL_DATA " << num_faces << "\n";
    
    // Step 3: Write surfdata as CELL_DATA
    // out << "VECTORS face_normals float\n";
    // for (auto fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit) {
    //   auto it = face_normals.find(*fit);
    //   if (it != face_normals.end()) {
    //     const auto& n = it->second;
    //     out << n.x() << " " << n.y() << " " << n.z() << "\n";
    //   } else {
    //     out << "0.0 0.0 0.0\n";  // or some fallback value
    //   }
    // }

    // SCALARS

    out << "SCALARS face_area float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int iface=0;iface<num_faces;iface++) {     
      out << static_cast<float>(surfdata_a[iface].area) << "\n";
    }

    out << "SCALARS pressure float 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<num_faces;iface++) {             
      out << static_cast<float>(surfdata_a[iface].pressure) << "\n";
    }
    out << "SCALARS temperature float 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<ntotalfaces;iface++) {      
      out << static_cast<float>(surfdata_a[iface].temperature) << "\n";
    }
    out << "SCALARS gradT float 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<ntotalfaces;iface++) {      
      out << static_cast<float>(surfdata_a[iface].dTdn) << "\n";
    }

    out << "SCALARS calculated float 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<num_faces;iface++) { 
      if (surfdata_a[iface].pointfound) {
        out << 1.0 << "\n";
      }
      else
      {
        out << 0.0 << "\n";
      }

    }

    out.close();
    amrex::Print() << "----------------------------------\n";
    amrex::Print() << "Surface mesh plotted to " << filename << "\n";
    amrex::Print() << "----------------------------------\n";
  }
  //////////////////////////////////////////////////////////////////

private:
    // Taylor expansion around IB point (only up to QLS)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void extrapolate(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& stateNormal, Real dgp, Real dim) {
      Real sgn_dgp = -dgp ; // negative sign as taylor expansion is around IB point, IM and GP are in opposite directions
      for (int kk=0; kk<= cls_t::QLS; kk++) {
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

      // AMREX_ASSERT_WITH_MESSAGE( sumfluid >= 2,"Less than 2 interpolation points are fluid points");
      if (sumfluid < 2) {
        amrex::Print() << i << " " << j << " " << k << " " << sumfluid << std::endl;
        amrex::Warning("Less than 2 interpolation points are fluid points");
      }

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

  // @brief Change coordinates of velocity
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
  void global2local( int iip, Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& primsNormal, 
                    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm, const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1, const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2) {

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

  ////////////////////////////////////////////////////////////////
  /// \brief computeIB calculates the primitive array at IB (surface)
  /// based on values interpolated on the normal 
  /// WARN !! at present does first order
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void computeIB(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& primsNormal, const cls_t* cls) {
    
    // slip velocity (in local coordinates)
    primsNormal(1,cls_t::QU) = 0.0_rt; // un
    primsNormal(1,cls_t::QV) = primsNormal(2,cls_t::QV); // ut1
    primsNormal(1,cls_t::QW) = primsNormal(2,cls_t::QW); // ut2


    Real Yw[NUM_SPECIES]={0.0};

    // zerograd pressure and T (adiabatic)   
    primsNormal(1,cls_t::QPRES) = primsNormal(2,cls_t::QPRES); 
    primsNormal(1,cls_t::QT)    = primsNormal(2,cls_t::QT);

#if NUM_SPECIES > 1    
    Real sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yw[n]   =  primsNormal(2,cls_t::QFS+n);
      sumY += sumY;
    }
    for (int n = 0; n < NUM_SPECIES; ++n) { 
      primsNormal(1,cls_t::QFS+n)   =  Yw[n]/sumY;
    }
#endif                      

  }
  ////////////////////////////////////////////////////////////////
  /// \brief reads STL greometry from input
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

      // create face to surfdata map //  SNM
      // auto map = boost::make_assoc_property_map(face2state);
      //  for(face_descriptor f : faces(geom))
      //  {
      //    surfdata data;
      //    put(map, f, data);
      //    // std::cout << "face plane" << f->plane() << "\n";
      //  }
      // SNM

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