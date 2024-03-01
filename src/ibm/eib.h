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
/// \param NIMPS Number of image points (integer) 
///
template <int NIMPS>
struct gpData_t {
  gpData_t() {
    // index cube for image point interpolation
    indexCube.push_back(Array1D<int, 0, IDIM>{0, 0, 0});
    indexCube.push_back(Array1D<int, 0, IDIM>{0, 1, 0});
    indexCube.push_back(Array1D<int, 0, IDIM>{1, 1, 0});
    indexCube.push_back(Array1D<int, 0, IDIM>{1, 0, 0});
    indexCube.push_back(Array1D<int, 0, IDIM>{0, 0, 1});
    indexCube.push_back(Array1D<int, 0, IDIM>{0, 1, 1});
    indexCube.push_back(Array1D<int, 0, IDIM>{1, 1, 1});
    indexCube.push_back(Array1D<int, 0, IDIM>{1, 0, 1});
  }
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
  Gpu::ManagedVector<Array2D<Real, 0, NIMPS - 1, 0, IDIM>> imp_xyz;
  Gpu::ManagedVector<Array2D<int, 0, NIMPS - 1, 0, IDIM>> imp_ijk;
  // for imp1 [(i,j,k), (i+1,j,k), (i,j+1,k), (i,j,k+1),
  //  ... ]
  Gpu::ManagedVector<Array2D<Real, 0, NIMPS - 1, 0, 7>> imp_ipweights;
  // Anti-clockwise order k=0 plane first, then k=1 plane.
  Gpu::ManagedVector<Array1D<int, 0, IDIM>> indexCube;
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
  Real cim = cim_tparm::num / cim_tparm::den;

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
      di_a[i] =
        cim * sqrt(pow(dx_a[i][0], 2) + pow(dx_a[i][1], 2) + pow(dx_a[i][2], 2));
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

private:
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