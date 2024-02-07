#ifndef IBM_H_
#define IBM_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_Derive.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_StateDescriptor.H>
#include <prob.h>

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

using namespace amrex;

namespace IBM {

// CGAL types -------------------------------------------
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
// CGAL types -------------------------------------------

// surfdata holds the relevant surface data for each element.
struct surfdata {
  Array<Real, NPRIM> state = {0.0};
  Array<Real, AMREX_SPACEDIM> displace = {0.0};
};

// GP holds the relevant data for each ghost point.
// struct gp {

//     // attributes
//     Array<int , AMREX_SPACEDIM> idx={0};
//     Point_and_primitive_id closest; // closest surface point (ib point) and
//     face Real disGP; Vector<Array<Real,AMREX_SPACEDIM>> imps;
//     Vector<Array<int,AMREX_SPACEDIM>> impscell;
//     Vector<Array<Real,8>> ipweights;
//     // for imp1 [(i,j,k), (i+1,j,k), (i,j+1,k), (i,j,k+1),
//     //  ... ]
// };

struct gpData_t {
  gpData_t() {
    // index cube for image point interpolation
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{0, 0, 0});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{0, 1, 0});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{1, 1, 0});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{1, 0, 0});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{0, 0, 1});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{0, 1, 1});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{1, 1, 1});
    indexCube.push_back(Array1D<int, 0, AMREX_SPACEDIM - 1>{1, 0, 1});
  }

  // CPU only attributes
  int ngps;
  Vector<Point_and_primitive_id>
      closest_cgal;  // closest surface point (ib point) and face

  // GPU/CPU attributes
  Gpu::ManagedVector<Array1D<int, 0, AMREX_SPACEDIM - 1>> gp_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, AMREX_SPACEDIM - 1>> normal, tangent1,
      tangent2, ib_xyz;
  Gpu::ManagedVector<Real> disGP, disIM;
  Gpu::ManagedVector<int> geomIdx;

  Gpu::ManagedVector<Array2D<Real, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1>>
      imp_xyz;
  Gpu::ManagedVector<Array2D<int, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1>> imp_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, NIMPS - 1, 0, 7>> imp_ipweights;
  // for imp1 [(i,j,k), (i+1,j,k), (i,j+1,k), (i,j,k+1),
  //  ... ]

  // Anti-clockwise order k=0 plane first, then k=1 plane.
  Gpu::ManagedVector<Array1D<int, 0, AMREX_SPACEDIM - 1>> indexCube;
  // = { {0,0,0}, {0,1,0}, {1,1,0}, {1,0,0}, {0,0,1}, {0,1,1}, {1,1,1},
  // {1,0,1}}; // index cube for image point interpolation

  // void resize (int n) {
  //   gp_ijk.resize(n); ib_xyz.resize(n);
  //   disGP.resize(n); disIM.resize(n);
  //   imp_xyz.resize(n); imp_ijk.resize(n); imp_ipweights.resize(n);
  // }
};

// IBFab holds the solid point and ghost point boolean arrays
class IBFab : public BaseFab<bool> {
 public:
  // attributes //
  // Vector<gp> gpArray; // problem with this approach is that Vector type is
  // not available on the GPUs.

  gpData_t gpData;  // a single structure which holds GP data in vectors which
                    // can be resized on CPU or GPU.

  // AsyncArray<gp>* gpArrayAsync(1); // designed to be used for temporary
  // transfers, in routines, not for storage.

  // int ngps = 0;
  ///////////////

  // constructors and destructors //
  // using Box
  explicit IBFab(const Box& b, int ncomp, bool alloc = true,
                 bool shared = false, Arena* ar = nullptr);
  // using IBFab
  explicit IBFab(const IBFab& rhs, MakeType make_type, int scomp, int ncomp);

  // IBFab (IBFab&& rhs) noexcept = default;
  // IBFab& operator= (IBFab&&) noexcept = default;

  // IBFab (const IBM::IBFab&) {};
  // IBFab& operator= (const IBFab&) = delete;

  ~IBFab();
  ///////////////////////////////////////////////////////////////////////////

  // methods //
  /////////////////////////////
};

// IBMultiFab holds an array of IBFab on a level
class IBMultiFab : public FabArray<IBFab> {
 public:
  // constructor from BoxArray and DistributionMapping
  explicit IBMultiFab(
      const BoxArray& bxs, const DistributionMapping& dm, const int nvar,
      const int ngrow, const MFInfo& info = MFInfo(),
      const FabFactory<IBFab>& factory = DefaultFabFactory<IBFab>());

  // move constructor
  IBMultiFab(IBMultiFab&& rhs) noexcept;

  // destructor
  ~IBMultiFab();

  void copytoRealMF(MultiFab& mf, int ibcomp, int mfcomp);

 protected:
};

// IB is the main class. It holds an array of IBMultiFab, one for each AMR
// level; and it also holds the geometry
class IB {
 public:
  // attributes //
  Vector<Real> disIM;  // image point distance per level

  int MAX_LEVEL =
      0;  // need to save this as Amr protects this -- set in ib.init(...)
  int NGHOST_IB =
      0;  // number of ghost points in IBMultiFab -- set in ib.init(...)
  const int NVAR_IB = 2;  // number of variables in IBMultiFab
  // const int NIMPS = 2; // number of image points per ghost point
  int NGEOM = 1;  // number of geometries

  // pointer to Amr class instance
  Amr* pamr;

  // vector of refinement ratio per level in each direction
  Vector<IntVect> ref_ratio;

  // vector of cell sizes per level in each direction
  Vector<GpuArray<Real, AMREX_SPACEDIM>> cellSizes;

  // Immersed boundary MultiFab array
  Vector<IBMultiFab*> ibMFa;

  // Level set MultiFab array
  Vector<MultiFab> lsMFa;

  // IB explicit geometry
  // Polyhedron geom;
  Vector<Polyhedron> VGeom;

  // AABB tree
  // Tree* treePtr;
  Vector<Tree*> VtreePtr;

  // Instead of std::map you may use std::unordered_map, boost::unordered_map
  // or CGAL::Unique_hash_map
  // CGAL::Unique_hash_map<face_descriptor,Vector> fnormals;
  // boost::unordered_map<vertex_descriptor,Vector> vnormals;
  // face element normals
  // std::map<face_descriptor,Vector_CGAL> fnormals;
  Vector<std::map<face_descriptor, Vector_CGAL>> Vfnormals;

  // face state data map
  // std::map<face_descriptor,surfdata> face2state;
  // std::map stores information in a binary tree, it has log(N) complexity for
  // key-value pair insertion and value retrieval for a key. could also fit
  // normals into this -- however might need to modify compute_normals routine

  // face displacement map
  // std::map<face_descriptor,Vector_CGAL> fdisplace;
  // Vector or Real[3] don't work

  // in out testing function
  Vector<inside_t*> VInOutFunc;

  explicit IB();
  ~IB();

  // methods //
  void init(Amr* pointer_amr, const int nghost);

  void static compute_plane_equations(Polyhedron::Facet& f);

  void buildMFs(const BoxArray& bxa, const DistributionMapping& dm, int lev);

  void destroyMFs(int lev);

  void readGeom();

  void computeMarkers(int lev);

  void initialiseGPs(int lev);

  void computeGPs(int lev, MultiFab& consmf, MultiFab& primsmf,
                  const PROB::ProbClosures& closures);

  void moveGeom();

  void computeSurf(int lev);

 private:
  void computeIPweights(
      Array2D<Real, 0, NIMPS - 1, 0, 7>& weights,
      Array2D<Real, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
      Array2D<int, 0, NIMPS - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      GpuArray<Real, AMREX_SPACEDIM>& dxyz, const Array4<bool> markersFab,
      auto const idxCube);

  // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  // void extrapolateGP(Array2D<Real,0,NIMPS+1,0,NPRIM-1>& primStateNormal, Real
  // dgp, Real dim);

  // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void ComputeGPState(long ii,const
  // IBM::IBFab& ibFab);
};

// declare main IB class instance
inline IB ib;

}  // namespace IBM
#endif