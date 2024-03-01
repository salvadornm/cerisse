#ifndef IBM_H_
#define IBM_H_





using namespace amrex;

namespace IBM {

// surfdata holds the relevant surface data for each element.
// struct surfdata {
//   Array<Real, NPRIM> state = {0.0};
//   Array<Real, AMREX_SPACEDIM> displace = {0.0};
// };

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


// 
/// \brief IB is the main class. It holds an array of IBMultiFab, one for each AMR level; and it also holds the geometry
///
/// \param NIMPS Number of image points 
///
template <int NIMPS, typename cls_t>
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



  // methods //
  void init(Amr* pointer_amr, const int nghost);

  void static compute_plane_equations(Polyhedron::Facet& f);

  void buildMFs(const BoxArray& bxa, const DistributionMapping& dm, int lev);

  void destroyMFs(int lev);

  void readGeom();

  void computeMarkers(int lev);

  void initialiseGPs(int lev);

  void computeGPs(int lev, MultiFab& consmf, MultiFab& primsmf,
                  const cls_t& closures);

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

}  // namespace IBM
#endif