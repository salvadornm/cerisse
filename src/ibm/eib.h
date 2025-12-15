#ifndef eib_H_
#define eib_H_

#include <ratio>
#include <limits>

#include <IBMultiFab.h>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_Derive.H>

#include <eib_cgal.h>

// ----------------------------------------------------------------------------
// Optional memory shrink helper
// Define EIB_ENABLE_SHRINK (e.g. add -DEIB_ENABLE_SHRINK to CXXFLAGS) to
// enable conditional capacity shrinking for vectors that may experience large
// transient peaks (ghost/image points, per-fab face indices). This helps
// reduce long-term resident memory without adding overhead on every step.
// A shrink is attempted only if capacity > size * EIB_SHRINK_RATIO.
// Override EIB_SHRINK_RATIO (default 4) at compile time to tune aggressiveness.

#ifdef EIB_ENABLE_SHRINK
#ifndef EIB_SHRINK_RATIO
#define EIB_SHRINK_RATIO 4
#endif
template <class Vec>
inline void eib_shrink_if_excess(Vec& v) {
  if (v.capacity() > v.size() * static_cast<size_t>(EIB_SHRINK_RATIO)) {
    Vec tmp(v); // copy into right-sized temporary then swap
    tmp.swap(v);
  }
}
#define EIB_SHRINK(v) eib_shrink_if_excess(v)
#else
#define EIB_SHRINK(v) ((void)0)
#endif
// ----------------------------------------------------------------------------

// index dimension
static constexpr int IDIM = AMREX_SPACEDIM - 1;     

// Box extra width for ghost point search, no more than cls_t::NGHOST - 1
static constexpr int GP_BOX_EXTRA = 0;

// Minimum number of valid fluid points required in the interpolation stencil
//
// For the first image point, if the number of available interpolation points is
// less than INTERP_THRESHOLD, an error will be raised and the corresponding IB
// information will be reported.
//
// For all subsequent image points, if the number of interpolation points is
// insufficient, the image point will be discarded and the order of the
// extrapolation will be reduced accordingly.
// In 2D (4 points total), we require 2.
// In 3D (8 points total), we require 3.
#if (AMREX_SPACEDIM == 2)
static constexpr int INTERP_THRESHOLD = 2;
#else
static constexpr int INTERP_THRESHOLD = 3;
#endif

// Number of attempts for the first image point placement
static constexpr int N_ATTEMPTS = 3;  
// Factor to extend the first image point placement along normal
static constexpr Real IMP_FACTOR[N_ATTEMPTS] = {1.0, 1.5, 2.0};  

//----------------------------------------------------------------------------
// \brief Class to store ghost point arrays
// \param eorder_tparm Number of image points (integer) 
//
template <int eorder_tparm, int iorder_tparm>
struct gpData_t {
  // CPU only attributes
  gpData_t() : ngps(0) {}
  int ngps;    

  // ideal number of interpolation points for each image point
  static constexpr int  N_InterP = IntPow<iorder_tparm + 1, AMREX_SPACEDIM>::value;

  // closest surface point (ib point) and face ID
  //Vector<Point_and_primitive_id> closest_cgal;  

  // GPU/CPU attributes
  // Ghost point data
  Gpu::ManagedVector<Array1D< int, 0, IDIM>> gp_ijk;                        // Ghost point indices
  Gpu::ManagedVector<Array1D<Real, 0, IDIM>> ib_xyz;                        // IB point coordinates
  Gpu::ManagedVector<Real> disGP;                                           // Distance from IB point to ghost point
    
  // Surface identification
  Gpu::ManagedVector<int> geomIdx;                                          // Geometry index
  Gpu::ManagedVector<int> elemIdx;                                          // face/edge element index

  // Image point data arrays
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>> imp_xyz;  // Image point coordinates
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>> imp_ijk;  // Image point indices
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm - 1>> disIM;             // Distance from IB point to image points
      
  // Interpolation data for image points
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm - 1>> imp_ninterp;       // Actual number of interpolation points used for each image point
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>> imp_ipweights;

  // Helper function to resize all vectors
  void resize(int n) {
      ngps = n;
      
      gp_ijk.resize(n);
      ib_xyz.resize(n);
      disGP.resize(n);
      
      geomIdx.resize(n);
      elemIdx.resize(n);

      imp_xyz.resize(n);
      imp_ijk.resize(n);
      disIM.resize(n);
      
      imp_ninterp.resize(n);
      imp_ip_ijk.resize(n);
      imp_ipweights.resize(n);
  }

  // Clear and free memory
  void clear() {
      ngps = 0;
      
      gp_ijk.clear();      gp_ijk.shrink_to_fit();
      ib_xyz.clear();      ib_xyz.shrink_to_fit();
      disGP.clear();       disGP.shrink_to_fit();
      
      geomIdx.clear();     geomIdx.shrink_to_fit();
      elemIdx.clear();     elemIdx.shrink_to_fit();

      imp_xyz.clear();     imp_xyz.shrink_to_fit();
      imp_ijk.clear();     imp_ijk.shrink_to_fit();
      disIM.clear();       disIM.shrink_to_fit();
      
      imp_ninterp.clear(); imp_ninterp.shrink_to_fit();
      imp_ip_ijk.clear();  imp_ip_ijk.shrink_to_fit();
      imp_ipweights.clear(); imp_ipweights.shrink_to_fit();
  }

};

//----------------------------------------------------------------------------
// \brief Class to store surface data (SoA structure)
// Per-face surface data container used for reconstruction and output
// \param eorder_tparm_surf Number of image points used for surfdata reconstruction
//
template <int eorder_tparm_surf, int iorder_tparm_surf>
struct surfData_t{
  // CPU only attributes
  surfData_t() : nfaces(0) {}
  int nfaces;

  // ideal number of interpolation points for each image point
  static constexpr int  N_InterP = IntPow<iorder_tparm_surf + 1, AMREX_SPACEDIM>::value;

  // Surface identification
  Gpu::ManagedVector<int> elemIdx;                                                                    // Global face index across all geometries
  Gpu::ManagedVector<int> geomIdx;                                                                    // Geometry index
  
  // Indexing info
  Gpu::ManagedVector<int> ifab;      // Local FAB index on this MPI rank
  Gpu::ManagedVector<int> lev;       // AMR level
  Gpu::ManagedVector<int> rank;      // Owning MPI rank
  Gpu::ManagedVector<int> pointfound;// Whether this face has been located (int for GPU compatibility)

  // Surface fields (per face)
  Gpu::ManagedVector<Real> pressure;         // reconstructed local surface pressure
  Gpu::ManagedVector<Real> tau1;             // reconstructed local surface shear stress 1
  Gpu::ManagedVector<Real> tau2;             // reconstructed local surface shear stress 2
  Gpu::ManagedVector<Real> temperature;      // reconstructed temperature
  Gpu::ManagedVector<Real> dTdn;             // reconstructed grad(T)·n

  // Image point data (per face)
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1, 0, IDIM>> imp_xyz;                       // Physical-space coordinates of image points placed along the outward normal
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm_surf - 1, 0, IDIM>> imp_ijk;                       // Index of the “bottom-left” grid cell associated with each image point
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm_surf - 1>> disIM;                                  // Normal distances from the surface (IB point) to each image point
  
  // Interpolation data for image points (per face)
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm_surf - 1>> imp_ninterp;                            // Actual number of interpolation points used for each image point
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm_surf - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;   // Indices of the 8-point interpolation stencil for each image point
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1, 0, N_InterP - 1>> imp_ipweights;         // Trilinear interpolation weights for the 8-point stencil of each image point
  
  // Helper function to resize all vectors
  // Note: resize() initializes new elements to 0 / default constructor.
  // If you need specific default values (e.g. -1 for indices), set them manually after resize.
  void resize(int n) {
      int old_n = nfaces;
      nfaces = n;
      
      elemIdx.resize(n);
      geomIdx.resize(n);
      ifab.resize(n);
      lev.resize(n);
      rank.resize(n);
      pointfound.resize(n);
      
      pressure.resize(n);
      tau1.resize(n);
      tau2.resize(n);
      temperature.resize(n);
      dTdn.resize(n);

      imp_xyz.resize(n);
      imp_ijk.resize(n);
      disIM.resize(n);
      imp_ninterp.resize(n);
      imp_ip_ijk.resize(n);
      imp_ipweights.resize(n);

      // Initialize new elements with specific defaults if n > old_n
      if (n > old_n) {
          // Use parallel for or standard fill for initialization
          // Here we use simple loops for safety, assuming this runs on CPU during setup
          for (int i = old_n; i < n; ++i) {
              ifab[i] = -1;
              lev[i]  = -1;
              rank[i] = -99;
              pointfound[i] = 0; // false
          }
      }
  }

  // Clear and free memory
  void clear() {
      nfaces = 0;
      
      // clear() only sets size to 0 but keeps capacity.
      // shrink_to_fit() forces memory deallocation.
      elemIdx.clear(); elemIdx.shrink_to_fit();
      geomIdx.clear(); geomIdx.shrink_to_fit();
      ifab.clear();    ifab.shrink_to_fit();
      lev.clear();     lev.shrink_to_fit();
      rank.clear();    rank.shrink_to_fit();
      pointfound.clear(); pointfound.shrink_to_fit();
      
      pressure.clear(); pressure.shrink_to_fit();
      tau1.clear();     tau1.shrink_to_fit();
      tau2.clear();     tau2.shrink_to_fit();
      temperature.clear(); temperature.shrink_to_fit();
      dTdn.clear();     dTdn.shrink_to_fit();

      imp_xyz.clear();  imp_xyz.shrink_to_fit();
      imp_ijk.clear();  imp_ijk.shrink_to_fit();
      disIM.clear();    disIM.shrink_to_fit();
      imp_ninterp.clear(); imp_ninterp.shrink_to_fit();
      imp_ip_ijk.clear();  imp_ip_ijk.shrink_to_fit();
      imp_ipweights.clear(); imp_ipweights.shrink_to_fit();
  }

};

//============================================================================
///---------------------------- main class -----------------------------------

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
  static constexpr int  iorder_tparm = param::interp_order; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm = param::extrap_order; // number of image points used for ghost point extrapolation
  static constexpr Real cim = param::alpha;
  
  static constexpr int  iorder_tparm_surf = param::interp_order_surf; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm_surf = param::extrap_order_surf; // number of image points used for surfdata reconstruction
  static constexpr Real cim_surf  = param::alpha_surf;

  // ideal number of interpolation points for each image point(ghost point extrapolation and surface reconstruction)
  static constexpr int  N_InterP      = IntPow<iorder_tparm + 1, AMREX_SPACEDIM>::value;
  static constexpr int  N_InterP_surf = IntPow<iorder_tparm_surf + 1, AMREX_SPACEDIM>::value;

  // Enumeration used by functions check_interpolation_stencil to determine behavior when a check fails.
  enum class CheckMode {
      Silent,      // Do not output anything, just return status
      Warn,        // Output a warning message, return status
      Abort        // Abort execution immediately on failure
  };

  using GPDATA = gpData_t<eorder_tparm, iorder_tparm>;
  using SURFDATA = surfData_t<eorder_tparm_surf, iorder_tparm_surf>;

  // MultiFabs pointer to Amr class instance
  Amr* amr_p;                                             
  // Immersed boundary MultiFab array (uint8_t multifab array)
  Vector<IBMultiFab<uint8_t, GPDATA>*> bmf_a;  

  // parameters for cell size and refinement ratio
  Vector<IntVect> rratio_a;                                 // vector of refinement ratio per level in each direction
  Vector<GpuArray<Real, AMREX_SPACEDIM>> dx_a;              // vector of cell sizes per level in each direction
  Vector<Real> diag_a;                                      // vector of cell diagonal length per level
  Vector<Real> di_a;                                        // image point distance per level

  // geometry related data
  int ngeom = 0;                                            // number of geometries
  Vector<GeomType> geom_a;                                  // IB explicit geometry
  Vector<Tree*> tree_pa;                                    // AABB tree per geometry
  Gpu::ManagedVector<LocalFrame> LocalFrame_a;              // local orthonormal frame matrix (flattened)
  Gpu::ManagedVector<SurfElem> SurfElem_a;                  // surface element area and coordinates (flattened)
  Gpu::ManagedVector<int> geom_offsets;                     // Start index for each geometry in flattened arrays
  Vector<std::map<PrimitiveID, int>> IdxMap_a;              // face/edge element index per geometry
  Vector<inside_t*> inout_fa;                               // in out testing function per geometry

  // surface related data
  int ntotalfaces = 0;                                      // number of faces/edges across all geometries
  Gpu::ManagedVector<SURFDATA> surfdata_a;                  // surface/edge data (independent of number of geometries)
  Vector<LevelFaceCSR> faces_per_level;                     // faces/edges integers per fab and per level 

  // Destructor to release allocated memory
  ~eib_t() noexcept
  {
    // Release per-level IBMultiFab pointers if any remain
    for (auto*& p : bmf_a) {
      if (p) { delete p; p = nullptr; }
    }
    bmf_a.clear();

    // Release CGAL AABB trees
    for (auto*& t : tree_pa) {
      if (t) { delete t; t = nullptr; }
    }
    tree_pa.clear();

    // Release inside/outside testers
    for (auto*& f : inout_fa) {
        if (f) { delete f; f = nullptr; }
    }
    inout_fa.clear();

    // Clear containers holding geometry and associated data
    geom_a.clear();
    LocalFrame_a.clear();
    SurfElem_a.clear();
    geom_offsets.clear();
    IdxMap_a.clear();

    surfdata_a.clear();
    faces_per_level.clear();  
  }

  /**
   * \brief Initializes the Immersed Boundary (IB) method structures and geometry.
   *
   * This function sets up the AMR pointer, resizes internal data structures based on the maximum AMR level,
   * computes grid metrics (cell sizes, diagonals) for all levels, and loads the IB geometry from files.
   *
   * \param pointer_amr Pointer to the main Amr class instance.
   */
  void init(Amr* pointer_amr)
  {
    amr_p = pointer_amr;
    rratio_a = amr_p->refRatio();
    int lmax = amr_p->maxLevel();
    bmf_a.resize(lmax + 1);

    dx_a.resize(lmax + 1);
    dx_a[0] = amr_p->Geom(0).CellSizeArray();
    for (int i = 1; i <= lmax; i++) {
      for (int j = 0; j < AMREX_SPACEDIM; j++) {
        dx_a[i][j] = dx_a[i - 1][j] / rratio_a[i - 1][j];
      }
    }

    di_a.resize(lmax + 1);
    diag_a.resize(lmax + 1);
    for (int i = 0; i <= lmax; i++) {
      diag_a[i] = std::sqrt(
        AMREX_D_TERM( std::pow(dx_a[i][0], 2),
                      + std::pow(dx_a[i][1], 2),
                      + std::pow(dx_a[i][2], 2)) );
      di_a[i] = cim * diag_a[i];
    }

    // read geometry from file
    read_geom();
  };

  // create IBMultiFab at a level and store pointers to it
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    bmf_a[lev] =
      new IBMultiFab<uint8_t, GPDATA>(bxa, dm, 2, cls_t::NGHOST);
      // lsMFa[lev].define(bxa, dm, 1, NGHOST_IB);
  }

  // destroy IBMultiFab at a level
  void destroy_mf(int lev)
  {
    if (lev < static_cast<int>(bmf_a.size())) { 
      delete bmf_a[lev]; 
      bmf_a[lev] = nullptr;  // Prevent dangling pointer
    }
  }

  /**
   * \brief Computes the solid/fluid markers and identifies ghost points for the Immersed Boundary Method.
   *
   * This function iterates over the grid points per level to determine if they are inside (solid) or outside (fluid)
   * the immersed boundary geometry using CGAL. It populates the ibMarkers where:
   * - Component 0 indicates if a point is solid (true) or fluid (false).
   * - Component 1 indicates if a solid point is a ghost point (neighbor to a fluid point).
   *
   * \param lev The current AMR level.
   */
  void computeMarkers(int lev)
  {
    auto& mfab = *bmf_a[lev];
    // assuming same number of ghost points in all directions
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      auto& ibFab = mfab.get(mfi);
      const Box& bx = mfi.tilebox();
      const auto& ibMarkers = mfab.array(mfi); // boolean array

      // Ensure gpData containers do not accumulate across rebuilds
      // (post_regrid/post_restart may be called multiple times over the run).
      // Clear per-fab ghost-point data before recomputing markers.
      ibFab.gpData.clear();

      // compute sld markers - cannot use ParallelFor - CGAL call causes problems
      amrex::LoopOnCpu(amrex::grow(bx, cls_t::NGHOST), [&](int i, int j, int k) {
        // NOTE: ibMarkers are accessed on CPU here always. This relies on amrex.the_arena_is_managed=1 option in the inputs. 
        // TODO: remove this and transfer bool for all i for given j,k to GPU in async arrays.

        // initialise to false
        ibMarkers(i, j, k, 0) = static_cast<uint8_t>(0);
        ibMarkers(i, j, k, 1) = static_cast<uint8_t>(0);

        Real x = prob_lo[0] + (0.5_rt + Real(i)) * dx_a[lev][0];
        Real y = prob_lo[1] + (0.5_rt + Real(j)) * dx_a[lev][1];
#if (AMREX_SPACEDIM == 2)
        Point gridpoint(x, y);
#else
        Real z = prob_lo[2] + (0.5_rt + Real(k)) * dx_a[lev][2];
        Point gridpoint(x, y, z);
#endif

        for (int ii = 0; ii < ngeom; ii++) {
          
          inside_t& inside = *inout_fa[ii];
          CGAL::Bounded_side result = inside(gridpoint);
          IB_WarnOnBoundary(ii, lev, i, j, k, result, gridpoint);

          // if point inside any IB geometry, mark as solid, move on to
          // another point. This minimises the number of inout testing
          // (expensive) calls.
          if (int(result) == int(CGAL::ON_BOUNDED_SIDE)) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      }); // end LoopOnCpu for sld markers

      // compute ghost markers
      // Note: Although this loop does not depend on CGAL, we intentionally keep it on the CPU 
      // to avoid the memory transfers between Host and Device.
      // 1. The previous step (solid marker computation) runs on CPU due to CGAL dependencies.
      // 2. If moved this step (Ghost identification) to the GPU...
      // 3. The subsequent step (initialiseGPs) must run on CPU again for CGAL closest-point searches.
      // TODO: move to GPU only if CGAL dependency is removed or ported.
      ibFab.gpData.ngps = 0;
      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA), [&](int i, int j, int k) {
        bool ghost = false;
        if (ibMarkers(i, j, k, 0)) {
          // check neighbors in x directions
          ghost = ghost || (!ibMarkers(i-1, j,   k,   0));
          ghost = ghost || (!ibMarkers(i+1, j,   k,   0));
          // check neighbors in y directions
          ghost = ghost || (!ibMarkers(i,   j-1, k,   0));
          ghost = ghost || (!ibMarkers(i,   j+1, k,   0));
#if (AMREX_SPACEDIM == 3)
          // check neighbors in z directions
          ghost = ghost || (!ibMarkers(i,   j,   k-1, 0));
          ghost = ghost || (!ibMarkers(i,   j,   k+1, 0));
#endif
          ibFab.gpData.ngps += ghost;

          if (ghost) {
            // store GP index
            //ibFab.gpData.gp_ijk.push_back(make_vec<int>(i, j, k));
            ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
          } else {
            ibMarkers(i, j, k, 1) = static_cast<uint8_t>(0);
          }
        } // end if solid
      }); // end LoopOnCpu for ghost markers

  // After constructing the ghost point lists for this FAB, optionally
  // reclaim excess capacity from previous high-water allocations.
  // Since we are using clear() at the beginning of the loop (implied or explicit),
  // shrink_to_fit() inside clear() handles memory reclamation.
  // The EIB_SHRINK macro is no longer needed here.
    } // end MFIter
  } // end computeMarkers

  /**
   * \brief Initialises geometric and interpolation data for Ghost Points (GPs).
   *
   * This function iterates over all ghost points identified in the `computeMarkers` step.
   * For each ghost point, it performs the following operations:
   * 1. Identifies the specific geometry (body) the ghost point belongs to.
   * 2. Finds the closest point on the surface (IB point) using CGAL AABB trees.
   * 3. Computes the normal distance from the ghost point to the surface.
   * 4. Constructs a local orthonormal frame (normal, tangent1, tangent2) at the IB point.
   * 5. Projects "Image Points" (IMs) into the fluid domain along the surface normal.
   * 6. Computes trilinear interpolation weights and indices for these Image Points.
   *
   * All computed data is stored in the `gpData` structure for use in boundary condition reconstruction.
   *
   * \param lev The current AMR level index.
   */
  void initialiseGPs(int lev) {
    auto& mfab = *bmf_a[lev];
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      auto& ibFab = mfab.get(mfi);
      auto& gpData = ibFab.gpData;

      // Reserve memory to avoid reallocations during push_back
      if (gpData.ngps > 0) {
          gpData.gp_ijk.reserve(gpData.ngps);
          gpData.disGP.reserve(gpData.ngps);
          gpData.ib_xyz.reserve(gpData.ngps);

          gpData.geomIdx.reserve(gpData.ngps);
          gpData.elemIdx.reserve(gpData.ngps);

          gpData.imp_xyz.reserve(gpData.ngps);
          gpData.imp_ijk.reserve(gpData.ngps);
          gpData.disIM.reserve(gpData.ngps);

          gpData.imp_ninterp.reserve(gpData.ngps);
          gpData.imp_ipweights.reserve(gpData.ngps);
          gpData.imp_ip_ijk.reserve(gpData.ngps);       
      }

      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);  // uint8_t array

      // we need a CPU loop here (cannot be GPU loop) as CGAL tree seach for
      // closest element to a point needs to be called. instead of looping through
      // previously indexed gps, we loop through the whole ghost point field as it
      // is available on GPU and CPU at all times. Unlike the gp indexes, which
      // are only stored on GPU memory. Array1D<int,0,AMREX_SPACEDIM-1>& idx =
      // ibFab.gpData.gp_ijk[ii];

      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA), [&](int i, int j, int k) {
        // for each ghost point
        if (ibMarkers(i, j, k, 1)) {
          Real x = prob_lo[0] + (0.5_rt + i) * dx_a[lev][0];
          Real y = prob_lo[1] + (0.5_rt + j) * dx_a[lev][1];
#if (AMREX_SPACEDIM == 2)
          Point gp(x, y);
#else
          Real z = prob_lo[2] + (0.5_rt + k) * dx_a[lev][2];
          Point gp(x, y, z);
#endif
          gpData.gp_ijk.push_back(make_vec<int>(i, j, k));
          
          // find and store geometry index for this GP.
          // Since this is a ghost point, it must be a solid point, 
          // so ibMarkers(i,j,k,0) stores (geometry_index + 1).
          int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;

          AMREX_ASSERT_WITH_MESSAGE(geomIdx >= 0 && geomIdx < ngeom, 
                  "Invalid geometry index in initialiseGPs");
          gpData.geomIdx.push_back(geomIdx);

          // closest surface/edge point and surface/edge 
          Point_and_primitive_id closest_elem =
              tree_pa[geomIdx]->closest_point_and_primitive(gp);

          // map PrimitiveID to integer index and store that.
          PrimitiveID elm = closest_elem.second;
          int f_idx = IdxMap_a[geomIdx].at(elm);
          gpData.elemIdx.push_back(f_idx);

          // This closest point (cp) is between the face plane and the gp
          Point cp = closest_elem.first;

          // IB point -------------------------------------------
          Real disGP = sqrt(CGAL::squared_distance(gp, cp));
          AMREX_ASSERT_WITH_MESSAGE(
              disGP < diag_a[lev],  "Ghost point and IB point distance larger than mesh diagonal");
          gpData.disGP.push_back(disGP);

          // ib_xyz
          gpData.ib_xyz.push_back(make_vec<Real>(cp));
          
          // local frame at f_idx 
          const auto& localframe = LocalFrame_a[f_idx];

          // IM points -------------------------------------------
          Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
          Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;
          Array1D<Real, 0, eorder_tparm - 1> disIM;
          Array1D< int, 0, eorder_tparm - 1> imp_ninterp;
          
          // find the closest bottom left point of the image point
          // 
          //  demo In 2D, same idea in 3D.
          //     i,j+1 (2) ---------------------     i+1,j+1 (3)
          //     |                                  |
          //     |         P                        |
          //     |                                  |
          //     |                                  |
          //     i,j  (1) ----------------------      i+1,j  (4)

          for (int jj = 0; jj < eorder_tparm; jj++) {

            // Starting point for this image point:
            //   - for jj == 0: use the original control point cp
            //   - for jj > 0 : start from the previous image point position
            Point cp_start = (jj == 0) ? cp :
#if (AMREX_SPACEDIM == 2)
            Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1));
#else
            Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1), imp_xyz(jj - 1, 2));
#endif
          
            if (jj == 0) {
              // =======================================================
              // First image point: try multiple distances (based on N_ATTEMPTS and IMP_FACTOR) and take the best one
              //   - candidate 1: IMP_FACTOR[0] * di_a[lev]
              //   - candidate 2: IMP_FACTOR[1] * di_a[lev] (only if candidate 1 is not ideal)
              //   - candidate x: ......
              // select the first candidate that yields enough fluid points,
              // or the candidate with the maximum number of fluid points.
              // if no candidates has fluid points more than INTERP_THRESHOLD, abort
              // =======================================================
              search_optimal_image_point<eorder_tparm, iorder_tparm>(cp_start, localframe, lev, prob_lo, bxg, ibMarkers, gpData,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);

            } 
            else {
              // =======================================================
              // Subsequent image points:
              //   - move one more di_a[lev] along the normal from the previous image point
              //   - no optimization, no multiple attempts, never abort even if not enough interp points found
              //   - only check that the point stays inside the grown box
              //   - mark as invalid if it leaves the box, set number of interp points to zero.
              // =======================================================
              search_image_point<eorder_tparm, iorder_tparm>(jj, cp_start, localframe, lev, prob_lo, bxg, ibMarkers, gpData,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);

            } // end if 
          } // end loop on image points

          // store gpData
          // Note: We are using push_back here because we are building the list dynamically.
          // The resize() method added to gpData_t is intended for cases where we know the size in advance
          // or want to clear/reset. For this initialization loop, push_back is still appropriate
          // unless we pre-calculate the number of ghost points.
          gpData.imp_xyz.push_back(imp_xyz);
          gpData.imp_ijk.push_back(imp_ijk);
          gpData.disIM.push_back(disIM);
          gpData.imp_ninterp.push_back(imp_ninterp);

          // Interpolation points' (ips) weights for each image point
          Array2D<Real, 0, eorder_tparm - 1 , 0, N_InterP -1 > imp_ipweights;
          Array3D< int, 0, eorder_tparm - 1 , 0, N_InterP -1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;
          
          computeIPweights<eorder_tparm, iorder_tparm>(
                          imp_ipweights, imp_ip_ijk, imp_xyz, imp_ijk, imp_ninterp,
                          prob_lo, dx_a[lev], ibMarkers);
          
          // store
          gpData.imp_ipweights.push_back(imp_ipweights);
          gpData.imp_ip_ijk.push_back(imp_ip_ijk);

          } //end if (ibMarkers(i,j,k,1))
      });//end loop on bx
    
    } //end MFIter
  }

/**
 * \brief Reconstructs and applies immersed-boundary ghost-cell primitive states on the current FAB.
 *
 * For each ghost point on this MFIter tile and AMR level, this routine interpolates image-point data,
 * applies the wall model, extrapolates surface states back to the ghost cell, and writes the resulting
 * primitive variables into the `prims` array.
 *
 * \param mfi   MFIter pointing to the current FAB.
 * \param cons  Conservative variables (read-only).
 * \param prims Primitive variables to be updated at ghost cells.
 * \param cls   Pointer to the physics/closure class.
 * \param lev   Current AMR level index.
 */ 
void computeGPs(const MFIter& mfi,
                const Array4<Real>& cons,
                const Array4<Real>& prims,
                const cls_t* cls,
                int& lev)
{
  auto& mfab               = *bmf_a[lev];
  const auto& ibFab        = mfab.get(mfi);

  // Ghost-point related data (geometry + interpolation)
  auto const gp_ijk        = ibFab.gpData.gp_ijk.data();
  auto const imp_ipweights = ibFab.gpData.imp_ipweights.data();
  auto const imp_ip_ijk    = ibFab.gpData.imp_ip_ijk.data();
  auto const disGP         = ibFab.gpData.disGP.data();
  auto const disIM         = ibFab.gpData.disIM.data();
  auto const ib_xyz        = ibFab.gpData.ib_xyz.data();
  auto const imp_ninterp   = ibFab.gpData.imp_ninterp.data();
  auto const elemIdx       = ibFab.gpData.elemIdx.data();

  auto const* lf_ptr = LocalFrame_a.data();

  // Local copy of prims on grown box (GPU-friendly source field)
  const Box& bxg = mfi.growntilebox(cls->NGHOST);
  FArrayBox primf(bxg, cls_t::NPRIM, The_Async_Arena());
  Array4<Real> const& prims0 = primf.array();

  ParallelFor(bxg, cls_t::NPRIM,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
      prims0(i,j,k,n) = prims(i,j,k,n);
  });

  int const ngps = ibFab.gpData.ngps;

  ParallelFor(ngps, [=, copy=this] AMREX_GPU_DEVICE (int ii) noexcept
  {
    // --------------------------------------------------------------------
    // 1) Reconstruct local orthonormal frame from stored LocalFrame_a
    //     localframe(ii)(row, col): row = {0: n, 1: t1, (2: t2 in 3D)}
    // --------------------------------------------------------------------
    int elem_idx = elemIdx[ii];
    const auto& frame = lf_ptr[elem_idx];

#if (AMREX_SPACEDIM == 2)
    Array1D<Real, 0, AMREX_SPACEDIM - 1>  nvec = {   frame.normal[0],   frame.normal[1] };
    Array1D<Real, 0, AMREX_SPACEDIM - 1> t1vec = { frame.tangent1[0], frame.tangent1[1] };
    Array1D<Real, 0, AMREX_SPACEDIM - 1> t2vec = { 0.0, 0.0 }; 
#else
    Array1D<Real, 0, AMREX_SPACEDIM - 1>  nvec = {   frame.normal[0],   frame.normal[1],   frame.normal[2] };
    Array1D<Real, 0, AMREX_SPACEDIM - 1> t1vec = { frame.tangent1[0], frame.tangent1[1], frame.tangent1[2] };
    Array1D<Real, 0, AMREX_SPACEDIM - 1> t2vec = { frame.tangent2[0], frame.tangent2[1], frame.tangent2[2] };
#endif

    // --------------------------------------------------------------------
    // 2) Storage for primitive variables along the normal:
    //    index 0: ghost cell
    //    index 1: surface / reference point (convention dependent)
    //    index 2..(1+eorder_tparm): image points
    // --------------------------------------------------------------------
    Array2D<Real, 0, eorder_tparm + 1, 0, cls_t::NPRIM - 1> primsNormal;
    for (int p = 0; p <= eorder_tparm + 1; ++p) {
        for (int n = 0; n < cls_t::NPRIM; ++n) {
            primsNormal(p,n) = Real(0.0);
        }
    }

    // 3) Interpolate primitive variables at all image points from prims0
    copy->interpolateIMs<eorder_tparm, iorder_tparm>(imp_ip_ijk[ii], imp_ipweights[ii], prims0, primsNormal);

    // 4) Transform velocities at image points (> 1) to local frame
    for (int iip = 2; iip < 2 + eorder_tparm; ++iip) {
        copy->global2local<eorder_tparm>(iip, primsNormal, nvec, t1vec, t2vec);
    }

    // 5) Apply wall model at IB surface to set surface states (u, P, T, Y, ...)
    wallmodel::compute_surfIB(ib_xyz[ii], nvec, primsNormal, cls);

    // 6) Extrapolate from surface/image points back to ghost point along n
    copy->extrapolate<eorder_tparm>(primsNormal, imp_ninterp[ii], disGP[ii], disIM[ii]);

    // 7) Transform ghost-point velocity back to global coordinates
    int idx = 0;
    copy->local2global<eorder_tparm>(idx, primsNormal, nvec, t1vec, t2vec);

    // 8) Extract primitive variables at ghost point (slot 0)
    Real P = primsNormal(0, cls_t::QPRES);
    Real T = primsNormal(0, cls_t::QT);

    Real Y[NUM_SPECIES] = { Real(0.0) };
#if NUM_SPECIES > 1
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = primsNormal(0, cls_t::QFS + n);
    }
#endif

    Real ux = primsNormal(0, cls_t::QU);
    Real uy = primsNormal(0, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    Real uz = primsNormal(0, cls_t::QW);
#else
    Real uz = Real(0.0);
#endif

    // 9) Enforce thermodynamic consistency and fill complete primitive vector Q
    Real Q[cls_t::NPRIM];
    cls->ensurePTYfillq(P, T, Y, ux, uy, uz, Q);

    // 10) Write ghost-cell primitive variables back into prims
    int i = gp_ijk[ii](0);
    int j = gp_ijk[ii](1);
#if (AMREX_SPACEDIM == 3)
    int k = gp_ijk[ii](2);
#else
    int k = 0;  // in 2D, k-index is always 0
#endif

    for (int n = 0; n < cls_t::NPRIM; ++n) {
        prims(i,j,k,n) = Q[n];
    }
  }); // end ParallelFor over ghost points
}

/**
  * \brief Compute surface indices and interpolation data for all faces at given level
  *
  * Algorithm:
  *  1. Build spatial lookup (global_fab_idx -> local_fab_idx)
  *  2. For each face: compute mirror point, find owning FAB, compute interpolation weights
  *  3. Build CSR structure for GPU-friendly access
  *
  * \param lev AMR level
  */
void compute_surface_index(int lev) {

  // local rank of this process
  int myrank = amrex::ParallelDescriptor::MyProc();
  amrex::Print() << "Compute Surface Index at LEVEL " << lev << std::endl;

  auto& mfab = *bmf_a[lev];
  const int nfab_local = mfab.local_size();

  const BoxArray& ba = mfab.boxArray();
  const DistributionMapping& dm = mfab.DistributionMap();

  const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
  const auto& domain = amr_p->Geom(lev).Domain();

  // ========================================================================
  // Phase 0: Initialize surfdata_soa
  // ========================================================================
  // Clear previous data to ensure clean state (especially after regrid)
  // This also releases memory if capacity was large
  surfdata_soa.clear();
  
  // Resize to accommodate all faces (currently Replicated mode)
  // TODO: In future Distributed mode, this should only resize to n_local_faces
  surfdata_soa.resize(ntotalfaces);
  
  // Reset metadata
  // Note: resize() already handles basic initialization, but we set specific defaults here
  // We can use parallel for on CPU for faster initialization
  amrex::ParallelFor(ntotalfaces, [&](int i) {
      surfdata_soa.pointfound[i] = 0; // false
      surfdata_soa.rank[i] = -99;
      surfdata_soa.ifab[i] = -1;
      surfdata_soa.lev[i] = -1;
  });

  // ========================================================================
  // Phase 1: Build spatial lookup structures
  // ========================================================================
  // global to local fab index mapping
  Vector<int> global_to_local(ba.size(), -1);
  
  // fab arrays and boxes for fast access
  Vector<Array4<uint8_t const>> fab_markers(nfab_local);
  Vector<Box> fab_bx(nfab_local);
  Vector<Box> fab_bxg(nfab_local);
  
  for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {

      int gidx = mfi.index();
      int lidx = mfi.LocalIndex();
      global_to_local[gidx] = lidx;
      
      fab_markers[lidx] = mfab.const_array(mfi);
      fab_bx[lidx] = mfi.tilebox();
      fab_bxg[lidx] = mfi.growntilebox(cls_t::NGHOST);
  }

  // ========================================================================
  // Fast cell->FAB lookup using BoxArray
  // Returns local_fab index if the cell (i,j,k) AND its neighbor (i+1,j+1,k+1) 
  // are both contained in the same FAB (either valid or grown box).
  // Otherwise returns -1.
  auto find_local_fab_for_cell = [&](int i, int j, int k) -> int {
    
    IntVect iv_lo = AMREX_D_DECL(i, j, k);
    IntVect iv_hi = AMREX_D_DECL(i+1, j+1, k+1);

    // Check if points are inside global domain first
    if (!domain.contains(iv_lo) || !domain.contains(iv_hi)) 
        return -1;

    // Construct a small box covering both points
    Box query_box(iv_lo, iv_hi);

    // --------------------------------------------------------------------
    // Single Pass: Search grown box directly
    // Since grown box includes valid box, and we have FillPatch,
    // checking grown box is sufficient and covers stencil requirements.
    // --------------------------------------------------------------------
    for (int gidx = 0; gidx < static_cast<int>(ba.size()); ++gidx) {
        if (dm[gidx] != myrank) continue;
        
        int lidx = global_to_local[gidx];
        if (lidx < 0) continue;
        
        // Check if the FAB's grown box contains BOTH points
        if (fab_bxg[lidx].contains(query_box)) {
            return lidx;
        }
    }
    return -1;  // Not found
  };
  
  // ========================================================================
  // Phase 2: Process all faces - temporary storage for per-FAB lists
  // ========================================================================
  Vector<Vector<int>> faces_per_fab(nfab_local);
  int faces_found = 0;
  int faces_notfound = 0;

  // Loop over all faces (Replicated geometry)
  // Each core checks if it owns the image points for a face
  for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {
      
      // Get geometry data from pre-computed arrays
      const LocalFrame& localframe = LocalFrame_a[f_idx];
      const SurfElem& surfelem = SurfElem_a[f_idx];

      // Determine geometry index directly from SurfElem
      int geomIdx = surfelem.geomIdx;
 
      // Set metadata in surfdata_soa needed for search functions
      surfdata_soa.geomIdx[f_idx] = geomIdx;
      surfdata_soa.elemIdx[f_idx] = f_idx;

      // --------------------------------------------------------------------
      // Find owning FAB and Compute Image Points
      // --------------------------------------------------------------------

      // Loop over image points
      for (int jj = 0; jj < eorder_tparm_surf; ++jj) {

          Point cp_start;
          if (jj == 0) {
              cp_start = surfelem.centroid;
          } else {
#if (AMREX_SPACEDIM == 2)
              cp_start = Point(surfdata_soa.imp_xyz[f_idx](jj - 1, 0), surfdata_soa.imp_xyz[f_idx](jj - 1, 1));
#else
              cp_start = Point(surfdata_soa.imp_xyz[f_idx](jj - 1, 0), surfdata_soa.imp_xyz[f_idx](jj - 1, 1), surfdata_soa.imp_xyz[f_idx](jj - 1, 2));
#endif
          }

          if (jj == 0) {

          int local_fab = -1;
          for (int attempt = 0; attempt < N_ATTEMPTS; ++attempt) {

            // Compute best image point position for this attempt
            int best_fluid = 0;
            Array1D<Real, 0, AMREX_SPACEDIM - 1> best_xyz;
            Array1D< int, 0, AMREX_SPACEDIM - 1> best_ijk;
            
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                best_xyz[d] = cp_start[d] + IMP_FACTOR[attempt] * di_a[lev] * localframe.normal[d];
                best_ijk[d] = static_cast<int>(std::floor(
                    (best_xyz[d] - prob_lo[d]) / dx_a[lev][d] - 0.5
                ));
            }

            int mi = best_ijk[0];
            int mj = best_ijk[1];
#if (AMREX_SPACEDIM == 3)
            int mk = best_ijk[2];
#else
            int mk = 0;
#endif
            local_fab = find_local_fab_for_cell(mi, mj, mk);

            if (local_fab >= 0) {
              // Found valid FAB, check if image point is in fluid
              auto const ibMarkers = fab_markers[local_fab];
#if (AMREX_SPACEDIM == 2) 
              int n_fluid = valid_mirror(mi, mj,  0, ibMarkers);
#else
              int n_fluid = valid_mirror(mi, mj, mk, ibMarkers);
#endif
              if (n_fluid > best_fluid)
              {
                  best_fluid = n_fluid;
                  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                      surfdata_soa.imp_xyz[f_idx](jj, d) = best_xyz[d];
                      surfdata_soa.imp_ijk[f_idx](jj, d) = best_ijk[d];
                  }
                  surfdata_soa.disIM[f_idx](jj) = IMP_FACTOR[attempt] * di_a[lev];
                  surfdata_soa.imp_ninterp[f_idx](jj) = n_fluid;
                  surfdata_soa.pointfound[f_idx] = 1; // true
              }     
            } // end if local_fab
          } // end loop on attempts
          
          if (local_fab < 0) {
              // Mark as not found
              surfdata_soa.imp_ninterp[f_idx](jj) = 0; 
          } 
          } // end if jj==0
          else {
              // Subsequent image points:
              //   - move one more di_a[lev] along the normal from the previous image point
              //   - no optimization, no multiple attempts, never abort even if not enough interp points found
              //   - only check that the point stays inside the global domain
              //   - mark as invalid if it leaves the box, set number of interp points to zero.

              for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                surfdata_soa.imp_xyz[f_idx](jj, d) = cp_start[d] + di_a[lev] * localframe.normal[d];
                surfdata_soa.imp_ijk[f_idx](jj, d) = static_cast<int>(std::floor(
                    (surfdata_soa.imp_xyz[f_idx](jj, d) - prob_lo[d]) / dx_a[lev][d] - 0.5
                ));
            }

            int mi = surfdata_soa.imp_ijk[f_idx](jj, 0);
            int mj = surfdata_soa.imp_ijk[f_idx](jj, 1);
#if (AMREX_SPACEDIM == 3)
            int mk = surfdata_soa.imp_ijk[f_idx](jj, 2);
#else
            int mk = 0;
#endif

            int local_fab = find_local_fab_for_cell(mi, mj, mk);

            if (local_fab >= 0) {
              // Found valid FAB, check if image point is in fluid
              auto const ibMarkers = fab_markers[local_fab];
#if (AMREX_SPACEDIM == 2) 
              int n_fluid = valid_mirror(mi, mj,  0, ibMarkers);
#else
              int n_fluid = valid_mirror(mi, mj, mk, ibMarkers);

#endif
              surfdata_soa.disIM[f_idx](jj) = (jj + 1) * di_a[lev];
              surfdata_soa.imp_ninterp[f_idx](jj) = n_fluid;
              surfdata_soa.pointfound[f_idx] = 1; // true
            } else {
              // Mark as invalid
              surfdata_soa.imp_ninterp[f_idx](jj) = 0;
              surfdata_soa.pointfound[f_idx] = 0; // false
            }
          } // end else (jj > 0)
      } // end loop over image points
      
      // If point was found (at least the first one), finalize metadata
      if (surfdata_soa.pointfound[f_idx]) {
          // We need to find the local_fab again for the first point to store it
          // Or we could have stored it in the loop. 
          // For simplicity, re-evaluate for the first point (jj=0)
          int mi = surfdata_soa.imp_ijk[f_idx](0, 0);
          int mj = surfdata_soa.imp_ijk[f_idx](0, 1);
#if (AMREX_SPACEDIM == 3)
          int mk = surfdata_soa.imp_ijk[f_idx](0, 2);
#else
          int mk = 0;
#endif
          int local_fab = find_local_fab_for_cell(mi, mj, mk);
          
          if (local_fab >= 0) {
              surfdata_soa.ifab[f_idx] = local_fab;
              surfdata_soa.lev[f_idx] = lev;
              surfdata_soa.rank[f_idx] = myrank;
              
              faces_per_fab[local_fab].push_back(f_idx);
              faces_found++;
              
              // Compute Interpolation Weights
              auto const ibMarkers = fab_markers[local_fab];
              computeIPweights<eorder_tparm_surf, iorder_tparm_surf>(
                  surfdata_soa.imp_ipweights[f_idx], 
                  surfdata_soa.imp_ip_ijk[f_idx], 
                  surfdata_soa.imp_xyz[f_idx], 
                  surfdata_soa.imp_ijk[f_idx], 
                  surfdata_soa.imp_ninterp[f_idx],
                  prob_lo, dx_a[lev], ibMarkers);
          } else {
              // Should not happen if pointfound is true, but safety check
              surfdata_soa.pointfound[f_idx] = 0;
          }
      }

  }  // end loop over faces
              surf.imp_ninterp(jj) = 0; 















      // --------------------------------------------------------------------
      // Compute Interpolation Weights
      // --------------------------------------------------------------------
      computeIPweights<eorder_tparm_surf, iorder_tparm_surf>(
          surf.imp_ipweights, surf.imp_ip_ijk, surf.imp_xyz, surf.imp_ijk, surf.imp_ninterp,
          prob_lo, dx_a[lev], ibMarkers);

      // --------------------------------------------------------------------
      // Finalize
      // --------------------------------------------------------------------
      surf.ifab = local_fab;
      surf.lev = lev;
      surf.rank = myrank;
      surf.pointfound = true;

      faces_per_fab[local_fab].push_back(f_idx);
      faces_found++;

  }  // end loop over faces

  // ========================================================================
  // Phase 3: Build CSR structure
  // ========================================================================
  if (static_cast<int>(faces_per_level.size()) <= lev) {
      faces_per_level.resize(lev + 1);
  }
  
  auto& csr = faces_per_level[lev];
  csr.clear();
  csr.fab_offsets.resize(nfab_local + 1);
  
  csr.fab_offsets[0] = 0;
  for (int ifab = 0; ifab < nfab_local; ++ifab) {
      csr.fab_offsets[ifab + 1] = csr.fab_offsets[ifab] 
                                + static_cast<int>(faces_per_fab[ifab].size());
  }
  
  int total = csr.fab_offsets[nfab_local];
  csr.face_indices.resize(total);
  
  int idx = 0;
  for (int ifab = 0; ifab < nfab_local; ++ifab) {
      for (int f : faces_per_fab[ifab]) {
          csr.face_indices[idx++] = f;
      }
  }

  amrex::Print() << "  Level " << lev << ": " << faces_found 
                  << " faces found, " << faces_notfound << " not found\n";
}










///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief compute surface properties for each surface face
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_surface_props(MultiFab& stateprops,const cls_t* cls,int lev) {
  

  int myrank = amrex::ParallelDescriptor::MyProc(); 
  amrex::Print() << " Compute Surface Properties at LEVEL" << lev << std::endl;

  auto& mfab = *bmf_a[lev];
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

  // loop over mfi
  for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) { 
    
    //const int ifab= mfi.index();
    const int ifab = mfi.LocalIndex();
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
      
      const GeomType& mesh   = geom_a[ii];
                        
      const int nfaces_in_fab_lev = intfaces_in_fab[lev][ifab].size();
      //const int nfaces_in_fab_lev = nfaces_infab_inlevel[lev][ifab];

     // printf(" myrank =%d nfaces_in_fab(%d) = %d nfaces2=%d\n", myrank, ifab, nfaces_in_fab,nfaces2);

      for (int j = 0; j < nfaces_in_fab_lev; ++j){
        
        int iface = intfaces_in_fab[lev][ifab][j];
        
        auto& surf_dat = surfdata_a[iface];

        // extract arrays from surface data
        auto const norm      = surf_dat.norm;
        auto const ip_ijk    = surf_dat.ip_ijk;
        auto const ipweights = surf_dat.ipweights;
        auto const ib_xyz    = surf_dat.ib_xyz;     
        auto const o_dis     = surf_dat.o_dis;

        // only compute properties at the correct level and fab
        if ((surf_dat.lev == lev) && (surf_dat.ifab == ifab)) {
                      
          Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1> primsNormal={0.0};                

          
          // check point is correct           
          int i1 = ip_ijk(0,0,0);int j1 = ip_ijk(0,0,1);
#if (AMREX_SPACEDIM == 3)
          int k1 = ip_ijk(0,0,2);          
          if (! bxg.contains(i1, j1, k1) ){
#else
          if (! bxg.contains(i1, j1) ){
#endif
          
          printf(" ID=%d iface=%d, lev=%d, ifab = %d, j= %d \n", myrank,iface, lev, ifab, j);          
          amrex::Error("Interpolation point outside fab");
          }

          // interpolate primitive variables at the mirror points
          interpolateIMs(ip_ijk,ipweights,prims,primsNormal);
          // apply wall model
          wallmodel::compute_surfIB(ib_xyz,norm,primsNormal,cls);   

          // compute one-sided gradients dT/dn du/dn
          Real dTdn = (primsNormal(2,cls_t::QT) - primsNormal(1,cls_t::QT))*o_dis;

          // store values
          surf_dat.pressure       = primsNormal(1,cls_t::QPRES); 
          surf_dat.temperature    = primsNormal(1,cls_t::QT);          
          surf_dat.dTdn           = dTdn;                      

        }
        
      } //end loop faces          
    } //end loop geometry 
    //...........................................
  } // end looping mfi
}
} 


  // ////////////////////////////////////////////////////////////////
  // \brief put all surface data to ioproc for plotting
  //        This is called per level
  void gather_surfdata_to_rank0(int lev) {

    int myrank = amrex::ParallelDescriptor::MyProc();
    int nprocs = amrex::ParallelDescriptor::NProcs();
    
    // amrex::Print() << " Gathered surface data from nranks= " << nprocs << std::endl;
    
    // Step 1: Pack valid entries locally
    std::vector<surfData_t<iorder_tparm>> local_data;
    //for (const auto& dat : surfdata_a) {
    
    for (int j=0;j<ntotalfaces;j++) { 
        auto& dat = surfdata_a[j];

        bool tosend = dat.pointfound && (dat.lev == lev);
       
        if (tosend) {
            local_data.push_back(dat);
        }
    }

    int nlocal = local_data.size(); // amount of entries to send


    // Step 2: Gather counts
    std::vector<int> counts = amrex::ParallelDescriptor::Gather(nlocal, 0);


    // if (myrank == 0) {
    // for (int i = 0; i < nprocs; ++i)
    //     amrex::Print() << "Rank " << i << " sent " << counts[i] << " entries.\n";
    // }

    // Step 3: Calculate displacements and total size
    std::vector<int> displs;
    int total = 0;
    if (myrank == 0) {
        displs.resize(nprocs, 0);
        for (int i = 1; i < nprocs; ++i)
            displs[i] = displs[i - 1] + counts[i - 1];
        for (int c : counts) total += c;
    }

    using T = surfData_t<iorder_tparm>;    

    // Compute sizes in bytes
    size_t typesize = sizeof(T);
    std::vector<char> sendbuf(reinterpret_cast<char*>(local_data.data()),
                              reinterpret_cast<char*>(local_data.data()) + nlocal * typesize);
    std::vector<char> recvbuf;
    if (myrank == 0) recvbuf.resize(total * typesize);

    // Convert counts and displs to bytes
    std::vector<int> counts_bytes, displs_bytes;
    if (myrank == 0) {
      counts_bytes.resize(nprocs);
      displs_bytes.resize(nprocs);
      for (int i = 0; i < nprocs; ++i) {
        counts_bytes[i] = counts[i] * typesize;
        displs_bytes[i] = displs[i] * typesize;
      }
    }

    // Now use Gatherv on bytes
    amrex::ParallelDescriptor::Gatherv(
        sendbuf.data(), nlocal * typesize,
        recvbuf.data(), counts_bytes, displs_bytes, 0);
    
    // Reconstruct on rank 0
    if (myrank == 0) {
      AMREX_ASSERT(surfdata_a.size() == ntotalfaces);
      T* recv_ptr = reinterpret_cast<T*>(recvbuf.data());
      for (int i = 0; i < total; ++i) {
        const auto& dat = recv_ptr[i];
        surfdata_a[dat.iface] = dat;               
      }

      //std::cout << "Gathered " << total << " surface entries.\n";
    }

    //amrex::Print() << "Gathered surface data from all ranks.\n";
  }
  ////////////////////////////////////////////////////////////////
  //  \brief plot surface mesh to file
  //  \param igeom geometry index
  //  \param filename output file name
  //  \note uses CGAL Polygon_mesh_processing IO functions
  void plot_surface(const amrex::Real time,const int igeom, const std::string& filename) {
    const Polyhedron& mesh   = geom_a[igeom];
    const auto& face_normals = fnorm_a[igeom];

    std::ofstream out(filename);

    int myrank = amrex::ParallelDescriptor::MyProc();

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
    out << "VECTORS face_normals float\n";
    for (auto fd : faces(mesh)) {
      Array1D<Real, 0, AMREX_SPACEDIM - 1> n = {
                fnorm_a[igeom][fd][0], fnorm_a[igeom][fd][1],
                fnorm_a[igeom][fd][2]};
        out << n(0) << " " << n(1) << " " << n(2) << "\n";
    }

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
      

      //printf("PLOTTING iface=%d T=%f \n", iface, surfdata_a[iface].temperature);
      if (surfdata_a[iface].temperature <  100.0_rt) {
        printf("ERROR  iface=%d T=%f \n", iface, surfdata_a[iface].temperature);
        exit(0);
      }
      //if (iface==10) exit(0);
      

    }
    out << "SCALARS gradT float 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<ntotalfaces;iface++) {      
      out << static_cast<float>(surfdata_a[iface].dTdn) << "\n";
    }

    out << "SCALARS rank int 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<ntotalfaces;iface++) {      
      out << static_cast<float>(surfdata_a[iface].rank) << "\n";
    }

    out << "SCALARS level int 1\n";
    out << "LOOKUP_TABLE default\n";    
    for (int iface=0;iface<ntotalfaces;iface++) {      
      out << static_cast<float>(surfdata_a[iface].lev) << "\n";
    }


    // !!! important reset of point founds      
    for (int iface=0;iface<num_faces;iface++) { 
      if ( !(myrank == surfdata_a[iface].rank)) {
        surfdata_a[iface].pointfound = false;
      }
    }

    out.close();
    amrex::Print() << "----------------------------------\n";
    amrex::Print() << "Surface mesh plotted to " << filename << "\n";
    amrex::Print() << "----------------------------------\n";

  }

//////////////////////////////////////////////////////////////////



//============================================================================
///--------------------------- private functions -----------------------------
private:
  /*//////////////////////////////////////////////////////////////////
   * \brief Checks if an image point (mirror point) location is valid for interpolation.
   *
   * A mirror point is considered valid if its interpolation stencil (the surrounding 2x2 or 2x2x2 cell block)
   * contains enough fluid points. Specifically, it returns true if at least 3 of the surrounding grid points
   * are marked as fluid (value 0 in ibMarkers).
   *
   * \param i The i-index of the bottom-left (or reference) grid point of the stencil.
   * \param j The j-index of the bottom-left (or reference) grid point of the stencil.
   * \param k The k-index of the bottom-left (or reference) grid point of the stencil.
   * \param ibMarkers The marker array indicating solid (non-zero) or fluid (zero) status.
   * \return the number of fluid points in the stencil.
   *///////////////////////////////////////////////////////////////////
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int valid_mirror(int i, int j, int k, const Array4<const uint8_t>& ibMarkers) const {
    
    int fluid_count = 0;

#if (AMREX_SPACEDIM == 2)
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
          int ii = i + di;
          int jj = j + dj;
          // Check if point is fluid (marker == 0)
          if (ibMarkers(ii, jj, 0) == 0) {
              ++fluid_count;
          }
      }
    }
    // In 2D, we require at least 3 fluid points out of 4
    return fluid_count;
#else
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        for (int dk = 0; dk <= 1; ++dk) {
          int ii = i + di;
          int jj = j + dj;
          int kk = k + dk;
          // Check if point is fluid (marker == 0)
          if (ibMarkers(ii, jj, kk, 0) == 0) {
            ++fluid_count;
          }
        }
      }
    }
    // In 3D, we require at least 3 fluid points out of 8 (heuristic, can be adjusted)
    return fluid_count;
#endif
  }

  // Type trait to detect if a type is gpData_t (has gp_ijk member)
  template <typename T, typename = void>
  struct is_gpData_t : std::false_type {};

  template <typename T>
  struct is_gpData_t<T, std::void_t<decltype(std::declval<T>().gp_ijk)>> : std::true_type {};

  /*//////////////////////////////////////////////////////////////////
  * \brief Find the best location for the first image point.
  *
  * Starting from cp_start and moving along the outward normal in localframe,
  * this function tests candidate image-point locations and selects the one
  * with the largest number of valid fluid cells in its interpolation stencil.
  * It then stores the chosen image-point position, index, normal distance,
  * and stencil size.
  *
  * \tparam eorder_t    Number of image points (exterpolation order).
  * \tparam iorder_t    Number of interpolation points per image point.
  * \tparam IPDATA      Type of Data structure that holds image points.
  * \tparam GP_OR_SURF 1 = ghost point data structure, 0 = surface data structure.
  *
  * \param cp_start     Starting point physical coordinates where image point search begins.
  * \param localframe   Local frame containing the normal vector.
  * \param lev          AMR level.
  * \param prob_lo      Lower bounds of the problem domain.
  * \param bxg          Grown box of the current FAB.
  * \param ibMarkers    Marker array: non-zero = solid, zero = fluid.
  * \param ipData       Data structure that holds image points.
  * \ Outputs:
  * \param imp_xyz      Output: coordinates of image points.
  * \param imp_ijk      Output: grid indices of image points.
  * \param disIM        Output: normal distances of image points.
  * \param imp_ninterp  Output: number of fluid points in the stencil
  *                     for each image point.
  *///////////////////////////////////////////////////////////////////
  template <int eorder_t, int iorder_t, typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
  AMREX_FORCE_INLINE //AMREX_GPU_HOST_DEVICE
  void search_optimal_image_point(
      const Point& cp_start,
      const LocalFrame& localframe,
      int lev,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const Box& bxg,
      const Array4<uint8_t const>& ibMarkers,
      const IPDATA& ipData,
      // Outputs (passed by reference to update the arrays at index 0)
      Array2D<Real, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
      Array2D< int, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
      Array1D<Real, 0, eorder_t - 1>& disIM,
      Array1D< int, 0, eorder_t - 1>& imp_ninterp)
  {

    // =======================================================
    // First image point: try multiple distances (based on N_ATTEMPTS and IMP_FACTOR) and take the best one
    //   - candidate 1: IMP_FACTOR[0] * di_a[lev]
    //   - candidate 2: IMP_FACTOR[1] * di_a[lev] (only if candidate 1 is not ideal)
    //   - candidate x: ......
    // =======================================================
    int best_fluid = 0;
    Array1D<Real, 0, AMREX_SPACEDIM - 1> best_xyz;
    Array1D< int, 0, AMREX_SPACEDIM - 1> best_ijk;
    
    for (int attempt = 0; attempt < N_ATTEMPTS; ++attempt) {
    
      // Candidate attempt: IMP_FACTOR[attempt] * di along the outward normal
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          best_xyz(d) = cp_start[d] + IMP_FACTOR[attempt] * di_a[lev] * localframe.normal[d];
          best_ijk(d) = int(std::floor(
              (best_xyz(d) - prob_lo[d]) / dx_a[lev][d] - 0.5
          ));
      }
    
      // check stencil validity and count number of fluid points in stencil
      // First image point is the most critical one, so we abort if stencil is invalid immediately 
      // (as the first imp goes outside the box, the further imps will be even worse)
#if (AMREX_SPACEDIM == 2)
      bool in_box = check_interpolation_stencil<IPDATA>(best_ijk(0), best_ijk(1), 0, 
                                                bxg, lev,
                                                ipData,  
                                                (attempt == 0) ? CheckMode::Abort : CheckMode::Silent);
      int n_fluid = (in_box) ? valid_mirror(best_ijk(0), best_ijk(1), 0, ibMarkers) : -1;
#else
      bool in_box = check_interpolation_stencil<IPDATA>(best_ijk(0), best_ijk(1), best_ijk(2),
                                                bxg, lev,
                                                ipData,  
                                                (attempt == 0) ? CheckMode::Abort : CheckMode::Silent);
      int n_fluid = (in_box) ? valid_mirror(best_ijk(0), best_ijk(1), best_ijk(2), ibMarkers) : -1;
#endif

      // if this is the best so far, store it 
      if (n_fluid > best_fluid) {
        best_fluid = n_fluid;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            imp_xyz(0, d) = best_xyz(d);
            imp_ijk(0, d) = best_ijk(d);
        }

        // store number of fluid points found
        imp_ninterp(0) = n_fluid; 
        // store normal distance for the first image point
        disIM(0) = IMP_FACTOR[attempt] * di_a[lev];
      }

      // if ideal number of fluid points found, break
      if (best_fluid == N_InterP) {
        break;
      }
    } // end loop on attempt

    if (best_fluid < INTERP_THRESHOLD) {
        // Extract context info for detailed error message
        int current_geom;
        int current_elem;
        const char* point_label;
        Real p_x = 0.0, p_y = 0.0, p_z = 0.0;

        if constexpr (GP_OR_SURF == 1) {
             current_geom = ipData.geomIdx.back();
             current_elem = ipData.elemIdx.back();
             point_label = "Ghost Point";
             
             int gp_i = ipData.gp_ijk.back()[0];
             int gp_j = ipData.gp_ijk.back()[1];
             p_x = prob_lo[0] + (0.5_rt + gp_i) * dx_a[lev][0];
             p_y = prob_lo[1] + (0.5_rt + gp_j) * dx_a[lev][1];
#if (AMREX_SPACEDIM == 3)
             int gp_k = ipData.gp_ijk.back()[2];
             p_z = prob_lo[2] + (0.5_rt + gp_k) * dx_a[lev][2];
#endif
        } else {
             current_geom = ipData.geomIdx;
             current_elem = ipData.elemIdx;
             point_label = "Surface Point";
             
             p_x = cp_start[0];
             p_y = cp_start[1];
#if (AMREX_SPACEDIM == 3)
             p_z = cp_start[2];
#endif
        }
        
        Array1D<Real,0,AMREX_SPACEDIM-1> centroid;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            centroid(d) = SurfElem_a[current_elem].centroid[d];
        }

        // Use printf for better portability (works on GPU if needed) and avoid stringstream overhead
#if (AMREX_SPACEDIM == 3)
        std::printf("Not enough valid interpolation points found for the first image point!\n"
                    "  Level: %d\n"
                    "  %s: (%f, %f, %f)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Face Centroid:  (%f, %f, %f)\n"
                    "  Best Fluid Points Found: %d (Threshold: %d)\n",
                    lev, point_label, p_x, p_y, p_z, 
                    current_geom, current_elem, 
                    centroid[0], centroid[1], centroid[2],
                    best_fluid, INTERP_THRESHOLD);
#else
        std::printf("Not enough valid interpolation points found for the first image point!\n"
                    "  Level: %d\n"
                    "  %s: (%f, %f)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d \n"
                    "  Face Centroid:  (%f, %f)\n"
                    "  Best Fluid Points Found: %d (Threshold: %d)\n",
                    lev, point_label, p_x, p_y, 
                    current_geom, current_elem, 
                    centroid[0], centroid[1],
                    best_fluid, INTERP_THRESHOLD);
#endif
        std::fflush(stdout);
        amrex::Abort("Not enough valid interpolation points found for the first image point!");
    } // end check best fluid
  }

  // Helper function to compute image points without optimal searching
  template <int order_t, int iorder_t, typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
  AMREX_FORCE_INLINE
  void search_image_point(
      int jj,
      const Point& cp_start,
      const LocalFrame& localframe,
      int lev,
      const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
      const Box& bxg,
      const Array4<uint8_t const>& ibMarkers,
      const IPDATA& ipData,
      // Outputs (passed by reference to update arrays at index jj)
      Array2D<Real, 0, order_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
      Array2D< int, 0, order_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
      Array1D<Real, 0, order_t - 1>& disIM,
      Array1D< int, 0, order_t - 1>& imp_ninterp)
  {
    // Calculate position: move one more di_a[lev] along the normal
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        imp_xyz(jj, d) = cp_start[d] + di_a[lev] * localframe.normal[d];
        imp_ijk(jj, d) = int(std::floor(
            (imp_xyz(jj, d) - prob_lo[d]) / dx_a[lev][d] - 0.5
        ));
    }

    // Check stencil validity
#if (AMREX_SPACEDIM == 2)
    bool in_box = check_interpolation_stencil<IPDATA>(imp_ijk(jj, 0), imp_ijk(jj, 1), 0, 
                                              bxg, lev,
                                              ipData,  
                                              CheckMode::Silent);
    int fluid = (in_box) ? valid_mirror(imp_ijk(jj, 0), imp_ijk(jj, 1), 0, ibMarkers) : -1;
#else
    bool in_box = check_interpolation_stencil<IPDATA>(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2),
                                              bxg, lev,
                                              ipData,  
                                              CheckMode::Silent);   
    int fluid = (in_box) ? valid_mirror(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2), ibMarkers) : -1;
#endif

    // Accumulate normal distance along the ray
    disIM(jj) = (jj > 0) ? disIM(jj - 1) + di_a[lev] : di_a[lev];

    // Store number of fluid points found (0 if outside box)
    imp_ninterp(jj) = fluid;
  }

  /*//////////////////////////////////////////////////////////////////
    * \brief Computes interpolation weights for Image Points (IPs).
    *
    * This function calculates the trilinear (3D) or bilinear (2D) interpolation weights
    * (TODO : to be extended to higher order if needed according to N_InterP)
    * for a set of image points. It determines the stencil (surrounding grid points)
    * for each image point and computes the weights based on the relative position
    * within the cell.
    *
    * Key features:
    * - Handles both 2D and 3D cases via AMREX_SPACEDIM.
    * - Checks if stencil points are in the fluid or solid domain.
    * - Zeros out weights for solid points and renormalizes the remaining weights
    *   to ensure conservation.
    *
    * \tparam eorder_t Number of image points to process (template parameter).
    * \param[out] weights     Computed interpolation weights [image_idx][corner_idx].
    * \param[out] ip_ijk      Indices of the stencil points [image_idx][corner_idx][dim].
    * \param[in]  imp_xyz     Physical coordinates of the image points.
    * \param[in]  imp_ijk     Base grid index (bottom-left-back) of the cell containing the image point.
    * \param[in]  imp_ninterp Array indicating whether to interpolate each image point (1,2...) or not (0).
    * \param[in]  prob_lo     Physical coordinates of the domain lower bound.
    * \param[in]  dxyz        Grid spacing in each dimension.
    * \param[in]  ibFab       Marker array indicating fluid (0) or solid (1) state.
    *///////////////////////////////////////////////////////////////////
  template <int eorder_t, int iorder_t, int N_InterP = IntPow<iorder_t + 1, AMREX_SPACEDIM>::value>
  AMREX_FORCE_INLINE //AMREX_GPU_HOST_DEVICE
  void computeIPweights(
      Array2D<Real,0,eorder_t-1,0,N_InterP-1>&                     weights,
      Array3D< int,0,eorder_t-1,0,N_InterP-1,0,AMREX_SPACEDIM-1>&  ip_ijk,
      Array2D<Real,0,eorder_t-1,0,AMREX_SPACEDIM-1>&               imp_xyz,
      Array2D< int,0,eorder_t-1,0,AMREX_SPACEDIM-1>&               imp_ijk,
      Array1D< int,0,eorder_t-1>&                                  imp_ninterp,
      const GpuArray<Real, AMREX_SPACEDIM>&                        prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>&                        dxyz,
      Array4<uint8_t const>&                                       ibFab) const
  {
    // N_InterP Number of Interpolation Points (corners) per image point
    // NOW only works if N_InterP == 4 in 2D, or N_InterP == 8 in 3D (trilinear (3D) or bilinear (2D))

    // 1) Loop over all image points
    for (int iim = 0; iim < eorder_t; ++iim) {

      if (imp_ninterp(iim) < INTERP_THRESHOLD) {
        // If this image point is invalid, set weights to zero
        for (int corner = 0; corner < N_InterP; ++corner) {
      
          weights(iim, corner) = Real(0.0);
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = -99; // Invalid index
          }
        }
        continue; // Skip to next image point
      }

      int base_ijk[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {base_ijk[d] = imp_ijk(iim, d);}

      Real frac[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real lo = prob_lo[d] + Real(base_ijk[d] + 0.5_rt) * dxyz[d];
        frac[d] = (imp_xyz(iim, d) - lo) / dxyz[d];
      }

      int  sumfluid   = 0;
      Real sumweights = Real(0.0);

      // 3) Loop over all corners of the stencil cell
      for (int corner = 0; corner < N_InterP; ++corner) {

        int  ijk[AMREX_SPACEDIM];
        Real w = Real(1.0);

        // 3a) Determine corner offsets and weight contribution per dimension
        // Each 'corner' in [0 .. N_InterP-1] is interpreted as a binary code:
        //
        // 2D case (AMREX_SPACEDIM == 2, N_InterP = 4):
        //   corner | binary | (dx, dy)
        //   --------------------------
        //     0    |  00    | (0, 0)
        //     1    |  01    | (1, 0)
        //     2    |  10    | (0, 1)
        //     3    |  11    | (1, 1)
        //
        // 3D case (AMREX_SPACEDIM == 3, N_InterP = 8):
        //   corner | binary | (dx, dy, dz)
        //   --------------------------------
        //     0    |  000   | (0, 0, 0)
        //     1    |  001   | (1, 0, 0)
        //     2    |  010   | (0, 1, 0)
        //     3    |  011   | (1, 1, 0)
        //     4    |  100   | (0, 0, 1)
        //     5    |  101   | (1, 0, 1)
        //     6    |  110   | (0, 1, 1)
        //     7    |  111   | (1, 1, 1)
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int bit = (corner >> d) & 1; // 0 or 1 along this dimension

            ijk[d] = base_ijk[d] + bit;

            const Real fd = frac[d];
            w *= (bit ? fd : (Real(1.0) - fd));
        }

        // 3b) Store stencil indices
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = ijk[d];
        }

        // 4) Check if this stencil cell is fluid or solid
        int ii = ijk[0];
        int jj = ijk[1];
#if (AMREX_SPACEDIM == 3)
        int kk = ijk[2];
#else
        int kk = 0;   // in 2D, k-index is always 0 in Array4
#endif
        int fluid = !ibFab(ii, jj, kk, 0); // true if not solid cell

        weights(iim, corner) = w * Real(fluid);
        sumfluid   += fluid;
        sumweights += weights(iim, corner);
      }

      // keep the Assert as a sanity check for numerical issues.
      AMREX_ASSERT_WITH_MESSAGE(
          sumweights > Real(0.0),
          "computeIPweights: sum of raw weights is zero (unexpected numerical error).");

      // 5) Sanity check: ensure that the number of fluid stencil points matches imp_ninterp
      if (sumfluid != imp_ninterp(iim)) {
        amrex::Abort("computeIPweights: mismatch in fluid stencil count.");
      } 

      // 6) Renormalise weights so they sum to 1 over all N_InterP corners
      Real inv_sum = Real(1.0) / sumweights;
      Real check_sum = Real(0.0);

      for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) *= inv_sum;
          check_sum += weights(iim, corner);
      }

      AMREX_ASSERT_WITH_MESSAGE(
          std::abs(check_sum - Real(1.0)) < Real(1.0e-9),
          "Interpolation point weights do not sum to 1.0");

    } // end loop over image points
  }

  /*////////////////////////////////////////////////////////////////
  * \brief Interpolates primitive variables at image points using a given IP stencil and weights.
  *
  * For each image point (iim = 0..order_t-1), this routine accumulates contributions from its N_InterP
  * stencil corners into row (iim+2) of primsNormal, where rows 0 and 1 are reserved for the
  * ghost point and IB/surface state respectively.
  *
  * Row convention for primsNormal:
  *   0 : ghost point (GP)
  *   1 : IB/surface reference point
  *   2..(1+eorder_t) : image points along the normal
  *////////////////////////////////////////////////////////////////
  template <int eorder_t, int iorder_t, int N_InterP = IntPow<iorder_t + 1, AMREX_SPACEDIM>::value>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void interpolateIMs(
      const Array3D< int, 0, eorder_t - 1, 0, N_InterP - 1, 0, AMREX_SPACEDIM-1>&  imp_ip_ijk,
      const Array2D<Real, 0, eorder_t - 1, 0, N_InterP - 1>&                       imp_ipweights,
      const Array4<Real>&                                                          prims,
      Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM-1>&                           primsNormal) noexcept
  {
      // For each image point
      for (int iim = 0; iim < eorder_t; ++iim) {
          // For each interpolation point (corner) of its stencil
          for (int iip = 0; iip < N_InterP; ++iip) {

              const Real w = imp_ipweights(iim, iip);
              // If weight is zero (e.g. solid point or invalid marker), skip to avoid potential invalid memory access
              if (w == 0.0) continue;

              const int ii = imp_ip_ijk(iim, iip, 0);
              const int jj = imp_ip_ijk(iim, iip, 1);
          #if (AMREX_SPACEDIM == 3)
              const int kk = imp_ip_ijk(iim, iip, 2);
          #else
              // In 2D, k-index is always 0 in Array4
              const int kk = 0;
          #endif

              // Accumulate contributions for all primitive variables at row (iim+2)
              for (int n = 0; n < cls_t::NPRIM; ++n) {
                  primsNormal(iim + 2, n) += prims(ii, jj, kk, n) * w;
              } // end loop over primitive variables
          } // end loop over interpolation points (corners)
      } // end loop over image points
  }

   /*//////////////////////////////////////////////////////////////
   * \brief Extrapolates primitive variables from Image Points/Surface to the Ghost Point.
   *
   * This function solves a linear system (or uses simple linear interpolation) to determine
   * the state at the Ghost Point (primsNormal(0, :)) based on the states at the Surface (1)
   * and Image Points (2...).
   *
   * It assumes a polynomial profile P(d) = a + b*d + c*d^2 ... along the normal.
   *
   * \tparam eorder_t    Number of image points used (extrapolation order).
   * \param prims        Array of primitive variables along the normal.
   *                     Row 0: Ghost Point (Output)
   *                     Row 1: Surface Point (Input, from wall model)
   *                     Row 2..eorder_t+1: Image Points (Input, interpolated)
   * \param imp_ninterp  Number of interpolate points for each image point (-1 indicates outside valid range).
   * \param disGP        Distance from Surface to Ghost Point (> 0).
   * \param disIM        Array of distances from Surface to Image Points.
   *//////////////////////////////////////////////////////////////
  template <int eorder_t>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void extrapolate(Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM - 1>& prims, 
             const Array1D< int, 0, eorder_t - 1>& imp_ninterp,
             const Real disGP, const Array1D<Real, 0, eorder_t - 1>& disIM) const
  {
      // Determine effective order based on INTERP_THRESHOLD
      int eff_order = eorder_t;
      for (int k = 0; k < eorder_t; ++k) {
          if (imp_ninterp(k) < INTERP_THRESHOLD) {
              eff_order = k;
              break;
          }
      }

      // only extrapolate up to QLS (Last Species), skipping aux vars like QC, QG, QEINT.
      // Aux vars will be recomputed later via EOS (ensurePTYfillq).
      for (int n = 0; n <= cls_t::QLS; ++n) {
          
        
          // ----------------------------------------------------------------
          // CASE 2: Quadratic Extrapolation (eff_order >= 2)
          // Uses Surface Point (1), IM1 (2), IM2 (3)
          // ----------------------------------------------------------------
          if (eff_order >= 2) {
              if constexpr (eorder_t >= 2) {
                  Real u0 = prims(1, n);
                  Real u1 = prims(2, n);
                  Real u2 = prims(3, n);
                  
                  Real x1 = disIM(0);
                  Real x2 = disIM(1);
                  AMREX_ASSERT(x1 > 0 && x2 > 0);
                  
                  Real x = -disGP;
                  
                  Real L0 = (x - x1) * (x - x2) / (x1 * x2);
                  Real L1 = x * (x - x2) / (x1 * (x1 - x2));
                  Real L2 = x * (x - x1) / (x2 * (x2 - x1));
                  
                  prims(0, n) = u0 * L0 + u1 * L1 + u2 * L2;
              } else {
                  // Fallback if template eorder_t < 2 but eff_order >= 2 (impossible)
                  Real val_surf = prims(1, n);
                  Real val_im1  = prims(2, n);
                  Real d_im1    = disIM(0);
                  Real slope = (val_im1 - val_surf) / d_im1;
                  prims(0, n) = val_surf - slope * disGP;
              }
          }
          // ----------------------------------------------------------------
          // CASE 1: Linear Extrapolation (eff_order == 1)
          // Uses Surface Point (1) and First Image Point (2)
          // ----------------------------------------------------------------
          else if (eff_order == 1) {
              Real val_surf = prims(1, n);
              Real val_im1  = prims(2, n);
              Real d_im1    = disIM(0);

              Real slope = (val_im1 - val_surf) / d_im1;
              prims(0, n) = val_surf - slope * disGP;
          }
          // ----------------------------------------------------------------
          // CASE 0: Constant Extrapolation (eff_order == 0)
          // Uses Surface Point (1) only
          // ----------------------------------------------------------------
          else {
              prims(0, n) = prims(1, n);
          }
      }
  }

  /*//////////////////////////////////////////////////////////////
   * \brief Transforms velocity from global Cartesian coordinates to local surface-aligned coordinates.
   *
   * The local frame is defined by an orthonormal basis {n, t1, t2}. This routine performs a pure
   * rotation such that:
   *   (u_n, u_t1, u_t2)^T = R * (u_x, u_y, u_z)^T,
   * where R = [ n^T ; t1^T ; t2^T ].
   *
   * \param iip          Row index in primsNormal to be transformed.
   * \param primsNormal  Primitive-variable matrix storing velocities in-place.
   * \param norm         Surface unit normal vector.
   * \param tan1         First tangent vector.
   * \param tan2         Second tangent vector (only used in 3D).
   *////////////////////////////////////////////////////////////////
  template <int eorder_t>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void global2local(
      int iip,
      Array2D<Real,0,eorder_t+1,0,cls_t::NPRIM-1>& primsNormal,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2)
  {
    const Real ux = primsNormal(iip, cls_t::QU);
    const Real uy = primsNormal(iip, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    const Real uz = primsNormal(iip, cls_t::QW);
#else
    const Real uz = Real(0.0);
#endif

    // Normal component
    primsNormal(iip, cls_t::QU) =
        ux * norm(0) + uy * norm(1)
#if (AMREX_SPACEDIM == 3)
        + uz * norm(2)
#endif
        ;

    // Tangential component 1
    primsNormal(iip, cls_t::QV) =
        ux * tan1(0) + uy * tan1(1)
#if (AMREX_SPACEDIM == 3)
        + uz * tan1(2)
#endif
        ;

#if (AMREX_SPACEDIM == 3)
    // Tangential component 2 (only meaningful in 3D)
    primsNormal(iip, cls_t::QW) =
        ux * tan2(0) + uy * tan2(1) + uz * tan2(2);
#endif
  }

  /*////////////////////////////////////////////////////////////////
   * \brief Transforms velocity from local surface-aligned coordinates back to global Cartesian coordinates.
   *
   * This applies the inverse rotation of global2local:
   *   (u_x, u_y, u_z)^T = R^T * (u_n, u_t1, u_t2)^T,
   * where R = [ n^T ; t1^T ; t2^T ].
   *
   * \param jj           Row index in primsNormal to be transformed.
   * \param primsNormal  Primitive-variable matrix storing velocities in-place.
   * \param norm         Surface unit normal vector.
   * \param tan1         First tangent vector.
   * \param tan2         Second tangent vector (only used in 3D).
   *////////////////////////////////////////////////////////////////
  template <int eorder_t>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void local2global(
      int jj,
      Array2D<Real,0,eorder_t+1,0,cls_t::NPRIM-1>& primsNormal,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2)
  {
    const Real un  = primsNormal(jj, cls_t::QU);
    const Real ut1 = primsNormal(jj, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    const Real ut2 = primsNormal(jj, cls_t::QW);
#else
    const Real ut2 = Real(0.0);
#endif

    // Global x-component
    primsNormal(jj, cls_t::QU) =
        un  * norm(0) + ut1 * tan1(0)
#if (AMREX_SPACEDIM == 3)
      + ut2 * tan2(0)
#endif
      ;

    // Global y-component
    primsNormal(jj, cls_t::QV) =
        un  * norm(1) + ut1 * tan1(1)
#if (AMREX_SPACEDIM == 3)
      + ut2 * tan2(1)
#endif
      ;

#if (AMREX_SPACEDIM == 3)
    // Global z-component (only in 3D)
    primsNormal(jj, cls_t::QW) =
        un * norm(2) + ut1 * tan1(2) + ut2 * tan2(2);
#endif
  }

  /*////////////////////////////////////////////////////////////////
  * \brief Checks if the interpolation stencil is within the valid box.
  *
  * \param i, j, k   Indices of the bottom-left corner of the interpolation stencil.
  * \param bx        The box to check against.
  * \param lev       Current AMR level (for logging).
  * \param ipData    Reference to the ghost point data structure (used to extract context like geomIdx).
  * \param mode      Action to take on failure (Silent, Warn, or Abort).
  * 
  * \return true if the stencil is valid (inside bxg), false otherwise.
  /*///////////////////////////////////////////////////////////////
  template <typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool check_interpolation_stencil(int i, int j, int k, 
                                  const amrex::Box& bx, 
                                  int lev,
                                  const IPDATA& ipData, 
                                  CheckMode mode = CheckMode::Silent) const
  {
    // Check if the full stencil (2x2 in 2D, 2x2x2 in 3D) is contained in the box
    // We check the bottom-left (i,j,k) and top-right (i+1,j+1,k+1) corners.
    bool is_valid = bx.contains(amrex::IntVect(AMREX_D_DECL(i, j, k))) && 
                    bx.contains(amrex::IntVect(AMREX_D_DECL(i+1, j+1, k+1)));

    if (!is_valid) {
        if (mode == CheckMode::Silent) {
            return false;
        }

        // Extract context info from ipData. 
        int current_geom;
        int current_elem;
        
        if constexpr (GP_OR_SURF == 1) {
             current_geom = ipData.geomIdx.back();
             current_elem = ipData.elemIdx.back();
        } else {
             current_geom = ipData.geomIdx;
             current_elem = ipData.elemIdx;
        }

        Array1D<Real,0,AMREX_SPACEDIM-1> centroid;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            centroid(d) = SurfElem_a[current_elem].centroid[d];
        }

#if (AMREX_SPACEDIM == 3)
        std::printf("Interpolation stencil out of box bounds!\n"
                    "  Level: %d\n"
                    "  Stencil Base: (%d, %d, %d)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Face Centroid:  (%f, %f, %f)\n",
                    lev, i, j, k,
                    current_geom, current_elem,
                    centroid[0], centroid[1], centroid[2]);
#else
        std::printf("Interpolation stencil out of box bounds!\n"
                    "  Level: %d\n"
                    "  Stencil Base: (%d, %d)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Face Centroid:  (%f, %f)\n",
                    lev, i, j,
                    current_geom, current_elem,
                    centroid[0], centroid[1]);
#endif
        std::fflush(stdout);

        if (mode == CheckMode::Warn) {
            amrex::Warning("Interpolation stencil out of box bounds!");
        } else if (mode == CheckMode::Abort) {
            amrex::Abort("Interpolation stencil out of box bounds!");
        }
    }

    return is_valid;
  }

  /*////////////////////////////////////////////////////////////////
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
  */////////////////////////////////////////////////////////////////

  /*////////////////////////////////////////////////////////////////
   * \brief Reads geometry files and initializes CGAL data structures.
   *
   * This function performs the following tasks:
   * 1. Reads geometry filenames from the input file (parameter: `ib.filename`).
   * 2. Clears and resizes internal containers for geometry data.
   * 3. Iterates through each geometry file:
   *    - For 2D: Reads a polygon, constructs an AABB tree for distance queries, and sets up an in/out testing functor.
   *    - For 3D: Reads a polygon mesh, ensures it is triangulated and outward-oriented, constructs an AABB tree, and sets up an in/out testing functor.
   * 4. Computes and stores geometric properties (centroids, normals, areas/lengths) for all faces/edges into flattened arrays (`SurfElem_a`, `LocalFrame_a`) for efficient GPU access.
   * 5. Builds a mapping (`IdxMap_a`) from CGAL primitive IDs to linear indices.
   *
   * \note This function handles both 2D (Polygon) and 3D (Polyhedron) geometries based on `AMREX_SPACEDIM`.
   */////////////////////////////////////////////////////////////////
  void read_geom()
  {
    ParmParse pp;
    Vector<std::string> files_a;
    bool plot_surf = false;

    pp.getarr("ib.filename", files_a);
    pp.query("ib.plot_surf", plot_surf);

    // Basic validation
    if (files_a.empty()) {
      amrex::Warning("ib.filename is empty: no geometry files provided");
    }

    // Safely delete any previously allocated pointers
    for (auto* p : tree_pa)   { delete p; }
    for (auto* p : inout_fa)  { delete p; }
    tree_pa.clear();
    inout_fa.clear();

    // Resize all containers to match new geometry count
    this->ngeom = static_cast<int>(files_a.size());
    this->geom_a.resize(this->ngeom);
    this->tree_pa.resize(this->ngeom);
    this->LocalFrame_a.clear();
    this->SurfElem_a.clear();
    this->geom_offsets.resize(this->ngeom);
    this->IdxMap_a.resize(this->ngeom);
    this->inout_fa.resize(this->ngeom);
    this->ntotalfaces = 0;

    // Determine minimum dx across all levels for setting polygon tolerance
    Real min_dx = std::numeric_limits<Real>::max();
    for (int lev = 0; lev < dx_a.size(); ++lev) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            min_dx = std::min(min_dx, dx_a[lev][d]);
        }
    }
    min_dx *= Real(0.5);

    for (int i = 0; i < ngeom; i++) {
      Print() << "----------------------------------" << std::endl;

#if (AMREX_SPACEDIM == 2)

      // Read 2D polygon from file
      if (!read_polygon_2d(files_a[i], geom_a[i], min_dx)) {
        amrex::Abort(std::string("Invalid 2D geometry filename: ") + files_a[i]);
      }
      Print() << "Geometry (i=" << i << ") " << files_a[i] << " read" << std::endl;
      Print() << "Number of vertices in polygon: "<< geom_a[i].size() << "\n";
      
      // constructs AABB tree and computes internal KD-tree
      // data structure to accelerate distance queries
      // Note: geom_a[i] is Polygon2D, which holds a persistent polygon 
      // and a vector of segments which we use its iterators directly.
      tree_pa[i] = new Tree(geom_a[i].edges_begin(), geom_a[i].edges_end());
      tree_pa[i]->build();
      tree_pa[i]->accelerate_distance_queries();
      Print() << "AABB tree constructed" << std::endl;

      // in/out functor: set polygon pointer
      inout_fa[i] = new inside_t(geom_a[i]);
      Print() << "2D in/out testing functor constructed for polygon " << files_a[i] << "\n";

      // "face count" in 2D can be understood as the number of edges
      ntotalfaces += geom_a[i].size();

#elif (AMREX_SPACEDIM == 3)

      // Read 3D Polyhedron mesh from file
      if (!PMP::IO::read_polygon_mesh(files_a[i], geom_a[i])) {
        amrex::Abort(std::string("Invalid geometry filename: ") + files_a[i]);
      }
      Print() << "Geometry (i=" << i << ") " << files_a[i] << " read" << std::endl;
      Print() << "Is geometry only made of triangles? "
              << geom_a[i].is_pure_triangle() << std::endl;
      Print() << "Number of facets " << geom_a[i].size_of_facets() << std::endl;

      if (!geom_a[i].is_pure_triangle()) {
        amrex::Abort(std::string("Geometry is not pure triangle: ") + files_a[i] +
                      ". Please provide triangulated mesh.");
      }

      // Ensure outward orientation for consistent normal vectors
      if (!PMP::is_outward_oriented(geom_a[i])) {
          PMP::reverse_face_orientations(geom_a[i]);
          Print() << "Info: Reversed face orientations to be outward for " 
                  << files_a[i] << "\n";
      }

      // constructs AABB tree and computes internal KD-tree
      // data structure to accelerate distance queries
      tree_pa[i] =
        new Tree(faces(geom_a[i]).first, faces(geom_a[i]).second, geom_a[i]);
      tree_pa[i]->build(); 
      tree_pa[i]->accelerate_distance_queries();
      Print() << "AABB tree constructed" << std::endl;

      // plane class also computes orthogonal direction to the face. However, the
      // orthogonal vector is not normalised.
      std::for_each(geom_a[i].facets_begin(), geom_a[i].facets_end(),
                    compute_plane_equations);
      Print() << "Plane equations per face computed" << std::endl;

      // make inside/outside function for current geometry only
      //inout_fa[i] = new inside_t(geom_a[i]);
      inout_fa[i] = new inside_t(*tree_pa[i]);
      Print() << "In out testing function constructed for geometry " << files_a[i] << "\n";
      
      //compute total face count across all geometries
      ntotalfaces += geom_a[i].size_of_facets();
    
#endif

    // Store offset for current geometry
    this->geom_offsets[i] = static_cast<int>(this->LocalFrame_a.size());

    // Build geometry cache:
    // This step flattens the geometry data (faces/edges) into linear arrays (SurfElem_a, LocalFrame_a)
    // and creates a mapping (IdxMap_a) from CGAL's internal IDs to these linear indices.
    // This allows for efficient O(1) access to geometric properties on the GPU using a simple integer index.
    build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a, IdxMap_a[i], this->geom_offsets[i], i);
    } // end loop over geometries

    // Sanity check: verify that the cache grew by the expected amount
    if ((static_cast<int>(SurfElem_a.size()) != ntotalfaces) || (static_cast<int>(LocalFrame_a.size()) != ntotalfaces)) {
        amrex::Print() << "Error: Mismatch in geometry cache size for geom " << "\n"
                       << "  Expected total faces     : " << ntotalfaces << "\n"
                       << "  Actual SurfElem_a size   : " << SurfElem_a.size() << "\n"
                       << "  Actual LocalFrame_a size : " << LocalFrame_a.size() << "\n";
        amrex::Abort("build_geometry_cache failed to add the correct number of elements.");
    }

    // Check for geometry consistency (no intersections, no containment)
    check_ibm_geometry_consistency(ngeom, geom_a.data(), inout_fa.data(), files_a.data());

    Print() << "----------------------------------" << std::endl;
    Print() << "----------------------------------" << std::endl;
    
    if (plot_surf) {

      Print() << "Total number of faces across all geometries: " << ntotalfaces << std::endl;
      // Initialize surfdata container to this total; safe defaults
      surfdata_a.resize(ntotalfaces);
    }
  } // end read_geom


}; // end class eib_t
#endif // EIB_H_