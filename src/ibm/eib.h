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

static constexpr int  NIP = (1 << AMREX_SPACEDIM);  // number of interpolation points for each image point
static constexpr int IDIM = AMREX_SPACEDIM - 1;     // index dimension

///----------------------------------------------------------------------------
/// \brief Class to store ghost point arrays
/// \param eorder_tparm Number of image points (integer) 
///
template <int eorder_tparm>
struct gpData_t {
  // CPU only attributes
  gpData_t() : ngps(0) {}
  int ngps;           

  // closest surface point (ib point) and face ID
  //Vector<Point_and_primitive_id> closest_cgal;  

  // GPU/CPU attributes
  // Ghost point data
  Gpu::ManagedVector<Array1D< int, 0, IDIM>> gp_ijk;                // Ghost point indices
  Gpu::ManagedVector<Array1D<Real, 0, IDIM>> ib_xyz;                // IB point coordinates
  Gpu::ManagedVector<Real> disGP;                                   // Distance from IB point to ghost point
  Gpu::ManagedVector<int> geomIdx;                                  // Geometry index
  Gpu::ManagedVector<int> faceIdx;                                  // face/edge element index
  //Gpu::ManagedVector<Array2D<Real, 0, IDIM, 0, IDIM>> localframe;   // Local orthonormal frame matrix: columes = {normal, tangent1(, tangent2)}
  //Gpu::ManagedVector<Array1D<Real, 0, IDIM>> normal;                // Surface normal vectors
  //Gpu::ManagedVector<Array1D<Real, 0, IDIM>> tangent1;              // First tangent vectors
  //Gpu::ManagedVector<Array1D<Real, 0, IDIM>> tangent2;              // Second tangent vectors
  
  // Image point data arrays
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm - 1>> disIM;     // Distance from IB point to image points
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>> imp_xyz;  // Image point coordinates
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>> imp_ijk;  // Image point indices
  
  // Interpolation data for image points
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm - 1, 0, NIP - 1, 0, IDIM>> imp_ip_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, NIP - 1>> imp_ipweights;

};

///----------------------------------------------------------------------------
/// \brief Class to store surface data
/// Per-face surface data container used for reconstruction and output
/// \param iorder_tparm Number of image points used for surfdata reconstruction
///
template <int iorder_tparm>
struct surfData_t{
  // CPU only attributes        
  surfData_t() 
    : ifab(-1),
      lev(-1),
      iface(-1),  
      igeom(-1),
      rank(-99),
      pressure(0.0),
      tau1(0.0),
      tau2(0.0),
      temperature(0.0),
      dTdn(0.0),
      pointfound(false)
  {}

  // Image point data (per face)
  Array2D<Real, 0, iorder_tparm - 1, 0, IDIM> imp_xyz;              // Physical-space coordinates of image points placed along the outward normal
  Array2D< int, 0, iorder_tparm - 1, 0, IDIM> imp_ijk;              // Index of the “bottom-left” grid cell associated with each image point
  Array2D<Real, 0, iorder_tparm - 1, 0, NIP - 1> ipweights;         // Trilinear interpolation weights for the 8-point stencil of each image point
  Array3D< int, 0, iorder_tparm - 1, 0, NIP - 1, 0, IDIM> ip_ijk;   // Indices of the 8-point interpolation stencil for each image point
  Array1D<Real, 0, iorder_tparm - 1> o_dis;                         // Normal distances from the surface (IB point) to each image point (per j)

  // Surface geometry (per face)
  // IB point (face reference point) coordinates, aera and Local orthonormal frame matrix: columes = {normal, tangent1, tangent2}
  // Array1D<Real, 0, IDIM> ib_xyz;
  // Array2D<Real, 0, IDIM, 0, IDIM> localframe;
  // amrex::Real aera;    // face area

  // ifab: local FAB index on this MPI rank (MFIter::LocalIndex)
  // lev: AMR level; iface: global face index across all geometries
  // igeom: geometry index; rank: owning MPI rank at time of packing
  int  ifab,lev,iface,igeom,rank; 

  // Surface fields (per face)
  amrex::Real temperature, dTdn;             // reconstructed temperature and grad(T)·n
  amrex::Real pressure, tau1, tau2;          // reconstructed local surface force
  //Array1D<Real, 0, IDIM> local_velocity    // surface local velocity

  // Whether this face has been located/mapped to a FAB in this pass
  bool pointfound;

};


///main class ----------------------------------------------------------------------

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
  static constexpr int  iorder_tparm = param::interp_order; // number of image points used for surfdata reconstruction
  static constexpr int  eorder_tparm = param::extrap_order; // number of image points used for ghost point extrapolation
  static constexpr Real cim          = param::alpha;

  // MultiFabs
  Amr* amr_p;                                               // pointer to Amr class instance
  Vector<IBMultiFab<uint8_t, gpData_t<eorder_tparm>>*> bmf_a;  // Immersed boundary MultiFab array (bool multifab array)

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
  Vector<std::map<PrimitiveID, int>> idxmap_a;              // face/edge element index per geometry
  Vector<inside_t*> inout_fa;                               // in out testing function per geometry

  // surface related data
  int ntotalfaces = 0;                          // number of faces/edges across all geometries
  Vector<surfData_t<iorder_tparm>> surfdata_a;  // surface/edge data (independent of number of geometries)
  Vector<Vector<Vector<int>>> intfaces_in_fab;  // faces integers per fab and per level [lev][ifab][face_index]

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
    idxmap_a.clear();

    surfdata_a.clear();
    intfaces_in_fab.clear();  
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

  // create IBMultiFabs at a level and store pointers to it
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    bmf_a[lev] =
      new IBMultiFab<uint8_t, gpData_t<eorder_tparm>>(bxa, dm, 2, cls_t::NGHOST);
      // lsMFa[lev].define(bxa, dm, 1, NGHOST_IB);
  }

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
   * the immersed boundary geometry using CGAL. It populates the `ibMarkers` MultiFab where:
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
      {
        auto& gpData = ibFab.gpData;
        gpData.ngps = 0;
        gpData.gp_ijk.clear();
        gpData.disGP.clear();
        gpData.disIM.clear();
        gpData.imp_ipweights.clear();
        gpData.imp_ip_ijk.clear();
        gpData.ib_xyz.clear();
        gpData.geomIdx.clear();
        gpData.faceIdx.clear();
        gpData.imp_xyz.clear();
        gpData.imp_ijk.clear();
        //gpData.closest_cgal.clear();
        //gpData.localframe.clear();
        //gpData.normal.clear();
        //gpData.tangent1.clear();
        //gpData.tangent2.clear();
      }

      // compute sld markers (including ghost points) - cannot use ParallelFor - CGAL call causes problems
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
      // TODO: move to GPU.
      ibFab.gpData.ngps = 0;
      amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
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
          ibMarkers(i, j, k, 1) = ghost;
          ibFab.gpData.ngps += ghost;

          if (ghost) {
            // store GP index
            ibFab.gpData.gp_ijk.push_back(make_vec<int>(i, j, k));
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
          } else {
            ibMarkers(i, j, k, 1) = static_cast<uint8_t>(0);
          }
        } // end if solid
      }); // end LoopOnCpu for ghost markers

#ifdef EIB_ENABLE_SHRINK
  // After constructing the ghost point lists for this FAB, optionally
  // reclaim excess capacity from previous high-water allocations.
  {
    auto& gpData = ibFab.gpData;
    EIB_SHRINK(gpData.gp_ijk);
    EIB_SHRINK(gpData.disGP);
    EIB_SHRINK(gpData.disIM);
    EIB_SHRINK(gpData.imp_ipweights);
    EIB_SHRINK(gpData.imp_ip_ijk);
    //EIB_SHRINK(gpData.localframe);
    //EIB_SHRINK(gpData.normal);
    //EIB_SHRINK(gpData.tangent1);
    //EIB_SHRINK(gpData.tangent2);
    EIB_SHRINK(gpData.ib_xyz);
    EIB_SHRINK(gpData.geomIdx);
    EIB_SHRINK(gpData.faceIdx);
    //EIB_SHRINK(gpData.closest_cgal);
    EIB_SHRINK(gpData.imp_xyz);
    EIB_SHRINK(gpData.imp_ijk);
  }
#endif
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
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    const Box& bx = mfi.tilebox();
    auto const ibMarkers = mfab.array(mfi);  // boolean array

    // we need a CPU loop here (cannot be GPU loop) as CGAL tree seach for
    // closest element to a point needs to be called. instead of looping through
    // previously indexed gps, we loop through the whole ghost point field as it
    // is available on GPU and CPU at all times. Unlike the gp indexes, which
    // are only stored on GPU memory. Array1D<int,0,AMREX_SPACEDIM-1>& idx =
    // ibFab.gpData.gp_ijk[ii];

    amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
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
        
        // find and store geometery index for this GP.
        // Since this is a ghost point, it must be a solid point, 
        // so ibMarkers(i,j,k,0) stores (geometry_index + 1).
        int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;

        AMREX_ASSERT_WITH_MESSAGE(geomIdx >= 0 && geomIdx < ngeom, 
                "Invalid geometry index in initialiseGPs");
        ibFab.gpData.geomIdx.push_back(geomIdx);

        // closest surface/edge point and surface/edge --------------------------
        Point_and_primitive_id closest_elem =
            tree_pa[geomIdx]->closest_point_and_primitive(gp);

        // map PrimitiveID to integer index and store that.
        PrimitiveID elm = closest_elem.second;
        int f_idx = idxmap_a[geomIdx].at(elm);
        ibFab.gpData.faceIdx.push_back(f_idx);

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
        const auto& frame = LocalFrame_a[f_idx];

        // IM points -------------------------------------------
        Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
        Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;
        Array1D<Real, 0, eorder_tparm - 1> disIM;
        
        // find image point and the bottom left point closest to the image
        // point
        //  In 2D, same idea in 3D.
        //     i,j+1 (2) ---------------------     i+1,j+1 (3)
        //     |                                  |
        //     |         P                        |
        //     |                                  |
        //     |                                  |
        //     i,j  (1) ----------------------      i+1,j  (4)

        for (int jj = 0; jj <= eorder_tparm - 1; jj++) {
          for (int kk = 0; kk < AMREX_SPACEDIM; kk++) {
            imp_xyz(jj, kk) = cp[kk] + Real(jj + 1) * di_a[lev] * frame.normal[kk];
            imp_ijk(jj, kk) = static_cast<int>(
              std::floor((imp_xyz(jj, kk) - prob_lo[kk]) / dx_a[lev][kk] - 0.5_rt));
          }

          AMREX_ASSERT_WITH_MESSAGE(
              bxg.contains(AMREX_D_TERM(imp_ijk(jj, 0),
                                        imp_ijk(jj, 1),
                                        imp_ijk(jj, 2))) &&
              bxg.contains(AMREX_D_TERM(imp_ijk(jj, 0) + 1,
                                        imp_ijk(jj, 1) + 1,
                                        imp_ijk(jj, 2) + 1)),
              "Interpolation point outside fab");

          disIM(jj) = Real(jj + 1) * di_a[lev];
        }

        // store
        gpData.disIM.push_back(disIM);
        gpData.imp_xyz.push_back(imp_xyz);
        gpData.imp_ijk.push_back(imp_ijk);

        // Interpolation points' (ips) weights for each image point
        Array2D<Real, 0, eorder_tparm - 1 , 0, NIP -1 > ipweights;
        Array3D< int, 0, eorder_tparm - 1 , 0, NIP -1, 0, AMREX_SPACEDIM - 1> ip_ijk;
        
        computeIPweights<eorder_tparm>(
                ipweights, ip_ijk, imp_xyz, imp_ijk,
                prob_lo, dx_a[lev], ibMarkers, lev);
        
        // *store*
        gpData.imp_ipweights.push_back(ipweights);
        gpData.imp_ip_ijk.push_back(ip_ijk);

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
                int& lev){
  
  auto& mfab               = *bmf_a[lev];
  const auto& ibFab        = mfab.get(mfi);

  // Ghost-point related data (geometry + interpolation)
  auto const gp_ijk        = ibFab.gpData.gp_ijk.data();
  auto const imp_ijk       = ibFab.gpData.imp_ijk.data();
  auto const imp_ipweights = ibFab.gpData.imp_ipweights.data();
  auto const imp_ip_ijk    = ibFab.gpData.imp_ip_ijk.data();
  auto const disGP         = ibFab.gpData.disGP.data();
  auto const disIM         = ibFab.gpData.disIM.data();
  auto const localframe    = ibFab.gpData.localframe.data();
  auto const ib_xyz        = ibFab.gpData.ib_xyz.data();

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
    // 1) Reconstruct local orthonormal frame from stored 2D/3D localframe
    //     localframe(ii)(row, col): row = {0: n, 1: t1, (2: t2 in 3D)}
    // --------------------------------------------------------------------
    Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        nvec(d)  = localframe[ii](0,d);
        t1vec(d) = localframe[ii](1,d);
#if (AMREX_SPACEDIM == 3)
        t2vec(d) = localframe[ii](2,d);
#else
        // In 2D, t2 is not used; keep as zero for interface compatibility.
        t2vec(d) = Real(0.0);
#endif
      }

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
      copy->interpolateIMs<eorder_tparm>(imp_ip_ijk[ii], imp_ipweights[ii], prims0, primsNormal);

      // 4) Transform velocities at image points (> 1) to local frame
      for (int iip = 2; iip < 2 + eorder_tparm; ++iip) {
          copy->global2local<eorder_tparm>(iip, primsNormal, nvec, t1vec, t2vec);
      }

      // 5) Apply wall model at IB surface to set surface states (u, P, T, Y, ...)
      wallmodel::compute_surfIB(ib_xyz[ii], nvec, primsNormal, cls);

      // 6) Extrapolate from surface/image points back to ghost point along n
      copy->extrapolate(primsNormal, disGP[ii], disIM[ii]);

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
      int i = gp_ijk;
      int j = gp_ijk;
#if (AMREX_SPACEDIM == 3)
      int k = gp_ijk;
#else
      int k = 0;  // in 2D, k-index is always 0
#endif

      for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(i,j,k,n) = Q[n];
      }
  }); // end ParallelFor over ghost points
  }



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief compute surface indexes and store them (core, fab, lev)
// is independet of geomenries, it will store faces 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_surface_index(int lev) {

  int myrank = amrex::ParallelDescriptor::MyProc(); 
  amrex::Print() << " Compute Surface Index at LEVEL " << lev << std::endl;

  auto& mfab = *bmf_a[lev];
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

  int faces_notfound = 0;
  int iface = -1;    // face counter

  // Initialization policy:
  //  - If empty: allocate and initialize all faces.
  //  - If size mismatch (e.g., geometry changed): reinitialize to safe defaults.
  if (surfdata_a.empty()) {

    surfdata_a.resize(ntotalfaces);
    for (int f = 0; f < ntotalfaces; ++f) {
      surfdata_a[f].pointfound = false;
      surfdata_a[f].rank = -99;
    }
  } else if (static_cast<int>(surfdata_a.size()) != ntotalfaces) {

     
  }

  // resize vectors depending number of fabs to store face information
  int nfab   = mfab.local_size(); 
  int maxlev = amr_p->maxLevel();
  
  if (intfaces_in_fab.size() <= lev)
    intfaces_in_fab.resize(lev + 1); // Add space for current level if needed
  // If level exists, clear each fab’s vector before filling
  intfaces_in_fab[lev].resize(nfab);
  for (int ifab = 0; ifab < nfab; ++ifab) {
    intfaces_in_fab[lev][ifab].clear(); // Important if reusing
  }  

  // Initialize the number of faces in each fab at this level  
  //for (int ifab = 0; ifab < nfab; ++ifab)  nfaces_infab_inlevel[lev][ifab] = 0;
  
  for (int ii = 0; ii < ngeom; ii++) { 
    const GeomType& mesh  = geom_a[ii];

    // loop over the faces of geometry, iface is global counter of faces
    for (auto fd : faces(mesh)) {
      iface++;

      // create a surface point if first level (otherwise get the surfdata)

      surfData_t<iorder_tparm> surf_dat; 
      
      if (lev == 0) {
        surf_dat.pointfound = false;
        surf_dat.rank = -99;
      } 
      else {
        if (iface < static_cast<int>(surfdata_a.size())) {
        auto& surf_dat = surfdata_a[iface]; 
        } 
        else 
        {
        amrex::Abort("iface out of bounds for surfdata_a during lev > 0");
        }
      }

      // Use pre-computed data from flattened arrays
      const auto& frame = LocalFrame_a[iface];
      const auto& elem  = SurfElem_a[iface];

      Array1D<Real, 0, AMREX_SPACEDIM - 1> norm;
      Array1D<Real, 0, AMREX_SPACEDIM - 1> ib_xyz;
      for(int d=0; d<AMREX_SPACEDIM; ++d) {
          norm(d) = frame.normal[d];
          ib_xyz(d) = elem.centroid[d];
      }

      // Reconstruct CGAL Point for compatibility with existing code below
#if (AMREX_SPACEDIM == 3)
      Point face_center(elem.centroid[0], elem.centroid[1], elem.centroid[2]);
#else
      Point face_center(elem.centroid[0], elem.centroid[1]);
#endif
      
      // find closest  point i,j,k  to the face_center
      const int i1 = int((face_center.x() - prob_lo[0]) / dx_a[lev][0] - 0.5_rt);
      const int j1 = int((face_center.y() - prob_lo[1]) / dx_a[lev][1] - 0.5_rt);
      const int k1 = int((face_center.z() - prob_lo[2]) / dx_a[lev][2] - 0.5_rt);
        
      const bool is_inside = amr_p->Geom(lev).Domain().contains(i1, j1, k1);
        
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

          // Use pre-computed area
          double area = elem.size;
                   
          // mirror point 
          int i =  imp_ijk(jj, 0);int j =  imp_ijk(jj, 1);int k = imp_ijk(jj, 2);
                            

          // locate the mirror points, loop over fabs
          bool face_present_core = false; 
          for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {             
    
            //const int ifab= mfi.index();
            const int ifab = mfi.LocalIndex();  // Use this instead
            auto& ibFab = mfab.get(mfi);                
            const Box& bx = mfi.tilebox();
            const Box& bxg = mfi.growntilebox(cls_t::NGHOST);

            // mirror point located
            if (bx.contains(i, j, k)) {

              auto const ibMarkers = mfab.array(mfi);
             

              // if point not valid as mirror, move mirror points    
              if (!valid_mirror(i, j, k, ibMarkers)) {
                                
                int jj = 0;int counter = 0;
                do {  
                  Real o_dis = 0.0;            
                  for (int kk = 0; kk < AMREX_SPACEDIM; kk++) {
                    imp_xyz(jj, kk) += di_a[lev] * norm(kk);
                    imp_ijk(jj, kk) = floor((imp_xyz(jj, kk) - prob_lo[kk]) / dx_a[lev][kk] - 0.5_rt);        
                    o_dis += (face_center[kk] - imp_xyz(jj, kk))*(face_center[kk] - imp_xyz(jj, kk));          
                  }                  
                  o_dis = 1.0/sqrt(o_dis);
                  i =  imp_ijk(jj, 0); j =  imp_ijk(jj, 1);k = imp_ijk(jj, 2);
                  counter++;
                } while ((!valid_mirror(i, j, k, ibMarkers)) && (counter < 10) && (bxg.contains(i, j, k)) );
                if (counter ==5){
                  amrex::Error(" Valid mirror point not found after 5 attempts \n");
                }
                if (!bxg.contains(i, j, k)) {
                  amrex::Error(" Valid mirror point not found in the box (including ghost) !! \n");
                }
              }

              // Interpolation points' (ips) weights for each image point
              Array2D<Real, 0, eorder_tparm - 1, 0, 7> ipweights;
              Array3D< int, 0, eorder_tparm - 1, 0, 7, 0, AMREX_SPACEDIM - 1> ip_ijk;
              computeIPweights<eorder_tparm>(
                    ipweights, ip_ijk, imp_xyz, imp_ijk,
                    prob_lo, dx_a[lev], ibMarkers, lev);

              face_present_core = true; 
             
              // store face information in the surf_dat
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
              surf_dat.rank  = myrank; // store rank of the processor where the face is located  
                                                          
              // if point not found previously in that level update variable
              if (!surf_dat.pointfound)
              {               
                surf_dat.pointfound = true;                                 
              }                
              intfaces_in_fab[lev][ifab].push_back(iface);
              // store or update data
              surfdata_a[iface] = surf_dat;  

              //nfaces_infab_inlevel[lev][ifab]++; // increment number of faces in fab at this level
                                            
            } // end contains
                          
          } // end loop over fabs
        if (!face_present_core) faces_notfound++;
        }
        else // face outside domain
        {          
          faces_notfound++;
        }
    }  // end loop faces

  } // end loop geometries


  // amrex::Print() << " Faces not found  (outside domain/lost) " << faces_notfound << " out of " << ntotalfaces << std::endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief function to check if a mirror points is valid
// returns true if at least 3 out of 8 surrounding points are fluid
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool valid_mirror(int i, int j, int k, const amrex::Array4<const bool>& ibMarkers) {
  
  int fluid_count = 0;

  for (int di = 0; di <= 1; ++di) {
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk) {
        int ii = i + di;
        int jj = j + dj;
        int kk = k + dk;
          if (!ibMarkers(ii, jj, kk, 0)) {
            ++fluid_count;
          }
        }
      }
    }
    return fluid_count >= 3;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief compute surface properties for each surface face
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_surface_props(MultiFab& stateprops,const cls_t* cls,int lev) {
  

  int myrank = amrex::ParallelDescriptor::MyProc(); 

  amrex::Print() << " Compute Surface Properties at LEVEL" << lev << std::endl;

  // loop over mfi
  auto& mfab = *bmf_a[lev];
  GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();
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
      
      const Polyhedron& mesh   = geom_a[ii];
                        
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
          int i1 = ip_ijk(0,0,0);int j1 = ip_ijk(0,0,1); int k1 = ip_ijk(0,0,2);          
          if (! bxg.contains(i1, j1, k1) ){
          
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

  /*//////////////////////////////////////////////////////////////////////////
    * \brief Computes interpolation weights for Image Points (IPs).
    *
    * This function calculates the trilinear (3D) or bilinear (2D) interpolation weights
    * for a set of image points. It determines the stencil (surrounding grid points)
    * for each image point and computes the weights based on the relative position
    * within the cell.
    *
    * Key features:
    * - Handles both 2D and 3D cases via AMREX_SPACEDIM.
    * - Checks if stencil points are in the fluid or solid domain using `ibFab`.
    * - Zeros out weights for solid points and renormalizes the remaining weights
    *   to ensure conservation.
    *
    * \tparam order_t Number of image points to process (template parameter).
    * \param[out] weights   Computed interpolation weights [image_idx][corner_idx].
    * \param[out] ip_ijk    Indices of the stencil points [image_idx][corner_idx][dim].
    * \param[in]  imp_xyz   Physical coordinates of the image points.
    * \param[in]  imp_ijk   Base grid index (bottom-left-back) of the cell containing the image point.
    * \param[in]  prob_lo   Physical coordinates of the domain lower bound.
    * \param[in]  dxyz      Grid spacing in each dimension.
    * \param[in]  ibFab     Marker array indicating fluid (0) or solid (1) state.
    *//////////////////////////////////////////////////////////////////////////////
  template <int order_t>
  AMREX_FORCE_INLINE
  void computeIPweights(
      Array2D<Real,0,order_t-1,0,NIP-1>&                    weights,
      Array3D< int,0,order_t-1,0,NIP-1,0,AMREX_SPACEDIM-1>& ip_ijk,
      Array2D<Real,0,order_t-1,0,AMREX_SPACEDIM-1>&         imp_xyz,
      Array2D< int,0,order_t-1,0,AMREX_SPACEDIM-1>&         imp_ijk,
      const GpuArray<Real, AMREX_SPACEDIM>&                 prob_lo,
      const GpuArray<Real, AMREX_SPACEDIM>&                 dxyz,
      const Array4<uint8_t>&                                ibFab,
      int                                                   lev)
  {
      // 1) Loop over all image points
    for (int iim = 0; iim < order_t; ++iim) {

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
      for (int corner = 0; corner < NIP; ++corner) {

        int  ijk[AMREX_SPACEDIM];
        Real w = Real(1.0);

        // 3a) Determine corner offsets and weight contribution per dimension
        // Each 'corner' in [0 .. NIP-1] is interpreted as a binary code:
        //
        // 2D case (AMREX_SPACEDIM == 2, NIP = 4):
        //   corner | binary | (dx, dy)
        //   --------------------------
        //     0    |  00    | (0, 0)
        //     1    |  01    | (1, 0)
        //     2    |  10    | (0, 1)
        //     3    |  11    | (1, 1)
        //
        // 3D case (AMREX_SPACEDIM == 3, NIP = 8):
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
        const int fluid = !ibFab(ii, jj, kk, 0); // true if not IB cell

        weights(iim, corner) = w * Real(fluid);
        sumfluid   += fluid;
        sumweights += weights(iim, corner);
      }

      // 5) Basic sanity: at least 2 fluid points, and non-zero total weight
      if (sumfluid < 2) {
      amrex::Print()
      << "Warning: Less than 2 interpolation points are fluid points\n"
      << "  level (lev)       = " << lev << "\n"
      << "  image index (iim) = " << iim << "\n"
      << "  base cell index   = ("
      << base_ijk[0] << ", "
      << base_ijk[1]
#if (AMREX_SPACEDIM == 3)
      << ", " << base_ijk[2]
#endif
      << ")\n"
      << "  image point xyz   = ("
      << imp_xyz(iim,0) << ", "
      << imp_xyz(iim,1)
#if (AMREX_SPACEDIM == 3)
      << ", " << imp_xyz(iim,2)
#endif
      << ")\n"
      << "  sumfluid          = " << sumfluid << "\n";

    amrex::Warning("Less than 2 interpolation points are fluid points");
    }

      AMREX_ASSERT_WITH_MESSAGE(
          sumweights > Real(0.0),
          "computeIPweights: sum of raw weights is zero.");

      // 6) Renormalise weights so they sum to 1 over all NIP corners
      Real inv_sum = Real(1.0) / sumweights;
      Real check_sum = Real(0.0);

      for (int corner = 0; corner < NIP; ++corner) {
          weights(iim, corner) *= inv_sum;
          check_sum += weights(iim, corner);
      }

      AMREX_ASSERT_WITH_MESSAGE(
          std::abs(check_sum - Real(1.0)) < Real(1.0e-9),
          "Interpolation point weights do not sum to 1.0");
    } // end loop over image points
  }

  /*//////////////////////////////////////////////////////////////////////////
  * \brief Interpolates primitive variables at image points using a given IP stencil and weights.
  *
  * For each image point (iim = 0..order_t-1), this routine accumulates contributions from its NIP
  * stencil corners into row (iim+2) of primsNormal, where rows 0 and 1 are reserved for the
  * ghost point and IB/surface state respectively.
  *
  * Row convention for primsNormal:
  *   0 : ghost point (GP)
  *   1 : IB/surface reference point
  *   2..(1+order_t) : image points along the normal
  *////////////////////////////////////////////////////////////////////////
  template <int order_t>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void interpolateIMs(
      const Array3D< int, 0, order_t-1, 0, NIP-1, 0, AMREX_SPACEDIM-1>& imp_ip_ijk,
      const Array2D<Real, 0, order_t-1, 0, NIP-1>&                      ipweights,
      const Array4<Real>&                                               prims,
      Array2D<Real, 0, order_t+1, 0, cls_t::NPRIM-1>&                   primsNormal) noexcept
  {
      // For each image point
      for (int iim = 0; iim < order_t; ++iim) {
          // For each interpolation point (corner) of its stencil
          for (int iip = 0; iip < NIP; ++iip) {
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
                  primsNormal(iim + 2, n) += prims(ii, jj, kk, n) * ipweights(iim, iip);
              } // end loop over primitive variables
          } // end loop over interpolation points (corners)
      } // end loop over image points
  }

  /*////////////////////////////////////////////////////////////////
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
  template <int order_t>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void global2local(
      int iip,
      Array2D<Real,0,order_t+1,0,cls_t::NPRIM-1>& primsNormal,
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

  /**
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
   */
  template <int order_t>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void local2global(
      int jj,
      Array2D<Real,0,order_t+1,0,cls_t::NPRIM-1>& primsNormal,
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
    this->idxmap_a.resize(this->ngeom);
    this->inout_fa.resize(this->ngeom);
    this->ntotalfaces = 0;

    for (int i = 0; i < ngeom; i++) {
      Print() << "----------------------------------" << std::endl;

#if (AMREX_SPACEDIM == 2)

      // Determine minimum dx across all levels for setting polygon tolerance
      Real min_dx = std::numeric_limits<Real>::max();
      for (int lev = 0; lev < dx_a.size(); ++lev) {
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              min_dx = std::min(min_dx, dx_a[lev][d]);
          }
      }
      min_dx *= Real(0.5);

      // Read 2D polygon from file
      if (!read_polygon_2d(files_a[i], geom_a[i], min_dx)) {
        amrex::Abort(std::string("Invalid 2D geometry filename: ") + files_a[i]);
      }
      Print() << "Geometry (i=" << i << ") " << files_a[i] << " read" << std::endl;
      Print() << "Number of vertices in polygon: "<< geom_a[i].size() << "\n";
      
      // constructs AABB tree and computes internal KD-tree
      // data structure to accelerate distance queries
      tree_pa[i] = 
        new Tree(geom_a[i].edges_begin(), geom_a[i].edges_end());
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
      tree_pa[i]->accelerate_distance_queries();
      Print() << "AABB tree constructed" << std::endl;

      // plane class also computes orthogonal direction to the face. However, the
      // orthogonal vector is not normalised.
      std::for_each(geom_a[i].facets_begin(), geom_a[i].facets_end(),
                    compute_plane_equations);
      Print() << "Plane equations per face computed" << std::endl;

      // make inside/outside function for current geometry only
      inout_fa[i] = new inside_t(geom_a[i]);
      Print() << "In out testing function constructed for geometry " << files_a[i] << "\n";
      
      //compute total face count across all geometries
      ntotalfaces += geom_a[i].size_of_facets();
    
#endif

    // Store offset for current geometry
    this->geom_offsets[i] = static_cast<int>(this->LocalFrame_a.size());

    // build local frame, surface element area and index mapping
    build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a, idxmap_a[i], this->geom_offsets[i]);
    } // end loop over geometries

    // Check for geometry consistency (no intersections, no containment)
    check_ibm_geometry_consistency(ngeom, geom_a.data(), inout_fa.data(), files_a.data());

    Print() << "----------------------------------" << std::endl;
    Print() << "----------------------------------" << std::endl;
    
    if (plot_surf) {

      Print() << "Total number of faces across all geometries: " << ntotalfaces << std::endl;
      // Initialize surfdata container to this total; safe defaults
      surfdata_a.resize(ntotalfaces);
      for (int f = 0; f < ntotalfaces; ++f) {
        surfdata_a[f].iface = f; 
      }
      intfaces_in_fab.clear(); // face indexing will be built later per level
    }
  } // end read_geom

}; // end class eib_t
#endif // EIB_H_