#ifndef IBM_SOLVER_H_
#define IBM_SOLVER_H_

#include <AMReX_Scan.H>
#include <AMReX_Reduce.H>

#include "ibm_containers.h"

#include <iomanip>  // std::setprecision, std::fixed

//===================================================================================
///-------------------------------- main class --------------------------------------
///
/// \brief ibm_solver_t is explicit geometry (triangulation based) immersed boundary method
/// class. It holds an array of IBMultiFab, one for each AMR level; and it also holds
/// the geometry
///
template <typename wallmodel, typename param, typename cls_t>
class ibm_solver_t
{
public:
  // constant factor for image point
  static constexpr int  iorder_tparm = param::interp_order; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm = param::extrap_order; // number of image points used for ghost point extrapolation
  static constexpr Real cim = param::alpha;
  
  static constexpr int  iorder_tparm_surf = param::interp_order_surf; // number of weighted points used for image point construction
  static constexpr int  eorder_tparm_surf = param::extrap_order_surf; // number of image points used for surface reconstruction
  static constexpr Real cim_surf  = param::alpha_surf;

  // number of ghost layers needed for IB method
  static constexpr int  ghost_layers = param::ghost_layers;
  static_assert(ghost_layers <= cls_t::NGHOST,
              "IBM ghost_layers exceeds cls_t::NGHOST");
  
  // true: interior of closed geometry is solid; false: interior is fluid
  static constexpr bool interior_is_solid = param::interior_is_solid; 

  // ideal number of interpolation points for each image point(ghost point extrapolation and surface reconstruction)
  static constexpr int  N_InterP      = ipow(iorder_tparm + 1, AMREX_SPACEDIM);
  static constexpr int  N_InterP_surf = ipow(iorder_tparm_surf + 1, AMREX_SPACEDIM);

  using GPDATA = gpData_t<eorder_tparm, iorder_tparm>;
  using SURFIMP = surfImp_t<eorder_tparm_surf, iorder_tparm_surf>;
  using GPSTORE = GPStore<eorder_tparm, iorder_tparm>;
  using GPSTOREVIEW = GPStoreView<eorder_tparm, iorder_tparm>;

  // MultiFabs pointer to Amr class instance
  Amr* amr_p;                                             
  // Immersed boundary MultiFab array (uint8_t multifab array)
  Vector<IBMultiFab<uint8_t, GPDATA>*> bmf_a;  
  // Level-wide flattened ghost-point storage (CSR indexed by local FAB)
  Vector<GPSTORE> gpstore_a;

  // parameters for cell size and refinement ratio
  Vector<IntVect> rratio_a;                                 // vector of refinement ratio per level in each direction
  Vector<GpuArray<Real, AMREX_SPACEDIM>> dx_a;              // vector of cell sizes per level in each direction
  Vector<Real> diag_a;                                      // vector of cell diagonal length per level
  Vector<Real> di_a;                                        // image point distance per level for ghost point extrapolation
  Vector<Real> di_a_surf;                                   // image point distance per level for surface reconstruction

  // geometry related data
  int ngeom = 0;                                            // number of geometries
  Vector<GeomType> geom_a;                                  // IB explicit geometry
#ifdef AMREX_USE_CGAL
  Vector<Tree> tree_a;                                      // CGAL AABB tree per geometry
  Vector<PrimitiveIndexMap> idxmap_a;                       // CGAL primitive-to-index map per geometry
#else
  Vector<BVH> bvh_a;                                        // BVH per geometry (replaces CGAL AABB tree)
#endif
  Vector<inside_t*> inout_fa;                               // in out testing function per geometry
  Vector<Bbox> bbox_a;                                      // bounding box per geometry (world-frame, for fast rejection)
  Vector<Bbox> bbox_body_a;                                 // bounding box per geometry (body-frame, constant after init)
  Vector<RigidTransform> transform_a;                       // rigid-body transform per geometry (body→world)
  Vector<std::string> geom_names;                           // geometry names stripped from input filenames

  Gpu::ManagedVector<LocalFrame> LocalFrame_a;              // local orthonormal frame matrix (flattened, body-frame)
  Gpu::ManagedVector<SurfElem> SurfElem_a;                  // surface element area and coordinates (flattened, body-frame)
  Gpu::ManagedVector<int> geom_offsets;                     // Start index for each geometry in flattened arrays
 
  // surface related data
  int ntotalfaces = 0;                                      // number of faces/edges across all geometries
  SURFIMP surfimp_soa;                                      // face/edge image point information (SoA structure)
  surfPhys_t surfphys_soa;                                  // face/edge physical data and identification (SoA structure)
  Vector<FaceCSR> faces_per_level;                          // faces integers per fab and per level [lev] (CSR format)

  /** 
   * \brief Destructor to release allocated memory
   */
  ~ibm_solver_t() noexcept
  {
    // Release per-level IBMultiFab pointers if any remain
    for (auto*& p : bmf_a) {
      if (p) { delete p; p = nullptr; }
    }

    // Release inside/outside testers
    for (auto*& f : inout_fa) {
        if (f) { delete f; f = nullptr; }
    }
  }

  /**
   * \brief Explicitly release all GPU-managed memory before amrex::Finalize().
   *
   * Must be called while AMReX arenas are still alive. The global inline
   * IBM::ib outlives amrex::Finalize(), so its implicit destructor would
   * free arena memory after the arena is destroyed (static destruction
   * order problem).  Calling cleanup() first leaves the destructor with
   * nothing to free.
   */
  void cleanup() noexcept
  {
    // Delete raw-pointer members first (their internals use arena memory)
    for (auto*& p : bmf_a)    { if (p) { delete p; p = nullptr; } }
    for (auto*& f : inout_fa) { if (f) { delete f; f = nullptr; } }

    // Swap-with-empty idiom: guarantees capacity→0 and arena memory freed NOW,
    // so the post-Finalize destructor finds nothing to deallocate.
    { decltype(bmf_a)        tmp; tmp.swap(bmf_a);        }
    { decltype(inout_fa)     tmp; tmp.swap(inout_fa);     }
    { decltype(LocalFrame_a) tmp; tmp.swap(LocalFrame_a); }
    { decltype(SurfElem_a)   tmp; tmp.swap(SurfElem_a);   }
    { decltype(geom_offsets) tmp; tmp.swap(geom_offsets);  }

    // Compound types: destroying elements calls their ManagedVector destructors
    { decltype(geom_a) tmp; tmp.swap(geom_a); }
#ifdef AMREX_USE_CGAL
    { decltype(tree_a)   tmp; tmp.swap(tree_a);   }
    { decltype(idxmap_a) tmp; tmp.swap(idxmap_a); }
#else
    { decltype(bvh_a) tmp; tmp.swap(bvh_a); }
#endif
    { decltype(bbox_a)      tmp; tmp.swap(bbox_a);      }
    { decltype(bbox_body_a) tmp; tmp.swap(bbox_body_a); }
    { decltype(transform_a) tmp; tmp.swap(transform_a); }
    { decltype(geom_names)  tmp; tmp.swap(geom_names);  }
    { decltype(gpstore_a)  tmp; tmp.swap(gpstore_a);  }

    // These structs have their own clear() with shrink_to_fit()
    surfimp_soa.clear();
    surfphys_soa.clear();

    { decltype(faces_per_level) tmp; tmp.swap(faces_per_level); }

    amr_p = nullptr;
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
    gpstore_a.resize(lmax + 1);
    faces_per_level.resize(lmax + 1);

    dx_a.resize(lmax + 1);
    dx_a[0] = amr_p->Geom(0).CellSizeArray();
    for (int i = 1; i <= lmax; i++) {
      for (int j = 0; j < AMREX_SPACEDIM; j++) {
        dx_a[i][j] = dx_a[i - 1][j] / rratio_a[i - 1][j];
      }
    }

    di_a.resize(lmax + 1);
    di_a_surf.resize(lmax + 1);
    diag_a.resize(lmax + 1);

    for (int i = 0; i <= lmax; i++) {
      diag_a[i] = std::sqrt(
        AMREX_D_TERM( std::pow(dx_a[i][0], 2),
                      + std::pow(dx_a[i][1], 2),
                      + std::pow(dx_a[i][2], 2)) );

      di_a[i] = cim * diag_a[i];
      di_a_surf[i] = cim_surf * diag_a[i];
    }

    // read geometry from file
    Real t_geom = amrex::second();
    read_geom();
    t_geom = amrex::second() - t_geom;
    amrex::Print() << "  IBM init: read_geom       = " << t_geom << " s\n";
  }

  /**
   * \brief create IBMultiFab at a level and store pointers to it
   * \param bxa BoxArray for the level
   * \param dm DistributionMapping for the level
   * \param lev The AMR level
   */
  void build_mf(const BoxArray& bxa, const DistributionMapping& dm, int lev)
  {
    // Prevent memory leak if MF already exists
    if (lev < 0 || lev >= static_cast<int>(bmf_a.size())) {
        amrex::Abort("ibm_solver_t::build_mf: lev is out of bounds");
    }

    // Prevent memory leak if MF already exists
    if (bmf_a[lev] != nullptr) {
        destroy_mf(lev);
    }

    // Use default MFInfo (configured to use The_Managed_Arena in IBMultiFab.h)
    bmf_a[lev] = new IBMultiFab<uint8_t, GPDATA>(bxa, dm, 2, cls_t::NGHOST);
  }

  /**
   * \brief destroy IBMultiFab at a level
   * \param lev The AMR level of the IBMultiFab to be destroyed.
   */
  void destroy_mf(int lev)
  {
    if (lev >= 0 && lev < static_cast<int>(bmf_a.size())) { 
      // safe delete: delete nullptr is valid in C++
      delete bmf_a[lev]; 
      bmf_a[lev] = nullptr;  // Prevent dangling pointer
    }
    if (lev >= 0 && lev < static_cast<int>(gpstore_a.size())) {
      gpstore_a[lev].clear();
    }
  }

  /**
   * \brief Rebuild all geometry-level data from current vertex positions.
   *
   * After externally modifying vertex positions in geom_a (e.g. for moving
   * geometry), call this to reconstruct BVH trees, inside/outside testers,
   * bounding boxes, and the surface cache (LocalFrame_a, SurfElem_a).
   * Must be called BEFORE rebuildIBM().
   */
  /**
   * \brief Lightweight rigid-body transform update (replaces full BVH rebuild).
   *
   * For rigid-body FSI, the geometry shape never changes — only its position
   * and orientation.  Instead of moving all vertices and rebuilding the BVH
   * tree, InsideTester, LocalFrame, and SurfElem every time step, we store
   * the geometry in its reference (body) frame and only update the transform.
   *
   * All query functions (computeMarkers, initialiseGPs, computeAllGPs, etc.)
   * apply the inverse transform to query points before using the body-frame
   * BVH, then forward-transform the results back to the world frame.
   * The results are mathematically identical to a full rebuild.
   *
   * \param geomIdx  Index of the geometry to update.
   * \param T        New rigid-body transform (body → world).
   */
  void updateRigidTransform(int geomIdx, const RigidTransform& T)
  {
      transform_a[geomIdx] = T;
      // Update world-frame bounding box from body-frame bbox + new transform
      bbox_a[geomIdx] = T.transform_bbox(bbox_body_a[geomIdx]);
  }

  /**
   * \brief Full geometry rebuild (for deformable bodies or re-initialization).
   *
   * This is the expensive path that rebuilds BVH, InsideTester, LocalFrame,
   * and SurfElem from scratch.  For rigid-body FSI, use updateRigidTransform()
   * instead — it is O(1) per geometry.
   */
  void rebuildGeometryData()
  {
    LocalFrame_a.clear();
    SurfElem_a.clear();

    for (int i = 0; i < ngeom; i++) {
#ifdef AMREX_USE_CGAL
      tree_a[i].clear();
#if (AMREX_SPACEDIM == 2)
      tree_a[i].insert(geom_a[i].edges_begin(), geom_a[i].edges_end());
#elif (AMREX_SPACEDIM == 3)
      tree_a[i].insert(faces(geom_a[i]).first, faces(geom_a[i]).second, geom_a[i]);
#endif
      tree_a[i].build();

      delete inout_fa[i];
      inout_fa[i] = new inside_t(geom_a[i]);
#else
      bvh_a[i].build(geom_a[i]);

      delete inout_fa[i];
      inout_fa[i] = new inside_t(geom_a[i], bvh_a[i]);
#endif

#if defined(AMREX_USE_CGAL) && (AMREX_SPACEDIM == 3)
      bbox_body_a[i] = PMP::bbox(geom_a[i]);
#else
      bbox_body_a[i] = geom_a[i].bbox();
#endif
      bbox_a[i] = transform_a[i].transform_bbox(bbox_body_a[i]);

      geom_offsets[i] = static_cast<int>(LocalFrame_a.size());
#ifdef AMREX_USE_CGAL
      build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                           idxmap_a[i], geom_offsets[i], i);
#elif defined(AMREX_USE_GPU)
      build_geometry_cache_gpu(geom_a[i], SurfElem_a, LocalFrame_a, i);
#else
      build_geometry_cache(geom_a[i], SurfElem_a, LocalFrame_a,
                           geom_offsets[i], i);
#endif
    }
    geom_offsets[ngeom] = static_cast<int>(LocalFrame_a.size());

    if (!interior_is_solid) {
      convert_inout(LocalFrame_a);
    }
  }

  /**
   * \brief Computes the solid/fluid markers and identifies ghost points for the Immersed Boundary Method.
   *
   * This function iterates over the grid points per level to determine if they are inside (solid) or outside (fluid)
   * the immersed boundary geometry using BVH-accelerated inside/outside tests. It populates the ibMarkers where:
   * - Component 0 indicates if a point is solid (true) or fluid (false).
   * - Component 1 indicates if a solid point is a ghost point (neighbor to a fluid point).
   *
   * \param lev The current AMR level.
   */
  void computeMarkers(int lev)
  {
    BL_PROFILE("IBM::computeMarkers");
    auto& mfab = *bmf_a[lev];
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

    // Collect per-FAB ghost-point counts for CSR allocation
    const int nfabs_local = mfab.local_size();
    Vector<int> gp_counts(nfabs_local, 0);
    int ifab_local = 0;

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      auto& ibFab = mfab.get(mfi);
      const Box& bx = mfi.tilebox();
      const auto& ibMarkers = mfab.array(mfi);

      // Clear legacy per-fab ghost-point data
      ibFab.gpData.clear();

#ifdef AMREX_USE_GPU
      // ================================================================
      // GPU path: ParallelFor for solid + ghost markers
      // ================================================================
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
          "ngeom exceeds MAX_NGEOM; increase MAX_NGEOM in ibm_containers.h");

      // Pre-extract GPU-friendly views (trivially-copyable POD structs)
      GpuArray<InsideTesterView, MAX_NGEOM> inout_views;
      GpuArray<AABB, MAX_NGEOM> bboxes;
      GpuArray<RigidTransform, MAX_NGEOM> transforms;
      for (int ii = 0; ii < ngeom; ii++) {
          inout_views[ii] = inout_fa[ii]->view();
          bboxes[ii] = bbox_a[ii];           // world-frame bbox
          transforms[ii] = transform_a[ii];  // body→world transform
      }

      const int ngeom_local = ngeom;
      const auto dx_lev = dx_a[lev];

      // --- Step 1: Compute solid markers (component 0) ---
      // Query points are in world frame; we fast-reject with world-frame bbox,
      // then inverse-transform to body frame for the BVH inside/outside test.
      const Box& bxg = amrex::grow(bx, cls_t::NGHOST);
      amrex::ParallelFor(bxg,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        ibMarkers(i, j, k, 0) = static_cast<uint8_t>(0);
        ibMarkers(i, j, k, 1) = static_cast<uint8_t>(0);

        Point gridpoint = make_grid_point(prob_lo, dx_lev, i, j, k);

        for (int ii = 0; ii < ngeom_local; ii++) {
          // Fast rejection with world-frame bounding box
          if (!bbox_contains(bboxes[ii], gridpoint)) {
              continue;
          }

          // Transform query point to body frame for BVH query
          Point gp_body = transforms[ii].to_body(gridpoint);
          BoundedSide result = inout_views[ii](gp_body);

          if (result == BoundedSide::Inside || result == BoundedSide::OnBoundary) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      }); // end ParallelFor for solid markers

      // --- Step 2: Compute ghost markers + count ghost points (fused) ---
      constexpr int gl = ghost_layers;
      const Box& bxgp = amrex::grow(bx, GP_BOX_EXTRA);
      {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bxgp, reduce_data,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
          if (!ibMarkers(i, j, k, 0)) return { 0 };

          bool ghost = false;
          for (int d = 1; d <= gl; d++) {
            ghost = ghost || (!ibMarkers(i-d, j,   k,   0));
            ghost = ghost || (!ibMarkers(i+d, j,   k,   0));
            ghost = ghost || (!ibMarkers(i,   j-d, k,   0));
            ghost = ghost || (!ibMarkers(i,   j+d, k,   0));
#if (AMREX_SPACEDIM == 3)
            ghost = ghost || (!ibMarkers(i,   j,   k-d, 0));
            ghost = ghost || (!ibMarkers(i,   j,   k+d, 0));
#endif
            if (ghost) break;
          }

          if (ghost) {
            ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
            return { 1 };
          }
          return { 0 };
        });
        ReduceTuple hv = reduce_data.value(reduce_op);
        int ngps_fab = amrex::get<0>(hv);
        ibFab.gpData.ngps = ngps_fab;
        gp_counts[ifab_local] = ngps_fab;
      }

#else
      // ================================================================
      // CPU path: original LoopOnCpu implementation
      // ================================================================

      // compute solid markers
      amrex::LoopOnCpu(amrex::grow(bx, cls_t::NGHOST), [&](int i, int j, int k) {
        ibMarkers(i, j, k, 0) = static_cast<uint8_t>(0);
        ibMarkers(i, j, k, 1) = static_cast<uint8_t>(0);

        Point gridpoint = make_grid_point(prob_lo, dx_a[lev], i, j, k);

        for (int ii = 0; ii < ngeom; ii++) {
          if (!bbox_contains(bbox_a[ii], gridpoint)) {
              continue;
          }

          // Transform to body frame for inside/outside test
          Point gp_body = transform_a[ii].to_body(gridpoint);
          inside_t& inside = *inout_fa[ii];
          BoundedSide result = inside(gp_body);
          IB_WarnOnBoundary(ii, lev, i, j, k, result, gridpoint);

          if (result == BoundedSide::Inside || result == BoundedSide::OnBoundary) {
            ibMarkers(i, j, k, 0) = static_cast<uint8_t>(ii + 1);
            break;
          }
        }
      }); // end LoopOnCpu for solid markers

      // compute ghost markers
      int ngps_fab = 0;
      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA), [&](int i, int j, int k) {
        bool ghost = false;
        if (ibMarkers(i, j, k, 0)) {
          for(int d = 1; d <= ghost_layers; d++) {
            ghost = ghost || (!ibMarkers(i-d, j,   k,   0));
            ghost = ghost || (!ibMarkers(i+d, j,   k,   0));
            ghost = ghost || (!ibMarkers(i,   j-d, k,   0));
            ghost = ghost || (!ibMarkers(i,   j+d, k,   0));
#if (AMREX_SPACEDIM == 3)
            ghost = ghost || (!ibMarkers(i,   j,   k-d, 0));
            ghost = ghost || (!ibMarkers(i,   j,   k+d, 0));
#endif
            if (ghost) break;
          }
          ngps_fab += ghost;

          if (ghost) {
            ibMarkers(i, j, k, 1) = ibMarkers(i, j, k, 0);
          } else {
            ibMarkers(i, j, k, 1) = static_cast<uint8_t>(0);
          }
        }
      }); // end LoopOnCpu for ghost markers
      ibFab.gpData.ngps = ngps_fab;
      gp_counts[ifab_local] = ngps_fab;

#endif // AMREX_USE_GPU

    } // end MFIter

    // Pre-allocate level-wide flattened GPStore from per-FAB counts (CSR)
    gpstore_a[lev].allocate(gp_counts);

  } // end computeMarkers

  /**
   * \brief Initialises geometric and interpolation data for Ghost Points (GPs).
   *
   * This function iterates over all ghost points identified in the `computeMarkers` step.
   * For each ghost point, it performs the following operations:
   * 1. Identifies the specific geometry (body) the ghost point belongs to.
   * 2. Finds the closest point on the surface (IB point) using BVH trees.
   * 3. Computes the normal distance from the ghost point to the surface.
   * 4. Constructs a local orthonormal frame (normal, tangent1, tangent2) at the IB point.
   * 5. Projects "Image Points" (IMs) into the fluid domain along the surface normal.
   * 6. Computes trilinear interpolation weights and indices for these Image Points.
   *
   * All computed data is stored in the level-wide `GPStore` (CSR layout) for use in boundary condition reconstruction.
   *
   * \param lev The current AMR level index.
   */
  void initialiseGPs(int lev) {
    BL_PROFILE("IBM::initialiseGPs");
    auto& mfab = *bmf_a[lev];
    auto& gpstore = gpstore_a[lev];
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = amr_p->Geom(lev).ProbLoArray();

#ifdef AMREX_USE_GPU
    // ================================================================
    // GPU path: ParallelFor with atomic counter for GP index
    // ================================================================

    // Pre-build GPU-capturable BVH4QueryViews for each geometry
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngeom <= MAX_NGEOM,
        "ngeom exceeds MAX_NGEOM; increase MAX_NGEOM in ibm_containers.h");
    GpuArray<BVH4QueryView, MAX_NGEOM> bqv;
    GpuArray<int, MAX_NGEOM + 1> geom_off;
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    for (int ii = 0; ii < ngeom; ii++) {
        bqv[ii]        = bvh_a[ii].query_view(geom_a[ii]);
        geom_off[ii]   = geom_offsets[ii];
        transforms[ii] = transform_a[ii];
    }
    geom_off[ngeom] = geom_offsets[ngeom];

    const auto dx_lev       = dx_a[lev];
    const Real di_lev       = di_a[lev];
    const Real diag_lev     = diag_a[lev];
    const int  ngeom_local  = ngeom;
    const auto* lf_ptr      = LocalFrame_a.data();

    auto gpview = gpstore.view();  // GPU-capturable POD with writable _w pointers

    int ifab_local = 0;
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      auto& ibFab = mfab.get(mfi);

      const int gp_offset = gpstore.fab_offsets[ifab_local];
      const int ngps_fab  = gpstore.fab_offsets[ifab_local + 1] - gp_offset;

      if (ngps_fab == 0) {
          ibFab.gpData.ngps = 0;
          continue;
      }

      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);

      // Atomic counter for assigning GP flat indices
      Gpu::DeviceScalar<int> d_gp_count(0);
      int* p_gp_count = d_gp_count.dataPtr();

      const Box bxgp = amrex::grow(bx, GP_BOX_EXTRA);
      amrex::ParallelFor(bxgp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (!ibMarkers(i, j, k, 1)) return;

        // Atomically claim a slot in the flat GP array
        int local_idx = Gpu::Atomic::Add(p_gp_count, 1);
        const int gidx = gp_offset + local_idx;

        Point gp = make_grid_point(prob_lo, dx_lev, i, j, k);
        gpview.gp_ijk_w[gidx] = make_vec<int>(i, j, k);

        int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;
        AMREX_ASSERT(geomIdx >= 0 && geomIdx < ngeom_local);
        gpview.geomIdx_w[gidx] = geomIdx;

        // Transform query point to body frame, then BVH closest-point query
        const auto& T = transforms[geomIdx];
        Point gp_body = T.to_body(gp);
        ClosestPointResult closest_elem = bqv[geomIdx].closest_point_query(gp_body);
        int f_idx = closest_elem.prim_id + geom_off[geomIdx];
        gpview.elemIdx_w[gidx] = f_idx;

        // Transform closest point back to world frame
        Point cp = T.to_world(closest_elem.point);

        // Distance (invariant under rigid transform, but computed in world frame)
        Real disGP_val = std::sqrt(point_distance_sq(gp, cp));
        AMREX_ASSERT(disGP_val < diag_lev);
        gpview.disGP_w[gidx] = disGP_val;
        gpview.ib_xyz_w[gidx] = make_vec<Real>(cp);

        // Rotate body-frame LocalFrame to world frame
        const LocalFrame& lf_body = lf_ptr[f_idx];
        LocalFrame localframe;
        T.rotate_to_world(lf_body.normal,   localframe.normal);
        T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
        T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif

        // Image points
        Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
        Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;
        Array1D<Real, 0, eorder_tparm - 1> disIM;
        Array1D< int, 0, eorder_tparm - 1> imp_ninterp;

        for (int jj = 0; jj < eorder_tparm; jj++) {
          Point cp_start;
          if (jj == 0) {
            cp_start = cp;
          } else {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) cp_start[d] = imp_xyz(jj - 1, d);
          }

          if (jj == 0) {
            search_optimal_image_point<eorder_tparm, iorder_tparm>(
                cp_start, localframe,
                lev, prob_lo, dx_lev, di_lev,
                bxg, ibMarkers,
                gpview, gidx,
                imp_xyz, imp_ijk, disIM, imp_ninterp);
          } else {
            search_image_point<eorder_tparm, iorder_tparm>(
                jj, cp_start, localframe,
                lev, prob_lo, dx_lev, di_lev,
                bxg, ibMarkers,
                gpview, gidx,
                imp_xyz, imp_ijk, disIM, imp_ninterp);
          }
        }

        gpview.imp_xyz_w[gidx]       = imp_xyz;
        gpview.imp_ijk_w[gidx]       = imp_ijk;
        gpview.disIM_w[gidx]         = disIM;
        gpview.imp_ninterp_w[gidx]   = imp_ninterp;

        Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1> imp_ipweights;
        Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;

        computeIPweights<eorder_tparm, iorder_tparm, GPSTOREVIEW>(
            imp_ipweights, imp_ip_ijk,
            imp_xyz, imp_ijk, imp_ninterp,
            prob_lo, dx_lev, ibMarkers);

        gpview.imp_ipweights_w[gidx] = imp_ipweights;
        gpview.imp_ip_ijk_w[gidx]    = imp_ip_ijk;

        ibMarkers(i, j, k, 1) = static_cast<uint8_t>(imp_ninterp(0));
      }); // end ParallelFor

      // Verify count matches
      int h_gp_count = d_gp_count.dataValue();
      if (h_gp_count != ngps_fab) {
        amrex::Abort("Error in initialiseGPs (GPU): mismatch in ghost point count");
      }
      ibFab.gpData.ngps = ngps_fab;
    } // end MFIter

#else
    // ================================================================
    // CPU path: original LoopOnCpu implementation
    // ================================================================
    int ifab_local = 0;
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab_local) {
      auto& ibFab = mfab.get(mfi);

      // Get the pre-allocated range in the flattened GPStore for this FAB
      const int gp_offset = gpstore.fab_offsets[ifab_local];
      const int ngps_fab  = gpstore.fab_offsets[ifab_local + 1] - gp_offset;

      if (ngps_fab == 0) {
          ibFab.gpData.ngps = 0;
          continue;
      }

      const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
      const Box& bx = mfi.tilebox();
      auto const ibMarkers = mfab.array(mfi);  // uint8_t array

      // CPU loop: Write directly into gpstore at pre-allocated offsets (no push_back).
      int gp_count = 0; 
      amrex::LoopOnCpu(amrex::grow(bx, GP_BOX_EXTRA), [&](int i, int j, int k) {
        // for each ghost point
        if (ibMarkers(i, j, k, 1)) {
          const int gidx = gp_offset + gp_count;  // global flat index in GPStore
          gp_count++;

          Point gp = make_grid_point(prob_lo, dx_a[lev], i, j, k);
          gpstore.gp_ijk[gidx] = make_vec<int>(i, j, k);
          
          int geomIdx = static_cast<int>(ibMarkers(i, j, k, 0)) - 1;
          AMREX_ASSERT_WITH_MESSAGE(geomIdx >= 0 && geomIdx < ngeom, 
                  "Invalid geometry index in initialiseGPs");
          gpstore.geomIdx[gidx] = geomIdx;

#ifdef AMREX_USE_CGAL
          ClosestPointResult closest_elem =
              cgal_closest_point_query(tree_a[geomIdx], idxmap_a[geomIdx],
                                       gp, this->geom_offsets[geomIdx]);
#else
          ClosestPointResult closest_elem =
              bvh_a[geomIdx].closest_point_query(gp, geom_a[geomIdx]);
#endif

          int f_idx = closest_elem.prim_id + this->geom_offsets[geomIdx];
          gpstore.elemIdx[gidx] = f_idx;

          Point cp = closest_elem.point;

          Real disGP = std::sqrt(point_distance_sq(gp, cp));

          AMREX_ASSERT_WITH_MESSAGE(
              disGP < diag_a[lev],  "Ghost point and IB point distance larger than mesh diagonal");
          gpstore.disGP[gidx] = disGP;

          gpstore.ib_xyz[gidx] = make_vec<Real>(cp);
          
          const auto& localframe = LocalFrame_a[f_idx];

          Array2D<Real, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_xyz;
          Array2D< int, 0, eorder_tparm - 1, 0, AMREX_SPACEDIM - 1> imp_ijk;
          Array1D<Real, 0, eorder_tparm - 1> disIM;
          Array1D< int, 0, eorder_tparm - 1> imp_ninterp;

          for (int jj = 0; jj < eorder_tparm; jj++) {
            Point cp_start;
            if (jj == 0) {
              cp_start = cp;
            } else {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
              cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1));
#else
              cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1), imp_xyz(jj - 1, 2));
#endif
#else
              for (int d = 0; d < AMREX_SPACEDIM; ++d) cp_start[d] = imp_xyz(jj - 1, d);
#endif
            }
          
            if (jj == 0) {
              search_optimal_image_point<eorder_tparm, iorder_tparm>(
                                      cp_start, localframe, 
                                      lev, prob_lo, dx_a[lev], di_a[lev],
                                      bxg, ibMarkers, 
                                      gpstore, gidx,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);
            } 
            else {
              search_image_point<eorder_tparm, iorder_tparm>(
                                      jj, cp_start, localframe, 
                                      lev, prob_lo, dx_a[lev], di_a[lev],
                                      bxg, ibMarkers, 
                                      gpstore, gidx,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);
            }
          }

          gpstore.imp_xyz[gidx] = imp_xyz;
          gpstore.imp_ijk[gidx] = imp_ijk;
          gpstore.disIM[gidx] = disIM;
          gpstore.imp_ninterp[gidx] = imp_ninterp;

          Array2D<Real, 0, eorder_tparm - 1 , 0, N_InterP -1 > imp_ipweights;
          Array3D< int, 0, eorder_tparm - 1 , 0, N_InterP -1, 0, AMREX_SPACEDIM - 1> imp_ip_ijk;
          
          computeIPweights<eorder_tparm, iorder_tparm, GPSTORE>(
              imp_ipweights, 
              imp_ip_ijk, 
              imp_xyz, 
              imp_ijk, 
              imp_ninterp,
              prob_lo, dx_a[lev], ibMarkers);
          
          gpstore.imp_ipweights[gidx] = imp_ipweights;
          gpstore.imp_ip_ijk[gidx] = imp_ip_ijk;

          ibMarkers(i, j, k, 1) = static_cast<uint8_t>(imp_ninterp(0)); 

          } //end if (ibMarkers(i,j,k,1))
      });//end loop on bx
    
      if(gp_count != ngps_fab) {
        amrex::Abort("Error in initialiseGPs: mismatch in ghost point count");
      }
      ibFab.gpData.ngps = ngps_fab;
    } //end MFIter

#endif // AMREX_USE_GPU

    // NOTE: gpstore.shrink() disabled — shrink_to_fit() on PODVector
    // corrupts data on some platforms (observed with AMReX ManagedVector on CPU).

    // Print GP reconstruction quality diagnostics
    // reportGPDiagnostics(lev);
  }

  // ========================================================================
  // GP reconstruction quality diagnostics
  // ========================================================================

  /**
   * \brief Report statistics on ghost-point reconstruction quality.
   *
   * Scans the GPStore after initialiseGPs and prints:
   *   - Total GP count
   *   - Histogram of first image-point fluid stencil count (imp_ninterp[0])
   *   - Effective extrapolation order distribution (0th / 1st / 2nd)
   *   - disGP statistics (min, mean, max) in units of dx diagonal
   *   - Weight concentration: max single-point weight among all GPs
   *
   * \param lev  AMR level index.
   */
  void reportGPDiagnostics(int lev)
  {
    auto& gpstore = gpstore_a[lev];
    const int ngps = gpstore.total_ngps;
    const int istep = amr_p->levelSteps(0);

    if (ngps == 0) {
      amrex::Print() << "[GP-Diag] Step " << istep << " Level " << lev << ": 0 ghost points\n";
      return;
    }

    // Ensure GPU data is accessible on host
    Gpu::streamSynchronize();

    // Grid metadata for coordinate conversion
    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto& dx = dx_a[lev];

    constexpr int IDEAL = N_InterP;
    int hist[IDEAL + 1] = {};

    // Cross-tabulation: severity × effective order
    // severity bins: ideal (n_fluid==IDEAL), mild (IDEAL/2 < n < IDEAL),
    //               severe (n <= IDEAL/2)
    int sev_ideal = 0, sev_mild = 0, sev_severe = 0;
    int sev_mild_order[3] = {}, sev_severe_order[3] = {};  // indexed by eff_order

    int order_0 = 0, order_1 = 0, order_2 = 0;

    Real disGP_min = std::numeric_limits<Real>::max();
    Real disGP_max = Real(0.0);
    Real disGP_sum = Real(0.0);

    Real max_single_weight = Real(0.0);
    Real sum_weight_deficit = Real(0.0);
    int  n_low_stencil = 0;

    for (int ii = 0; ii < ngps; ++ii) {
      int n0 = gpstore.imp_ninterp[ii](0);
      int idx = amrex::max(0, amrex::min(IDEAL, n0));
      hist[idx]++;

      // Effective order
      int eff = eorder_tparm;
      for (int k = 0; k < eorder_tparm; ++k) {
        if (gpstore.imp_ninterp[ii](k) < INTERP_THRESHOLD_GP) {
          eff = k; break;
        }
      }
      if      (eff == 0) order_0++;
      else if (eff == 1) order_1++;
      else               order_2++;

      // Severity classification
      if (n0 == IDEAL) {
        sev_ideal++;
      } else if (n0 > IDEAL / 2) {
        sev_mild++;
        sev_mild_order[amrex::min(eff, 2)]++;
      } else {
        sev_severe++;
        sev_severe_order[amrex::min(eff, 2)]++;
      }

      // disGP
      Real d = gpstore.disGP[ii];
      disGP_min = amrex::min(disGP_min, d);
      disGP_max = amrex::max(disGP_max, d);
      disGP_sum += d;

      // Weight concentration
      Real wmax = Real(0.0);
      for (int c = 0; c < N_InterP; ++c) {
        wmax = amrex::max(wmax, gpstore.imp_ipweights[ii](0, c));
      }
      max_single_weight = amrex::max(max_single_weight, wmax);

      if (n0 < IDEAL) {
        n_low_stencil++;
        sum_weight_deficit += wmax;
      }
    }

    Real diag = diag_a[lev];
    Real inv_diag = (diag > Real(1e-30)) ? Real(1.0) / diag : Real(1.0);

    // --- Console output ---
    amrex::Print()
      << "\n[GP-Diag] Step " << istep << " Level " << lev
      << "  |  eorder=" << eorder_tparm
      << "  iorder=" << iorder_tparm
      << "  ghost_layers=" << ghost_layers
      << "  alpha=" << cim
      << "\n";

    amrex::Print()
      << "[GP-Diag]   Total GPs: " << ngps << "\n";

    // Fluid stencil histogram (layered)
    amrex::Print() << "[GP-Diag]   Stencil quality (ideal=" << IDEAL << "):\n";
    amrex::Print() << "[GP-Diag]     IDEAL  (" << IDEAL << "/" << IDEAL << "): "
                   << sev_ideal << " GPs (" << std::fixed << std::setprecision(1)
                   << Real(100.0) * Real(sev_ideal) / Real(ngps) << "%)\n";
    if (sev_mild > 0) {
      amrex::Print() << "[GP-Diag]     MILD   (>" << IDEAL/2 << "/" << IDEAL << "): "
                     << sev_mild << " GPs (" << std::setprecision(1)
                     << Real(100.0) * Real(sev_mild) / Real(ngps) << "%)";
      // Detail per n_fluid
      for (int k = IDEAL - 1; k > IDEAL / 2; --k) {
        if (hist[k] > 0) amrex::Print() << "  [" << k << "/" << IDEAL << "]=" << hist[k];
      }
      amrex::Print() << "\n";
    }
    if (sev_severe > 0) {
      amrex::Print() << "[GP-Diag]     SEVERE (<=" << IDEAL/2 << "/" << IDEAL << "): "
                     << sev_severe << " GPs (" << std::setprecision(1)
                     << Real(100.0) * Real(sev_severe) / Real(ngps) << "%)";
      for (int k = IDEAL / 2; k >= 0; --k) {
        if (hist[k] > 0) amrex::Print() << "  [" << k << "/" << IDEAL << "]=" << hist[k];
      }
      amrex::Print() << "\n";
    }

    // Cross-tabulation: severity × order
    if (sev_mild + sev_severe > 0) {
      amrex::Print() << "[GP-Diag]   Severity x Order cross-tab:\n";
      if (sev_mild > 0) {
        amrex::Print() << "[GP-Diag]     MILD  : ";
        if (sev_mild_order[0]) amrex::Print() << "0th=" << sev_mild_order[0] << " ";
        if (sev_mild_order[1]) amrex::Print() << "1st=" << sev_mild_order[1] << " ";
        if (sev_mild_order[2]) amrex::Print() << "2nd=" << sev_mild_order[2] << " ";
        amrex::Print() << "\n";
      }
      if (sev_severe > 0) {
        amrex::Print() << "[GP-Diag]     SEVERE: ";
        if (sev_severe_order[0]) amrex::Print() << "0th=" << sev_severe_order[0] << " ";
        if (sev_severe_order[1]) amrex::Print() << "1st=" << sev_severe_order[1] << " ";
        if (sev_severe_order[2]) amrex::Print() << "2nd=" << sev_severe_order[2] << " ";
        amrex::Print() << "\n";
      }
      // Weight stats
      Real avg_wmax = sum_weight_deficit / Real(n_low_stencil);
      amrex::Print()
        << "[GP-Diag]   Weight: avg_max=" << std::setprecision(3) << avg_wmax
        << "  worst_max=" << std::setprecision(3) << max_single_weight
        << "\n";
    }

    // Effective order distribution
    amrex::Print() << "[GP-Diag]   Effective order:";
    if (order_0 > 0) amrex::Print() << "  0th=" << order_0;
    if (order_1 > 0) amrex::Print() << "  1st=" << order_1;
    if (order_2 > 0) amrex::Print() << "  2nd=" << order_2;
    amrex::Print() << "\n";

    // disGP statistics
    amrex::Print()
      << "[GP-Diag]   disGP (diag units): min="
      << std::setprecision(4) << disGP_min * inv_diag
      << "  mean=" << std::setprecision(4) << (disGP_sum / Real(ngps)) * inv_diag
      << "  max=" << std::setprecision(4) << disGP_max * inv_diag
      << "  (diag=" << std::setprecision(6) << diag << ")\n";

    amrex::Print() << "\n";
  }

  // ========================================================================
  // Ghost‐point reconstruction — single‐kernel, level‐wide
  // ========================================================================

  /**
   * \brief Level-wide single-kernel ghost-point reconstruction.
   *
   * Prerequisite: \c prims_mf must already be filled with primitives
   * (via cons2prims) for every local FAB on this rank.
   *
   * Builds device-accessible Array4 pointers (one per local FAB) and
   * launches a single ParallelFor over all ghost points at this level.
   *
   * \param prims_mf  MultiFab of primitive variables (same BA/DM as bmf_a[lev]).
   *                  Ghost cells will be overwritten with IB-reconstructed values.
   * \param cls       Pointer to the physics/closure class (device-accessible).
   * \param lev       Current AMR level index.
   */
  void computeAllGPs(MultiFab& prims_mf,
                     const cls_t* cls,
                     int lev)
  {
    BL_PROFILE("IBM::computeAllGPs");
    auto& gpstore = gpstore_a[lev];
    if (gpstore.total_ngps == 0) return;

    auto gpview = gpstore.view();
    auto const* lf_ptr = LocalFrame_a.data();
    const int nfabs_local = gpview.nfabs;

    // Capture per-geometry transforms for rotating body-frame LocalFrame to world
    AMREX_ALWAYS_ASSERT(ngeom <= MAX_NGEOM);
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    GpuArray<int, MAX_NGEOM + 1> geom_off;
    for (int ii = 0; ii < ngeom; ii++) {
        transforms[ii] = transform_a[ii];
        geom_off[ii]   = geom_offsets[ii];
    }
    geom_off[ngeom] = geom_offsets[ngeom];
    const int ngeom_local = ngeom;

    // Build device array of Array4 pointers (one per local FAB).
    // No snapshot needed: GPs write only to ghost cells (ibMarkers==1)
    // while image-point interpolation reads only from fluid cells (ibMarkers==0).
    auto& mfab = *bmf_a[lev];

    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
      }
      Gpu::copyAsync(Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      Gpu::streamSynchronize();
    }

    auto* prims_arr = d_prims.data();

    // Single kernel launch over all ghost points on this level
    auto* copy = this;
    const int total_ngps = gpview.total_ngps;
    const int num_lframes = static_cast<int>(LocalFrame_a.size());

    ParallelFor(total_ngps, [=] AMREX_GPU_DEVICE (int ii) noexcept
    {
      // Determine which FAB this GP belongs to
      int ifab = gpview.gp_fab[ii];
      auto prims = prims_arr[ifab];  // read image points & write ghost cells

      // 1) Reconstruct local orthonormal frame — rotate from body to world frame
      int elem_idx = gpview.elemIdx[ii];
      const auto& frame = lf_ptr[elem_idx];

      // Find which geometry this GP belongs to (for transform lookup)
      int geomIdx = gpview.geomIdx[ii];
      const auto& xform = transforms[geomIdx];

      int type_solid_bc = 0;

      // Rotate body-frame vectors to world frame
      amrex::Real n_w[AMREX_SPACEDIM], t1_w[AMREX_SPACEDIM], t2_w[AMREX_SPACEDIM];
      xform.rotate_to_world(frame.normal,   n_w);
      xform.rotate_to_world(frame.tangent1, t1_w);
#if (AMREX_SPACEDIM == 3)
      xform.rotate_to_world(frame.tangent2, t2_w);
#endif

      Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          nvec(d)  = n_w[d];
          t1vec(d) = t1_w[d];
#if (AMREX_SPACEDIM == 3)
          t2vec(d) = t2_w[d];
#else
          t2vec(d) = Real(0.0);
#endif
      }

      // 2) Storage for primitive variables along the normal
      Array2D<Real, 0, eorder_tparm + 1, 0, cls_t::NPRIM - 1> primsNormal;
      for (int p = 0; p <= eorder_tparm + 1; ++p) {
          for (int n = 0; n < cls_t::NPRIM; ++n) {
              primsNormal(p,n) = Real(0.0);
          }
      }

      // 3) Interpolate primitive variables at all image points
      copy->template interpolateIMs<eorder_tparm, iorder_tparm>(
          gpview.imp_ip_ijk[ii], gpview.imp_ipweights[ii], prims, primsNormal);

      // 4) Transform velocities at image points to local frame
      for (int iip = 2; iip < 2 + eorder_tparm; ++iip) {
          copy->template global2local<eorder_tparm>(iip, primsNormal, nvec, t1vec, t2vec);
      }

      // 5) Apply wall model at IB surface
      ibm_detail::dispatch_compute_surfIB<eorder_tparm, wallmodel>(
          gpview.ib_xyz[ii], nvec, t1vec, t2vec, primsNormal, type_solid_bc, cls);

      // 6) Extrapolate from surface/image points back to ghost point
      copy->template extrapolate<eorder_tparm>(
          primsNormal, gpview.imp_ninterp[ii], gpview.disGP[ii], gpview.disIM[ii]);

      // 7) Transform ghost-point velocity back to global coordinates
      int idx = 0;
      copy->template local2global<eorder_tparm>(idx, primsNormal, nvec, t1vec, t2vec);

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

      // 9) Enforce thermodynamic consistency
      Real Q[cls_t::NPRIM];
      cls->ensurePTYfillq(P, T, Y, ux, uy, uz, Q);

      // 10) Write ghost-cell primitive variables back into prims
      int i = gpview.gp_ijk[ii](0);
      int j = gpview.gp_ijk[ii](1);
#if (AMREX_SPACEDIM == 3)
      int k = gpview.gp_ijk[ii](2);
#else
      int k = 0;
#endif

      for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(i,j,k,n) = Q[n];
      }
    }); // end ParallelFor over all ghost points
  }

  /**
   * \brief Fix conservative state for cells freshly exposed by geometry movement.
   *
   * When the immersed body moves, some cells transition from solid to fluid.
   * Their conservative values are stale (whatever was there when they were solid)
   * and may produce NaN when converted to primitive variables.
   *
   * Strategy:
   *   Pass 0: zero-initialise density of freshly-exposed cells so pass 2 can
   *           reliably distinguish "already fixed" from "stale solid data."
   *   Pass 1: average from immediate fluid neighbours (4/6-connected).
   *   Pass 2: wider 2-ring search for cells that had no valid pass-1 donors.
   *
   * \param old_markers  Snapshot of markers BEFORE geometry rebuild (comp 0 only, same BA/DM).
   * \param state_mf     Conservative state MultiFab to repair (valid + ghost cells).
   * \param lev          AMR level.
   */
  void fixExposedCells(const FabArray<BaseFab<uint8_t>>& old_markers,
                       MultiFab& state_mf,
                       int lev)
  {
    BL_PROFILE("IBM::fixExposedCells");
    auto& mfab = *bmf_a[lev];
    const int ncons = cls_t::NCONS;

    // Pass 0: zero-initialise freshly-exposed cells so that pass 2 can
    // distinguish "fixed by pass 1" (URHO > 0) from "still needs fixing."
    // Without this, stale solid-cell data (which may have URHO > 1e-10)
    // would cause pass 2 to skip cells that were never handled by pass 1.
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& old_mk = old_markers.const_array(mfi);
      auto const& new_mk = mfab.const_array(mfi);
      auto const& state  = state_mf.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (old_mk(i,j,k,0) != 0 && new_mk(i,j,k,0) == 0) {
          for (int n = 0; n < ncons; ++n)
            state(i,j,k,n) = Real(0.0);
        }
      });
    }

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& old_mk = old_markers.const_array(mfi);   // old comp-0 markers
      auto const& new_mk = mfab.const_array(mfi);          // new markers (comp 0)
      auto const& state  = state_mf.array(mfi);

      // Pass 1: fill freshly-exposed cells from valid fluid neighbours
      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        // Only process cells that changed from solid → fluid
        if (old_mk(i,j,k,0) == 0) return;   // was already fluid
        if (new_mk(i,j,k,0) != 0) return;   // still solid

        // Average conservative state from immediate fluid neighbours
        // (cells that are fluid in BOTH old and new markers → reliable data)
        Real sum[cls_t::NCONS] = {};
        int count = 0;

        // 6-connected (2D: 4-connected) neighbourhood
        constexpr int offsets[][3] = {
          {-1,0,0},{1,0,0},{0,-1,0},{0,1,0}
#if (AMREX_SPACEDIM == 3)
          ,{0,0,-1},{0,0,1}
#endif
        };

        for (const auto& off : offsets) {
          int ii = i + off[0], jj = j + off[1];
          AMREX_D_TERM(;, ;, int kk = k + off[2];)
#if (AMREX_SPACEDIM == 2)
          int kk = 0;
#endif
          // Donor must be fluid in old state (has valid conservative data)
          // and also fluid in new state (not about to be covered)
          if (old_mk(ii,jj,kk,0) == 0 && new_mk(ii,jj,kk,0) == 0) {
            for (int n = 0; n < ncons; ++n)
              sum[n] += state(ii,jj,kk,n);
            count++;
          }
        }

        if (count > 0) {
          Real inv = Real(1.0) / count;
          for (int n = 0; n < ncons; ++n)
            state(i,j,k,n) = sum[n] * inv;
        }
        // If count==0, all neighbours are also newly exposed or solid.
        // Leave for Pass 2 (wider stencil) below.
      });
    }

    // Pass 2: sweep again for any remaining unfixed cells (all neighbours were also exposed).
    // Use a 2-ring search. This is very rare for typical motions.
    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& old_mk = old_markers.const_array(mfi);
      auto const& new_mk = mfab.const_array(mfi);
      auto const& state  = state_mf.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (old_mk(i,j,k,0) == 0) return;
        if (new_mk(i,j,k,0) != 0) return;

        // Check if this cell was already fixed in pass 1
        // (density should be reasonable if fixed)
        if (state(i,j,k, cls_t::URHO) > Real(1.0e-10)) return;

        // Wider search: 2-ring
        Real sum[cls_t::NCONS] = {};
        int count = 0;
        for (int dj = -2; dj <= 2; ++dj) {
          for (int di = -2; di <= 2; ++di) {
#if (AMREX_SPACEDIM == 3)
            for (int dk = -2; dk <= 2; ++dk) {
#else
            { int dk = 0;
#endif
              if (di == 0 && dj == 0 && dk == 0) continue;
              int ii = i+di, jj = j+dj, kk = k+dk;
              // Any cell that is now fluid with reasonable density
              if (new_mk(ii,jj,kk,0) == 0 && state(ii,jj,kk, cls_t::URHO) > Real(1.0e-10)) {
                for (int n = 0; n < ncons; ++n)
                  sum[n] += state(ii,jj,kk,n);
                count++;
              }
            }
          }
        }
        if (count > 0) {
          Real inv = Real(1.0) / count;
          for (int n = 0; n < ncons; ++n)
            state(i,j,k,n) = sum[n] * inv;
        }
      });
    }
  }

  /**
    * \brief Compute surface indices and interpolation data for all faces/edges at given level
    *
    * Algorithm:
    *  1. Build spatial lookup (global_fab_idx -> local_fab_idx)
    *  2. For each face: compute mirror point, find owning FAB, compute interpolation weights
    *  3. Build CSR structure for GPU-friendly access
    *
    * \param lev AMR level
    */
  void computeSurfIndices(int lev)
  {
    BL_PROFILE("IBM::computeSurfIndices");

    int myrank = amrex::ParallelDescriptor::MyProc();
    amrex::Print()  << "Compute Surface Index at LEVEL " << lev  << std::endl;

    auto& mfab = *bmf_a[lev];

    const BoxArray& ba_global     = mfab.boxArray();
    const DistributionMapping& dm = mfab.DistributionMap();

    const int nfab_local  = mfab.local_size();
    const int nfab_global = ba_global.size();

    const auto prob_lo = amr_p->Geom(lev).ProbLoArray();
    const auto& domain = amr_p->Geom(lev).Domain();

    // ========================================================================
    // Phase 0: Initialize surfimp_soa / surfphys_soa
    // ========================================================================
    if (lev == amr_p->finestLevel()) {
      if ((surfimp_soa.elemIdx.size() != ntotalfaces) || (surfphys_soa.elemIdx.size() != ntotalfaces)) {
        surfimp_soa.resize(ntotalfaces);
        surfphys_soa.resize(ntotalfaces);
      }
      surfphys_soa.reset();
    }

    // ========================================================================
    // Phase 1: Build per-FAB lookup structures
    // ========================================================================
    Vector<Array4<uint8_t const>> fab_markers(nfab_local);
    Vector<Box> fab_bxg(nfab_local);
    Vector<int> global_to_local_fab(nfab_global, -1);

    for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi) {
        int lidx = mfi.LocalIndex();
        int gidx = mfi.index();
        fab_markers[lidx] = mfab.const_array(mfi);
        fab_bxg[lidx] = mfi.growntilebox(cls_t::NGHOST);
        global_to_local_fab[gidx] = lidx;
    }

    // ========================================================================
    // Phase 2: Fast locate — assign each face to a FAB (lightweight)
    //
    // This pass only does centroid→cell→FAB lookup and writes metadata.
    // The heavy image-point computation is deferred to Phase 3.
    // ========================================================================
    // Per-face FAB assignment: >=0 means locally owned, -1 means not local
    std::vector<int> face_local_fab(ntotalfaces, -1);

    int faces_found = 0;
    int faces_out_domain = 0;
    int faces_out_level  = 0;
    int faces_out_rank   = 0;
    int faces_in_finer   = 0;

    // Precompute inverse dx for faster centroid→cell conversion
    amrex::GpuArray<Real, AMREX_SPACEDIM> inv_dx;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        inv_dx[d] = Real(1.0) / dx_a[lev][d];
    }

    {
      std::vector<std::pair<int, Box>> isects;

      for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {

        // Skip faces already assigned to a finer level
        if (surfphys_soa.elemfound[f_idx] && surfphys_soa.lev[f_idx] > lev) {
            faces_in_finer++;
            continue;
        }

        // Convert body-frame centroid to world frame, then to cell index
        const SurfElem& surfelem = SurfElem_a[f_idx];
        int gIdx = getGeomIdx(f_idx);
        const auto& T = transform_a[gIdx];
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(surfelem.centroid[0], surfelem.centroid[1]);
#else
        Point c_body(surfelem.centroid[0], surfelem.centroid[1], surfelem.centroid[2]);
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = surfelem.centroid[dd];
#endif
        Point c_world = T.to_world(c_body);

        IntVect iv;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv[d] = static_cast<int>(std::floor(
                    (c_world[d] - prob_lo[d]) * inv_dx[d]));
        }

        if (!domain.contains(iv)) { faces_out_domain++; continue; }

        ba_global.intersections(Box(iv, iv), isects, true, IntVect::TheZeroVector());
        if (isects.empty()) { faces_out_level++; continue; }

        int gidx       = isects[0].first;
        int owner_rank = dm[gidx];
        int local_fab  = global_to_local_fab[gidx];

        surfphys_soa.ifab[f_idx] = local_fab;
        surfphys_soa.rank[f_idx] = owner_rank;
        surfphys_soa.lev[f_idx]  = lev;

        if (owner_rank != myrank) {
            surfphys_soa.elemfound[f_idx] = -1;
            faces_out_rank++;
            continue;
        }

        surfphys_soa.elemfound[f_idx] = 1;
        face_local_fab[f_idx] = local_fab;
        faces_found++;
      }
    } // isects freed here

    // ========================================================================
    // Phase 3: Group locally-owned faces by FAB for cache locality
    //
    // Faces in the same FAB access the same ibMarkers array.  Processing
    // them together keeps marker data in cache and avoids thrashing when
    // calling search_optimal_image_point / computeIPweights.
    // ========================================================================
    std::vector<std::vector<int>> fab_faces(nfab_local);
    for (int f_idx = 0; f_idx < ntotalfaces; ++f_idx) {
        if (face_local_fab[f_idx] >= 0) {
            fab_faces[face_local_fab[f_idx]].push_back(f_idx);
        }
    }

    // ========================================================================
    // Phase 4: Compute image points and interpolation weights per FAB
    //
    // Each face writes exclusively to its own f_idx slot in surfimp_soa /
    // surfphys_soa — no cross-face data dependencies.  The inner helpers
    // (search_optimal_image_point, computeIPweights) are stateless and
    // thread-safe, so the outer FAB loop is safe to parallelize with OpenMP.
    // ========================================================================
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int lfab = 0; lfab < nfab_local; ++lfab) {
      const auto& flist = fab_faces[lfab];
      if (flist.empty()) continue;

      // Pin the marker array and box for this FAB — all faces in flist
      // will read from this same data, maximizing cache reuse.
      auto const ibMarkers = fab_markers[lfab];
      auto const bxg       = fab_bxg[lfab];

      for (int f_idx : flist) {
        surfimp_soa.elemIdx[f_idx]  = f_idx;
        surfphys_soa.elemIdx[f_idx] = f_idx;

        // Rotate body-frame LocalFrame to world frame
        const LocalFrame& lf_body = LocalFrame_a[f_idx];
        const SurfElem&   surfelem = SurfElem_a[f_idx];
        int gIdx = getGeomIdx(f_idx);
        const auto& T = transform_a[gIdx];

        LocalFrame localframe;
        T.rotate_to_world(lf_body.normal,   localframe.normal);
        T.rotate_to_world(lf_body.tangent1, localframe.tangent1);
#if (AMREX_SPACEDIM == 3)
        T.rotate_to_world(lf_body.tangent2, localframe.tangent2);
#endif

        // Transform body-frame centroid to world frame
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(surfelem.centroid[0], surfelem.centroid[1]);
#else
        Point c_body(surfelem.centroid[0], surfelem.centroid[1], surfelem.centroid[2]);
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = surfelem.centroid[dd];
#endif
        Point surf_centroid = T.to_world(c_body);

        // Temporary storage for this face's image point data
        Array2D<Real, 0, eorder_tparm_surf - 1, 0, IDIM> imp_xyz;
        Array2D< int, 0, eorder_tparm_surf - 1, 0, IDIM> imp_ijk;
        Array1D<Real, 0, eorder_tparm_surf - 1> disIM;
        Array1D< int, 0, eorder_tparm_surf - 1> imp_ninterp;

        for (int jj = 0; jj < eorder_tparm_surf; jj++) {
          Point cp_start;
          if (jj == 0) {
            cp_start = surf_centroid;
          } else {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
            cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1));
#else
            cp_start = Point(imp_xyz(jj - 1, 0), imp_xyz(jj - 1, 1), imp_xyz(jj - 1, 2));
#endif
#else
            for (int d = 0; d < AMREX_SPACEDIM; ++d) cp_start[d] = imp_xyz(jj - 1, d);
#endif
          }

          if (jj == 0) {
            search_optimal_image_point<eorder_tparm_surf, iorder_tparm_surf>(
                                          cp_start, localframe,
                                          lev, prob_lo, dx_a[lev], di_a_surf[lev],
                                          bxg, ibMarkers,
                                          surfimp_soa, f_idx,
                                          imp_xyz, imp_ijk, disIM, imp_ninterp);
          }
          else {
            search_image_point<eorder_tparm_surf, iorder_tparm_surf>(
                                      jj, cp_start, localframe,
                                      lev, prob_lo, dx_a[lev], di_a_surf[lev],
                                      bxg, ibMarkers,
                                      surfimp_soa, f_idx,
                                      imp_xyz, imp_ijk, disIM, imp_ninterp);
          }
        } // end loop on image points

        surfimp_soa.imp_xyz[f_idx]     = imp_xyz;
        surfimp_soa.imp_ijk[f_idx]     = imp_ijk;
        surfimp_soa.disIM[f_idx]       = disIM;
        surfimp_soa.imp_ninterp[f_idx] = imp_ninterp;
        surfphys_soa.ip_quality[f_idx] = imp_ninterp(0);

        Array3D< int, 0, eorder_tparm_surf - 1, 0, N_InterP_surf - 1, 0, IDIM> imp_ip_ijk;
        Array2D<Real, 0, eorder_tparm_surf - 1, 0, N_InterP_surf - 1> imp_ipweights;

        computeIPweights<eorder_tparm_surf, iorder_tparm_surf, SURFIMP>(
            imp_ipweights, imp_ip_ijk,
            imp_xyz, imp_ijk, imp_ninterp,
            prob_lo, dx_a[lev], ibMarkers);

        surfimp_soa.imp_ip_ijk[f_idx]    = imp_ip_ijk;
        surfimp_soa.imp_ipweights[f_idx] = imp_ipweights;
      } // end loop over faces in this FAB
    } // end loop over FABs

    // Build CSR at the coarsest level (surface is built from finest to coarsest)
    if (lev == 0) buildCSR();

  }

  /**
   * \brief Computes surface properties (pressure, temperature, gradients) for each surface face.
   *
   * This function iterates over all surface element owned by the current process. For each face:
   * 1. Interpolates primitive variables at image points using the pre-computed weights.
   * 2. Applies the wall model to determine surface state (e.g., no-slip, adiabatic/isothermal).
   * 3. Computes gradients (e.g., dT/dn) at the surface.
   * 4. Stores the results back into the Surface Data SoA structure.
   *
   * \param stateprops MultiFab containing the fluid state properties.
   * \param cls        Pointer to the physics/closure class.
   * \param lev        Current AMR level.
   */
  void computeSURFs(MultiFab& prims_mf, const cls_t* cls, int lev) 
  {
    BL_PROFILE("IBM::computeSURFs");

    // Skip if no CSR data for this level
    if (lev >= static_cast<int>(faces_per_level.size())) return;
    auto& csr = faces_per_level[lev];
    const int nfaces_local = static_cast<int>(csr.face_indices.size());
    if (nfaces_local == 0) return;

    amrex::Print() << "computeSURFs at LEVEL " << lev
                   << " (" << nfaces_local << " faces)" << std::endl;

    auto& mfab = *bmf_a[lev];
    const int nfabs_local = mfab.local_size();

    // Build device array of Array4<Real> pointers (one per local FAB)
    Gpu::DeviceVector<Array4<Real>> d_prims(nfabs_local);
    {
      Vector<Array4<Real>> h_prims(nfabs_local);
      int ifab = 0;
      for (MFIter mfi(mfab, false); mfi.isValid(); ++mfi, ++ifab) {
        h_prims[ifab] = prims_mf.array(mfi);
      }
      Gpu::copyAsync(Gpu::hostToDevice, h_prims.begin(), h_prims.end(), d_prims.begin());
      Gpu::streamSynchronize();
    }

    auto* prims_arr = d_prims.data();

    // Device-accessible SoA pointers
    const int* d_face_indices  = csr.face_indices.data();

    auto* sp_pressure    = surfphys_soa.pressure.data();
    auto* sp_temperature = surfphys_soa.temperature.data();
    auto* sp_dTdn        = surfphys_soa.dTdn.data();
    auto* sp_tau1        = surfphys_soa.tau1.data();
    auto* sp_tau2        = surfphys_soa.tau2.data();
    const int* sp_ifab   = surfphys_soa.ifab.data();

    const auto* si_imp_ip_ijk    = surfimp_soa.imp_ip_ijk.data();
    const auto* si_imp_ipweights = surfimp_soa.imp_ipweights.data();
    const auto* si_disIM         = surfimp_soa.disIM.data();

    auto const* lf_ptr = LocalFrame_a.data();
    auto const* se_ptr = SurfElem_a.data();

    // Capture transforms and geom_offsets for body→world rotation
    AMREX_ALWAYS_ASSERT(ngeom <= MAX_NGEOM);
    GpuArray<RigidTransform, MAX_NGEOM> transforms;
    GpuArray<int, MAX_NGEOM + 1> geom_off;
    for (int ii = 0; ii < ngeom; ii++) {
        transforms[ii] = transform_a[ii];
        geom_off[ii]   = geom_offsets[ii];
    }
    geom_off[ngeom] = geom_offsets[ngeom];
    const int ngeom_local = ngeom;

    auto* copy = this;

    ParallelFor(nfaces_local, [=] AMREX_GPU_DEVICE (int ii) noexcept
    {
      // Global face index from CSR
      int f_idx = d_face_indices[ii];
      int ifab = sp_ifab[f_idx];
      auto prims = prims_arr[ifab];

      // Find geometry index for this face (to get the right transform)
      int gIdx = 0;
      for (int g = 0; g < ngeom_local; ++g) {
          if (f_idx >= geom_off[g] && f_idx < geom_off[g + 1]) { gIdx = g; break; }
      }
      const auto& T = transforms[gIdx];

      // 1) Reconstruct local orthonormal frame — rotate body→world
      const auto& frame = lf_ptr[f_idx];

      amrex::Real n_w[AMREX_SPACEDIM], t1_w[AMREX_SPACEDIM], t2_w[AMREX_SPACEDIM];
      T.rotate_to_world(frame.normal,   n_w);
      T.rotate_to_world(frame.tangent1, t1_w);
#if (AMREX_SPACEDIM == 3)
      T.rotate_to_world(frame.tangent2, t2_w);
#endif

      Array1D<Real, 0, AMREX_SPACEDIM - 1> nvec, t1vec, t2vec;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          nvec(d)  = n_w[d];
          t1vec(d) = t1_w[d];
#if (AMREX_SPACEDIM == 3)
          t2vec(d) = t2_w[d];
#else
          t2vec(d) = Real(0.0);
#endif
      }

      // Surface centroid coordinates — transform body→world
      Array1D<Real, 0, AMREX_SPACEDIM - 1> xyz;
      {
#ifdef AMREX_USE_CGAL
#if (AMREX_SPACEDIM == 2)
        Point c_body(se_ptr[f_idx].centroid[0], se_ptr[f_idx].centroid[1]);
#else
        Point c_body(se_ptr[f_idx].centroid[0], se_ptr[f_idx].centroid[1], se_ptr[f_idx].centroid[2]);
#endif
        Point c_world = T.to_world(c_body);
        xyz(0) = c_world.x(); xyz(1) = c_world.y();
#if (AMREX_SPACEDIM == 3)
        xyz(2) = c_world.z();
#endif
#else
        Point c_body;
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) c_body[dd] = se_ptr[f_idx].centroid[dd];
        Point c_world = T.to_world(c_body);
        for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) xyz(dd) = c_world[dd];
#endif
      }

      // 2) Zero-initialize primsNormal
      Array2D<Real, 0, eorder_tparm_surf + 1, 0, cls_t::NPRIM - 1> primsNormal;
      for (int p = 0; p <= eorder_tparm_surf + 1; ++p) {
          for (int n = 0; n < cls_t::NPRIM; ++n) {
              primsNormal(p,n) = Real(0.0);
          }
      }

      // 3) Interpolate primitive variables at image points
      copy->template interpolateIMs<eorder_tparm_surf, iorder_tparm_surf>(
          si_imp_ip_ijk[f_idx], si_imp_ipweights[f_idx], prims, primsNormal);

      // 4) Transform velocities at image points to local frame
      for (int iip = 2; iip < 2 + eorder_tparm_surf; ++iip) {
          copy->template global2local<eorder_tparm_surf>(iip, primsNormal, nvec, t1vec, t2vec);
      }

      // 5) Apply wall model at IB surface (fills primsNormal slot 1)
      int type_solid_bc = 0;
      ibm_detail::dispatch_compute_surfIB<eorder_tparm_surf, wallmodel>(
          xyz, nvec, t1vec, t2vec, primsNormal, type_solid_bc, cls);

      // 6) Extract surface quantities
      Real P_surf = primsNormal(1, cls_t::QPRES);
      Real T_surf = primsNormal(1, cls_t::QT);

      // One-sided gradient: dT/dn = (T_image - T_surface) / distance
      Real dis = si_disIM[f_idx](0);
      Real inv_dis = (dis > Real(0.0)) ? Real(1.0) / dis : Real(0.0);
      Real dTdn = (primsNormal(2, cls_t::QT) - T_surf) * inv_dis;

      // Wall shear stress: tau = mu * du_tangential / dn
      // (in local frame: QU=normal, QV=tangent1, QW=tangent2)
      Real mu = cls->visc(T_surf);
      Real tau1 = mu * (primsNormal(2, cls_t::QV) - primsNormal(1, cls_t::QV)) * inv_dis;
#if (AMREX_SPACEDIM == 3)
      Real tau2 = mu * (primsNormal(2, cls_t::QW) - primsNormal(1, cls_t::QW)) * inv_dis;
#else
      Real tau2 = Real(0.0);
#endif

      // 7) Store results
      sp_pressure[f_idx]    = P_surf;
      sp_temperature[f_idx] = T_surf;
      sp_dTdn[f_idx]        = dTdn;
      sp_tau1[f_idx]        = tau1;
      sp_tau2[f_idx]        = tau2;
    }); // end ParallelFor

    // Ensure GPU writes are visible to CPU before gatherSurfData / plotSURF
    Gpu::streamSynchronize();
  }


//============================================================================
///--------------------------- private functions -----------------------------
private:

  // Interpolation, extrapolation, coordinate-transform helpers
  #include "ibm_solver_interp.h"

  // Geometry I/O, VTK output, MPI gather, CSR builder
  #include "ibm_solver_io.h"

}; // end class ibm_solver_t

#endif // IBM_SOLVER_H_
