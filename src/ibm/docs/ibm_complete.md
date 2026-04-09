# IBM Module — Complete Reference

This document consolidates all IBM documentation for `src/ibm` (CERISSE project, BVH backend, March 2026).

---

## 1. Code Layout

The IBM code is modularized around the custom BVH backend:

| File | Role |
| --- | --- |
| `bvh_defs.h` | Pure data definitions: `Point`, `Vec`, `BoundedSide`, `AABB`, `ClosestPointResult`, `BVHNode`, `TriMesh`/`Polygon2D`, `LocalFrame`, `SurfElem` |
| `bvh.h` | BVH algorithms: Morton encoding, primitive closest-point utilities, BVH build, BVH closest-point queries on CPU/GPU |
| `eib_bvh.h` | IBM/BVH bridge: `InsideTester` + `InsideTesterView`, geometry readers (STL/OFF/polygon), cache builders, geometry consistency checks |
| `eib_data.h` | IBM solver SoA data containers and constants: `gpData_t`, `surfImp_t`, `surfPhys_t`, `FaceCSR`, interpolation thresholds and factors |
| `eib.h` | Main solver class template `eib_t<wallmodel,param,cls_t>` and runtime kernels: markers, ghost-point setup, GP reconstruction, surface reconstruction/output |
| `eib_io.h` | Geometry-loading stage called by `eib_t::init`: reads `ib.filename`, builds BVH, initializes in/out testers, builds flattened surface caches |
| `eib_interp.h` | Interpolation/image-point search helpers used by GP and SURF reconstruction |
| `IBMultiFab.h` | IBM-specific FAB/MultiFab containers combining marker fields and per-FAB ghost-point SoA payload |
| `ib_walltypes.h` | Wall model templates used by `compute_surfIB` dispatch |
| `eib_geometry_bvh.h` | Compatibility forwarding header → `eib_bvh.h` |
| `eib_bvh_geom.h` | Compatibility forwarding header → `eib_bvh.h` |
| `bvh_spatial_index.h` | Compatibility forwarding header → `bvh.h` |

---

## 2. CGAL → BVH Type Mapping

One-to-one type mappings between the original CGAL-based geometry stack and the custom BVH implementation.

| Functionality | CGAL Type | Custom BVH Type | Notes |
| --- | --- | --- | --- |
| Point representation | `Point = Kernel::Point_2 / Point_3` | `Point = GpuArray<Real, AMREX_SPACEDIM>` | Same semantic role; custom type is GPU-friendly fixed-size storage. |
| Vector representation | `Vector_CGAL = Kernel::Vector_2 / Vector_3` | `Vec = GpuArray<Real, AMREX_SPACEDIM>` | Same geometric meaning (direction/difference vectors). |
| Inside/outside enum | `CGAL::Bounded_side` | `enum class BoundedSide` | Preserves the same three states: inside, boundary, outside. |
| Bounding box type | `Bbox = CGAL::Bbox_2 / Bbox_3` | `AABB` (alias `Bbox = AABB`) | Equivalent pruning/containment role for spatial queries. |
| 2D geometry container | `GeomType = Polygon2D` (CGAL wrapper) | `GeomType = Polygon2D` (custom) | Both represent polygon geometry in 2D. |
| 3D geometry container | `GeomType = Polyhedron` (`CGAL::Polyhedron_3`) | `GeomType = TriMesh` | Both represent surface meshes; custom stack uses triangle-only explicit indexing. |
| Spatial acceleration structure | `Tree = CGAL::AABB_tree<Traits>` | `BVH` | Both accelerate nearest-point and ray-intersection traversal. |
| Acceleration internals | `Primitive + Traits + Tree internals` | `BVHNode + BVH::nodes` | CGAL hides node internals; custom stack exposes explicit node arrays. |
| Closest-point query result | `Point_and_primitive_id` | `ClosestPointResult` | Both return closest point and primitive id; custom result also stores distance directly. |
| Primitive identifier | `PrimitiveID = Tree::Primitive_id` | `int prim_id` | Same purpose: identifies edge/face primitive hit by query. |
| Element descriptor | `elm_descriptor` (iterator or `face_descriptor`) | Direct integer indexing | Custom stack replaces handle/descriptor abstractions with linear ids. |
| Inside tester functor | `inside_t = Side_of_triangle_mesh` (or 2D wrapper) | `inside_t = InsideTester` | Unified callable interface retained: `operator()(Point)`. Uses BVH-backed queries in 3D and BVH-accelerated culling in 2D. |
| Primitive-id map | `PrimitiveIndexMap = std::map<PrimitiveID, int>` | Not required | Removed because custom pipeline uses direct integer primitive indexing. |

---

## 3. Initialization Flow

### Entry Point

```
IBM::ib.init(&amr)
  └── eib_t::init()           [eib.h]
        └── read_geom()       [eib_io.h]
  (later, in AMR regrid hooks)
        └── computeMarkers()  [eib.h]
        └── initialiseGPs()   [eib.h]
```

### Stage A — `eib_t::init` (Level Metrics)

- Stores `amr_p` and refinement ratios.
- Allocates per-level containers (`bmf_a`, `faces_per_level`).
- Computes per-level cell sizes Δx, cell diagonal h, and image-point distances:

  d_GP = α · h,  d_surf = α_surf · h

- Calls `read_geom()`.

### Stage B — `read_geom` (Geometry Load + BVH Build)

Geometry filenames are read from `ib.filename`. For each geometry index i:

1. **Read geometry**
   - 2D: `read_polygon_2d(files_a[i], geom_a[i], min_dx)` — subdivides edges so all segment lengths ≤ min_dx/2
   - 3D: `read_mesh(files_a[i], geom_a[i])` — reads STL/OFF/VTK
2. **Build BVH**: `bvh_a[i].build(geom_a[i])` — Morton-sorted median split
3. **Construct inside/outside tester**:
   - Both 2D and 3D: `inside_t(geom_a[i], bvh_a[i])`
4. **Store bbox**: `bbox_a[i] = geom_a[i].bbox()`
5. **Append flattened per-element cache**: `SurfElem_a`, `LocalFrame_a`, `geom_offsets[i]`

After all geometries: validate cache sizes, run geometry consistency checks, optionally flip normals if `interior_is_solid = false`.

### Stage C — `computeMarkers` (Solid/Ghost Classification)

Fills per-cell marker components (comp 0: fluid/solid id, comp 1: ghost-cell id):

1. Solid pass: `InsideTesterView` in `ParallelFor` with bbox rejection — GPU path available.
2. Ghost pass: flags interface-adjacent solids.
3. Prefix-sum compaction extracts `gp_ijk` list.

### Stage D — `initialiseGPs` (Per-Ghost Geometry Metadata)

For each ghost point:

1. Determine owning geometry index.
2. BVH closest-point query → IB foot point and primitive id.
3. Map primitive id to global face/edge index via `geom_offsets`.
4. Compute distance to boundary and local frame.
5. Place image points and build interpolation stencils/weights.
6. Store into `gpData_t` SoA arrays.

### Data State After Initialization

| Object | Description |
| --- | --- |
| `geom_a` | Geometry containers (`TriMesh` or `Polygon2D`) |
| `bvh_a` | Per-geometry BVH trees |
| `inout_fa` | Inside/outside testers (`InsideTester*`) |
| `bbox_a` | Coarse geometry bounding boxes |
| `SurfElem_a` | Flattened surface element descriptors |
| `LocalFrame_a` | Flattened normal/tangent basis per element |
| `geom_offsets` | Start index per geometry in flattened arrays |
| `IBMultiFab` (per level) | Marker field + per-FAB `gpData_t` payload |

---

## 4. Notes on 2D vs 3D

- **3D**: `InsideTesterView` uses BVH-guided ray/closest-point logic on triangle faces.
- **2D**: `InsideTesterView` uses BVH-accelerated polygon edge culling for both boundary checks and ray-cast inside/outside tests. This was added to handle subdivided polygons (e.g., 40 → 736 vertices after `read_polygon_2d`) which caused GPU hangs with the previous O(N) linear scan.
- Both dimensions share the same dataflow in `eib_t`, with geometry-specific logic selected by `AMREX_SPACEDIM` at compile time.

---

## 5. Primary Runtime Path (after init)

```
RHS call
  └── computeGPs()          reconstruct ghost-point flow state
  └── computeSURFs()        reconstruct surface quantities
  └── plotSURF()            write surface output (drag, heat flux)
```

Surface post-processing integrates over `SurfElem_a` using precomputed `LocalFrame_a`. The polygon subdivision in Stage B ensures sufficient resolution for drag and heat flux integration.
