# IBM BVH Type Mapping

This document lists one-to-one type mappings between the original CGAL-based IBM geometry stack and the custom BVH-based implementation, focusing on equivalent functionality.

Current source placement:
- Data types: `src/ibm/bvh_defs.h`
- BVH algorithms: `src/ibm/bvh.h`
- IBM/BVH bridge and geometry readers: `src/ibm/eib_bvh.h`
- Solver data containers: `src/ibm/eib_data.h`
- Geometry loading flow: `src/ibm/eib_io.h`

| Functionality | CGAL Type | Custom BVH Type | Mapping Notes |
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
| Inside tester functor | `inside_t = Side_of_triangle_mesh` (or 2D wrapper) | `inside_t = InsideTester` | Unified callable interface retained: `operator()(Point)`. Current path uses BVH-backed queries in 3D and BVH-accelerated culling for 2D polygon checks. |
| Primitive-id map | `PrimitiveIndexMap = std::map<PrimitiveID, int>` | Not required | Removed because custom pipeline uses direct integer primitive indexing. |
