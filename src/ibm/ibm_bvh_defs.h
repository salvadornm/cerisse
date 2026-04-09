#ifndef IBM_BVH_DEFS_H_
#define IBM_BVH_DEFS_H_

// ============================================================================
// ibm_bvh_defs.h — BVH and geometry data structure definitions
//
// Contains all pure data types used by the BVH acceleration module:
//   1. Primitive types    : Point, Vec, BoundedSide, AABB, ClosestPointResult
//   2. BVH node           : BVHNode
//   3. Geometry containers: TriMesh (3D) or Polygon2D (2D)
//   4. Surface/frame types: LocalFrame, SurfElem
// ============================================================================

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>

#include <limits>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <cstdint>
#include <functional>

// NOTE: Do NOT use 'using namespace amrex;' in headers — it pollutes every
// translation unit that includes this file.  All AMReX types are qualified
// explicitly with amrex:: below.

// ============================================================================
// 1. PRIMITIVE TYPES
// ============================================================================

/// Point represented as GpuArray<Real, DIM>
using Point = amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>;

/// Construct a Point from cell-centred grid coordinates (i,j,k).
/// Accepts any array-like type for prob_lo and dx (e.g. GpuArray, Real*).
template <typename ArrayLike>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Point make_grid_point (const ArrayLike& prob_lo, const ArrayLike& dx,
                       int i, int j, int k) {
    Point p;
    p[0] = prob_lo[0] + (amrex::Real(0.5) + amrex::Real(i)) * dx[0];
    p[1] = prob_lo[1] + (amrex::Real(0.5) + amrex::Real(j)) * dx[1];
#if (AMREX_SPACEDIM == 3)
    p[2] = prob_lo[2] + (amrex::Real(0.5) + amrex::Real(k)) * dx[2];
#else
    amrex::ignore_unused(k);
#endif
    return p;
}

/// Vector (same storage as Point, different semantics)
using Vec   = amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>;

/// Enumeration replacing CGAL::Bounded_side
enum class BoundedSide : int {
    Inside    =  1,   // ON_BOUNDED_SIDE
    OnBoundary = 0,   // ON_BOUNDARY
    Outside   = -1    // ON_UNBOUNDED_SIDE
};

/// Axis-Aligned Bounding Box
struct AABB {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> lo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> hi;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    AABB () {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lo[d] =  std::numeric_limits<amrex::Real>::max();
            hi[d] = -std::numeric_limits<amrex::Real>::max();
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    AABB (const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& lo_,
          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& hi_) : lo(lo_), hi(hi_) {}

    /// Expand to include a point
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void expand (const Point& p) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lo[d] = amrex::min(lo[d], p[d]);
            hi[d] = amrex::max(hi[d], p[d]);
        }
    }

    /// Merge with another AABB
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void merge (const AABB& other) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lo[d] = amrex::min(lo[d], other.lo[d]);
            hi[d] = amrex::max(hi[d], other.hi[d]);
        }
    }

    /// Check if a point is inside the box (inclusive)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool contains (const Point& p) const {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (p[d] < lo[d] || p[d] > hi[d]) return false;
        }
        return true;
    }

    /// Check overlap with another AABB
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool overlaps (const AABB& other) const {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (lo[d] > other.hi[d] || hi[d] < other.lo[d]) return false;
        }
        return true;
    }

    /// Surface area (for SAH heuristic; in 2D returns perimeter)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real surface_area () const {
#if (AMREX_SPACEDIM == 2)
        amrex::Real dx = hi[0] - lo[0];
        amrex::Real dy = hi[1] - lo[1];
        return amrex::Real(2.0) * (dx + dy);
#else
        amrex::Real dx = hi[0] - lo[0];
        amrex::Real dy = hi[1] - lo[1];
        amrex::Real dz = hi[2] - lo[2];
        return amrex::Real(2.0) * (dx*dy + dy*dz + dz*dx);
#endif
    }

    /// Centroid of the box
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Point centroid () const {
        Point c;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            c[d] = amrex::Real(0.5) * (lo[d] + hi[d]);
        }
        return c;
    }
};

/// Closest-point query result (replaces CGAL Point_and_primitive_id)
struct ClosestPointResult {
    Point  point;        ///< Closest point on the surface
    int    prim_id;      ///< Face/edge index (linear)
    amrex::Real distance;     ///< Distance from query to closest point

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ClosestPointResult ()
        : prim_id(-1), distance(std::numeric_limits<amrex::Real>::max())
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) point[d] = amrex::Real(0.0);
    }
};

/// Alias kept for drop-in compatibility
using Bbox = AABB;

// ============================================================================
// 2. BVH NODE
// ============================================================================

struct BVHNode {
    AABB box;
    int  left;     ///< Left child index  (-1 for leaf)
    int  right;    ///< Right child index (-1 for leaf)
    int  prim_id;  ///< Primitive index (only valid for leaf, -1 otherwise)

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_leaf () const { return left == -1 && right == -1; }
};

// ============================================================================
// 3. GEOMETRY CONTAINERS
// ============================================================================

#if (AMREX_SPACEDIM == 3)

/// \brief Simple triangle mesh — replaces CGAL::Polyhedron_3
///
/// Stores vertices as flat Real array and faces as int[3] indices.
/// All data is GPU-friendly (contiguous arrays, no pointers-of-pointers).
struct TriMesh {
    amrex::Gpu::ManagedVector<Point>                   vertices;
    amrex::Gpu::ManagedVector<amrex::GpuArray<int, 3>> faces;     ///< Each face has 3 vertex indices

    int num_vertices () const { return static_cast<int>(vertices.size()); }
    int num_faces    () const { return static_cast<int>(faces.size()); }

    /// Bulk-assign from std::vectors (avoids per-element push_back on ManagedVector)
    void assign (const std::vector<Point>& v, const std::vector<amrex::GpuArray<int,3>>& f) {
        vertices.resize(v.size());
        std::copy(v.begin(), v.end(), vertices.begin());
        faces.resize(f.size());
        std::copy(f.begin(), f.end(), faces.begin());
    }

    int size_of_facets   () const { return num_faces(); }
    int size_of_vertices () const { return num_vertices(); }
    bool is_pure_triangle () const { return true; }
    bool empty () const { return vertices.empty() || faces.empty(); }

    AABB bbox () const {
        AABB box;
        for (const auto& v : vertices) box.expand(v);
        return box;
    }

    AABB face_bbox (int f) const {
        AABB box;
        for (int k = 0; k < 3; ++k) box.expand(vertices[faces[f][k]]);
        return box;
    }

    /// Alias expected by BVH build
    AABB face_aabb (int f) const { return face_bbox(f); }

    // --- Host-only methods (access std::vector, not callable from device) ---

    Point face_centroid (int f) const {
        Point c;
        for (int d = 0; d < 3; ++d) {
            c[d] = (vertices[faces[f][0]][d]
                  + vertices[faces[f][1]][d]
                  + vertices[faces[f][2]][d]) / amrex::Real(3.0);
        }
        return c;
    }

    Vec face_normal_raw (int f) const {
        const auto& p0 = vertices[faces[f][0]];
        const auto& p1 = vertices[faces[f][1]];
        const auto& p2 = vertices[faces[f][2]];
        Vec v1 = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        Vec v2 = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        return Vec{v1[1]*v2[2] - v1[2]*v2[1],
                   v1[2]*v2[0] - v1[0]*v2[2],
                   v1[0]*v2[1] - v1[1]*v2[0]};
    }

    amrex::Real face_area (int f) const {
        Vec n = face_normal_raw(f);
        amrex::Real len2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        return amrex::Real(0.5) * std::sqrt(len2);
    }

private:
    /// Hash functor for undirected edges stored as ordered pairs.
    struct EdgeHash {
        std::size_t operator()(const std::pair<int,int>& e) const {
            std::size_t h1 = std::hash<int>{}(e.first);
            std::size_t h2 = std::hash<int>{}(e.second);
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };

public:
    /// Check if the mesh is closed (watertight).
    /// Every edge must appear in exactly two faces with opposite orientations.
    /// Uses unordered_map for O(n) average-case performance.
    bool is_closed () const {
        std::unordered_map<std::pair<int,int>, int, EdgeHash> edge_count;
        edge_count.reserve(faces.size() * 3);  // At most 3 edges per face

        for (const auto& f : faces) {
            for (int k = 0; k < 3; ++k) {
                int a = f[k], b = f[(k+1)%3];
                auto key = std::make_pair(std::min(a,b), std::max(a,b));
                edge_count[key]++;
            }
        }
        for (const auto& [key, cnt] : edge_count) {
            if (cnt != 2) return false;
        }
        return true;
    }

    /// Ensure outward-facing normals using signed volume test.
    /// Precondition: mesh must be closed (watertight).
    /// If the signed volume is negative, all face winding orders are flipped.
    /// \param skip_closed_check  Pass true when the caller has already verified
    ///        is_closed(), avoiding a redundant O(n) edge traversal.
    void ensure_outward_orientation (bool skip_closed_check = false) {
        if (!skip_closed_check) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(is_closed(),
                "ensure_outward_orientation: mesh must be closed (watertight)");
        }

        amrex::Real vol = amrex::Real(0.0);
        for (const auto& f : faces) {
            const auto& a = vertices[f[0]];
            const auto& b = vertices[f[1]];
            const auto& c = vertices[f[2]];
            vol += a[0]*(b[1]*c[2] - b[2]*c[1])
                 + a[1]*(b[2]*c[0] - b[0]*c[2])
                 + a[2]*(b[0]*c[1] - b[1]*c[0]);
        }
        vol /= amrex::Real(6.0);
        if (vol < amrex::Real(0.0)) {
            for (auto& f : faces) std::swap(f[1], f[2]);
            amrex::Print() << "Info: Reversed face orientations to be outward\n";
        }
    }
};

using GeomType = TriMesh;

#elif (AMREX_SPACEDIM == 2)

/// \brief 2D Polygon — replaces CGAL Polygon2D wrapper
///
/// Stores vertices as a vector of Point (GpuArray<Real,2>).
/// Edges are implicitly (vertex[i] -> vertex[(i+1)%n]).
struct Polygon2D {
    amrex::Gpu::ManagedVector<Point> verts;

    void clear ()                        { verts.clear(); }
    void push_back (const Point& p)      { verts.push_back(p); }
    void assign (const std::vector<Point>& src) {
        verts.resize(src.size());
        std::copy(src.begin(), src.end(), verts.begin());
    }
    std::size_t size () const            { return verts.size(); }
    const Point& vertex (std::size_t i) const { return verts[i]; }

    auto vertices_begin () const { return verts.cbegin(); }
    auto vertices_end   () const { return verts.cend(); }

    Point edge_source (int i) const { return verts[i]; }
    Point edge_target (int i) const { return verts[(i + 1) % verts.size()]; }

    amrex::Real area () const {
        amrex::Real a = amrex::Real(0.0);
        int n = static_cast<int>(verts.size());
        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            a += verts[i][0] * verts[j][1] - verts[j][0] * verts[i][1];
        }
        return a * amrex::Real(0.5);
    }

    bool is_clockwise_oriented () const { return area() < amrex::Real(0.0); }

    void reverse_orientation () { std::reverse(verts.begin(), verts.end()); }

    /// Check if the polygon is simple (no self-intersecting edges).
    /// NOTE: O(n^2) brute-force — acceptable for one-time validation at load,
    ///       but not suitable for per-step use on large polygons.
    bool is_simple () const {
        int n = static_cast<int>(verts.size());
        if (n < 3) return false;
        for (int i = 0; i < n; ++i) {
            int i1 = (i + 1) % n;
            for (int j = i + 2; j < n; ++j) {
                if (i == 0 && j == n - 1) continue;  // adjacent edges share a vertex
                int j1 = (j + 1) % n;
                if (segments_intersect(verts[i], verts[i1], verts[j], verts[j1]))
                    return false;
            }
        }
        return true;
    }

    AABB bbox () const {
        AABB box;
        for (const auto& v : verts) box.expand(v);
        return box;
    }

    AABB edge_aabb (int i) const {
        int j = (i + 1) % static_cast<int>(verts.size());
        AABB box;
        box.expand(verts[i]);
        box.expand(verts[j]);
        return box;
    }

    BoundedSide bounded_side (const Point& p) const {
        int n = static_cast<int>(verts.size());
        int crossings = 0;
        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            const amrex::Real yi = verts[i][1], yj = verts[j][1];
            const amrex::Real xi = verts[i][0], xj = verts[j][0];
            if (point_on_segment(p, verts[i], verts[j])) {
                return BoundedSide::OnBoundary;
            }
            if ((yi <= p[1] && yj > p[1]) || (yj <= p[1] && yi > p[1])) {
                amrex::Real t = (p[1] - yi) / (yj - yi);
                amrex::Real x_cross = xi + t * (xj - xi);
                if (p[0] < x_cross) crossings++;
            }
        }
        return (crossings % 2 == 1) ? BoundedSide::Inside : BoundedSide::Outside;
    }

    void finalize () {}

private:
    static bool segments_intersect (const Point& a0, const Point& a1,
                                    const Point& b0, const Point& b1) {
        amrex::Real d1 = cross2d(b0, b1, a0);
        amrex::Real d2 = cross2d(b0, b1, a1);
        amrex::Real d3 = cross2d(a0, a1, b0);
        amrex::Real d4 = cross2d(a0, a1, b1);
        if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
            ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
            return true;
        return false;
    }

    static amrex::Real cross2d (const Point& a, const Point& b, const Point& c) {
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
    }

    static bool point_on_segment (const Point& p, const Point& a, const Point& b) {
        amrex::Real cross = (b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]);
        amrex::Real len2  = (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]);
        if (len2 < amrex::Real(1e-30)) return false;
        amrex::Real eps = amrex::Real(1e-10) * len2;
        if (cross * cross > eps) return false;
        amrex::Real dx = b[0] - a[0], dy = b[1] - a[1];
        amrex::Real t = ((p[0]-a[0])*dx + (p[1]-a[1])*dy) / len2;
        return t >= amrex::Real(-1e-10) && t <= amrex::Real(1.0 + 1e-10);
    }
};

using Polygon  = Polygon2D;
using GeomType = Polygon2D;

#endif // AMREX_SPACEDIM

// ============================================================================
// 4. SURFACE / FRAME TYPES
// ============================================================================

struct LocalFrame {
    amrex::Real normal[AMREX_SPACEDIM];
    amrex::Real tangent1[AMREX_SPACEDIM];
#if (AMREX_SPACEDIM == 3)
    amrex::Real tangent2[AMREX_SPACEDIM];
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    LocalFrame () {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            normal[i]   = amrex::Real(0.0);
            tangent1[i] = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
            tangent2[i] = amrex::Real(0.0);
#endif
        }
    }

#if (AMREX_SPACEDIM == 3)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    LocalFrame (const amrex::Real n[3], const amrex::Real t1[3], const amrex::Real t2[3]) {
        for (int i = 0; i < 3; ++i) {
            normal[i]   = n[i];
            tangent1[i] = t1[i];
            tangent2[i] = t2[i];
        }
    }
#else
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    LocalFrame (const amrex::Real n[2], const amrex::Real t1[2]) {
        for (int i = 0; i < 2; ++i) {
            normal[i]   = n[i];
            tangent1[i] = t1[i];
        }
    }
#endif
};

struct SurfElem {
    int  geomIdx;
    amrex::Real centroid[AMREX_SPACEDIM];
    amrex::Real measure;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfElem () : geomIdx(-1), measure(amrex::Real(0.0)) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) centroid[i] = amrex::Real(0.0);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfElem (const amrex::Real* c, amrex::Real s, int g_idx) : geomIdx(g_idx), measure(s) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) centroid[i] = c[i];
    }
};

#endif // IBM_BVH_DEFS_H_
