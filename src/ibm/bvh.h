#ifndef BVH_H_
#define BVH_H_

// ============================================================================
// bvh.h — BVH algorithmic code
//
// Implements all algorithmic (functional) BVH logic:
//   1. Morton encoding  : namespace morton (morton_code, expand_bits_*)
//   2. Proximity queries: closest_point_on_segment/triangle, distance helpers
//   3. BVH class        : build() (CPU, Morton-sorted median split),
//                         closest_point_query() (CPU/GPU stack traversal)
// ============================================================================

#include "bvh_types.h"

#include <numeric>
#include <vector>
#include <algorithm>
#include <cstdint>

// ============================================================================
// 1. MORTON ENCODING
// ============================================================================

namespace morton {

/// Expand a 10-bit integer to 30 bits for 3D Morton code
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
uint32_t expand_bits_3d (uint32_t v) {
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

/// Expand a 16-bit integer to 32 bits for 2D Morton code
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
uint32_t expand_bits_2d (uint32_t v) {
    v = (v | (v << 8)) & 0x00FF00FFu;
    v = (v | (v << 4)) & 0x0F0F0F0Fu;
    v = (v | (v << 2)) & 0x33333333u;
    v = (v | (v << 1)) & 0x55555555u;
    return v;
}

/// Compute Morton code for a normalised point [0,1]^DIM
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
uint32_t morton_code (const Point& p, const AABB& scene_box) {
    Point norm;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real range = scene_box.hi[d] - scene_box.lo[d];
        norm[d] = (range > Real(1e-30))
                ? (p[d] - scene_box.lo[d]) / range
                : Real(0.5);
        norm[d] = amrex::max(Real(0.0), amrex::min(Real(1.0), norm[d]));
    }
#if (AMREX_SPACEDIM == 2)
    uint32_t xi = static_cast<uint32_t>(norm[0] * Real(65535.0));
    uint32_t yi = static_cast<uint32_t>(norm[1] * Real(65535.0));
    return (expand_bits_2d(xi) << 1) | expand_bits_2d(yi);
#else
    uint32_t xi = static_cast<uint32_t>(norm[0] * Real(1023.0));
    uint32_t yi = static_cast<uint32_t>(norm[1] * Real(1023.0));
    uint32_t zi = static_cast<uint32_t>(norm[2] * Real(1023.0));
    return (expand_bits_3d(xi) << 2) | (expand_bits_3d(yi) << 1) | expand_bits_3d(zi);
#endif
}

} // namespace morton

// ============================================================================
// 2. PROXIMITY PRIMITIVES
// ============================================================================

/// Closest point on a line segment (a, b) to query point p
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Point closest_point_on_segment (const Point& p, const Point& a, const Point& b) {
    Vec ab, ap;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        ab[d] = b[d] - a[d];
        ap[d] = p[d] - a[d];
    }
    Real dot_ab_ab = Real(0.0), dot_ap_ab = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dot_ab_ab += ab[d] * ab[d];
        dot_ap_ab += ap[d] * ab[d];
    }
    Real t = (dot_ab_ab > Real(1e-30)) ? dot_ap_ab / dot_ab_ab : Real(0.0);
    t = amrex::max(Real(0.0), amrex::min(Real(1.0), t));
    Point cp;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) cp[d] = a[d] + t * ab[d];
    return cp;
}

#if (AMREX_SPACEDIM == 3)
/// Closest point on a triangle (v0, v1, v2) to query point p.
/// Uses the Voronoi-region approach (Ericson, "Real-Time Collision Detection", Ch.5).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Point closest_point_on_triangle (const Point& p,
                                 const Point& v0,
                                 const Point& v1,
                                 const Point& v2)
{
    auto sub = [](const Point& a, const Point& b) -> Vec {
        return Vec{a[0]-b[0], a[1]-b[1], a[2]-b[2]};
    };
    auto dot = [](const Vec& a, const Vec& b) -> Real {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    Vec ab = sub(v1, v0);
    Vec ac = sub(v2, v0);
    Vec ap = sub(p,  v0);
    Real d1 = dot(ab, ap);
    Real d2 = dot(ac, ap);
    if (d1 <= Real(0.0) && d2 <= Real(0.0)) return v0;
    Vec bp = sub(p, v1);
    Real d3 = dot(ab, bp);
    Real d4 = dot(ac, bp);
    if (d3 >= Real(0.0) && d4 <= d3) return v1;
    Real vc = d1*d4 - d3*d2;
    if (vc <= Real(0.0) && d1 >= Real(0.0) && d3 <= Real(0.0)) {
        Real v_val = d1 / (d1 - d3);
        return Point{v0[0]+v_val*ab[0], v0[1]+v_val*ab[1], v0[2]+v_val*ab[2]};
    }
    Vec cp_ = sub(p, v2);
    Real d5 = dot(ab, cp_);
    Real d6 = dot(ac, cp_);
    if (d6 >= Real(0.0) && d5 <= d6) return v2;
    Real vb = d5*d2 - d1*d6;
    if (vb <= Real(0.0) && d2 >= Real(0.0) && d6 <= Real(0.0)) {
        Real w = d2 / (d2 - d6);
        return Point{v0[0]+w*ac[0], v0[1]+w*ac[1], v0[2]+w*ac[2]};
    }
    Real va = d3*d6 - d5*d4;
    if (va <= Real(0.0) && (d4-d3) >= Real(0.0) && (d5-d6) >= Real(0.0)) {
        Real w = (d4-d3) / ((d4-d3) + (d5-d6));
        return Point{v1[0]+w*(v2[0]-v1[0]), v1[1]+w*(v2[1]-v1[1]), v1[2]+w*(v2[2]-v1[2])};
    }
    Real denom = Real(1.0) / (va + vb + vc);
    Real v_val = vb * denom;
    Real w     = vc * denom;
    return Point{v0[0]+ab[0]*v_val+ac[0]*w,
                 v0[1]+ab[1]*v_val+ac[1]*w,
                 v0[2]+ab[2]*v_val+ac[2]*w};
}
#endif // 3D

/// Squared distance between two points
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real point_distance_sq (const Point& a, const Point& b) {
    Real d2 = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real dd = a[d] - b[d];
        d2 += dd * dd;
    }
    return d2;
}

/// Minimum squared distance from point to AABB (for BVH traversal pruning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real point_aabb_distance_sq (const Point& p, const AABB& box) {
    Real d2 = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real v = p[d];
        if (v < box.lo[d]) { Real dd = box.lo[d] - v; d2 += dd*dd; }
        if (v > box.hi[d]) { Real dd = v - box.hi[d]; d2 += dd*dd; }
    }
    return d2;
}

// ============================================================================
// 3. BVHQueryView — GPU-capturable POD for BVH closest-point queries
// ============================================================================

/// \brief Trivially-copyable view into BVH + geometry data for GPU kernels.
struct BVHQueryView {
    const BVHNode* node_arr = nullptr;
    int root = -1;

#if (AMREX_SPACEDIM == 3)
    const Point*               verts     = nullptr;
    const GpuArray<int, 3>*    faces_arr = nullptr;
#elif (AMREX_SPACEDIM == 2)
    const Point* verts   = nullptr;
    int          n_verts = 0;
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query) const
    {
        ClosestPointResult best;
        best.distance = std::numeric_limits<Real>::max();
#if (AMREX_SPACEDIM == 3)
        if (root < 0 || verts == nullptr || faces_arr == nullptr || node_arr == nullptr)
            return best;
#else
        if (root < 0 || verts == nullptr || node_arr == nullptr || n_verts < 2)
            return best;
#endif
        constexpr int MAX_STACK = 128;
        int stack[MAX_STACK];
        int top = 0;
        stack[top++] = root;
        Real best_dist2 = std::numeric_limits<Real>::max();

        while (top > 0) {
            int idx = stack[--top];
            const BVHNode& nd = node_arr[idx];
            Real box_dist2 = point_aabb_distance_sq(query, nd.box);
            if (box_dist2 >= best_dist2) continue;
            if (nd.is_leaf()) {
#if (AMREX_SPACEDIM == 3)
                int fid = nd.prim_id;
                const auto& f = faces_arr[fid];
                Point cp = closest_point_on_triangle(query,
                    verts[f[0]], verts[f[1]], verts[f[2]]);
#else
                int eid = nd.prim_id;
                Point a = verts[eid];
                Point b = verts[(eid + 1) % n_verts];
                Point cp = closest_point_on_segment(query, a, b);
#endif
                Real d2 = point_distance_sq(query, cp);
                if (d2 < best_dist2) {
                    best_dist2    = d2;
                    best.point    = cp;
                    best.prim_id  = nd.prim_id;
                    best.distance = std::sqrt(d2);
                }
            } else {
                // Inline near-first push for both children
                int left  = nd.left;
                int right = nd.right;
                if (left >= 0 && right >= 0) {
                    Real d_left  = point_aabb_distance_sq(query, node_arr[left].box);
                    Real d_right = point_aabb_distance_sq(query, node_arr[right].box);
                    if (d_left < d_right) {
                        if (d_right < best_dist2 && top < MAX_STACK) stack[top++] = right;
                        if (d_left  < best_dist2 && top < MAX_STACK) stack[top++] = left;
                    } else {
                        if (d_left  < best_dist2 && top < MAX_STACK) stack[top++] = left;
                        if (d_right < best_dist2 && top < MAX_STACK) stack[top++] = right;
                    }
                } else {
                    auto push_if = [&](int child) {
                        if (child >= 0 && top < MAX_STACK) {
                            Real d = point_aabb_distance_sq(query, node_arr[child].box);
                            if (d < best_dist2) stack[top++] = child;
                        }
                    };
                    push_if(left);
                    push_if(right);
                }
            }
        }
        return best;
    }
};

// ============================================================================
// 4. BVH CLASS
// ============================================================================

/// \brief Binary BVH tree for closest-point queries.
///
/// Built on CPU with Morton code sorting. Traversed on CPU or GPU.
struct BVH {
    Gpu::ManagedVector<BVHNode> nodes;
    int root = -1;       ///< Index of root node in nodes[]
    int num_prims = 0;   ///< Number of leaf primitives

    bool empty () const { return nodes.empty(); }

    // ----- Build methods (CPU) ------------------------------------------

    /// Build from a flat array of primitive bounding boxes.
    void build (const std::vector<AABB>& prim_boxes) {
        int n = static_cast<int>(prim_boxes.size());
        num_prims = n;
        if (n == 0) { root = -1; return; }
        if (n == 1) {
            nodes.resize(1);
            nodes[0].box     = prim_boxes[0];
            nodes[0].left    = -1;
            nodes[0].right   = -1;
            nodes[0].prim_id = 0;
            root = 0;
            return;
        }
        AABB scene;
        for (auto& b : prim_boxes) scene.merge(b);

        std::vector<uint32_t> codes(n);
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        for (int i = 0; i < n; ++i) {
            codes[i] = morton::morton_code(prim_boxes[i].centroid(), scene);
        }
        std::sort(indices.begin(), indices.end(),
                  [&](int a, int b) { return codes[a] < codes[b]; });

        nodes.resize(2 * n - 1);
        int next_node = 0;
        root = build_recursive(indices.data(), 0, n, prim_boxes, next_node);
        nodes.resize(next_node);
    }

#if (AMREX_SPACEDIM == 3)
    void build (const TriMesh& mesh) {
        int nf = mesh.num_faces();
        std::vector<AABB> boxes(nf);
        for (int i = 0; i < nf; ++i) boxes[i] = mesh.face_aabb(i);
        build(boxes);
    }
#elif (AMREX_SPACEDIM == 2)
    void build (const Polygon2D& poly) {
        int ne = static_cast<int>(poly.size());
        std::vector<AABB> boxes(ne);
        for (int i = 0; i < ne; ++i) boxes[i] = poly.edge_aabb(i);
        build(boxes);
    }
#endif

    // ----- Query methods (CPU/GPU) --------------------------------------

#if (AMREX_SPACEDIM == 3)
    /// Find closest point on a TriMesh using BVH stack traversal.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (
        const Point& query,
        const Point* verts,
        const GpuArray<int,3>* faces_arr,
        const BVHNode* node_arr) const
    {
        ClosestPointResult best;
        best.distance = std::numeric_limits<Real>::max();
        if (root < 0 || verts == nullptr || faces_arr == nullptr || node_arr == nullptr) {
            return best;
        }
        constexpr int MAX_STACK = 128;
        int stack[MAX_STACK];
        int top = 0;
        stack[top++] = root;
        Real best_dist2 = std::numeric_limits<Real>::max();

        while (top > 0) {
            int idx = stack[--top];
            const BVHNode& nd = node_arr[idx];
            Real box_dist2 = point_aabb_distance_sq(query, nd.box);
            if (box_dist2 >= best_dist2) continue;
            if (nd.is_leaf()) {
                int fid = nd.prim_id;
                const auto& f = faces_arr[fid];
                Point cp = closest_point_on_triangle(query,
                    verts[f[0]], verts[f[1]], verts[f[2]]);
                Real d2 = point_distance_sq(query, cp);
                if (d2 < best_dist2) {
                    best_dist2    = d2;
                    best.point    = cp;
                    best.prim_id  = fid;
                    best.distance = std::sqrt(d2);
                }
            } else {
                push_children_near_first<MAX_STACK>(query, node_arr, nd.left, nd.right,
                                                    best_dist2, stack, top);
            }
        }
        return best;
    }
#endif

#if (AMREX_SPACEDIM == 2)
    /// Find closest point on polygon edges using BVH stack traversal.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (
        const Point& query,
        const Point* verts,
        int n_verts,
        const BVHNode* node_arr) const
    {
        ClosestPointResult best;
        best.distance = std::numeric_limits<Real>::max();
        if (root < 0 || verts == nullptr || node_arr == nullptr || n_verts < 2) {
            return best;
        }
        constexpr int MAX_STACK = 128;
        int stack[MAX_STACK];
        int top = 0;
        stack[top++] = root;
        Real best_dist2 = std::numeric_limits<Real>::max();

        while (top > 0) {
            int idx = stack[--top];
            const BVHNode& nd = node_arr[idx];
            Real box_dist2 = point_aabb_distance_sq(query, nd.box);
            if (box_dist2 >= best_dist2) continue;
            if (nd.is_leaf()) {
                int eid = nd.prim_id;
                Point a = verts[eid];
                Point b = verts[(eid + 1) % n_verts];
                Point cp = closest_point_on_segment(query, a, b);
                Real d2 = point_distance_sq(query, cp);
                if (d2 < best_dist2) {
                    best_dist2    = d2;
                    best.point    = cp;
                    best.prim_id  = eid;
                    best.distance = std::sqrt(d2);
                }
            } else {
                push_children_near_first<MAX_STACK>(query, node_arr, nd.left, nd.right,
                                                    best_dist2, stack, top);
            }
        }
        return best;
    }
#endif

    // ---- Convenience overloads that accept geometry objects directly ----

#if (AMREX_SPACEDIM == 3)
    AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query, const TriMesh& mesh) const {
        return closest_point_query(query,
            mesh.vertices.data(),
            mesh.faces.data(),
            nodes.data());
    }
#elif (AMREX_SPACEDIM == 2)
    AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query, const Polygon2D& poly) const {
        return closest_point_query(query,
            poly.verts.data(),
            static_cast<int>(poly.verts.size()),
            nodes.data());
    }
#endif

    // ---- BVHQueryView factory methods -----------------------------------

#if (AMREX_SPACEDIM == 3)
    BVHQueryView query_view (const TriMesh& mesh) const {
        BVHQueryView v;
        v.node_arr  = nodes.data();
        v.root      = root;
        v.verts     = mesh.vertices.data();
        v.faces_arr = mesh.faces.data();
        return v;
    }
#elif (AMREX_SPACEDIM == 2)
    BVHQueryView query_view (const Polygon2D& poly) const {
        BVHQueryView v;
        v.node_arr = nodes.data();
        v.root     = root;
        v.verts    = poly.verts.data();
        v.n_verts  = static_cast<int>(poly.verts.size());
        return v;
    }
#endif

private:
    template <int MAX_STACK>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void push_if_candidate (int child,
                                   const Point& query,
                                   const BVHNode* node_arr,
                                   Real best_dist2,
                                   int* stack,
                                   int& top) {
        if (child < 0) return;
        if (top >= MAX_STACK) return;
        const Real d = point_aabb_distance_sq(query, node_arr[child].box);
        if (d < best_dist2) stack[top++] = child;
    }

    template <int MAX_STACK>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void push_children_near_first (const Point& query,
                                          const BVHNode* node_arr,
                                          int left,
                                          int right,
                                          Real best_dist2,
                                          int* stack,
                                          int& top) {
        if (left < 0 || right < 0) {
            push_if_candidate<MAX_STACK>(left, query, node_arr, best_dist2, stack, top);
            push_if_candidate<MAX_STACK>(right, query, node_arr, best_dist2, stack, top);
            return;
        }

        const Real d_left  = point_aabb_distance_sq(query, node_arr[left].box);
        const Real d_right = point_aabb_distance_sq(query, node_arr[right].box);

        // LIFO stack: push farther child first so the nearer child is popped first.
        if (d_left < d_right) {
            if (d_right < best_dist2 && top < MAX_STACK) stack[top++] = right;
            if (d_left  < best_dist2 && top < MAX_STACK) stack[top++] = left;
        } else {
            if (d_left  < best_dist2 && top < MAX_STACK) stack[top++] = left;
            if (d_right < best_dist2 && top < MAX_STACK) stack[top++] = right;
        }
    }

    int build_recursive (const int* sorted, int lo, int hi,
                         const std::vector<AABB>& prim_boxes, int& next_node) {
        int idx = next_node++;
        if (hi - lo == 1) {
            nodes[idx].box     = prim_boxes[sorted[lo]];
            nodes[idx].left    = -1;
            nodes[idx].right   = -1;
            nodes[idx].prim_id = sorted[lo];
            return idx;
        }
        int mid = (lo + hi) / 2;
        nodes[idx].prim_id = -1;
        nodes[idx].left  = build_recursive(sorted, lo, mid, prim_boxes, next_node);
        nodes[idx].right = build_recursive(sorted, mid, hi, prim_boxes, next_node);
        nodes[idx].box = nodes[nodes[idx].left].box;
        nodes[idx].box.merge(nodes[nodes[idx].right].box);
        return idx;
    }
};

#endif // BVH_H_
