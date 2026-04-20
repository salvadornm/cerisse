#ifndef IBM_BVH_H_
#define IBM_BVH_H_

// ============================================================================
// ibm_bvh.h — BVH algorithmic code
//
// Implements all algorithmic (functional) BVH logic:
//   1. Morton encoding  : namespace morton (morton_code, expand_bits_*)
//   2. Proximity queries: closest_point_on_segment/triangle, distance helpers
//   3. BVH class        : build() (CPU, Morton-sorted median split),
//                         closest_point_query() (CPU/GPU stack traversal)
// ============================================================================

#include "ibm_bvh_defs.h"

#include <numeric>
#include <vector>
#include <algorithm>
#include <cstdint>

#if defined(AMREX_USE_CUDA)
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#endif

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
// 3. LBVH — GPU-parallel BVH construction (Karras 2012)
// ============================================================================

namespace lbvh {

/// Count leading zeros (portable host/device)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz32 (uint32_t x) {
    if (x == 0) return 32;
#if defined(__CUDA_ARCH__)
    return __clz(x);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_clz(x);
#else
    int n = 0;
    if (x <= 0x0000FFFFu) { n += 16; x <<= 16; }
    if (x <= 0x00FFFFFFu) { n +=  8; x <<=  8; }
    if (x <= 0x0FFFFFFFu) { n +=  4; x <<=  4; }
    if (x <= 0x3FFFFFFFu) { n +=  2; x <<=  2; }
    if (x <= 0x7FFFFFFFu) { n +=  1; }
    return n;
#endif
}

/// Longest common prefix length between sorted Morton codes i and j.
/// Returns -1 for out-of-range indices.  Extends with index comparison
/// when codes are identical (Karras 2012, §4).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int delta (const uint32_t* codes, int n, int i, int j) {
    if (j < 0 || j >= n) return -1;
    uint32_t xi = codes[i], xj = codes[j];
    if (xi != xj) return clz32(xi ^ xj);
    return 32 + clz32(static_cast<uint32_t>(i ^ j));
}

} // namespace lbvh

// ============================================================================
// 4. BVH4QueryView — GPU-capturable POD for 4-wide closest-point queries
// ============================================================================

/// \brief Trivially-copyable view into BVH4 + geometry data for GPU kernels.
///
/// Traverses the collapsed 4-wide tree.  Each internal node test evaluates
/// up to 4 child AABBs, sorts by distance, processes leaves inline, and
/// pushes internal children far-to-near (LIFO → near popped first).
struct BVH4QueryView {
    const BVH4Node* node4_arr = nullptr;
    int root4 = -1;

#if (AMREX_SPACEDIM == 3)
    const Point*            verts     = nullptr;
    const GpuArray<int, 3>* faces_arr = nullptr;
#elif (AMREX_SPACEDIM == 2)
    const Point* verts   = nullptr;
    int          n_verts = 0;
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query) const
    {
        ClosestPointResult best;
        best.distance = std::numeric_limits<Real>::max();
        if (root4 < 0 || node4_arr == nullptr || verts == nullptr) return best;

        // BVH4 halves tree depth → smaller stack suffices
        constexpr int MAX_STACK = 64;
        int stack[MAX_STACK];
        int top = 0;
        stack[top++] = root4;
        Real best_dist2 = std::numeric_limits<Real>::max();

        while (top > 0) {
            int idx = stack[--top];
            const BVH4Node& nd = node4_arr[idx];

            // Compute distance to each child AABB
            Real dist[4];
            for (int c = 0; c < nd.n_children; ++c)
                dist[c] = point_aabb_distance_sq(query, nd.child_box[c]);

            // Insertion-sort children by distance (ascending, max 4 elements)
            int order[4] = {0, 1, 2, 3};
            for (int c = 1; c < nd.n_children; ++c) {
                int key = order[c];
                Real key_d = dist[key];
                int j = c - 1;
                while (j >= 0 && dist[order[j]] > key_d) {
                    order[j + 1] = order[j];
                    --j;
                }
                order[j + 1] = key;
            }

            // Pass 1: process leaf children near-to-far (tightens bound early)
            for (int c = 0; c < nd.n_children; ++c) {
                int ci = order[c];
                if (dist[ci] >= best_dist2) break;
                if (!nd.child_is_leaf(ci)) continue;
                int pid = nd.child_prim(ci);
#if (AMREX_SPACEDIM == 3)
                const auto& f = faces_arr[pid];
                Point cp = closest_point_on_triangle(
                    query, verts[f[0]], verts[f[1]], verts[f[2]]);
#else
                Point cp = closest_point_on_segment(
                    query, verts[pid], verts[(pid + 1) % n_verts]);
#endif
                Real d2 = point_distance_sq(query, cp);
                if (d2 < best_dist2) {
                    best_dist2    = d2;
                    best.point    = cp;
                    best.prim_id  = pid;
                    best.distance = std::sqrt(d2);
                }
            }

            // Pass 2: push internal children far-to-near (LIFO → near popped first)
            for (int c = nd.n_children - 1; c >= 0; --c) {
                int ci = order[c];
                if (dist[ci] >= best_dist2) continue;
                if (nd.child_is_leaf(ci)) continue;
                if (top < MAX_STACK) stack[top++] = nd.child_idx[ci];
            }
        }
        return best;
    }
};

/// Legacy binary BVH query view — kept for InsideTester ray-casting
/// which still traverses the binary tree.
struct BVHQueryView {
    const BVHNode* node_arr = nullptr;
    int root = -1;

#if (AMREX_SPACEDIM == 3)
    const Point*            verts     = nullptr;
    const GpuArray<int, 3>* faces_arr = nullptr;
#elif (AMREX_SPACEDIM == 2)
    const Point* verts   = nullptr;
    int          n_verts = 0;
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query) const
    {
        ClosestPointResult best;
        best.distance = std::numeric_limits<Real>::max();
        if (root < 0 || node_arr == nullptr || verts == nullptr) return best;

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
                Point cp = closest_point_on_segment(query,
                    verts[eid], verts[(eid + 1) % n_verts]);
#endif
                Real d2 = point_distance_sq(query, cp);
                if (d2 < best_dist2) {
                    best_dist2    = d2;
                    best.point    = cp;
                    best.prim_id  = nd.prim_id;
                    best.distance = std::sqrt(d2);
                }
            } else {
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
                }
            }
        }
        return best;
    }
};

// ============================================================================
// 5. BVH CLASS — build (CPU recursive + GPU LBVH) + BVH4 collapse
// ============================================================================

/// \brief Binary BVH tree with optional 4-wide collapsed form.
///
/// Build: CPU recursive (Morton-sorted median split) or GPU LBVH (Karras 2012).
/// Query: through BVH4QueryView (collapsed 4-wide) for closest-point,
///        or BVHQueryView (binary) for InsideTester ray-casting.
struct BVH {
    Gpu::ManagedVector<BVHNode>  nodes;       ///< Binary tree nodes
    Gpu::ManagedVector<BVH4Node> nodes4;      ///< Collapsed 4-wide nodes
    int root      = -1;   ///< Root index in nodes[]
    int root4     = -1;   ///< Root index in nodes4[]
    int num_prims = 0;    ///< Number of leaf primitives

    bool empty () const { return nodes.empty(); }

    // ----- Build entry points -------------------------------------------

    /// Build from primitive bounding boxes.
    /// Uses GPU LBVH when available, else CPU recursive.
    /// Always collapses to BVH4 after binary build.
    void build (const std::vector<AABB>& prim_boxes) {
        int n = static_cast<int>(prim_boxes.size());
        num_prims = n;
        if (n == 0) { root = -1; root4 = -1; return; }
        if (n == 1) {
            nodes.resize(1);
            nodes[0].box     = prim_boxes[0];
            nodes[0].left    = -1;
            nodes[0].right   = -1;
            nodes[0].prim_id = 0;
            root = 0;
            collapse_to_bvh4();
            return;
        }

#ifdef AMREX_USE_GPU
        build_lbvh(prim_boxes);
#else
        build_cpu(prim_boxes);
#endif
        collapse_to_bvh4();
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

    // ----- BVH4QueryView factory (primary query path) -------------------

#if (AMREX_SPACEDIM == 3)
    BVH4QueryView query_view (const TriMesh& mesh) const {
        BVH4QueryView v;
        v.node4_arr = nodes4.data();
        v.root4     = root4;
        v.verts     = mesh.vertices.data();
        v.faces_arr = mesh.faces.data();
        return v;
    }
#elif (AMREX_SPACEDIM == 2)
    BVH4QueryView query_view (const Polygon2D& poly) const {
        BVH4QueryView v;
        v.node4_arr = nodes4.data();
        v.root4     = root4;
        v.verts     = poly.verts.data();
        v.n_verts   = static_cast<int>(poly.verts.size());
        return v;
    }
#endif

    // ----- Legacy binary BVHQueryView (for InsideTester) -----------------

#if (AMREX_SPACEDIM == 3)
    BVHQueryView binary_query_view (const TriMesh& mesh) const {
        BVHQueryView v;
        v.node_arr  = nodes.data();
        v.root      = root;
        v.verts     = mesh.vertices.data();
        v.faces_arr = mesh.faces.data();
        return v;
    }
#elif (AMREX_SPACEDIM == 2)
    BVHQueryView binary_query_view (const Polygon2D& poly) const {
        BVHQueryView v;
        v.node_arr = nodes.data();
        v.root     = root;
        v.verts    = poly.verts.data();
        v.n_verts  = static_cast<int>(poly.verts.size());
        return v;
    }
#endif

    // ----- Convenience CPU query (uses BVH4 internally) -----------------

#if (AMREX_SPACEDIM == 3)
    AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query, const TriMesh& mesh) const {
        auto v = query_view(mesh);
        return v.closest_point_query(query);
    }
#elif (AMREX_SPACEDIM == 2)
    AMREX_FORCE_INLINE
    ClosestPointResult closest_point_query (const Point& query, const Polygon2D& poly) const {
        auto v = query_view(poly);
        return v.closest_point_query(query);
    }
#endif

    // ----- CPU recursive build (Morton-sorted median split) -------------

    void build_cpu (const std::vector<AABB>& prim_boxes) {
        int n = static_cast<int>(prim_boxes.size());
        AABB scene;
        for (auto& b : prim_boxes) scene.merge(b);

        std::vector<uint32_t> codes(n);
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);
        for (int i = 0; i < n; ++i)
            codes[i] = morton::morton_code(prim_boxes[i].centroid(), scene);
        std::sort(indices.begin(), indices.end(),
                  [&](int a, int b) { return codes[a] < codes[b]; });

        nodes.resize(2 * n - 1);
        int next_node = 0;
        root = build_recursive(indices.data(), 0, n, prim_boxes, next_node);
        nodes.resize(next_node);
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

    // ----- GPU LBVH build (Karras 2012) ---------------------------------
    //
    // Layout: internal nodes at [0, n-2], leaf nodes at [n-1, 2n-2].
    // Total 2n-1 BVHNode entries (same as CPU build).
    // Root is always internal node 0.

    void build_lbvh (const std::vector<AABB>& prim_boxes) {
        using namespace amrex;
        const int n = static_cast<int>(prim_boxes.size());

        // --- 1. Scene bounding box ---
        AABB scene;
        for (auto& b : prim_boxes) scene.merge(b);

        // --- 2. Upload primitive boxes to device, compute Morton codes ---
        // Use DeviceVector for GPU-only temporaries (no managed memory page faults).
        Gpu::DeviceVector<AABB>     d_boxes(n);
        Gpu::DeviceVector<uint32_t> d_codes(n);
        Gpu::DeviceVector<int>      d_indices(n);

        // Gpu::htod_memcpy(d_boxes.data(), prim_boxes.data(), n * sizeof(AABB));  // Jiaye original line

        Gpu::copy(Gpu::hostToDevice, prim_boxes.begin(), prim_boxes.end(), d_boxes.begin());  // new line

        auto* box_ptr = d_boxes.data();
        auto* cod_ptr = d_codes.data();
        auto* idx_ptr = d_indices.data();
        const AABB scene_box = scene;

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
            cod_ptr[i] = morton::morton_code(box_ptr[i].centroid(), scene_box);
            idx_ptr[i] = i;
        });
        Gpu::streamSynchronize();

        // --- 3. Sort by Morton code ---
#if defined(AMREX_USE_CUDA)
        {
            thrust::device_ptr<uint32_t> t_codes(d_codes.data());
            thrust::device_ptr<int>      t_indices(d_indices.data());
            thrust::sort_by_key(thrust::device, t_codes, t_codes + n, t_indices);
            Gpu::streamSynchronize();
        }
#else
        {
            std::vector<uint32_t> h_codes(n);
            std::vector<int>      h_idx(n);
            //Gpu::dtoh_memcpy(h_codes.data(), d_codes.data(), n * sizeof(uint32_t)); // original line 
            Gpu::copy(Gpu::deviceToHost, d_codes.begin(), d_codes.end(), h_codes.begin()); // new

            //Gpu::dtoh_memcpy(h_idx.data(),   d_indices.data(), n * sizeof(int));  // original
            Gpu::copy(Gpu::deviceToHost,d_indices.begin(), d_indices.end(),h_idx.begin()); // new


            std::vector<int> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&](int a, int b) { return h_codes[a] < h_codes[b]; });
            std::vector<uint32_t> sorted_codes(n);
            std::vector<int>      sorted_idx(n);
            for (int i = 0; i < n; ++i) {
                sorted_codes[i] = h_codes[order[i]];
                sorted_idx[i]   = h_idx[order[i]];
            }
            //Gpu::htod_memcpy(d_codes.data(),   sorted_codes.data(), n * sizeof(uint32_t)); //  original
            //Gpu::htod_memcpy(d_indices.data(), sorted_idx.data(),   n * sizeof(int)); // original
            Gpu::copy(Gpu::hostToDevice, sorted_codes.begin(), sorted_codes.end(), d_codes.begin());
            Gpu::copy(Gpu::hostToDevice,sorted_idx.begin(), sorted_idx.end(),d_indices.begin());
        }
#endif

        // --- 4. Allocate tree nodes (managed — needed on host for BVH4 collapse) ---
        const int n_total = 2 * n - 1;
        nodes.resize(n_total);

        // Auxiliary arrays: device-only (freed automatically after build)
        Gpu::DeviceVector<int> d_parent(n_total, -1);
        Gpu::DeviceVector<int> d_counter(n - 1, 0);

        auto* nd_ptr     = nodes.data();
        auto* parent_ptr = d_parent.data();
        const auto* code_ptr = d_codes.data();
        const auto* sidx_ptr = d_indices.data();

        // --- 5. Initialize leaf nodes [n-1 .. 2n-2] ---
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
            const int leaf = (n - 1) + i;
            nd_ptr[leaf].box     = box_ptr[sidx_ptr[i]];
            nd_ptr[leaf].left    = -1;
            nd_ptr[leaf].right   = -1;
            nd_ptr[leaf].prim_id = sidx_ptr[i];
        });

        // --- 6. Karras internal node construction ---
        ParallelFor(n - 1, [=] AMREX_GPU_DEVICE (int i) noexcept {
            int d_fwd = lbvh::delta(code_ptr, n, i, i + 1);
            int d_bwd = lbvh::delta(code_ptr, n, i, i - 1);
            int d = (d_fwd > d_bwd) ? 1 : -1;

            int delta_min = lbvh::delta(code_ptr, n, i, i - d);
            int l_max = 2;
            while (lbvh::delta(code_ptr, n, i, i + l_max * d) > delta_min)
                l_max *= 2;

            int l = 0;
            for (int t = l_max / 2; t >= 1; t /= 2) {
                if (lbvh::delta(code_ptr, n, i, i + (l + t) * d) > delta_min)
                    l += t;
            }
            int j = i + l * d;

            int delta_node = lbvh::delta(code_ptr, n, i, j);
            int s = 0;
            int range_len = (i < j ? j : i) - (i < j ? i : j);
            for (int t = (range_len + 1) / 2; t >= 1; t = (t == 1) ? 0 : (t + 1) / 2) {
                if (lbvh::delta(code_ptr, n, i, i + (s + t) * d) > delta_node)
                    s += t;
                if (t == 1) break;
            }
            int gamma = i + s * d + amrex::min(d, 0);

            int left_child  = (amrex::min(i, j) == gamma)     ? (n - 1) + gamma       : gamma;
            int right_child = (amrex::max(i, j) == gamma + 1) ? (n - 1) + (gamma + 1) : gamma + 1;

            nd_ptr[i].left    = left_child;
            nd_ptr[i].right   = right_child;
            nd_ptr[i].prim_id = -1;
            parent_ptr[left_child]  = i;
            parent_ptr[right_child] = i;
        });
        Gpu::streamSynchronize();

        // --- 7. Bottom-up AABB propagation ---
        auto* cnt_ptr = d_counter.data();
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
            int current = parent_ptr[(n - 1) + i];
            while (current >= 0) {
                int old = Gpu::Atomic::Add(&cnt_ptr[current], 1);
                if (old == 0) return;
                nd_ptr[current].box = AABB();
                nd_ptr[current].box.merge(nd_ptr[nd_ptr[current].left].box);
                nd_ptr[current].box.merge(nd_ptr[nd_ptr[current].right].box);
                current = parent_ptr[current];
            }
        });
        Gpu::streamSynchronize();

        root = 0;
    }

    // ----- BVH4 collapse (CPU, runs once after build) -------------------
    //
    // Converts the binary tree into a 4-wide tree by merging pairs of
    // levels.  Each BVH4 node gathers up to 4 children (the grandchildren
    // of a binary node).  This halves tree depth and reduces stack usage.

    void collapse_to_bvh4 () {
        if (nodes.empty()) { root4 = -1; return; }

        // Copy binary nodes to host vector — avoids managed-memory page
        // faults during the recursive traversal on CPU.
        std::vector<BVHNode> h_nodes(nodes.begin(), nodes.end());

        std::vector<BVH4Node> tmp;
        tmp.reserve(num_prims);
        root4 = collapse_node(root, tmp, h_nodes);

        nodes4.resize(tmp.size());
        std::copy(tmp.begin(), tmp.end(), nodes4.begin());
    }

    /// Recursively collapse a binary subtree into BVH4 nodes.
    /// Returns a BVH4 node index, or (LEAF_FLAG | prim_id) for a leaf.
    int collapse_node (int bi, std::vector<BVH4Node>& out,
                       const std::vector<BVHNode>& h_nodes) {
        const auto& nd = h_nodes[bi];
        if (nd.is_leaf())
            return BVH4Node::LEAF_FLAG | nd.prim_id;

        // Gather up to 4 children by expanding internal children one level
        struct Child { int bin_idx; AABB box; };
        std::vector<Child> children;
        children.reserve(4);

        auto expand = [&](int child_bi) {
            const auto& c = h_nodes[child_bi];
            if (c.is_leaf() || static_cast<int>(children.size()) >= 3) {
                children.push_back({child_bi, c.box});
            } else {
                children.push_back({c.left,  h_nodes[c.left].box});
                children.push_back({c.right, h_nodes[c.right].box});
            }
        };

        expand(nd.left);
        expand(nd.right);

        // Allocate a BVH4 node
        int bvh4_idx = static_cast<int>(out.size());
        out.push_back(BVH4Node{});
        BVH4Node& n4 = out.back();
        n4.n_children = static_cast<int>(children.size());

        for (int c = 0; c < n4.n_children; ++c) {
            n4.child_box[c] = children[c].box;
            n4.child_idx[c] = collapse_node(children[c].bin_idx, out, h_nodes);
        }
        // Zero-fill unused slots
        for (int c = n4.n_children; c < 4; ++c) {
            n4.child_box[c] = AABB();
            n4.child_idx[c] = -1;
        }
        return bvh4_idx;
    }
};

#endif // IBM_BVH_H_
