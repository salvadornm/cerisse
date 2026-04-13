#ifndef IBM_BACKEND_BVH_H_
#define IBM_BACKEND_BVH_H_

// ============================================================================
// ibm_backend_bvh.h — BVH backend implementation
//
// Implements the interface between the IBM solver and BVH geometry:
//   1. InsideTesterView / InsideTester : GPU-compatible inside/outside tests
//   2. Geometry I/O                   : read_stl, read_off, read_polygon_2d
//   3. Geometry processing            : build_geometry_cache, convert_inout,
//                                       check_ibm_geometry_consistency
//   4. Utility functions              : make_vec, bbox_contains, IB_WarnOnBoundary
// ============================================================================

#include "ibm_bvh.h"

#include <AMReX_Print.H>

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <unordered_map>
#include <utility>
#include <cctype>

// ============================================================================
// 1. InsideTester — GPU-compatible inside/outside testing
//
// Two-layer design:
//   InsideTesterView — lightweight POD (raw pointers + ints), trivially copyable,
//                      safe to capture by value in GPU ParallelFor lambdas.
//   InsideTester     — owns managed-memory geometry; provides view() for GPU.
// ============================================================================

#if (AMREX_SPACEDIM == 3)

/// Lightweight, trivially-copyable device view for 3D inside/outside tests.
///
/// Uses a 1-ray fast path with lazy fallback to additional directions only
/// when the primary ray hits a degenerate configuration (near-parallel face,
/// vertex/edge grazing).  Most points resolve with a single BVH traversal;
/// pathological cases fall back to 2-ray then 3-ray majority voting.
struct InsideTesterView {
    const Point*           verts_     = nullptr;
    const GpuArray<int,3>* faces_     = nullptr;
    int                    n_faces_   = 0;
    const BVHNode*         bvh_nodes_ = nullptr;
    int                    bvh_root_  = -1;

    AMREX_GPU_HOST_DEVICE
    BoundedSide operator() (const Point& p) const {
        if (bvh_root_ < 0 || verts_ == nullptr || faces_ == nullptr || bvh_nodes_ == nullptr || n_faces_ <= 0) {
            return BoundedSide::Outside;
        }

        // Three pre-set ray directions (linearly independent, irrational oblique
        // third direction minimises alignment with mesh edges/faces).
        const Vec dirs[3] = {
            {Real(1.0), Real(0.0), Real(0.0)},
            {Real(0.0), Real(1.0), Real(0.0)},
            {Real(0.30151134457776363), Real(0.30151134457776363), Real(0.90453403373329089)}
        };

        // Fast path: single ray.  If no degenerate hit is detected the
        // result is authoritative and we return immediately (1 traversal).
        bool  degenerate = false;
        int   crossings  = ray_cast_count(p, dirs[0], degenerate);
        if (!degenerate) {
            return (crossings % 2 == 1) ? BoundedSide::Inside
                                        : BoundedSide::Outside;
        }

        // Fallback: the first ray hit a degenerate face.  Cast a second ray
        // in an independent direction.  If it is clean we trust it.
        bool  degen2 = false;
        int   cross2 = ray_cast_count(p, dirs[1], degen2);
        if (!degen2) {
            return (cross2 % 2 == 1) ? BoundedSide::Inside
                                     : BoundedSide::Outside;
        }

        // Both axis-aligned rays were degenerate — use the oblique third
        // direction and take its answer (best-effort).
        bool  degen3 = false;
        int   cross3 = ray_cast_count(p, dirs[2], degen3);
        return (cross3 % 2 == 1) ? BoundedSide::Inside
                                 : BoundedSide::Outside;
    }

private:
    /// Ray-cast with degeneracy detection.
    /// Sets `degenerate = true` if any triangle is near-parallel to the ray
    /// (|det| < degen_eps) so the caller can retry with a different direction.
    AMREX_GPU_HOST_DEVICE
    int ray_cast_count (const Point& p, const Vec& dir, bool& degenerate) const {
        degenerate = false;
        if (bvh_root_ < 0 || verts_ == nullptr || faces_ == nullptr || bvh_nodes_ == nullptr || n_faces_ <= 0) {
            return 0;
        }
        int crossings = 0;
        const Real eps = IBM_EPS::RAYCAST;
        const Real degen_eps = IBM_EPS::RAYCAST * Real(100.0);   // near-parallel threshold
        constexpr int MAX_STACK = 128;
        int stack[MAX_STACK];
        int top = 0;
        stack[top++] = bvh_root_;
        while (top > 0) {
            int idx = stack[--top];
            const BVHNode& nd = bvh_nodes_[idx];
            if (!ray_aabb_intersect(p, dir, nd.box)) continue;
            if (nd.is_leaf()) {
                int fid = nd.prim_id;
                const auto& f = faces_[fid];
                const Point& v0 = verts_[f[0]];
                const Point& v1 = verts_[f[1]];
                const Point& v2 = verts_[f[2]];
                Vec e1 = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
                Vec e2 = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};
                Vec h = {dir[1]*e2[2]-dir[2]*e2[1],
                         dir[2]*e2[0]-dir[0]*e2[2],
                         dir[0]*e2[1]-dir[1]*e2[0]};
                Real a_val = e1[0]*h[0] + e1[1]*h[1] + e1[2]*h[2];

                // Near-parallel: flag degenerate so caller can try another dir
                if (a_val > -degen_eps && a_val < degen_eps) {
                    degenerate = true;
                    continue;
                }

                Real inv_a = Real(1.0) / a_val;
                Vec s = {p[0]-v0[0], p[1]-v0[1], p[2]-v0[2]};
                Real u = inv_a * (s[0]*h[0] + s[1]*h[1] + s[2]*h[2]);
                if (u < -eps || u > Real(1.0)+eps) continue;
                Vec q = {s[1]*e1[2]-s[2]*e1[1],
                         s[2]*e1[0]-s[0]*e1[2],
                         s[0]*e1[1]-s[1]*e1[0]};
                Real v_val = inv_a * (dir[0]*q[0] + dir[1]*q[1] + dir[2]*q[2]);
                if (v_val < -eps || u + v_val > Real(1.0)+eps) continue;
                Real t = inv_a * (e2[0]*q[0] + e2[1]*q[1] + e2[2]*q[2]);
                if (t > eps) crossings++;
            } else {
                if (top < MAX_STACK) stack[top++] = nd.left;
                if (top < MAX_STACK) stack[top++] = nd.right;
            }
        }
        return crossings;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static bool ray_aabb_intersect (const Point& p, const Vec& dir, const AABB& box) {
        Real tmin = Real(0.0);
        Real tmax = std::numeric_limits<Real>::max();
        for (int d = 0; d < 3; ++d) {
            if (std::abs(dir[d]) < IBM_EPS::GEOM) {
                if (p[d] < box.lo[d] || p[d] > box.hi[d]) return false;
            } else {
                Real inv_d = Real(1.0) / dir[d];
                Real t1 = (box.lo[d] - p[d]) * inv_d;
                Real t2 = (box.hi[d] - p[d]) * inv_d;
                if (t1 > t2) { Real tmp = t1; t1 = t2; t2 = tmp; }
                tmin = amrex::max(tmin, t1);
                tmax = amrex::min(tmax, t2);
                if (tmin > tmax) return false;
            }
        }
        return true;
    }
};

/// Inside/outside tester that owns managed-memory copies of geometry data.
struct InsideTester {
    Gpu::ManagedVector<Point>           managed_verts_;
    Gpu::ManagedVector<GpuArray<int,3>> managed_faces_;
    InsideTesterView view_;

    InsideTester () = default;
    InsideTester (const TriMesh& mesh, const BVH& bvh) { setup(mesh, bvh); }

    void setup (const TriMesh& mesh, const BVH& bvh) {
        int nv = mesh.num_vertices();
        int nf = mesh.num_faces();
        managed_verts_.resize(nv);
        managed_faces_.resize(nf);
        for (int i = 0; i < nv; ++i) managed_verts_[i] = mesh.vertices[i];
        for (int i = 0; i < nf; ++i) managed_faces_[i] = mesh.faces[i];
        view_.verts_     = managed_verts_.data();
        view_.faces_     = managed_faces_.data();
        view_.n_faces_   = nf;
        view_.bvh_nodes_ = bvh.nodes.data();
        view_.bvh_root_  = bvh.root;
    }

    BoundedSide operator() (const Point& p) const { return view_(p); }
    InsideTesterView view () const { return view_; }
};

using inside_t = InsideTester;

#elif (AMREX_SPACEDIM == 2)

/// Lightweight, trivially-copyable device view for 2D inside/outside tests.
struct InsideTesterView {
    const Point* verts_   = nullptr;
    int          n_verts_ = 0;
    const BVHNode* bvh_nodes_ = nullptr;
    int            bvh_root_  = -1;

    AMREX_GPU_HOST_DEVICE
    BoundedSide operator() (const Point& p) const {
        if (verts_ == nullptr || n_verts_ <= 0) {
            return BoundedSide::Outside;
        }

        // Fast path: use BVH if available — single +x ray cast.
        if (bvh_nodes_ != nullptr && bvh_root_ >= 0) {
            constexpr int MAX_STACK = 128;

            int crossings = 0;
            int stack[MAX_STACK];
            int top = 0;
            stack[top++] = bvh_root_;

            while (top > 0) {
                int idx = stack[--top];
                const BVHNode& nd = bvh_nodes_[idx];

                if (!ray_aabb_intersect_x(p, nd.box)) continue;

                if (nd.is_leaf()) {
                    int i = nd.prim_id;
                    int j = (i + 1) % n_verts_;
                    const Real yi = verts_[i][1], yj = verts_[j][1];
                    const Real xi = verts_[i][0], xj = verts_[j][0];

                    if ((yi <= p[1] && yj > p[1]) || (yj <= p[1] && yi > p[1])) {
                        Real t = (p[1] - yi) / (yj - yi);
                        Real x_cross = xi + t * (xj - xi);
                        if (p[0] < x_cross) crossings++;
                    }
                } else {
                    if (top < MAX_STACK) stack[top++] = nd.left;
                    if (top < MAX_STACK) stack[top++] = nd.right;
                }
            }

            return (crossings % 2 == 1) ? BoundedSide::Inside : BoundedSide::Outside;
        }

        // Fallback path: linear edge scan (no BVH).
        int crossings = 0;
        for (int i = 0; i < n_verts_; ++i) {
            int j = (i + 1) % n_verts_;
            const Real yi = verts_[i][1], yj = verts_[j][1];
            const Real xi = verts_[i][0], xj = verts_[j][0];
            if ((yi <= p[1] && yj > p[1]) || (yj <= p[1] && yi > p[1])) {
                Real t = (p[1] - yi) / (yj - yi);
                Real x_cross = xi + t * (xj - xi);
                if (p[0] < x_cross) crossings++;
            }
        }
        return (crossings % 2 == 1) ? BoundedSide::Inside : BoundedSide::Outside;
    }

private:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static bool ray_aabb_intersect_x (const Point& p, const AABB& box) {
        // +x ray from p intersects box if y overlaps and box has some x >= p.x
        return (p[1] >= box.lo[1] && p[1] <= box.hi[1] && box.hi[0] >= p[0]);
    }
};

/// Inside/outside tester that owns managed-memory copies of polygon vertices.
struct InsideTester {
    Gpu::ManagedVector<Point> managed_verts_;
    InsideTesterView view_;

    InsideTester () = default;
    explicit InsideTester (const Polygon2D& p) { setup(p); }
    InsideTester (const Polygon2D& p, const BVH& bvh) { setup(p, bvh); }

    void setup (const Polygon2D& p) {
        int n = static_cast<int>(p.verts.size());
        managed_verts_.resize(n);
        for (int i = 0; i < n; ++i) managed_verts_[i] = p.verts[i];
        view_.verts_   = managed_verts_.data();
        view_.n_verts_ = n;
        view_.bvh_nodes_ = nullptr;
        view_.bvh_root_  = -1;
    }

    void setup (const Polygon2D& p, const BVH& bvh) {
        setup(p);
        view_.bvh_nodes_ = bvh.nodes.data();
        view_.bvh_root_  = bvh.root;
    }

    BoundedSide operator() (const Point& p) const { return view_(p); }
    InsideTesterView view () const { return view_; }
};

using inside_t = InsideTester;

#endif // AMREX_SPACEDIM

// ============================================================================
// 2. GEOMETRY I/O
// ============================================================================

#if (AMREX_SPACEDIM == 3)

// Forward declarations
inline bool read_stl_ascii  (const std::string& filename, TriMesh& mesh);
inline bool read_stl_binary (const std::string& filename, TriMesh& mesh);

/// Vertex deduplication tolerance: coordinates are rounded to 1/STL_INV_TOL
/// (i.e. 1e-7 length units) before hashing to merge coincident STL vertices.
static constexpr Real STL_INV_TOL = Real(1e7);

/// Vertex key for spatial deduplication
struct VertexKey {
    int64_t x, y, z;
    bool operator== (const VertexKey& o) const { return x == o.x && y == o.y && z == o.z; }
};
struct VertexKeyHasher {
    std::size_t operator()(const VertexKey& k) const {
        std::size_t h = 0;
        h ^= std::hash<int64_t>{}(k.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int64_t>{}(k.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int64_t>{}(k.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

/// Read a binary or ASCII STL file into a TriMesh. Returns true on success.
inline bool read_stl (const std::string& filename, TriMesh& mesh) {
    mesh.vertices.clear();
    mesh.faces.clear();
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.good()) {
        amrex::Print() << "read_stl: Cannot open file: " << filename << "\n";
        return false;
    }
    char header[80];
    ifs.read(header, 80);
    if (!ifs.good()) {
        amrex::Print() << "read_stl: File too short: " << filename << "\n";
        return false;
    }
    bool is_ascii = false;
    {
        ifs.seekg(0);
        std::string first_line;
        std::getline(ifs, first_line);
        ifs.seekg(0);
        if (first_line.find("solid") != std::string::npos) {
            std::string second_line;
            std::getline(ifs, second_line);
            while (std::getline(ifs, second_line)) {
                auto fp = second_line.find_first_not_of(" \t\r\n");
                if (fp != std::string::npos) {
                    if (second_line.find("facet") != std::string::npos ||
                        second_line.find("endsolid") != std::string::npos)
                        is_ascii = true;
                    break;
                }
            }
            ifs.seekg(0);
        }
    }
    return is_ascii ? read_stl_ascii(filename, mesh) : read_stl_binary(filename, mesh);
}

inline bool read_stl_ascii (const std::string& filename, TriMesh& mesh) {
    std::ifstream ifs(filename);
    if (!ifs.good()) return false;
    std::vector<Point> verts;
    std::vector<GpuArray<int,3>> fcs;
    std::unordered_map<VertexKey, int, VertexKeyHasher> vertex_map;
    auto get_or_add = [&](Real x, Real y, Real z) -> int {
        VertexKey key{static_cast<int64_t>(std::round(x * STL_INV_TOL)),
                      static_cast<int64_t>(std::round(y * STL_INV_TOL)),
                      static_cast<int64_t>(std::round(z * STL_INV_TOL))};
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) return it->second;
        int idx = static_cast<int>(verts.size());
        verts.push_back(Point{x, y, z});
        vertex_map[key] = idx;
        return idx;
    };
    std::string line;
    GpuArray<int,3> tri;
    int vi = 0;
    while (std::getline(ifs, line)) {
        auto pos = line.find_first_not_of(" \t");
        if (pos == std::string::npos) continue;
        std::string trimmed = line.substr(pos);
        if (trimmed.compare(0, 6, "vertex") == 0) {
            std::istringstream iss(trimmed.substr(6));
            Real x, y, z;
            if (!(iss >> x >> y >> z)) continue;
            tri[vi] = get_or_add(x, y, z);
            if (++vi == 3) { fcs.push_back(tri); vi = 0; }
        }
    }
    mesh.assign(verts, fcs);
    amrex::Print() << "read_stl (ASCII): " << mesh.num_vertices() << " vertices, "
                   << mesh.num_faces() << " faces from " << filename << "\n";
    return !mesh.faces.empty();
}

inline bool read_stl_binary (const std::string& filename, TriMesh& mesh) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.good()) return false;
    ifs.seekg(80);
    uint32_t num_triangles = 0;
    ifs.read(reinterpret_cast<char*>(&num_triangles), 4);
    if (!ifs.good() || num_triangles == 0) {
        amrex::Print() << "read_stl_binary: Invalid triangle count in " << filename << "\n";
        return false;
    }
    std::vector<Point> verts;
    std::vector<GpuArray<int,3>> fcs;
    std::unordered_map<VertexKey, int, VertexKeyHasher> vertex_map;
    vertex_map.reserve(static_cast<std::size_t>(num_triangles) * 2u);
    verts.reserve(static_cast<std::size_t>(num_triangles) * 3u);
    fcs.reserve(num_triangles);
    auto get_or_add = [&](float x, float y, float z) -> int {
        Real rx = static_cast<Real>(x), ry = static_cast<Real>(y), rz = static_cast<Real>(z);
        VertexKey key{static_cast<int64_t>(std::round(rx * STL_INV_TOL)),
                      static_cast<int64_t>(std::round(ry * STL_INV_TOL)),
                      static_cast<int64_t>(std::round(rz * STL_INV_TOL))};
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) return it->second;
        int idx = static_cast<int>(verts.size());
        verts.push_back(Point{rx, ry, rz});
        vertex_map[key] = idx;
        return idx;
    };
    for (uint32_t t = 0; t < num_triangles; ++t) {
        float buf[12];
        ifs.read(reinterpret_cast<char*>(buf), 48);
        uint16_t attr;
        ifs.read(reinterpret_cast<char*>(&attr), 2);
        if (!ifs.good()) {
            amrex::Print() << "read_stl_binary: Premature EOF at triangle " << t
                           << " in " << filename << "\n";
            return false;
        }
        GpuArray<int,3> tri;
        tri[0] = get_or_add(buf[3], buf[4],  buf[5]);
        tri[1] = get_or_add(buf[6], buf[7],  buf[8]);
        tri[2] = get_or_add(buf[9], buf[10], buf[11]);
        fcs.push_back(tri);
    }
    mesh.assign(verts, fcs);
    amrex::Print() << "read_stl (binary): " << mesh.num_vertices() << " vertices, "
                   << mesh.num_faces() << " faces from " << filename << "\n";
    return !mesh.faces.empty();
}

inline bool read_off (const std::string& filename, TriMesh& mesh) {
    std::ifstream ifs(filename);
    if (!ifs.good()) {
        amrex::Print() << "read_off: Cannot open file: " << filename << "\n";
        return false;
    }
    std::string line;
    std::getline(ifs, line);
    if (line.find("OFF") == std::string::npos) {
        amrex::Print() << "read_off: Missing OFF header in " << filename << "\n";
        return false;
    }
    int nv = 0, nf = 0, ne = 0;
    while (std::getline(ifs, line)) {
        auto fp = line.find_first_not_of(" \t\r\n");
        if (fp == std::string::npos || line[fp] == '#') continue;
        std::istringstream iss(line);
        if (iss >> nv >> nf >> ne) break;
    }
    std::vector<Point> verts;
    std::vector<GpuArray<int,3>> fcs;
    verts.reserve(nv);
    for (int i = 0; i < nv; ++i) {
        Real x, y, z;
        ifs >> x >> y >> z;
        verts.push_back(Point{x, y, z});
    }
    fcs.reserve(nf);
    for (int i = 0; i < nf; ++i) {
        int k;
        ifs >> k;
        if (k < 3) continue;
        std::vector<int> v(k);
        for (int j = 0; j < k; ++j) ifs >> v[j];
        bool valid = true;
        for (int j = 0; j < k; ++j) {
            if (v[j] < 0 || v[j] >= nv) { valid = false; break; }
        }
        if (!valid) {
            amrex::Print() << "read_off: Invalid vertex index in face " << i
                           << " of " << filename << "\n";
            continue;
        }
        for (int j = 1; j < k - 1; ++j) {
            fcs.push_back(GpuArray<int,3>{v[0], v[j], v[j+1]});
        }
    }
    mesh.assign(verts, fcs);
    amrex::Print() << "read_off: " << mesh.num_vertices() << " vertices, "
                   << mesh.num_faces() << " faces from " << filename << "\n";
    return !mesh.faces.empty();
}

/// Auto-detect file format, read 3D mesh, ensure outward face orientation.
inline bool read_mesh (const std::string& filename, TriMesh& mesh) {
    auto dot = filename.rfind('.');
    bool ok = false;
    if (dot != std::string::npos) {
        std::string ext = filename.substr(dot);
        for (auto& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (ext == ".stl")       ok = read_stl(filename, mesh);
        else if (ext == ".off")  ok = read_off(filename, mesh);
    }
    if (!ok) {
        amrex::Print() << "read_mesh: Unknown extension, trying STL for " << filename << "\n";
        ok = read_stl(filename, mesh);
    }
    if (ok && mesh.is_closed()) mesh.ensure_outward_orientation(true);
    return ok;
}

#endif // AMREX_SPACEDIM == 3

#if (AMREX_SPACEDIM == 2)

/// Read a 2D polygon from a text file (one vertex per line: "x y" or "x, y").
inline bool read_polygon_2d (const std::string& filename, Polygon2D& poly, Real dx = -1.0) {
    poly.clear();
    std::ifstream in(filename);
    if (!in) {
        amrex::Print() << "read_polygon_2d: Cannot open file: " << filename << "\n";
        return false;
    }
    std::string line;
    int line_number = 0;
    std::vector<Point> raw_points;
    while (std::getline(in, line)) {
        ++line_number;
        auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        auto last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);
        if (line[0] == '#' || line.compare(0, 2, "//") == 0) continue;
        for (char& c : line) {
            if (c == ',' || c == ';' || c == '\t') c = ' ';
        }
        std::istringstream iss(line);
        Real x, y;
        if (!(iss >> x >> y)) {
            amrex::Print() << "Warning: Skipping invalid line " << line_number
                           << " in " << filename << "\n";
            continue;
        }
        if (!std::isfinite(x) || !std::isfinite(y)) {
            amrex::Print() << "read_polygon_2d: Invalid coordinate on line "
                           << line_number << " in " << filename << "\n";
            return false;
        }
        raw_points.push_back(Point{x, y});
    }
    if (raw_points.empty()) {
        amrex::Print() << "read_polygon_2d: No valid vertices in " << filename << "\n";
        return false;
    }
    Real min_x = raw_points[0][0], max_x = min_x;
    Real min_y = raw_points[0][1], max_y = min_y;
    for (auto& p : raw_points) {
        min_x = std::min(min_x, p[0]); max_x = std::max(max_x, p[0]);
        min_y = std::min(min_y, p[1]); max_y = std::max(max_y, p[1]);
    }
    Real bbox_diag = std::sqrt((max_x-min_x)*(max_x-min_x) + (max_y-min_y)*(max_y-min_y));
    Real scale_ref = (dx > 0) ? dx : (bbox_diag > IBM_EPS::RAYCAST ? bbox_diag : 1.0);
    Real eps_dedup = IBM_EPS::DEDUP * scale_ref;
    Real area_eps  = IBM_EPS::RAYCAST * scale_ref * scale_ref;
    std::vector<Point> clean;
    clean.reserve(raw_points.size());
    clean.push_back(raw_points[0]);
    int n_dups = 0;
    for (size_t i = 1; i < raw_points.size(); ++i) {
        if (point_distance_sq(clean.back(), raw_points[i]) > eps_dedup * eps_dedup) {
            clean.push_back(raw_points[i]);
        } else {
            n_dups++;
        }
    }
    if (clean.size() > 1) {
        if (point_distance_sq(clean.back(), clean.front()) <= eps_dedup*eps_dedup) {
            clean.pop_back(); n_dups++;
        }
    }
    if (n_dups > 0) {
        amrex::Print() << "Info: Removed " << n_dups << " duplicate vertices (tol="
                       << eps_dedup << ") in " << filename << "\n";
    }
    // Bulk-assign to Polygon2D (uses Gpu::ManagedVector — avoid per-element push_back with CUDA)
    std::vector<Point> staging(clean.begin(), clean.end());
    poly.assign(staging);
    if (poly.size() < 3) {
        amrex::Print() << "read_polygon_2d: Need at least 3 vertices (found "
                       << poly.size() << ") in " << filename << "\n";
        return false;
    }
    if (std::abs(poly.area()) <= area_eps) {
        amrex::Print() << "read_polygon_2d: Near-zero area in " << filename << "\n";
        return false;
    }
    if (!poly.is_simple()) {
        amrex::Print() << "read_polygon_2d: Self-intersecting polygon in " << filename << "\n";
        return false;
    }
    std::size_t original_n = poly.size();
    if (dx > Real(0)) {
        std::vector<Point> ref_pts;
        int n = static_cast<int>(poly.size());
        ref_pts.reserve(n * 2); // rough estimate
        for (int i = 0; i < n; ++i) {
            const Point& p0 = poly.vertex(i);
            const Point& p1 = poly.vertex((i + 1) % n);
            ref_pts.push_back(p0);
            Real vx = p1[0] - p0[0], vy = p1[1] - p0[1];
            Real len = std::sqrt(vx*vx + vy*vy);
            if (len > dx) {
                int nseg = static_cast<int>(std::ceil(len / dx));
                Real inv = Real(1.0) / Real(nseg);
                for (int k = 1; k < nseg; ++k) {
                    Real t = inv * Real(k);
                    ref_pts.push_back(Point{p0[0]+t*vx, p0[1]+t*vy});
                }
            }
        }
        poly.assign(ref_pts);
    }
    if (poly.is_clockwise_oriented()) {
        poly.reverse_orientation();
        amrex::Print() << "Info: Reversed polygon orientation to CCW in " << filename << "\n";
    }
    if (dx > Real(0)) {
        amrex::Print() << "Loaded and refined polygon: " << original_n << " -> "
                       << poly.size() << " vertices from " << filename << "\n";
    } else {
        amrex::Print() << "Loaded polygon with " << poly.size()
                       << " vertices from " << filename << "\n";
    }
    return true;
}

#endif // AMREX_SPACEDIM == 2

// ============================================================================
// 3. GEOMETRY PROCESSING
// ============================================================================

/// Build surface element data and local frames from geometry.
/// Builds into std::vector first, then bulk-appends to Gpu::ManagedVector
/// to avoid per-element push_back overhead with CUDA managed memory.
inline void build_geometry_cache (
    const GeomType& geom,
    Gpu::ManagedVector<SurfElem>& surfelem,
    Gpu::ManagedVector<LocalFrame>& localframe,
    int offset = 0,
    int geomIdx = -1)
{
    amrex::ignore_unused(offset);
    std::vector<SurfElem> se_tmp;
    std::vector<LocalFrame> lf_tmp;
#if (AMREX_SPACEDIM == 2)
    int n = static_cast<int>(geom.size());
    se_tmp.reserve(n);
    lf_tmp.reserve(n);
    for (int i = 0; i < n; ++i) {
        Point a = geom.vertex(i);
        Point b = geom.vertex((i + 1) % n);
        Real c[2] = { Real(0.5)*(a[0]+b[0]), Real(0.5)*(a[1]+b[1]) };
        Real dx_ = b[0]-a[0], dy_ = b[1]-a[1];
        Real len = std::sqrt(dx_*dx_ + dy_*dy_);
        se_tmp.push_back(SurfElem(c, len, geomIdx));
        Real inv_len = (len > IBM_EPS::GEOM) ? Real(1.0)/len : Real(0.0);
        Real tx = dx_ * inv_len, ty = dy_ * inv_len;
        Real n_arr[2]  = { ty, -tx };
        Real t1_arr[2] = { tx, ty };
        lf_tmp.push_back(LocalFrame(n_arr, t1_arr));
    }
#elif (AMREX_SPACEDIM == 3)
    int nf = geom.num_faces();
    se_tmp.reserve(nf);
    lf_tmp.reserve(nf);
    for (int fi = 0; fi < nf; ++fi) {
        const auto& f = geom.faces[fi];
        const Point& p0 = geom.vertices[f[0]];
        const Point& p1 = geom.vertices[f[1]];
        const Point& p2 = geom.vertices[f[2]];
        Real c[3];
        for (int d = 0; d < 3; ++d) c[d] = (p0[d]+p1[d]+p2[d]) / Real(3.0);
        Vec v1 = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        Vec v2 = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        Vec n_raw = {v1[1]*v2[2]-v1[2]*v2[1],
                     v1[2]*v2[0]-v1[0]*v2[2],
                     v1[0]*v2[1]-v1[1]*v2[0]};
        Real n_len2 = n_raw[0]*n_raw[0] + n_raw[1]*n_raw[1] + n_raw[2]*n_raw[2];
        Real area = Real(0.5) * std::sqrt(n_len2);
        if (area < IBM_EPS::GEOM) {
            // Degenerate (zero-area) triangle: collapse to longest edge midpoint.
            // This avoids skipping the face entirely which would leave a hole
            // in closest-point queries near degenerate regions.
            Real e0_2 = (p1[0]-p0[0])*(p1[0]-p0[0])+(p1[1]-p0[1])*(p1[1]-p0[1])+(p1[2]-p0[2])*(p1[2]-p0[2]);
            Real e1_2 = (p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]);
            Real e2_2 = (p0[0]-p2[0])*(p0[0]-p2[0])+(p0[1]-p2[1])*(p0[1]-p2[1])+(p0[2]-p2[2])*(p0[2]-p2[2]);
            if (e0_2 >= e1_2 && e0_2 >= e2_2) {
                for (int d = 0; d < 3; ++d) c[d] = Real(0.5)*(p0[d]+p1[d]);
            } else if (e1_2 >= e2_2) {
                for (int d = 0; d < 3; ++d) c[d] = Real(0.5)*(p1[d]+p2[d]);
            } else {
                for (int d = 0; d < 3; ++d) c[d] = Real(0.5)*(p2[d]+p0[d]);
            }
            area = IBM_EPS::GEOM;  // assign minimal area so it participates in queries
            amrex::Print() << "Warning: degenerate face " << fi
                           << " collapsed to edge midpoint\n";
        }
        se_tmp.push_back(SurfElem(c, area, geomIdx));
        Real inv_n = Real(1.0) / std::sqrt(n_len2);
        Vec n_unit = {n_raw[0]*inv_n, n_raw[1]*inv_n, n_raw[2]*inv_n};
        Real t1_len2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
        Vec t1;
        if (t1_len2 > IBM_EPS::GEOM) {
            Real inv_t = Real(1.0) / std::sqrt(t1_len2);
            t1 = {v1[0]*inv_t, v1[1]*inv_t, v1[2]*inv_t};
        } else {
            Vec cross_x = {Real(0.0), -n_unit[2], n_unit[1]};
            Real cx2 = cross_x[0]*cross_x[0]+cross_x[1]*cross_x[1]+cross_x[2]*cross_x[2];
            if (cx2 < IBM_EPS::RAYCAST) {
                cross_x = {n_unit[2], Real(0.0), -n_unit[0]};
                cx2 = cross_x[0]*cross_x[0]+cross_x[1]*cross_x[1]+cross_x[2]*cross_x[2];
            }
            Real inv_c = Real(1.0) / std::sqrt(cx2);
            t1 = {cross_x[0]*inv_c, cross_x[1]*inv_c, cross_x[2]*inv_c};
        }
        Vec t2 = {n_unit[1]*t1[2]-n_unit[2]*t1[1],
                  n_unit[2]*t1[0]-n_unit[0]*t1[2],
                  n_unit[0]*t1[1]-n_unit[1]*t1[0]};
        Real n_arr[3]  = {n_unit[0], n_unit[1], n_unit[2]};
        Real t1_arr[3] = {t1[0], t1[1], t1[2]};
        Real t2_arr[3] = {t2[0], t2[1], t2[2]};
        lf_tmp.push_back(LocalFrame(n_arr, t1_arr, t2_arr));
    }
#endif
    // Bulk-append to ManagedVector (single resize + copy instead of many push_backs)
    std::size_t old_se = surfelem.size();
    std::size_t old_lf = localframe.size();
    surfelem.resize(old_se + se_tmp.size());
    localframe.resize(old_lf + lf_tmp.size());
    std::copy(se_tmp.begin(), se_tmp.end(), surfelem.begin() + old_se);
    std::copy(lf_tmp.begin(), lf_tmp.end(), localframe.begin() + old_lf);
}

/// GPU-parallel build of surface element data and local frames.
/// Pre-allocates output arrays and fills them with a single ParallelFor.
/// Degenerate faces (area ≈ 0) are included with zero measure — this avoids
/// GPU-side compaction and keeps face indices aligned with the geometry.
/// Called from build_geometry_cache when GPU is available and geometry is
/// already in managed memory (BVH path only — CGAL path is always CPU).
#ifdef AMREX_USE_GPU
inline void build_geometry_cache_gpu (
    const GeomType& geom,
    Gpu::ManagedVector<SurfElem>& surfelem,
    Gpu::ManagedVector<LocalFrame>& localframe,
    int geomIdx = -1)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    const int n = static_cast<int>(geom.size());
    const std::size_t old_se = surfelem.size();
    const std::size_t old_lf = localframe.size();
    surfelem.resize(old_se + n);
    localframe.resize(old_lf + n);

    auto* se_ptr = surfelem.data()   + old_se;
    auto* lf_ptr = localframe.data() + old_lf;
    const auto* v_ptr = geom.verts.data();
    const int gIdx = geomIdx;
    const int nv   = n;

    ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
        Point a = v_ptr[i];
        Point b = v_ptr[(i + 1) % nv];
        Real cx = Real(0.5) * (a[0] + b[0]);
        Real cy = Real(0.5) * (a[1] + b[1]);
        Real dx_ = b[0] - a[0], dy_ = b[1] - a[1];
        Real len = std::sqrt(dx_ * dx_ + dy_ * dy_);
        Real inv = (len > IBM_EPS::GEOM) ? Real(1.0) / len : Real(0.0);
        Real tx = dx_ * inv, ty = dy_ * inv;
        Real c[2] = {cx, cy};
        se_ptr[i] = SurfElem(c, len, gIdx);
        Real n_arr[2]  = {ty, -tx};
        Real t1_arr[2] = {tx, ty};
        lf_ptr[i] = LocalFrame(n_arr, t1_arr);
    });
    Gpu::streamSynchronize();

#elif (AMREX_SPACEDIM == 3)
    const int nf = geom.num_faces();
    const std::size_t old_se = surfelem.size();
    const std::size_t old_lf = localframe.size();
    surfelem.resize(old_se + nf);
    localframe.resize(old_lf + nf);

    auto* se_ptr = surfelem.data()   + old_se;
    auto* lf_ptr = localframe.data() + old_lf;
    const auto* vert_ptr = geom.vertices.data();
    const auto* face_ptr = geom.faces.data();
    const int gIdx = geomIdx;

    ParallelFor(nf, [=] AMREX_GPU_DEVICE (int fi) noexcept {
        const auto& f = face_ptr[fi];
        const Point& p0 = vert_ptr[f[0]];
        const Point& p1 = vert_ptr[f[1]];
        const Point& p2 = vert_ptr[f[2]];

        Real c[3];
        for (int d = 0; d < 3; ++d) c[d] = (p0[d] + p1[d] + p2[d]) / Real(3.0);

        Vec v1 = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        Vec v2 = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        Vec n_raw = {v1[1]*v2[2] - v1[2]*v2[1],
                     v1[2]*v2[0] - v1[0]*v2[2],
                     v1[0]*v2[1] - v1[1]*v2[0]};
        Real n_len2 = n_raw[0]*n_raw[0] + n_raw[1]*n_raw[1] + n_raw[2]*n_raw[2];
        Real area = Real(0.5) * std::sqrt(n_len2);

        if (area <= Real(0.0)) {
            // Degenerate face: write zero-measure entry (preserves index alignment)
            se_ptr[fi] = SurfElem(c, Real(0.0), gIdx);
            lf_ptr[fi] = LocalFrame();
            return;
        }

        se_ptr[fi] = SurfElem(c, area, gIdx);

        Real inv_n = Real(1.0) / std::sqrt(n_len2);
        Vec n_unit = {n_raw[0]*inv_n, n_raw[1]*inv_n, n_raw[2]*inv_n};

        Real t1_len2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
        Vec t1;
        if (t1_len2 > IBM_EPS::GEOM) {
            Real inv_t = Real(1.0) / std::sqrt(t1_len2);
            t1 = {v1[0]*inv_t, v1[1]*inv_t, v1[2]*inv_t};
        } else {
            Vec cross_x = {Real(0.0), -n_unit[2], n_unit[1]};
            Real cx2 = cross_x[0]*cross_x[0] + cross_x[1]*cross_x[1] + cross_x[2]*cross_x[2];
            if (cx2 < IBM_EPS::RAYCAST) {
                cross_x = {n_unit[2], Real(0.0), -n_unit[0]};
                cx2 = cross_x[0]*cross_x[0] + cross_x[1]*cross_x[1] + cross_x[2]*cross_x[2];
            }
            Real inv_c = Real(1.0) / std::sqrt(cx2);
            t1 = {cross_x[0]*inv_c, cross_x[1]*inv_c, cross_x[2]*inv_c};
        }

        Vec t2 = {n_unit[1]*t1[2] - n_unit[2]*t1[1],
                  n_unit[2]*t1[0] - n_unit[0]*t1[2],
                  n_unit[0]*t1[1] - n_unit[1]*t1[0]};

        Real n_arr[3]  = {n_unit[0], n_unit[1], n_unit[2]};
        Real t1_arr[3] = {t1[0], t1[1], t1[2]};
        Real t2_arr[3] = {t2[0], t2[1], t2[2]};
        lf_ptr[fi] = LocalFrame(n_arr, t1_arr, t2_arr);
    });
    Gpu::streamSynchronize();
#endif
}
#endif // AMREX_USE_GPU

/// Flip normals for interior-is-fluid mode.
inline void convert_inout (Gpu::ManagedVector<LocalFrame>& localframe_a) {
    for (auto& lf : localframe_a) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lf.normal[d]   = -lf.normal[d];
            lf.tangent1[d] = -lf.tangent1[d];
        }
    }
}

/// Check that IBM geometries do not intersect or contain each other.
inline void check_ibm_geometry_consistency (
    int                  ngeom,
    const GeomType*      geom_a,
    inside_t* const*     inout_fa,
    const std::string*   files_a)
{
    for (int i = 0; i < ngeom; ++i) {
        for (int j = i + 1; j < ngeom; ++j) {
            AABB bi = geom_a[i].bbox();
            AABB bj = geom_a[j].bbox();
            if (!bi.overlaps(bj)) continue;
#if (AMREX_SPACEDIM == 2)
            int ni = static_cast<int>(geom_a[i].size());
            int nj = static_cast<int>(geom_a[j].size());
            bool hit = false;
            for (int ei = 0; ei < ni && !hit; ++ei) {
                Point a0 = geom_a[i].vertex(ei);
                Point a1 = geom_a[i].vertex((ei+1)%ni);
                for (int ej = 0; ej < nj; ++ej) {
                    Point b0 = geom_a[j].vertex(ej);
                    Point b1 = geom_a[j].vertex((ej+1)%nj);
                    auto cross2 = [](const Point& a, const Point& b, const Point& c) {
                        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
                    };
                    Real d1 = cross2(b0, b1, a0), d2 = cross2(b0, b1, a1);
                    Real d3 = cross2(a0, a1, b0), d4 = cross2(a0, a1, b1);
                    if (((d1>0&&d2<0)||(d1<0&&d2>0)) && ((d3>0&&d4<0)||(d3<0&&d4>0))) {
                        hit = true; break;
                    }
                }
            }
            if (hit) {
                amrex::Print() << "ERROR: IBM geometries intersect:\n"
                               << "  geom " << i << " : " << files_a[i] << "\n"
                               << "  geom " << j << " : " << files_a[j] << "\n";
                amrex::Abort("IBM geometries must not intersect.");
            }
#endif
            constexpr int K = 8;
            auto check_containment = [&](int inner, int outer) {
                const auto& g = geom_a[inner];
                inside_t& tester = *inout_fa[outer];
#if (AMREX_SPACEDIM == 2)
                int np = static_cast<int>(g.size());
                int step = std::max(1, np / K);
                int ins = 0, checked = 0;
                for (int k = 0; k < np && checked < K; k += step, ++checked) {
                    BoundedSide res = tester(g.vertex(k));
                    if (res == BoundedSide::OnBoundary)
                        amrex::Abort("IBM geometries touching (vertex on boundary).");
                    if (res == BoundedSide::Inside) ins++;
                }
#else
                int np = g.num_vertices();
                int step = std::max(1, np / K);
                int ins = 0, checked = 0;
                for (int k = 0; k < np && checked < K; k += step, ++checked) {
                    BoundedSide res = tester(g.vertices[k]);
                    if (res == BoundedSide::OnBoundary)
                        amrex::Abort("IBM geometries touching (vertex on boundary).");
                    if (res == BoundedSide::Inside) ins++;
                }
#endif
                if (checked > 0 && ins == checked) {
                    amrex::Print() << "WARNING: geom " << inner << " is inside geom " << outer << "\n";
                }
            };
            check_containment(i, j);
            check_containment(j, i);
        }
    }
}

// ============================================================================
// 4. UTILITY FUNCTIONS
// ============================================================================

/// Check if Point is inside AABB
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool bbox_contains (const AABB& bb, const Point& p) {
    return bb.contains(p);
}

/// Create an AMReX Array1D from scalars
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Array1D<T, 0, AMREX_SPACEDIM-1>
make_vec (T x, T y, T z) noexcept {
    Array1D<T, 0, AMREX_SPACEDIM-1> a;
    AMREX_D_TERM(a(0) = x;, a(1) = y;, a(2) = z;);
    return a;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Array1D<T, 0, AMREX_SPACEDIM-1>
make_vec (T x, T y) noexcept {
    return make_vec<T>(x, y, T(0));
}

template <typename T, typename PointT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Array1D<T, 0, AMREX_SPACEDIM-1>
make_vec (const PointT& p) noexcept {
    Array1D<T, 0, AMREX_SPACEDIM-1> a;
#if (AMREX_SPACEDIM == 2)
    a(0) = static_cast<T>(p[0]);
    a(1) = static_cast<T>(p[1]);
#elif (AMREX_SPACEDIM == 3)
    a(0) = static_cast<T>(p[0]);
    a(1) = static_cast<T>(p[1]);
    a(2) = static_cast<T>(p[2]);
#endif
    return a;
}

/// Warn if a grid point lies on the IB surface (CPU-only: uses amrex::Print)
AMREX_FORCE_INLINE
void IB_WarnOnBoundary (int ii, int level, int i, int j, int k,
                        const BoundedSide& result, const Point& gridpoint)
{
    if (result == BoundedSide::OnBoundary) {
        amrex::Print()
            << "Warning: Grid point on IB surface\n"
            << "  geom ii = " << ii << "\n"
            << "  level   = " << level << "\n"
            << "  cell    = (" << i << ", " << j
#if (AMREX_SPACEDIM == 3)
            << ", " << k
#endif
            << ")\n"
            << "  point   = ("
            << gridpoint[0] << ", " << gridpoint[1]
#if (AMREX_SPACEDIM == 3)
            << ", " << gridpoint[2]
#endif
            << ")\n";
    }
}

#endif // IBM_BACKEND_BVH_H_
