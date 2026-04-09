#ifndef IBM_BACKEND_CGAL_H_
#define IBM_BACKEND_CGAL_H_

// ============================================================================
// ibm_backend_cgal.h — CGAL backend implementation (functions and algorithms)
//
// Implements the interface between the IBM solver and CGAL geometry:
//   1. Geometry I/O           : read_polygon_2d, check_3d_mesh_validity
//   2. Geometry processing    : build_geometry_cache, convert_inout,
//                               check_ibm_geometry_consistency
//   3. Utility functions      : make_vec, bbox_contains, IB_WarnOnBoundary,
//                               make_grid_point, point_distance_sq,
//                               cgal_closest_point_query
//
// Type definitions are in ibm_cgal_defs.h.
// ============================================================================

#include "ibm_cgal_defs.h"

#include <AMReX_Print.H>
#include <AMReX_GpuContainers.H>

#include <fstream>
#include <sstream>

//============================================================================
// FUNCTION IMPLEMENTATIONS
//============================================================================

//---------------------------------------------------------------------------
// Helper: Check if Point is inside BBox
//---------------------------------------------------------------------------
AMREX_FORCE_INLINE bool bbox_contains(const Bbox& bb, const Point& p)
{
#if (AMREX_SPACEDIM == 2)
    return p.x() >= bb.xmin() && p.x() <= bb.xmax() &&
           p.y() >= bb.ymin() && p.y() <= bb.ymax();
#else
    return p.x() >= bb.xmin() && p.x() <= bb.xmax() &&
           p.y() >= bb.ymin() && p.y() <= bb.ymax() &&
           p.z() >= bb.zmin() && p.z() <= bb.zmax();
#endif
}

//---------------------------------------------------------------------------
// Compatibility: make_grid_point — construct a CGAL Point from grid indices
//---------------------------------------------------------------------------
template <typename ArrayLike>
AMREX_FORCE_INLINE
Point make_grid_point(const ArrayLike& prob_lo, const ArrayLike& dx,
                      int i, int j, int k)
{
    Real x = prob_lo[0] + (Real(0.5) + Real(i)) * dx[0];
    Real y = prob_lo[1] + (Real(0.5) + Real(j)) * dx[1];
#if (AMREX_SPACEDIM == 3)
    Real z = prob_lo[2] + (Real(0.5) + Real(k)) * dx[2];
    return Point(x, y, z);
#else
    amrex::ignore_unused(k);
    return Point(x, y);
#endif
}

//---------------------------------------------------------------------------
// Compatibility: point_distance_sq — squared distance between two CGAL Points
//---------------------------------------------------------------------------
AMREX_FORCE_INLINE
Real point_distance_sq(const Point& a, const Point& b)
{
    Real d2 = (a.x() - b.x()) * (a.x() - b.x())
            + (a.y() - b.y()) * (a.y() - b.y());
#if (AMREX_SPACEDIM == 3)
    d2 += (a.z() - b.z()) * (a.z() - b.z());
#endif
    return d2;
}

//---------------------------------------------------------------------------
// Compatibility: closest_point_query via CGAL tree — returns ClosestPointResult
//   with local prim_id (subtract offset so caller can add geom_offsets[geomIdx])
//---------------------------------------------------------------------------
inline ClosestPointResult cgal_closest_point_query(
    const Tree& tree, const PrimitiveIndexMap& idxmap,
    const Point& query, int offset)
{
    auto ppid = tree.closest_point_and_primitive(query);
    ClosestPointResult res;
    res.point   = ppid.first;
    res.prim_id = idxmap.at(ppid.second) - offset;
    return res;
}

#if (AMREX_SPACEDIM == 2)

//----------------------------------------------------------------------------
// 2D Polygon File Reader (inline bool read_polygon_2d)
//----------------------------------------------------------------------------
/// \brief Reads a 2D polygon from a text file.
///        (to be constructed: extend to other 2D mesh formats in the future.)
///
/// File Format:
///   - One vertex per line
///   - Each line contains two coordinates: "x y"  or  "x, y"
///   - Empty lines and lines starting with '#' or '//' are skipped
///   - Commas, semicolons, and tabs are treated as delimiters
///
/// Example file:
///   \code
///   # Rectangle
///   0.0 0.0
///   1.0, 0.0
///   1.0 1.0
///   0.0, 1.0
///   \endcode
///
/// Post-processing:
///   - Validates that the polygon has at least 3 vertices
///   - Validates that the polygon is simple (no self-intersections)
///   - Optionally refines edges so that their length does not exceed dx = min_dx/2
///   - Enforces counter-clockwise (CCW) orientation
///
/// \param filename Path to the input file.
/// \param poly     Output polygon (will be cleared first).
/// \param dx       Optional maximum edge length. If dx > 0, each polygon edge
///                 is uniformly subdivided so that every segment length is
///                 <= min_dx/2. If dx <= 0, no refinement is performed.
///
/// \return true on success, false if file cannot be opened, input is invalid,
///         or the polygon fails validation.
inline bool read_polygon_2d(const std::string& filename, Polygon2D& poly, Real dx = -1.0)
{
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

        // Trim leading and trailing whitespace
        auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            // Entire line is whitespace: skip
            continue;
        }
        auto last = line.find_last_not_of(" \t\r\n");
        line = line.substr(first, last - first + 1);

        // Skip comment lines
        if (line[0] == '#' || line.compare(0, 2, "//") == 0) {
            continue;
        }

        // Replace common delimiters with spaces
        for (char& c : line) {
            if (c == ',' || c == ';' || c == '\t') {
                c = ' ';
            }
        }

        // Parse coordinates
        std::istringstream iss(line);
        Real x, y;

        if (!(iss >> x >> y)) {
            // Invalid line - skip with warning (do not abort the whole file)
            amrex::Print() << "Warning: Skipping invalid line " << line_number
                           << " in " << filename << "\n";
            continue;
        }

        // Validate numeric values
        if (!std::isfinite(x) || !std::isfinite(y)) {
            amrex::Print() << "read_polygon_2d: Invalid coordinate (NaN/Inf) on line "
                           << line_number << " in " << filename << "\n";
            return false;
        }

        raw_points.push_back(Point(x, y));
    }

    if (raw_points.empty()) {
        amrex::Print() << "read_polygon_2d: No valid vertices found in " << filename << "\n";
        return false;
    }

    // Compute Geometry Scale & Tolerances
    Real min_x = raw_points[0].x(), max_x = min_x;
    Real min_y = raw_points[0].y(), max_y = min_y;
    for (const auto& p : raw_points) {
        min_x = std::min(min_x, p.x()); max_x = std::max(max_x, p.x());
        min_y = std::min(min_y, p.y()); max_y = std::max(max_y, p.y());
    }
    Real bbox_diag = std::sqrt(std::pow(max_x - min_x, 2) + std::pow(max_y - min_y, 2));

    // Use user-provided dx as reference if available, otherwise bounding box
    Real scale_ref = (dx > 0) ? dx : (bbox_diag > 1e-12 ? bbox_diag : 1.0);

    // Determine tolerances
    // - deduplication: very tight (micro-gaps)
    // - area: check for collapse
    Real eps_dedup = 1e-6 * scale_ref;
    Real area_eps  = 1e-12 * scale_ref * scale_ref;

    // Vertex Deduplication & Cleanup
    std::vector<Point> clean_points;
    clean_points.reserve(raw_points.size());
    clean_points.push_back(raw_points[0]);

    int n_dups = 0;
    for (size_t i = 1; i < raw_points.size(); ++i) {
        Real dx_pt = clean_points.back().x() - raw_points[i].x();
        Real dy_pt = clean_points.back().y() - raw_points[i].y();
        Real d2 = dx_pt*dx_pt + dy_pt*dy_pt;
        if (d2 > eps_dedup * eps_dedup) {
            clean_points.push_back(raw_points[i]);
        } else {
            n_dups++;
        }
    }

    // Remove "closing" vertex if it duplicates the start vertex
    if (clean_points.size() > 1) {
        Real dx_end = clean_points.back().x() - clean_points.front().x();
        Real dy_end = clean_points.back().y() - clean_points.front().y();
        Real d2 = dx_end*dx_end + dy_end*dy_end;
        if (d2 <= eps_dedup * eps_dedup) {
            clean_points.pop_back();
            n_dups++;
        }
    }

    if (n_dups > 0) {
        amrex::Print() << "Info: Removed " << n_dups << " duplicate/close vertices (tol="
                       << eps_dedup << ") in " << filename << "\n";
    }

    // Build Polygon
    for (const auto& p : clean_points) {
        poly.push_back(p);
    }

    // Validate vertex count
    if (poly.size() < 3) {
        amrex::Print() << "read_polygon_2d: Polygon must have at least 3 unique vertices (found "
                       << poly.size() << ") in: " << filename << "\n";
        return false;
    }

    // Validate degenerate area (collinear points etc.)
    if (std::abs(poly.area()) <= area_eps) {
        amrex::Print() << "read_polygon_2d: Polygon has near-zero area (" << poly.area()
                       << " <= " << area_eps << ") in: " << filename << "\n";
        return false;
    }

    // Validate simplicity (no self-intersections)
    if (!poly.is_simple()) {
        amrex::Print() << "read_polygon_2d: Polygon has self-intersections in: "
                       << filename << "\n";
        return false;
    }

    std::size_t original_n = poly.size();

    // Optional edge refinement:
    //   If dx > 0, each edge (p0 -> p1) is subdivided into nseg segments,
    //   where nseg = ceil(|p1 - p0| / dx). New points are inserted by
    //   uniform linear interpolation along the edge:
    //
    //       t_k = k / nseg,  k = 1, ..., nseg-1
    //
    //   so that each resulting segment has length approximately |p1-p0| / nseg,
    //   which is <= dx. This avoids creating very short "leftover" segments.
    if (dx > Real(0)) {
        Polygon2D refined;

        const std::size_t n = poly.size();
        for (std::size_t i = 0; i < n; ++i) {
            const Point& p0 = poly.vertex(i);
            const Point& p1 = poly.vertex((i + 1) % n);

            Real x0 = p0.x();
            Real y0 = p0.y();
            Real x1 = p1.x();
            Real y1 = p1.y();

            // Always add the current vertex (start of the edge)
            refined.push_back(p0);

            Real vx = x1 - x0;
            Real vy = y1 - y0;
            Real len = std::sqrt(vx*vx + vy*vy);

            if (len > dx) {
                // Split this edge into nseg segments so that each segment length <= dx.
                int nseg = static_cast<int>(std::ceil(len / dx));
                if (nseg < 1) nseg = 1;  // Safety guard (should not happen when len > dx)

                Real inv = Real(1.0) / Real(nseg);
                for (int k = 1; k < nseg; ++k) {
                    Real t = inv * Real(k);  // 0 < t < 1
                    Real xn = x0 + t * vx;
                    Real yn = y0 + t * vy;
                    refined.push_back(Point(xn, yn));
                }
            }
        }

        // Replace the original polygon with the refined one.
        poly = std::move(refined);
    }

    // Enforce counter-clockwise orientation for consistency.
    // CGAL (and many geometric algorithms) commonly expect CCW polygons.
    // Note: reverse_orientation() internally calls finalize()
    if (poly.is_clockwise_oriented()) {
        poly.reverse_orientation();
        amrex::Print() << "Info: Reversed polygon orientation to CCW in "
                       << filename << "\n";
    } else {
        // If not reversed, we still need to build edges
        poly.finalize();
    }

    if (dx > Real(0)) {
        amrex::Print() << "Successfully loaded and refined polygon from "
                       << original_n << " to " << poly.size()
                       << " vertices (target max edge length <= " << dx
                       << ") from " << filename << "\n";
    } else {
        amrex::Print() << "Successfully loaded polygon with " << poly.size()
                       << " vertices from " << filename << "\n";
    }

    return true;
}

#endif  // AMREX_SPACEDIM == 2

#if (AMREX_SPACEDIM == 3)

/// \brief Check 3D mesh validity (closedness, etc.)
inline bool check_3d_mesh_validity(const Polyhedron& geom, const std::string& name)
{
    // 1. Check if mesh is closed (watertight)
    if (!geom.is_closed()) {
        amrex::Print() << "ERROR: Mesh " << name << " is not closed (watertight).\n"
                       << "       IBM requires a closed surface for inside/outside test.\n";
        return false;
    }
    // NOTE: Self-intersection check (PMP::does_self_intersect) omitted for performance.
    //       Enable if mesh quality issues are suspected.
    return true;
}

/// \brief Compute plane equation for a facet
inline void compute_plane_equations(Polyhedron::Facet& f)
{
    Polyhedron::Halfedge_handle h = f.halfedge();
    f.plane() = Polyhedron::Plane_3(
        h->opposite()->vertex()->point(),
        h->vertex()->point(),
        h->next()->vertex()->point()
    );
}

#endif // AMREX_SPACEDIM == 3

// -----------------------------------------------------------------------------
// build_geometry_cache
// -----------------------------------------------------------------------------
/// \brief Computes all geometric data (SurfElem, LocalFrame, and ID mapping) in a single pass.
///
/// This function iterates over the geometry elements (edges in 2D, faces in 3D) exactly once.
/// It populates the surface element data (centroid, size), the local orthonormal frames,
/// and the mapping from CGAL PrimitiveID to integer index.
///
/// Using a single pass ensures that the integer indices in idxmap perfectly align
/// with the array indices in `surfelem` and `localframe`.
///
/// \param geom       The input geometry (Polygon2D or Polyhedron).
/// \param surfelem   Output vector for surface element properties (appended to).
/// \param localframe Output vector for local frames (appended to).
/// \param idxmap     Output map from PrimitiveID to integer index (cleared before use).
/// \param offset     Global index offset for this geometry (default 0).
/// \param geomIdx    Index of the geometry (body) being processed.
inline void build_geometry_cache(
    const GeomType& geom,
    amrex::Gpu::ManagedVector<SurfElem>& surfelem,
    amrex::Gpu::ManagedVector<LocalFrame>& localframe,
    PrimitiveIndexMap& idxmap,
    int offset = 0,
    int geomIdx = -1)
{
    // Clear map for this geometry
    idxmap.clear();

#if (AMREX_SPACEDIM == 2)
    // -----------------------------------------------------------------------
    // 2D Implementation (Polygon edges)
    // -----------------------------------------------------------------------
    std::size_t n_elems = geom.size();

    int local_idx = 0;
    for (auto eit = geom.edges_begin(); eit != geom.edges_end(); ++eit, ++local_idx) {
        // --- 1. ID Mapping ---
        PrimitiveID pid = eit;
        idxmap[pid] = offset + local_idx;

        // --- 2. SurfElem (Centroid & Length) ---
        Segment s = *eit;
        Point mid = CGAL::midpoint(s.source(), s.target());
        Real len = std::sqrt(s.squared_length());
        Real c[2] = { mid.x(), mid.y() };
        surfelem.push_back(SurfElem(c, len, geomIdx));

        // --- 3. LocalFrame (Normal & Tangent) ---
        Vector_CGAL t_vec = s.to_vector();
        Real len2 = t_vec.squared_length();

        // Handle zero-length edge gracefully (though read_polygon_2d should prevent this)
        Vector_CGAL t1;
        if (len2 > 0) {
            t1 = t_vec / std::sqrt(len2);
        } else {
            t1 = Vector_CGAL(1.0, 0.0);
        }

        // Outward normal: rotate tangent -90 degrees (tx, ty) -> (ty, -tx)
        Vector_CGAL n(t1.y(), -t1.x());

        Real n_arr[2]  = { n.x(), n.y() };
        Real t1_arr[2] = { t1.x(), t1.y() };
        localframe.push_back(LocalFrame(n_arr, t1_arr));
    }

#elif (AMREX_SPACEDIM == 3)
    // -----------------------------------------------------------------------
    // 3D Implementation (Polyhedron faces)
    // -----------------------------------------------------------------------
    std::size_t n_elems = geom.size_of_facets();

    int local_idx = 0;
    for (auto f = faces(geom).first; f != faces(geom).second; ++f, ++local_idx) {
        // --- 1. ID Mapping ---
        elm_descriptor fd = *f;
        PrimitiveID pid = fd;
        idxmap[pid] = offset + local_idx;

        // --- 2. SurfElem (Centroid & Area) ---
        auto h = fd->halfedge();
        Point p1 = h->vertex()->point();
        Point p2 = h->next()->vertex()->point();
        Point p3 = h->next()->next()->vertex()->point();

        Point cent = CGAL::centroid(p1, p2, p3);
        Vector_CGAL v1 = p2 - p1;
        Vector_CGAL v2 = p3 - p1;
        Real area = 0.5 * std::sqrt(CGAL::cross_product(v1, v2).squared_length());
        if (area <= Real(0)) {
            amrex::Abort("build_geometry_cache: detected non-positive face area");
        }

        Real c[3] = { cent.x(), cent.y(), cent.z() };
        surfelem.push_back(SurfElem(c, area, geomIdx));

        // --- 3. LocalFrame (Normal, Tangent1, Tangent2) ---
        Vector_CGAL n = PMP::compute_face_normal(*f, geom);
        Real n_len2 = n.squared_length();
        if (n_len2 > 0) n = n / std::sqrt(n_len2);

        Vector_CGAL t1_raw = p2 - p1;
        Real t1_len2 = t1_raw.squared_length();
        Vector_CGAL t1;
        if (t1_len2 > 0) {
            t1 = t1_raw / std::sqrt(t1_len2);
        } else {
            // Fallback for degenerate edge
            t1 = CGAL::cross_product(n, Vector_CGAL(1,0,0));
            if (t1.squared_length() < 1e-12) t1 = CGAL::cross_product(n, Vector_CGAL(0,1,0));
            t1 = t1 / std::sqrt(t1.squared_length());
        }
        Vector_CGAL t2 = CGAL::cross_product(n, t1);

        Real n_arr[3]  = { n.x(), n.y(), n.z() };
        Real t1_arr[3] = { t1.x(), t1.y(), t1.z() };
        Real t2_arr[3] = { t2.x(), t2.y(), t2.z() };
        localframe.push_back(LocalFrame(n_arr, t1_arr, t2_arr));
    }
#endif
}

// -----------------------------------------------------------------------------
// convert_inout
// -----------------------------------------------------------------------------
/// \brief Flip stored local frames so that normals point inward.
///
/// If switch the interior_is_solid, geometry interior is treated as fluid,
/// typically flip the normals.
///
/// Tangent handling:
/// - We also flip tangent1. In 3D this preserves a right-handed orthonormal
///   frame because t_2 = n*t_1 stays
///   invariant when both n and t_1 are negated.
/// - tangent2 (3D) is left unchanged.
inline void convert_inout(amrex::Gpu::ManagedVector<LocalFrame>& localframe_a)
{
    for (auto& lf : localframe_a) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            lf.normal[d]   = -lf.normal[d];
            lf.tangent1[d] = -lf.tangent1[d];
        }
    }
}

//============================================================================
// GEOMETRY CONSISTENCY CHECKING
//============================================================================

/// \brief Check that IBM geometries do not intersect or touch each other.
///
/// This function validates that multiple IBM geometries are properly separated.
///
/// Level 1: Surface intersection/touching check
///   - In 2D: Checks if polygon edges intersect or touch
///   - In 3D: Checks if mesh surfaces intersect or touch
///
/// Level 2: Volumetric containment check
///   - Tests if one geometry is completely inside another
///   - Uses representative vertices from each geometry
///   - Containment is reported as a WARNING (not fatal)
///
/// \param ngeom Number of geometries to check
/// \param geom_a Array of geometry objects (Polygon_2 or Polyhedron_3)
/// \param inout_fa Array of inside/outside testers
/// \param files_a Array of filenames (for error reporting)
inline void check_ibm_geometry_consistency(
    int                 ngeom,
    const GeomType*     geom_a,
    inside_t* const*    inout_fa,
    const std::string*  files_a)
{
    for (int i = 0; i < ngeom; ++i) {
        for (int j = i + 1; j < ngeom; ++j) {

#if (AMREX_SPACEDIM == 2)
            //----------------------------------------------------------------
            // 2D CHECKS
            //----------------------------------------------------------------

            // Skip empty polygons
            if (geom_a[i].size() < 3 || geom_a[j].size() < 3) {
                continue;
            }

            // Quick BBox rejection (Optimization)
            // If bounding boxes do not overlap, polygons cannot intersect, touch, or contain each other.
            if (!CGAL::do_overlap(geom_a[i].bbox(), geom_a[j].bbox())) {
                continue;
            }

            // 1) Check for edge-edge intersections or touching
            //    Test every edge of polygon i against every edge of polygon j
            bool surfaces_intersect = false;

            for (auto ei = geom_a[i].edges_begin();
                 ei != geom_a[i].edges_end() && !surfaces_intersect;
                 ++ei) {
                for (auto ej = geom_a[j].edges_begin();
                     ej != geom_a[j].edges_end();
                     ++ej) {

                    // Check if segments intersect (including touching)
                    if (CGAL::do_intersect(*ei, *ej)) {
                        surfaces_intersect = true;
                        break;
                    }
                }
            }

            if (surfaces_intersect) {
                amrex::Print() << "ERROR: IBM geometries intersect or touch:\n"
                              << "  geom " << i << " : " << files_a[i] << "\n"
                              << "  geom " << j << " : " << files_a[j] << "\n";
                amrex::Abort(
                    "IBM geometries must not intersect or touch each other.");
            }

            // 2) Check for containment (Level 2)
            //    Robust check using multiple sample points (K=8)
            constexpr int K_samples = 8;

            auto check_containment_2d = [&](int idx_inner, int idx_outer) {
                const auto& poly_in = geom_a[idx_inner];
                // Use non-const ref just in case the functor is not const-correct in all versions
                inside_t& tester = *inout_fa[idx_outer];

                std::size_t n_pts = poly_in.size();
                if (n_pts < 3) return;

                std::size_t step = std::max(std::size_t(1), n_pts / K_samples);

                int inside_count = 0;
                int checked_count = 0;

                // stride sampling over vertices using iterators for robustness
                auto vbegin = poly_in.vertices_begin();
                auto vend   = poly_in.vertices_end();

                for (std::size_t k = 0; k < n_pts && checked_count < K_samples; k += step) {
                    auto vit = vbegin;
                    std::advance(vit, static_cast<long>(k));

                    const Point& p = *vit;
                    BoundedSide res = tester(p);
                    checked_count++;

                    if (res == BoundedSide::OnBoundary) {
                        amrex::Print() << "ERROR: IBM geometry touching detected (vertex on boundary):\n"
                                       << "  geom " << idx_inner << " : " << files_a[idx_inner] << "\n"
                                       << "  touches geom " << idx_outer << " : " << files_a[idx_outer] << "\n";
                        amrex::Abort("IBM geometries must not intersect or touch each other.");
                    }
                    if (res == BoundedSide::Inside) {
                        inside_count++;
                    }
                }

                // If ALL sampled points are inside, we declare containment.
                if (checked_count > 0 && inside_count == checked_count) {
                    amrex::Print() << "WARNING: IBM geometry containment detected (Robust check):\n"
                                   << "  geom " << idx_inner << " : " << files_a[idx_inner]
                                   << " is strictly inside\n"
                                   << "  geom " << idx_outer << " : " << files_a[idx_outer] << "\n";
                    amrex::Print() << "Continuing, but results may be undefined depending on the IBM setup.\n";
                }
            };

            check_containment_2d(i, j); // Check i inside j
            check_containment_2d(j, i); // Check j inside i

#else  // AMREX_SPACEDIM == 3
            //----------------------------------------------------------------
            // 3D CHECKS
            //----------------------------------------------------------------

            // 1) Surface intersection/touching check
            //    PMP::do_intersect returns true if there is any non-empty
            //    intersection between the two meshes:
            //      - triangle penetration
            //      - face/edge/vertex touching
            if (PMP::do_intersect(geom_a[i], geom_a[j])) {
                amrex::Print() << "ERROR: IBM geometries intersect or touch:\n"
                              << "  geom " << i << " : " << files_a[i] << "\n"
                              << "  geom " << j << " : " << files_a[j] << "\n";
                amrex::Abort(
                    "IBM geometries must not intersect or touch each other.");
            }

            // 2) Containment check (Level 2)
            if (geom_a[i].empty() || geom_a[j].empty()) continue;

            constexpr int K_samples = 8;

            auto check_containment_3d = [&](int idx_inner, int idx_outer) {
                const auto& mesh_in = geom_a[idx_inner];
                inside_t& tester = *inout_fa[idx_outer];

                // Polyhedron_3 supports size_of_vertices()
                std::size_t n_pts = mesh_in.size_of_vertices();
                if (n_pts == 0) return;

                std::size_t step = std::max(std::size_t(1), n_pts / K_samples);

                int inside_count = 0;
                int checked_count = 0;

                auto vit = mesh_in.vertices_begin();
                auto vend = mesh_in.vertices_end();

                while (vit != vend && checked_count < K_samples) {
                    const Point& p = vit->point();
                    BoundedSide res = tester(p);
                    checked_count++;

                    if (res == BoundedSide::OnBoundary) {
                        amrex::Print() << "ERROR: IBM geometry touching detected (vertex on boundary):\n"
                                       << "  geom " << idx_inner << " : " << files_a[idx_inner] << "\n"
                                       << "  touches geom " << idx_outer << " : " << files_a[idx_outer] << "\n";
                        amrex::Abort("IBM geometries must not intersect or touch each other.");
                    }
                    if (res == BoundedSide::Inside) {
                        inside_count++;
                    }

                    // advance by 'step', explicit separate loop for safety
                    for (std::size_t s = 0; s < step && vit != vend; ++s) {
                        ++vit;
                    }
                }

                if (checked_count > 0 && inside_count == checked_count) {
                     amrex::Print() << "WARNING: IBM geometry containment detected (Robust check):\n"
                                   << "  geom " << idx_inner << " : " << files_a[idx_inner]
                                   << " is strictly inside\n"
                                   << "  geom " << idx_outer << " : " << files_a[idx_outer] << "\n";
                    amrex::Print() << "Continuing, but results may be undefined depending on the IBM setup.\n";
                }
            };

            check_containment_3d(i, j);
            check_containment_3d(j, i);
#endif  // AMREX_SPACEDIM
        }
    }
}

//============================================================================
/// \brief Create an AMReX Array1D holding vector-like data.
///
/// This utility supports:
///   1) make_vec<T>(x, y, z)   : three scalars of any type T (int, Real, etc.)
///   2) make_vec<T>(x, y)      : two scalars, z is automatically set to zero
///   3) make_vec<T>(point)    : 2D/3D Point converted to Array1D<T>
///
/// Notes:
///   - When AMREX_SPACEDIM == 2, only a(0) and a(1) are used.
///   - When AMREX_SPACEDIM == 3, a(0), a(1), and a(2) are used.
///   - For the Point overload, the coordinates are cast to type T.

//------------------------------------------------------------------------------
// 1) Generic scalar version: x, y, z
//------------------------------------------------------------------------------
template <typename T>
AMREX_FORCE_INLINE
Array1D<T, 0, AMREX_SPACEDIM-1>
make_vec(T x, T y, T z) noexcept
{
    Array1D<T, 0, AMREX_SPACEDIM-1> a;
    AMREX_D_TERM( a(0) = x;,
                  a(1) = y;,
                  a(2) = z; );
    return a;
}

//------------------------------------------------------------------------------
// 2) Generic scalar version: x, y (z is set to zero automatically)
//------------------------------------------------------------------------------
template <typename T>
AMREX_FORCE_INLINE
Array1D<T, 0, AMREX_SPACEDIM-1>
make_vec(T x, T y) noexcept
{
    return make_vec<T>(x, y, T(0));
}

//------------------------------------------------------------------------------
// 3) Point-based version: convert a 2D/3D Point to Array1D<T>
//------------------------------------------------------------------------------
template <typename T, typename PointT>
AMREX_FORCE_INLINE
Array1D<T, 0, AMREX_SPACEDIM-1>
make_vec(const PointT& p) noexcept
{
    Array1D<T, 0, AMREX_SPACEDIM-1> a;

#if (AMREX_SPACEDIM == 2)
    a(0) = static_cast<T>(p.x());
    a(1) = static_cast<T>(p.y());
#elif (AMREX_SPACEDIM == 3)
    a(0) = static_cast<T>(p.x());
    a(1) = static_cast<T>(p.y());
    a(2) = static_cast<T>(p.z());
#endif

    return a;
}

//============================================================================
/// \brief Warn if a grid point lies exactly on the IB surface boundary.
///
AMREX_FORCE_INLINE
void IB_WarnOnBoundary(int ii,
                       int level,
                       int i, int j, int k,
                       const BoundedSide& result,
                       const Point& gridpoint)
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
            << gridpoint.x() << ", "
            << gridpoint.y()
#if (AMREX_SPACEDIM == 3)
            << ", " << gridpoint.z()
#endif
            << ")\n";
    }
}

#endif  // IBM_BACKEND_CGAL_H_
