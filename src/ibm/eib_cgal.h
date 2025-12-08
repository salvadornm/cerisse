#ifndef EIB_CGAL_H_
#define EIB_CGAL_H_

// AMReX headers (must come first)
#include <AMReX.H>
#include <AMReX_Array.H>

// Standard library headers
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <cmath>

// Basic CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/enum.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#if (AMREX_SPACEDIM == 2)

    // CGAL headers for AABB tree in 2D
    #include <CGAL/AABB_traits_2.h>
    #include <CGAL/AABB_segment_primitive_2.h> 

    // CGAL header for inout testing in 2D
    #include <CGAL/Polygon_2.h>
    #include <CGAL/Polygon_2_algorithms.h>

    // CGAL header for CGAL::do_intersect in 2D
    #include <CGAL/intersections.h>

#elif (AMREX_SPACEDIM == 3)

    // CGAL headers for AABB tree in 3D
    #include <CGAL/AABB_traits_3.h>
    #include <CGAL/AABB_face_graph_triangle_primitive.h>

    // CGAL header for inout testing in 3D
    #include <CGAL/Polyhedron_3.h>
    #include <CGAL/Side_of_triangle_mesh.h>

    // CGAL headers for Surface normal computation
    #include <CGAL/Polygon_mesh_processing/compute_normal.h>
    #include <CGAL/Polygon_mesh_processing/orientation.h>

    // CGAL headers for face area measurement
    #include <CGAL/centroid.h>
    #include <CGAL/Polygon_mesh_processing/measure.h>

    // BGL traits for Polyhedron_3 with PMP
    #include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

    // CGAL header providing mesh–mesh intersection predicates
    #include <CGAL/Polygon_mesh_processing/intersection.h>
    #include <boost/property_map/property_map.hpp>

    namespace PMP = CGAL::Polygon_mesh_processing;

#else
    #error "AMREX_SPACEDIM must be 2 or 3 for CGAL configuration !"
#endif

//============================================================================
// COMMON TYPE DEFINITIONS (Both 2D and 3D)
//============================================================================

// Kernel and basic types
using Real   = amrex::Real;
using Kernel = CGAL::Simple_cartesian<Real>;
using FT     = Kernel::FT;

#if (AMREX_SPACEDIM == 2)

    // -----------------------------------------------------------
    // 2D version
    // -----------------------------------------------------------

    // Basic geometric primitives in 2D
    using Point         = Kernel::Point_2;
    using Vector_CGAL   = Kernel::Vector_2;
    using Segment       = Kernel::Segment_2;
    using Polygon       = CGAL::Polygon_2<Kernel>;
    using GeomType      = Polygon;

    // AABB tree primitives (tree built over polygon edges)
    using Edge_iterator = Polygon::Edge_const_iterator;
    using Primitive     = CGAL::AABB_segment_primitive_2<Kernel, Edge_iterator>;
    using Traits        = CGAL::AABB_traits_2<Kernel, Primitive>;
    using Tree          = CGAL::AABB_tree<Traits>;

    // Query result type: (closest point, primitive ID)
    using Point_and_primitive_id = Tree::Point_and_primitive_id;

    // Face descriptor for 2D (edge iterator)
    using elm_descriptor = Edge_iterator;
    using PrimitiveID = Tree::Primitive_id;

    //----------------------------------------------------------------------------
    // Inside/Outside Tester for 2D Polygons
    //----------------------------------------------------------------------------
    /// \brief Functor for testing if a point is inside a 2D polygon.
    ///
    /// This class wraps CGAL::bounded_side_2 and provides a unified interface
    /// compatible with the 3D Side_of_triangle_mesh class.
    ///
    /// Requirements:
    ///   - Polygon must be counter-clockwise oriented
    ///   - Polygon must be simple (no self-intersections)
    ///
    /// Usage:
    ///   \code
    ///   Polygon poly = ...;
    ///   inside_t tester(poly);
    ///   Point query(1.0, 2.0);
    ///   if (tester(query) == CGAL::ON_BOUNDED_SIDE) {
    ///       // Point is inside polygon
    ///   }
    ///   \endcode
    class inside_t {
    private:
        const Polygon* poly = nullptr;

    public:
        /// Default constructor (creates uninitialized tester)
        inside_t() = default;

        /// Construct from a polygon reference
        /// \param p Reference to a CGAL Polygon_2 (must outlive this object)
        explicit inside_t(const Polygon& p) : poly(&p) {
            AMREX_ASSERT_WITH_MESSAGE(poly != nullptr, "inside_t: null polygon pointer");
            
            if (poly->is_clockwise_oriented()) {
                amrex::Abort("inside_t: Polygon must be counter-clockwise oriented");
            }
        }

        /// Check if this tester is properly initialized
        bool is_valid() const { 
            return poly != nullptr; 
        }

        /// Test if a point is inside, outside, or on the boundary
        /// \param p Query point
        /// \return CGAL::ON_BOUNDED_SIDE, CGAL::ON_BOUNDARY, or CGAL::ON_UNBOUNDED_SIDE
        CGAL::Bounded_side operator()(const Point& p) const {
            AMREX_ASSERT_WITH_MESSAGE(poly != nullptr, "inside_t: Cannot query uninitialized tester");
            
            return CGAL::bounded_side_2(
                poly->vertices_begin(),
                poly->vertices_end(),
                p,
                Kernel()
            );
        }
    };

#elif (AMREX_SPACEDIM == 3)

    // -----------------------------------------------------------
    // 3D version
    // -----------------------------------------------------------

    // Basic geometric primitives in 3D
    using Point       = Kernel::Point_3;
    using Vector_CGAL = Kernel::Vector_3;
    using Segment     = Kernel::Segment_3;
    using Polyhedron  = CGAL::Polyhedron_3<Kernel>;
    using GeomType    = Polyhedron;

    // AABB tree over triangular faces
    using Primitive   = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using Traits      = CGAL::AABB_traits_3<Kernel, Primitive>;
    using Tree        = CGAL::AABB_tree<Traits>;

    // Query result type: (closest point, primitive ID)
    using Point_and_primitive_id = Tree::Point_and_primitive_id;

    // Typedefs for face indexing and inside/outside classification
    using elm_descriptor  = boost::graph_traits<Polyhedron>::face_descriptor;
    using PrimitiveID     = Tree::Primitive_id;
    using inside_t        = CGAL::Side_of_triangle_mesh<Polyhedron, Kernel>;

#endif // AMREX_SPACEDIM


//============================================================================
// FUNCTION IMPLEMENTATIONS
//============================================================================

#if (AMREX_SPACEDIM == 2)

//----------------------------------------------------------------------------
// 2D Polygon File Reader (inline bool read_polygon_2d)
//----------------------------------------------------------------------------
/// \brief Reads a 2D polygon from a text file.
///        (Can be extended to other 2D mesh formats in the future.)
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
inline bool read_polygon_2d(const std::string& filename, Polygon& poly, Real dx = -1.0)
{
    poly.clear();
    
    std::ifstream in(filename);
    if (!in) {
        amrex::Print() << "read_polygon_2d: Cannot open file: " << filename << "\n";
        return false;
    }
    
    std::string line;
    int line_number = 0;
    
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
        
        poly.push_back(Point(x, y));
    }
    
    // Validate vertex count
    if (poly.size() < 3) {
        amrex::Print() << "read_polygon_2d: Polygon must have at least 3 vertices (found "
                       << poly.size() << ") in: " << filename << "\n";
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
        Polygon refined;
        refined.clear();
        refined.reserve(original_n);  // Reserve at least the original number of vertices

        const std::size_t n = poly.size();
        for (std::size_t i = 0; i < n; ++i) {
            const Point& p0 = poly[i];
            const Point& p1 = poly[(i + 1) % n];

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
        // Refinement only adds points along existing edges, so
        // simplicity (no self-intersections) is preserved.
        poly.swap(refined);
    }

    // Enforce counter-clockwise orientation for consistency.
    // CGAL (and many geometric algorithms) commonly expect CCW polygons.
    if (poly.is_clockwise_oriented()) {
        poly.reverse_orientation();
        amrex::Print() << "Info: Reversed polygon orientation to CCW in " 
                       << filename << "\n";
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

#endif

//============================================================================
// Local Coordinate Frame Structure
//============================================================================
/// \brief Represents a local orthonormal coordinate system on the IB surface.
///
/// This structure stores the basis vectors for the local frame at a specific
/// point (usually the centroid of a face or edge).
///
/// Basis vectors:
///   - normal:   Outward unit normal vector
///   - tangent1: First unit tangent vector
///   - tangent2: Second unit tangent vector (always present for alignment)
///
/// In 2D: The frame is {n, t1, 0}, where t1 is the edge tangent.
/// In 3D: The frame is {n, t1, t2}, forming a right-handed system.
struct LocalFrame {
    Real normal[AMREX_SPACEDIM];
    Real tangent1[AMREX_SPACEDIM];
#if (AMREX_SPACEDIM == 3)
    Real tangent2[AMREX_SPACEDIM];
#endif

    // Default constructor
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    LocalFrame() {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            normal[i] = 0.0;
            tangent1[i] = 0.0;
#if (AMREX_SPACEDIM == 3)
            tangent2[i] = 0.0;
#endif
        }
    }

    // Constructor
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
#if (AMREX_SPACEDIM == 3)
    LocalFrame(const Real n[3], const Real t1[3], const Real t2[3]) {
        for (int i = 0; i < 3; ++i) {
            normal[i] = n[i];
            tangent1[i] = t1[i];
            tangent2[i] = t2[i];
        }
    }
#else
    LocalFrame(const Real n[2], const Real t1[2]) {
        for (int i = 0; i < 2; ++i) {
            normal[i] = n[i];
            tangent1[i] = t1[i];
        }
    }
#endif
};

//============================================================================
// Surface Element Data Structure (Centroid & Size)
//============================================================================
/// \brief Stores geometric properties of a surface element (face in 3D, edge in 2D).
struct SurfElem {
    Real centroid[AMREX_SPACEDIM];
    Real size;        // Area in 3D, Length in 2D

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfElem() : size(0.0) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) centroid[i] = 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfElem(const Real* c, Real s) : size(s) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) centroid[i] = c[i];
    }
};

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
/// \param geom       The input geometry (Polygon_2 or Polyhedron_3).
/// \param surfelem   Output vector for surface element properties (appended to).
/// \param localframe Output vector for local frames (appended to).
/// \param idxmap     Output map from PrimitiveID to integer index (cleared before use).
/// \param offset     Global index offset for this geometry (default 0).
inline void build_geometry_cache(
    const GeomType& geom,
    amrex::Gpu::ManagedVector<SurfElem>& surfelem,
    amrex::Gpu::ManagedVector<LocalFrame>& localframe,
    std::map<PrimitiveID, int>& idxmap,
    int offset = 0)
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
        surfelem.push_back(SurfElem(c, len));

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
    // surfelem.reserve(surfelem.size() + n_elems);
    // localframe.reserve(localframe.size() + n_elems);

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
        surfelem.push_back(SurfElem(c, area));

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

//============================================================================
// GEOMETRY CONSISTENCY CHECKING
//============================================================================

/// \brief Check that IBM geometries do not intersect, touch, or contain each other.
///
/// This function validates that multiple IBM geometries are properly separated:
///
/// Level 1: Surface intersection/touching check
///   - In 2D: Checks if polygon edges intersect or touch
///   - In 3D: Checks if mesh surfaces intersect or touch
///
/// Level 2: Volumetric containment check
///   - Tests if one geometry is completely inside another
///   - Uses representative vertices from each geometry
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

            // 2) Check for containment (one polygon inside another)
            //    Use representative vertices from each polygon
            
            // Representative vertex from polygon i
            auto vi = geom_a[i].vertices_begin();
            const Point& pi = *vi;

            // Representative vertex from polygon j
            auto vj = geom_a[j].vertices_begin();
            const Point& pj = *vj;

            inside_t& inside_i = *inout_fa[i];
            inside_t& inside_j = *inout_fa[j];

            CGAL::Bounded_side side_i_in_j = inside_j(pi);
            CGAL::Bounded_side side_j_in_i = inside_i(pj);

            // Check if polygon i is strictly inside polygon j
            if (side_i_in_j == CGAL::ON_BOUNDED_SIDE) {
                amrex::Print() << "ERROR: IBM geometry containment detected:\n"
                              << "  geom " << i << " : " << files_a[i]
                              << " is strictly inside\n"
                              << "  geom " << j << " : " << files_a[j] << "\n";
                amrex::Abort(
                    "IBM geometries must not contain one another "
                    "(no full inclusion).");
            }

            // Check if polygon j is strictly inside polygon i
            if (side_j_in_i == CGAL::ON_BOUNDED_SIDE) {
                amrex::Print() << "ERROR: IBM geometry containment detected:\n"
                              << "  geom " << j << " : " << files_a[j]
                              << " is strictly inside\n"
                              << "  geom " << i << " : " << files_a[i] << "\n";
                amrex::Abort(
                    "IBM geometries must not contain one another "
                    "(no full inclusion).");
            }

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

            // 2) Containment check (pure inclusion without surface intersection)
            //    At this point surfaces do not intersect. Test whether a
            //    representative vertex of one closed mesh lies strictly inside
            //    the volume of the other.
            
            if (geom_a[i].empty() || geom_a[j].empty()) {
                continue;
            }

            // Representative vertex from geometry i
            auto vi = geom_a[i].vertices_begin();
            const Point& pi = vi->point();

            // Representative vertex from geometry j
            auto vj = geom_a[j].vertices_begin();
            const Point& pj = vj->point();

            inside_t& inside_i = *inout_fa[i];
            inside_t& inside_j = *inout_fa[j];

            CGAL::Bounded_side side_i_in_j = inside_j(pi);
            CGAL::Bounded_side side_j_in_i = inside_i(pj);

            // Check if mesh i is strictly inside mesh j
            if (side_i_in_j == CGAL::ON_BOUNDED_SIDE) {
                amrex::Print() << "ERROR: IBM geometry containment detected:\n"
                              << "  geom " << i << " : " << files_a[i]
                              << " is strictly inside\n"
                              << "  geom " << j << " : " << files_a[j] << "\n";
                amrex::Abort(
                    "IBM geometries must not contain one another "
                    "(no full inclusion).");
            }

            // Check if mesh j is strictly inside mesh i
            if (side_j_in_i == CGAL::ON_BOUNDED_SIDE) {
                amrex::Print() << "ERROR: IBM geometry containment detected:\n"
                              << "  geom " << j << " : " << files_a[j]
                              << " is strictly inside\n"
                              << "  geom " << i << " : " << files_a[i] << "\n";
                amrex::Abort(
                    "IBM geometries must not contain one another "
                    "(no full inclusion).");
            }
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
                       const CGAL::Bounded_side& result,
                       const Point& gridpoint)
{
    if (result == CGAL::ON_BOUNDARY) {
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

#endif  // EIB_CGAL_H_