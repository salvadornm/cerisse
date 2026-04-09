#ifndef IBM_CGAL_DEFS_H_
#define IBM_CGAL_DEFS_H_

// ============================================================================
// ibm_cgal_defs.h — CGAL geometry type definitions for the IBM solver
//
// Contains all pure type definitions and data structures used by the
// CGAL geometry backend:
//   1. CGAL kernel and basic types : Kernel, FT, Point, Vector_CGAL, Segment
//   2. Geometry containers         : Polygon2D (2D) or Polyhedron (3D)
//   3. AABB tree types             : Primitive, Traits, Tree
//   4. Compatibility layer         : BoundedSide, ClosestPointResult
//   5. Inside/outside tester       : inside_t
//   6. Surface/frame types         : LocalFrame, SurfElem
//   7. ID mapping types            : PrimitiveIndexMap
// ============================================================================

// AMReX headers (must come first)
#include <AMReX.H>
#include <AMReX_Array.H>

// Standard library headers
#include <string>
#include <map>
#include <vector>
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
    #include <CGAL/Polygon_mesh_processing/bbox.h>
    #include <boost/property_map/property_map.hpp>

    namespace PMP = CGAL::Polygon_mesh_processing;

#else
    #error "AMREX_SPACEDIM must be 2 or 3 for CGAL configuration !"
#endif

// ============================================================================
// 1. CGAL KERNEL AND BASIC TYPES
// ============================================================================

using Real   = amrex::Real;
using Kernel = CGAL::Simple_cartesian<Real>;
using FT     = Kernel::FT;

// ============================================================================
// 2. DIMENSION-DEPENDENT TYPES
// ============================================================================

#if (AMREX_SPACEDIM == 2)

    // -----------------------------------------------------------
    // 2D version
    // -----------------------------------------------------------

    // Basic geometric primitives in 2D
    using Point         = Kernel::Point_2;
    using Vector_CGAL   = Kernel::Vector_2;
    using Segment       = Kernel::Segment_2;

    // Custom Polygon2D class wrapping CGAL::Polygon_2 and edge storage
    class Polygon2D {
    private:
        CGAL::Polygon_2<Kernel>  poly_;
        std::vector<Segment>     edges_;

    public:
        Polygon2D() = default;

        // Enable explicit copy semantics
        Polygon2D(const Polygon2D&) = default;
        Polygon2D& operator=(const Polygon2D&) = default;

        // Move semantics (for efficient poly = std::move(refined))
        Polygon2D(Polygon2D&&) = default;
        Polygon2D& operator=(Polygon2D&&) = default;

        // Vertex operations
        void clear()                     { poly_.clear(); edges_.clear(); }
        void push_back(const Point& p)   { poly_.push_back(p); }
        std::size_t size() const         { return poly_.size(); }
        const Point& vertex(std::size_t i) const { return poly_.vertex(i); }

        // Rebuild edges from vertices (must be called after adding all vertices)
        void finalize() {
            edges_.clear();
            edges_.reserve(poly_.size());
            for (auto it = poly_.edges_begin(); it != poly_.edges_end(); ++it) {
                edges_.push_back(*it);
            }
        }

        // Edge iterators (for AABB tree)
        auto edges_begin() const { return edges_.cbegin(); }
        auto edges_end()   const { return edges_.cend(); }

        // Vertex iterators (for check_ibm_geometry_consistency)
        auto vertices_begin() const { return poly_.vertices_begin(); }
        auto vertices_end()   const { return poly_.vertices_end(); }

        // Properties
        bool is_simple() const             { return poly_.is_simple(); }
        bool is_clockwise_oriented() const { return poly_.is_clockwise_oriented(); }
        void reverse_orientation()         { poly_.reverse_orientation(); finalize(); }
        Real area() const                  { return poly_.area(); }

        // Inside/outside test
        CGAL::Bounded_side bounded_side(const Point& p) const {
            return poly_.bounded_side(p);
        }

        // Bounding box access
        CGAL::Bbox_2 bbox() const { return poly_.bbox(); }
    };

    using Polygon   = Polygon2D;
    using GeomType  = Polygon2D;
    using SegmentIterator = std::vector<Segment>::const_iterator;

    // AABB tree primitives (tree built over polygon edges)
    using Primitive        = CGAL::AABB_segment_primitive_2<Kernel, SegmentIterator>;
    using Traits           = CGAL::AABB_traits_2<Kernel, Primitive>;
    using Tree             = CGAL::AABB_tree<Traits>;

    // Query result type: (closest point, primitive ID)
    using Point_and_primitive_id = Tree::Point_and_primitive_id;

    // Face descriptor for 2D (SegmentIterator)
    using elm_descriptor = SegmentIterator;
    using PrimitiveID    = Tree::Primitive_id;

    // Bounding box type in 2D
    using Bbox           = CGAL::Bbox_2;

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
    using inside_cgal_t   = CGAL::Side_of_triangle_mesh<Polyhedron, Kernel>;

    // Bounding box type in 3D
    using Bbox            = CGAL::Bbox_3;

#endif // AMREX_SPACEDIM

// ============================================================================
// 3. COMPATIBILITY LAYER — Unified BoundedSide enum matching BVH interface
// ============================================================================

enum class BoundedSide : int {
    Inside     =  1,   // CGAL::ON_BOUNDED_SIDE
    OnBoundary =  0,   // CGAL::ON_BOUNDARY
    Outside    = -1    // CGAL::ON_UNBOUNDED_SIDE
};

AMREX_FORCE_INLINE
BoundedSide cgal_to_bounded_side(CGAL::Bounded_side bs)
{
    if (bs == CGAL::ON_BOUNDED_SIDE)  return BoundedSide::Inside;
    if (bs == CGAL::ON_BOUNDARY)      return BoundedSide::OnBoundary;
    return BoundedSide::Outside;
}

// ============================================================================
// 4. INSIDE/OUTSIDE TESTER
// ============================================================================

#if (AMREX_SPACEDIM == 2)

    /// \brief Functor for testing if a point is inside a 2D polygon.
    ///
    /// This class wraps Polygon2D::bounded_side and provides a unified interface
    /// compatible with the 3D Side_of_triangle_mesh class.
    ///
    /// Requirements:
    ///   - Polygon must be counter-clockwise oriented
    ///   - Polygon must be simple (no self-intersections)
    class inside_t {
    private:
        const Polygon2D* poly_ = nullptr;

    public:
        /// Default constructor (creates uninitialized tester)
        inside_t() = default;

        /// Construct from a Polygon2D reference
        /// \param p Reference to a Polygon2D (must outlive this object)
        explicit inside_t(const Polygon2D& p) : poly_(&p) {
            AMREX_ASSERT_WITH_MESSAGE(poly_ != nullptr, "inside_t: null polygon pointer");

            if (poly_->is_clockwise_oriented()) {
                amrex::Abort("inside_t: Polygon must be counter-clockwise oriented");
            }
        }

        /// Check if this tester is properly initialized
        bool is_valid() const {
            return poly_ != nullptr;
        }

        /// Test if a point is inside, outside, or on the boundary
        BoundedSide operator()(const Point& p) const {
            AMREX_ASSERT_WITH_MESSAGE(poly_ != nullptr, "inside_t: Cannot query uninitialized tester");
            return cgal_to_bounded_side(poly_->bounded_side(p));
        }
    };

#elif (AMREX_SPACEDIM == 3)

    /// \brief Inside/Outside Tester wrapper for 3D meshes — returns BoundedSide
    class inside_t {
    private:
        inside_cgal_t impl_;
    public:
        inside_t() = delete;
        explicit inside_t(const Polyhedron& mesh)
            : impl_(mesh) {}

        BoundedSide operator()(const Point& p) const {
            return cgal_to_bounded_side(impl_(p));
        }
    };

#endif // AMREX_SPACEDIM

// ============================================================================
// 5. CLOSEST POINT RESULT (compatibility with BVH interface)
// ============================================================================

struct ClosestPointResult {
    Point point;
    int   prim_id;
};

// ============================================================================
// 6. SURFACE / FRAME TYPES
// ============================================================================

/// \brief Represents a local orthonormal coordinate system on the IB surface.
///
/// Basis vectors:
///   - normal:   Outward unit normal vector
///   - tangent1: First unit tangent vector
///   - tangent2: Second unit tangent vector (3D only)
///
/// In 2D: The frame is {n, t1}, where t1 is the edge tangent.
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

/// \brief Stores geometric properties of a surface element (face in 3D, edge in 2D).
struct SurfElem {
    int geomIdx;         // Standard int (4 bytes) to avoid padding issues
    Real centroid[AMREX_SPACEDIM];
    Real measure;        // Area in 3D, Length in 2D

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfElem() : geomIdx(-1), measure(0.0) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) centroid[i] = 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfElem(const Real* c, Real s, int g_idx) : geomIdx(g_idx), measure(s) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) centroid[i] = c[i];
    }
};

// ============================================================================
// 7. ID MAPPING TYPES
// ============================================================================

/// \brief Map type for associating geometric primitives with integer indices.
using PrimitiveIndexMap = std::map<PrimitiveID, int>;

#endif // IBM_CGAL_DEFS_H_
