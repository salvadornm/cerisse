#ifndef IBM_H_
#define IBM_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_FabArray.H>
#include <AMReX_IntVect.H>
#include <AMReX_Derive.H>
#include <AMReX_StateDescriptor.H>

// basic CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>


// CGAL headers for AABB tree for closest point
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


// CGAL headers for AABB tree for surface data
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
// #include <CGAL/Unique_hash_map.h>
// #include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>

using namespace amrex;

namespace IBM {

  // GP holds the relevant data for each ghost point.
  struct gp {

      // attributes
      Array<int , AMREX_SPACEDIM> idx;
      int norm[3]={0}; // closest element's normal
      // Immersed boundary point data
      /////////////////////////////////////////////

  };

  // IBFab holds the solid point and ghost point boolean arrays
  class IBFab: public BaseFab<bool> 
  {
    public:
      // attributes //
      Vector<gp> gpArray;
      int ngps = 0;
      ///////////////

      // constructors and destructors //
      // using Box
      explicit IBFab (const Box& b, int ncomp, bool alloc=true, 
                      bool shared=false, Arena* ar = nullptr);
      // using IBFab
      explicit IBFab (const IBFab& rhs, MakeType make_type, int scomp, 
                      int ncomp);
     ~IBFab ();
     ///////////////////////////////////////////////////////////////////////////

      // methods //
      /////////////////////////////
  };

  // IBMultiFab holds an array of IBFab on a level
  class IBMultiFab: public FabArray<IBFab> {
    public:
      // constructor from BoxArray and DistributionMapping
      explicit IBMultiFab ( const BoxArray& bxs, 
                            const DistributionMapping& dm, 
                            const int nvar, 
                            const int ngrow, 
                            const MFInfo& info = MFInfo(), 
                            const FabFactory<IBFab>& factory = DefaultFabFactory<IBFab>());
      // destructor
      ~IBMultiFab ();

      void copytoRealMF(MultiFab &mf, int ibcomp, int mfcomp);
    protected:

  };

  // IB is the main class. It holds an array of IBMultiFab, one for each AMR level; and it also holds the geometry
  class IB {
    public:
      // attributes //

      int max_level=0; // need to save this as Amr protects this

      // pointer to Amr class instance
      Amr*pamr;

      // vector of refinement ratio per level in each direction
      Vector<IntVect> ref_ratio;

      // vector of cell sizes per level in each direction
      Vector<GpuArray<Real,AMREX_SPACEDIM>> cellSizes;

      // IBMultiFab array
      Vector<IBMultiFab*> mfa; 


      // CGAL types
      // typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
      typedef CGAL::Simple_cartesian<double> K2;
      typedef K2::Point_3 Point;
      typedef CGAL::Polyhedron_3<K2>                                  Polyhedron;
      Polyhedron geom; // IB geometry 

      // typedef CGAL::Simple_cartesian<double> K2;
      typedef K2::FT FT;
      typedef K2::Segment_3 Segment;
      typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
      typedef CGAL::AABB_traits<K2, Primitive> Traits;
      typedef CGAL::AABB_tree<Traits> Tree;
      typedef Tree::Point_and_primitive_id Point_and_primitive_id;
      Tree* treePtr;

      typedef K2::Vector_3 Vector_CGAL;
      // typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
      typedef boost::graph_traits<Polyhedron>::face_descriptor      face_descriptor;

      std::map<face_descriptor,Vector_CGAL> fnormals;
      // can use this map to create surfdata

      // constructors and destructors //
      // explicit IB ( const Vector<BoxArray>& bxs, 
      //               const Vector<DistributionMapping>& dm, 
      //               const int nvar, const int nghost, const int max_level);
      explicit IB ();
      ~IB ();
      

      // methods //
      void setMaxLevel(int max_lev);

      void initialise (Amr* pointer_amr, const int nvar, const int nghost);

      void closestElem(Point& query);

      void buildIBMultiFab (const BoxArray& bxa, const DistributionMapping& dm, int lev ,int nvar,int nghost);

      void destroyIBMultiFab (int lev);

      void computeMarkers (int lev);

      void initialiseGPs(int lev);

      void readGeom(const std::string filename);

  };

  // declare main IB class instance
  inline IB ib;

}
#endif