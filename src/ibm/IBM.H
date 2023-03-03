#ifndef IBM_H_
#define IBM_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_FabArray.H>
#include <AMReX_IntVect.H>
#include <AMReX_Derive.H>
#include <AMReX_StateDescriptor.H>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

using namespace amrex;

namespace IBM {

  // TODO: make it a linked list!
  // https://stackoverflow.com/questions/21482700/how-to-delete-this-linked-list-from-memory

  // GP holds the relevant data for each ghost point.
  struct gp {

      // attributes
      int idx[3]={0}; // GP index
      int norm[3]={0}; // closest element's normal
      // Immersed boundary point data
      /////////////////////////////////////////////

  };

  // IBFab holds the solid point and ghost point boolean arrays
  class IBFab: public BaseFab<bool> 
  {
    public:
      // attributes //
      gp* gpArray;
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
      // allocate space for gpArray
      void allocateGPs(int numGPs);
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

      // void copytoRealMF(IBMultiFab &ibmf, MultiFab &mf, int ibcomp, int mfcomp);
      void copytoRealMF(MultiFab &mf, int ibcomp, int mfcomp);
    protected:

  };

  // IB is the main class. It holds an array of IBMultiFab, one for each AMR level; and it also holds the geometry
  class IB {
    public:
      // attributes //

      int max_lev=0; // need to save this as Amr protects this

      // pointer to Amr class instance
      Amr*pamr;

      // vector of refinement ratio per level in each direction
      Vector<IntVect> ref_ratio;

      //vector of cell sizes per level in each direction
      Vector<GpuArray<Real,AMREX_SPACEDIM>> cellSizes;

      // IBMultiFab array
      Vector<IBMultiFab*> mfa[2]; 

      // geometry class
      typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
      typedef K::Point_3                                             Point;
      typedef CGAL::Polyhedron_3<K>                                  Polyhedron;

      Polyhedron geom;

      // constructors and destructors //
      // explicit IB ( const Vector<BoxArray>& bxs, 
      //               const Vector<DistributionMapping>& dm, 
      //               const int nvar, const int nghost, const int max_level);
      explicit IB ();
      ~IB ();
      

      // methods //
      void initialise (Amr* pointer_amr, const int nvar, const int nghost, const int max_level, const Vector<GpuArray<Real,AMREX_SPACEDIM>> dx);
      // compute IB field markers
      void compute_markers ();

      // read geometry
      void read_geom(const std::string filename);

  };

  // declare main IB class instance
  inline IB ib;

}
#endif