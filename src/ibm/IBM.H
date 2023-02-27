#ifndef IBM_H_
#define IBM_H_

#include <AMReX_FabArray.H>
#include <IBgeom.H>
#include <cgal_test.H>

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
    protected:

  };

  // IB is the main class. It holds an array of IBMultiFab, one for each AMR level; and it also holds the geometry
  class IB {
    public:
      // attributes //
      int max_lev=0;
      // IBMultiFab array
      Vector<IBMultiFab*> mfa[2]; 
      // geometry class

      // constructors and destructors //
      explicit IB ( const Vector<BoxArray>& bxs, 
                    const Vector<DistributionMapping>& dm, 
                    const int nvar, const int nghost, const int max_level);
      ~IB ();

      // methods
      void compute_markers ();

      //compute

  };

}
#endif