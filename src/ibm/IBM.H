#ifndef IBM_H_
#define IBM_H_

#include <AMReX_FabArray.H>
#include <geometry.H>
#include <cgal_test.H>

using namespace amrex;

namespace IBM {

  // TODO: make it a linked list!
  // https://stackoverflow.com/questions/21482700/how-to-delete-this-linked-list-from-memory
  struct gp {

      // attributes
      int idx[3]={0}; // GP index
      int norm[3]={0}; // closest element's normal

      // Need constructors?
      // Need a destructor?
  };


  class IBFab: public BaseFab<bool> 
  {
    public:
      // attributes
      gp* gpArray;
      int ngps = 0;

      // constructors
      explicit IBFab (const Box& b, int ncomp, bool alloc=true, 
                      bool shared=false, Arena* ar = nullptr);

      explicit IBFab (const IBFab& rhs, MakeType make_type, int scomp, 
                      int ncomp);
      // explicit specifier means input arguments are not converted/typecasted

      // destructor
     ~IBFab ();


      // void addGP(int val);
      void allocateGPs(int numGPs);
  };


  // Technically FabArray only takes BaseFab type
  class IBMultiFab: public FabArray<IBFab> 
  {
  // One per simulation - an array of n levels. Each level holds a FabArray with IB data for all local fabs. The class also stores IB geometry.
    public:

      explicit IBMultiFab ( const BoxArray& bxs, 
                            const DistributionMapping& dm, 
                            int nvar, 
                            int ngrow, 
                            const MFInfo& info = MFInfo(), 
                            const FabFactory<IBFab>& factory = DefaultFabFactory<IBFab>());
        // amrex::Print() << "IBMultiFab constructed" << std::endl;

      // destructor
      ~IBMultiFab ();

    protected:

  };



}
#endif