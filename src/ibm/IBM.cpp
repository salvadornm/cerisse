#include <IBM.H>
#include <CGAL/Side_of_triangle_mesh.h>

using namespace amrex;
using namespace IBM;

IBFab::IBFab (const Box& b, int ncomp, bool alloc, bool shared, Arena* ar)    
              : BaseFab<bool>(b,ncomp,alloc,shared,ar) {}
IBFab::IBFab (const IBFab& rhs, MakeType make_type, int scomp, int ncomp) 
              : BaseFab<bool>(rhs,make_type,scomp,ncomp) {}
IBFab::~IBFab () { }
// delete[] gpArray; TODO 04/03/2023

void IBFab::allocateGPs(int numGPs) {
  delete[] gpArray;
  gpArray = new gp[numGPs]; 
  ngps = numGPs; }

IBMultiFab::IBMultiFab ( const BoxArray& bxs, const DistributionMapping& dm, 
                        const int nvar, const int ngrow, const MFInfo& info, 
                        const FabFactory<IBFab>& factory )  :
                        FabArray<IBFab>(bxs,dm,nvar,ngrow,info,factory) {}
IBMultiFab::~IBMultiFab () {}

// for a single level
// void IBMultiFab::copytoRealMF(IBMultiFab &ibmf, MultiFab &mf, int ibcomp, int mfcomp) {
void IBMultiFab::copytoRealMF(MultiFab &mf, int ibcomp, int mfcomp) {

    for (MFIter mfi(*this,false); mfi.isValid(); ++mfi) {

        // const Box& ibbox = mfi.fabbox(); // box with ghost points
        const Box& ibbox = mfi.validbox(); // box with ghost points

        IBM::IBFab &ibfab = this->get(mfi);
        FArrayBox &realfab = mf.get(mfi);
        Array4<bool> ibMarkers = ibfab.array(); // boolean array
        Array4<Real> realfield = realfab.array(); // real array

        // assert that box sizes is same - TODO

        // const int *lo = fab.loVect();
        // const int *hi = fab.hiVect();
        // amrex::Print() << lo[0] << " " << lo[1] << " " << lo[2] << std::endl;
        // amrex::Print() << hi[0] << " " << hi[1] << " " << hi[2] << std::endl;

        amrex::ParallelFor(ibbox,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (ibMarkers(i,j,k,ibcomp)){
            realfield(i,j,k,mfcomp) = 1.0; }
          else {
            realfield(i,j,k,mfcomp) = 0.0; }

          if (ibMarkers(i,j,k,ibcomp+1)){
            realfield(i,j,k,mfcomp+1) = 1.0; }
          else {
            realfield(i,j,k,mfcomp+1) = 0.0; }

          });

        }

    }

// constructor and destructor
IB::IB (){}
IB::~IB () {}

void IB::setMaxLevel(int max_lev) {
  IB::max_level = max_lev;
  IB::mfa->resize(max_lev + 1);
  };

// initialise IB
void IB::initialise(Amr* pointer_amr, const int nvar, const int nghost) {

  // store pointer to main Amr class object's instance
  IB::pamr = pointer_amr ;

  ref_ratio = pamr->refRatio();
  max_level = pamr->maxLevel();
  // pamr->max_level() is a protected member of AmrInfo

  cellSizes.resize(IB::max_level+1);
  cellSizes[0] = pamr->Geom(0).CellSizeArray();
  for (int i=1;i<=IB::max_level;i++) {
  for (int j=0;j<=AMREX_SPACEDIM-1;j++) {
    cellSizes[i][j] = cellSizes[0][j]/ref_ratio[i-1][j];
  }}
}

// create IBMultiFabs at a level and store pointers to it
void IB::buildIBMultiFab (const BoxArray& bxa, const DistributionMapping& dm, int lev ,int nvar,int nghost) {
  mfa->at(lev) = new IBMultiFab(bxa,dm,nvar,nghost);
}

void IB::destroyIBMultiFab (int lev) {
  if (!mfa->empty()) {
      delete mfa->at(lev);
      // mfa->erase(lev); Do not use, shortens vector by one
  }
}


void IB::compute_markers (int lev) {

  CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(IB::geom);

  IBMultiFab *mfab = mfa->at(lev);
  int nghost = mfab->nGrow(0); // assuming same number of ghost points in all directions
  for (MFIter mfi(*mfab,false); mfi.isValid(); ++mfi) {
      IBM::IBFab &fab = mfab->get(mfi);
      const int *lo = fab.loVect();
      const int *hi = fab.hiVect();

      fab.setVal(false); // initialise sld and ghs to false
      Array4<bool> ibMarkers = fab.array(); // boolean array

      // compute sld markers (including at ghost points)
      for (int k = lo[2]; k <= hi[2]; ++k) {
      for (int j = lo[1]; j <= hi[1]; ++j) {
      for (int i = lo[0]; i <= hi[0]; ++i) {
          Real x=(0.5 + i)*cellSizes[lev][0];
          Real y=(0.5 + j)*cellSizes[lev][1];
          Real z=(0.5 + k)*cellSizes[lev][2];
          Point gridpoint(x,y,z);
          CGAL::Bounded_side res = inside(gridpoint);

          if (res == CGAL::ON_BOUNDED_SIDE) { 
              // soild marker
              ibMarkers(i,j,k,0) = true;}
          else if (res == CGAL::ON_BOUNDARY) { 
            amrex::Print() << "Grid point on IB surface" << " " << std::endl;
            exit(0);
          }
      }}};

      // compute ghs markers
      for (int k = lo[2]+nghost; k <= hi[2]-nghost; ++k) {
      for (int j = lo[1]+nghost; j <= hi[1]-nghost; ++j) {
      for (int i = lo[0]+nghost; i <= hi[0]-nghost; ++i) {

      bool fluidneigh = false;
      if (ibMarkers(i,j,k,0)) {
      for (int l = -1; l<=1; ++l) {
        fluidneigh = fluidneigh || (!ibMarkers(i+l,j,k,0));
        fluidneigh = fluidneigh || (!ibMarkers(i,j+l,k,0));
        fluidneigh = fluidneigh || (!ibMarkers(i,j,k+l,0));
      }
      ibMarkers(i,j,k,1) = fluidneigh; 

      //save ghs index

      };


      }}};
      // amrex::Print() << "fab ngps = " << fab.ngps << std::endl;
      // amrex::Print() << "------------------- " << std::endl;
  }
}


void IB::read_geom(const std::string filename) {
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(!PMP::IO::read_polygon_mesh(filename, IB::geom) || CGAL::is_empty(IB::geom) || !CGAL::is_triangle_mesh(IB::geom))
  {
    std::cerr << "Invalid input." << std::endl;
    exit(1);
  }
    amrex::Print() << "----------------------------------" << std::endl;
    amrex::Print() << "Geometry " << filename << " read"<< std::endl;
    amrex::Print() << "----------------------------------" << std::endl;
}


















    // // accessing derived variable list
    // DeriveList *derive = &(lev->get_derive_lst());
    // std::list<DeriveRec> list = derive->dlist();
    // amrex::Print() << list.front().name() << std::endl;
 
    // // accessing descriptor list
    // DescriptorList const *desList= &lev->get_desc_lst();
    // StateDescriptor const *state_var = &desList->operator[](0);
    // amrex::Print() << state_var->name(0) << std::endl;







///////////////////////////////CHECK IBMULTIFABS///////////////////
// int temp=0 ;
// for (MFIter mfi(*mfa[0],false); mfi.isValid(); ++mfi) // without tiling
// {
//     IBM::IBFab &fab = (*mfa[0])[mfi];
//     const int *lo = fab.loVect();
//     const int* hi = fab.hiVect();

//     amrex::Print() << "Level " << max_level << std::endl;
//     amrex::Print() << lo[0] << " " << lo[1] << " " << lo[2] << std::endl;
//     amrex::Print() << hi[0] << " " << hi[1] << " " << hi[2] << std::endl;

    
//     amrex::Print() << "fab allocated " << fab.isAllocated() << std::endl;
//     fab.allocateGPs(temp);
//     amrex::Print() << "fab ngps = " << fab.ngps << std::endl;

//     amrex::Print() << "------------------- " << std::endl;
//     temp += 1;
// }

///////////////////////////////////////////NORMAL MultiFab//////////////////////
  // MultiFab data( get_new_data(State_Type).boxArray(), get_new_data(State_Type).DistributionMap(), NUM_STATE, 1, MFInfo(), Factory());

  // const amrex::Real cur_time = state[State_Type].curTime();
  // FillPatch(*this, data, data.nGrow(), cur_time, State_Type, Density, NUM_STATE, 0);

  //   for (MFIter mfi(data,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  //   {
  //       FArrayBox& fab = data[mfi];
  //       const int *lo = fab.loVect();
  //       const int* hi = fab.hiVect();

  //       amrex::Print() << "Level " << level << std::endl;
  //       amrex::Print() << lo[0] << " " << lo[1] << " " << lo[2] << std::endl;
  //       amrex::Print() << hi[0] << " " << hi[1] << " " << hi[2] << std::endl;
        
  //   }
////////////////////////////////////////////////////////////////////////////////