#include <IBM.H>

using namespace amrex;
using namespace IBM;

IBFab::IBFab (const Box& b, int ncomp, bool alloc, bool shared, Arena* ar)    
              : BaseFab<bool>(b,ncomp,alloc,shared,ar) {}
IBFab::IBFab (const IBFab& rhs, MakeType make_type, int scomp, int ncomp) 
              : BaseFab<bool>(rhs,make_type,scomp,ncomp) {}
IBFab::~IBFab () { delete[] gpArray; }

void IBFab::allocateGPs(int numGPs) {
  delete[] gpArray;
  gpArray = new gp[numGPs]; 
  ngps = numGPs;
}

IBMultiFab::IBMultiFab ( const BoxArray& bxs, const DistributionMapping& dm, 
                        const int nvar, const int ngrow, const MFInfo& info, 
                        const FabFactory<IBFab>& factory )  :
                        FabArray<IBFab>(bxs,dm,nvar,ngrow,info,factory) {}
IBMultiFab::~IBMultiFab () {}


// constructor without geometry init
IB::IB (const Vector<BoxArray>& bxs, 
        const Vector<DistributionMapping>& dm, 
        const int nvar, const int nghost, const int max_level) {

        max_lev = max_level;
        // create IBMultiFabs at each level and store pointers to it
        mfa->resize(max_level);
        for (int lev=0; lev<=max_level; lev++) {
          IBMultiFab temp(bxs[lev],dm[lev],nvar,nghost);
          mfa->at(lev) = &temp;
        }


    // Compute_solidmarkers 

  }

IB::~IB () { 
  // https://stackoverflow.com/questions/6353149/does-vectorerase-on-a-vector-of-object-pointers-destroy-the-object-itself 
  // For a vector of pointers, we must delete each pointer to object individually.
 for (int lev=0; lev<=max_lev; lev++) {
          delete[] mfa->at(lev);
        }
}

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




 // 1. Read/Create and save geometry
  // void IBMinit() {
  //   // Create geometry   (TODO: or read geometry)
  //   sphere::Surface_mesh* sm;
  //   sm=sphere::create();

  //   const BoxArray& bxs = amrex::AmrLevel::get_new_data(CNS::State_Type).boxArray();
  //   // const DistributionMapping& dm = get_new_data(CNS::State_Type).DistributionMap();

  //   // IBMultiFab ibmf(amrex::AmrLevel::get_new_data(CNS::State_Type).boxArray(),
  //   //                 amrex::AmrLevel::get_new_data(State_Type).DistributionMap(),CNS::NUM_STATE,1);

  //   exit(0);


  //   // Compute_solidmarkers 
  //   test_point_in_body();

  // };



//   !initialise solid per lev
//   do lev=1,max_levs
//     ! create logical fabs sld(1) and ghs(2)
//     call lmultifab_build(solido(lev)%lfab,mla%la(lev),2,nhalo)
//     ! number of fabs
//     nf = nfabs(solido(lev)%lfab)
//     ! ghost points
//     ng_p = solido(lev)%lfab%ng

//     if (nf >= 1) then
//       ! allocate ibms per fab
//       allocate(solido(lev)%ibm(nf))

//       do ii = 1,nf
//         ! assign and initialise auxiliar real and logical pointers  
//         log_p            => dataptr(solido(lev)%lfab,ii)
//         log_p (:,:,:,:)  = .false.

//         ! array limits  
//         lo = lwb(get_box(solido(lev)%lfab,ii))
//         hi = upb(get_box(solido(lev)%lfab,ii))
//         call amr_nxnynz(lo,hi,ng_p)
        
//         ! loop over bodies in solid and update markers

// #if gts   
//       do nb = 1, max_bodies_lev
//           call compute_solidmarkergts(nb,log_p(:,:,:,1),dx(lev))
//       end do
// #else
//           call compute_solidmarker(log_p(:,:,:,1),dx(lev),ii,lev)
// #endif
//       end do
//     endif
//   enddo


// 2. compute and save solid and ghost point field


// 3. compute ghost point state


// 4. compute surface data