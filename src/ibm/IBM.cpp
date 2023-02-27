#include <IBM.H>


using namespace amrex;
using namespace IBM;

// constructor
IBFab::IBFab (const Box& b, int ncomp, bool alloc, bool shared, Arena* ar)    
              : BaseFab<bool>(b,ncomp,alloc,shared,ar) {
              }

// move constructor
IBFab::IBFab (const IBFab& rhs, MakeType make_type, int scomp, int ncomp) 
              : BaseFab<bool>(rhs,make_type,scomp,ncomp) {
              }

// no copy constructor

// destructor
IBFab::~IBFab () {
  delete[] gpArray;
}

void IBFab::allocateGPs(int numGPs) {
  gpArray = new gp[numGPs]; 
  ngps = numGPs;
}


// constructor
IBMultiFab::IBMultiFab ( const BoxArray& bxs, const DistributionMapping& dm, 
                        int nvar, int ngrow, const MFInfo& info, 
                        const FabFactory<IBFab>& factory )  :
                        FabArray<IBFab>(bxs,dm,nvar,ngrow,info,factory) {}

// destructor
IBMultiFab::~IBMultiFab () {}


  // void IBMinit() {
  //   // Create geometry   (TODO: or read geometry)
  //   sphere::Surface_mesh* sm;
  //   sm=sphere::create();

  //   const BoxArray& bxs = amrex::AmrLevel::get_new_data(CNS::State_Type).boxArray();
  //   // const DistributionMapping& dm = get_new_data(CNS::State_Type).DistributionMap();

  //   // IBMultiFab ibmf(amrex::AmrLevel::get_new_data(CNS::State_Type).boxArray(),
  //   //                 amrex::AmrLevel::get_new_data(State_Type).DistributionMap(),CNS::NUM_STATE,0);

  //   exit(-1);


  //   // IBMultifab define
  //   // int ncomp = 3;
  //   // int ngrow = 0;
  //   // MultiFab mf(ba, dm, ncomp, ngrow);

  //   // Compute_solidmarkers 
  //   test_point_in_body();

  // };










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


  //   // IBMultifab define
  //   // int ncomp = 4;
  //   // int ngrow = 1;
  //   // MultiFab mf(ba, dm, ncomp, ngrow);

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