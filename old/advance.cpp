#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <CNS_hydro_K.h>
#include <prob.h>
#ifdef AMREX_USE_GPIBM
#include <IBM.h>
#endif
using namespace amrex;

// Since we do not want to use expensive cudaMemCopy, we are storing all our
// data on the GPU to begin with. Concurrency on GPU using streams, parallel
// computation and data transfer, is not useful then. Therefore, we can have all
// grid point computations, per fab, in a single MFIter loop (single stream).

//////////////////////////////////////////////////////////////////////////////

//   // Euler flux corrections (overwrite numflxmf) //
//   // Recompute fluxes on planes adjacent to physical boundaries (Order
//   reduction) if (flux_euler==1 && !(CentralKEEP::order_keep==2)) {
//     CentralKEEP::Flux_2nd_Order_KEEP(geom,primsmf,numflxmf);
//   }
//   // Order reduction near IBM

//   // Artificial dissipation (adding to numflxmf)
//   // JST artificial dissipation shock capturing
//   if (art_diss==1) {
//   // make multifab for spectral radius and sensor for artificial dissipation
//     MultiFab lambdamf; lambdamf.define(consmf.boxArray(),
//     consmf.DistributionMap(), AMREX_SPACEDIM, NGHOST); MultiFab sensormf;
//     sensormf.define(consmf.boxArray(), consmf.DistributionMap(),
//     AMREX_SPACEDIM, NGHOST); lambdamf = 0.0_rt; sensormf = 0.0_rt; for
//     (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//       const Box& bx      = mfi.tilebox();
//       const Box& bxnodal = mfi.grownnodaltilebox(-1,0);

//       auto const& statefab = consmf.array(mfi);
//       auto const& sensor   = sensormf.array(mfi);
//       auto const& lambda   = lambdamf.array(mfi);
//       auto const& prims    = primsmf.array(mfi);
//       AMREX_D_TERM(auto const& nfabfx = numflxmf[0].array(mfi);,
//                     auto const& nfabfy = numflxmf[1].array(mfi);,
//                     auto const& nfabfz = numflxmf[2].array(mfi););
//       amrex::ParallelFor(bx,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//           ComputeSensorLambda(i,j,k,primsfab,lambda,sensor,cls);
//         });

//         amrex::ParallelFor(bxnodal, ncons,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//         {
//           JSTflux(i,j,k,n,lambda,sensor,statefab,nfabfx,nfabfy,nfabfz,cls);
//         });
//     }
//   }
// }

//////////////////////////////////////////////////////////////////////////////

// Viscous Fluxes ////////////////////////////////////////////////////////////
// Gpu::streamSynchronize(); // ensure all rhs terms computed before assembly
// We have a separate MFIter loop here than the Euler fluxes and the source
// terms, so the work can be further parallised. As different MFIter loops can
// be in different GPU streams.

// Although conservative FD (finite difference) derivatives of viscous fluxes
// are not requried in the boundary layer, standard FD are likely sufficient.
// However, considering grid and flow discontinuities (coarse-interface
// flux-refluxing and viscous derivatives near shocks), conservative FD
// derivatives are preferred.
//   if (rhs_visc) {
//     Array<MultiFab,AMREX_SPACEDIM>& pntvflxmf = Vpntvflxmf[level];
// #if AMREX_USE_GPIBM
//     IBM::IBMultiFab& ibMultiFab = *IBM::ib.ibMFa[level];
// #endif
//     //for each fab
//     for (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//       const Box& bxpflx  = mfi.growntilebox(1);
//       auto const& prims  = primsmf.array(mfi);

//       AMREX_D_TERM(auto const& pfabfx = pntvflxmf[0].array(mfi);,
//                   auto const& pfabfy  = pntvflxmf[1].array(mfi);,
//                   auto const& pfabfz  = pntvflxmf[2].array(mfi););

//       AMREX_D_TERM(auto const& nfabfx = numflxmf[0].array(mfi);,
//                   auto const& nfabfy  = numflxmf[1].array(mfi);,
//                   auto const& nfabfz  = numflxmf[2].array(mfi););

//       // compute u,v,w,T derivatives and compute physical viscous fluxes
//       amrex::ParallelFor(bxpflx,[=] AMREX_GPU_DEVICE (int i, int j, int k)
//       noexcept {
//         ViscousFluxes(i, j, k, primsfab, pfabfx, pfabfy, pfabfz, dxinv, cls);
//       });

//       // Physical boundary viscous flux corrections
//       // (overwrite pfabfx, pfabfy, pfabfz)
//       // TODO:: generalise to wall boundary in x and z directions
//       const Box& bx  = mfi.tilebox();
//       //xlo
//       //xhi

//       //ylo
//       if(geom.Domain().smallEnd(1)==bx.smallEnd(1)) {
//         if ((*h_phys_bc).lo(1)==6) {
//           int jj = bx.smallEnd(1)-1;
//           IntVect small = {bxpflx.smallEnd(0), jj, bxpflx.smallEnd(2)};
//           IntVect big   = {bxpflx.bigEnd(0)  , jj, bxpflx.bigEnd(2)  };
//           Box bxboundary(small,big);

//           amrex::ParallelFor(bxboundary,[=] AMREX_GPU_DEVICE (int i, int j,
//           int k) noexcept {
//               ViscousWallFluxes(i, j, k, 0, primsfab, pfabfx, pfabfy, pfabfz,
//               dxinv, cls);
//           });
//         }
//       }
//       //yhi
//       if(geom.Domain().bigEnd(1)==bx.bigEnd(1)) {
//         if ((*h_phys_bc).hi(1)==6) {
//           int jj = bx.bigEnd(1) + 1;
//           IntVect small = {bxpflx.smallEnd(0), jj, bxpflx.smallEnd(2)};
//           IntVect big   = {bxpflx.bigEnd(0)  , jj, bxpflx.bigEnd(2)  };
//           Box bxboundary(small,big);

//           amrex::ParallelFor(bxboundary,[=] AMREX_GPU_DEVICE (int i, int j,
//           int k) noexcept {
//               ViscousWallFluxes(i, j, k, 1, primsfab, pfabfx, pfabfy, pfabfz,
//               dxinv, cls);
//           });
//         }
//       }
//       //zlo
//       //zhi

// #if AMREX_USE_GPIBM
//       //IBM GP viscous flux correction
//       auto const& ibFab = ibMultiFab.get(mfi);
//       auto const& markers = ibMultiFab.array(mfi);
//       auto const gp_ijk   = ibFab.gpData.gp_ijk.data();
//       amrex::ParallelFor(ibFab.gpData.ngps, [=] AMREX_GPU_DEVICE (int ii)
//       {
//         ViscousFluxGP(gp_ijk[ii](0),gp_ijk[ii](1),gp_ijk[ii](2),markers,primsfab,pfabfx,pfabfy,pfabfz,dxinv,cls);
//       });
// #endif

//       // compute numerical viscous fluxes (add to numflxmf)
//       const Box& bxnodal  = mfi.grownnodaltilebox(-1,0); // extent is 0,N+1
//       in all directions -- -1 means for all directions.
//       amrex::surroundingNodes(bx) does the same amrex::ParallelFor(bxnodal,
//       ncons,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
//         ViscousNumericalFluxes(i, j, k, n, pfabfx, pfabfy, pfabfz, nfabfx,
//         nfabfy, nfabfz);
//       });
//     }
//   }

// Re-fluxing ////////////////////////////////////////////////////////////////
// if constexpr (do_reflux) {
// if (fr_as_crse) {
//     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
//         const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ?
//         dx[0]*dx[2] : dx[0]*dx[1]); const Real scale = -dt*dA;
//         fr_as_crse->CrseInit(numflxmf[idim], idim, 0, 0, ncons, scale,
//         FluxRegister::ADD);
//     }
// }
// if (fr_as_fine) {
//     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
//         const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ?
//         dx[0]*dx[2] : dx[0]*dx[1]); const Real scale = dt*dA;
//         fr_as_fine->FineAdd(numflxmf[idim], idim, 0, 0, ncons, scale);
//     }
// }
// }
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// Add source term to RHS ////////////////////////////////////////////////////
// if (rhs_source) {
//   for (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi){
//     const Box& bx   = mfi.tilebox();
//     auto const& dsdtfab = dSdt.array(mfi);
//     auto const& statefab = consmf.array(mfi);

//     amrex::ParallelFor(bx,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//     { user_source(i,j,k,statefab,dsdtfab,lprobparm,cls,dx); });
//   }
// }
// //////////////////////////////////////////////////////////////////////////////

// Set solid point RHS to 0 //////////////////////////////////////////////////
// #if AMREX_USE_GPIBM
//   IBM::IBMultiFab& mfab = *IBM::ib.ibMFa[level];
//   for (MFIter mfi(consmf, TilingIfNotGPU()); mfi.isValid(); ++mfi){
//     const Box& bx   = mfi.tilebox();
//     auto const& dsdtfab = dSdt.array(mfi);
//     IBM::IBFab &fab = mfab.get(mfi);
//     Array4<bool> ibMarkers = fab.array();
//     amrex::ParallelFor(bx, ncons,
//     [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//     {
//     dsdtfab(i,j,k,n) = dsdtfab(i,j,k,n)*(1 - int(ibMarkers(i,j,k,0)));
//     });
//   }
// #endif
// }
//////////////////////////////////////////////////////////////////////////////
