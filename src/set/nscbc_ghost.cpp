// nscbc_ghost.cpp
#include <CNS.h>
#include <prob.h>
#include <nscbc.h>

using namespace amrex;

void CNS::compute_nscbc_ghost_rhs(amrex::MultiFab& S,
                                  amrex::MultiFab& G_rhs,
                                  amrex::Real /*dt*/)
{
  const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
  const CNS::NSCBCParm* nscbc_parm = CNS::d_nscbc_parm;

  const int ncons  = cls_h.NCONS;
  const int nghost = cls_h.NGHOST;

  const Box dom  = Geom().Domain();
  const Box gdom = amrex::grow(dom, nghost);
  const auto dxinv = Geom().InvCellSizeArray();

  G_rhs.setVal(Real(0.0), 0, ncons, nghost);

  for (MFIter mfi(S, false); mfi.isValid(); ++mfi) {
    auto const& state = S.array(mfi);
    auto const& rhs   = G_rhs.array(mfi);

    const Box gbx = mfi.growntilebox(nghost);

    FArrayBox primfab(gbx, cls_h.NPRIM, The_Async_Arena());
    auto const& q = primfab.array();

    cls_h.cons2prims(mfi, state, q);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

      if (nscbc_lo[dir] > 0) {
        Box slab = gdom;
        slab.setSmall(dir, dom.smallEnd(dir) - nghost);
        slab.setBig  (dir, dom.smallEnd(dir) - 1);

        Box b = gbx & slab;

        // cut corners  (really)
        for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
          if (tdir != dir) {
            b.grow(tdir, -1);
          }
        }


        if (b.ok()) {
          const int nscbc_type = nscbc_lo[dir];

          ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const IntVect iv(AMREX_D_DECL(i,j,k));

            nscbc::add_lodi_rhs_to_cons<PROB::ProbClosures>(
                iv,
                dir,
                +1,              // derivative points inward from low ghost side
                dxinv,
                cls_d,
                q,
                rhs,
                true,
                nscbc_type,
                *nscbc_parm);
          });
        }
      }

      if (nscbc_hi[dir] > 0) {
        Box slab = gdom;
        slab.setSmall(dir, dom.bigEnd(dir) + 1);
        slab.setBig  (dir, dom.bigEnd(dir) + nghost);

        Box b = gbx & slab;

        // cut corners  (really)
        for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
          if (tdir != dir) {
            b.grow(tdir, -1);
          }
        }


        if (b.ok()) {
          const int nscbc_type = nscbc_hi[dir];

          ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const IntVect iv(AMREX_D_DECL(i,j,k));

            nscbc::add_lodi_rhs_to_cons<PROB::ProbClosures>(
                iv,
                dir,
                -1,              // derivative points inward from high ghost side
                dxinv,
                cls_d,
                q,
                rhs,
                true,
                nscbc_type,
                *nscbc_parm);
          });
        }
      }
    }
  }
}

/// @brief  Initialize the NSCBC ghost state from the current solution state at a given time. 
// This is typically called after the initial conditions are set, and before the first time step is taken, 
// to ensure that the NSCBC ghost cells have consistent values with the interior solution.
/// @param time 
void CNS::init_nscbc_ghost_state(amrex::Real time)
{
  if (!use_nscbc || !nscbc_ghost_state) return;

  const int ncons  = d_prob_closures->NCONS;
  const int nghost = d_prob_closures->NGHOST;

  MultiFab Sinit(grids, dmap, ncons, nghost, MFInfo(), Factory());

  FillPatch(*this, Sinit, nghost, time, State_Type, 0, ncons);

  MultiFab::Copy(*nscbc_ghost_state, Sinit, 0, 0, ncons, nghost);
}