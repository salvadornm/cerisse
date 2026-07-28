#include <CNS.h>
#include <prob.h>
#include <nscbc.h>

using namespace amrex;

//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::compute_nscbc_face_rhs( MultiFab& cell_state, bool second_order)
  {
    const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
    const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
    const CNS::NSCBCParm* nscbc_parm = CNS::d_nscbc_parm;

    const int nprim  = cls_h.NPRIM;
    const int nghost = cls_h.NGHOST;

    const Box domain = Geom().Domain();
    const auto dxinv = Geom().InvCellSizeArray();

    /*
     * Build cell-centred primitive variables from the RK-stage conservative
     * solution. These are used for q0 and q1 in the one-sided normal
     * derivative.
     */
    MultiFab cell_prims(
        cell_state.boxArray(),
        cell_state.DistributionMap(),
        nprim,
        nghost,
        MFInfo().SetArena(The_Async_Arena()));

    for (MFIter mfi(cell_state, false); mfi.isValid(); ++mfi) {
        cls_h.cons2prims( mfi, cell_state.array(mfi), cell_prims.array(mfi));
    }

    /*
     * Convert the current RK-stage face state UBC to QBC.
     */
    update_nscbc_face_primitives();

    /*
     * Clear all face RHS arrays before writing boundary planes.
     */
    clear_nscbc_face_rhs();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (!nscbc_qbc[dir] || !nscbc_rhs_bc[dir]) {
            continue;
        }

        MultiFab const& qbc_mf = *nscbc_qbc[dir];
        MultiFab& rhs_mf       = *nscbc_rhs_bc[dir];

        const int lo_face = domain.smallEnd(dir);
        const int hi_face = domain.bigEnd(dir) + 1;

        const IntVect edir = IntVect::TheDimensionVector(dir);

        /*
         * nscbc_qbc[dir], nscbc_rhs_bc[dir], and cell_prims have the same
         * number of boxes and DistributionMapping. The first two are
         * face-centred, while cell_prims is cell-centred.
         */
        for (MFIter mfi(qbc_mf, false); mfi.isValid(); ++mfi) {

          const Box& face_valid_box = mfi.validbox();

          auto const& qbc = qbc_mf.const_array(mfi);
          auto const& q = cell_prims.const_array(mfi);
          auto const& rhs_bc = rhs_mf.array(mfi);
          // ==============================================================
          // Low physical boundary
          // ==============================================================
          if (nscbc_lo[dir] > 0 && face_valid_box.smallEnd(dir) <= lo_face && face_valid_box.bigEnd(dir)   >= lo_face) {

            Box boundary_box = face_valid_box;

            boundary_box.setSmall(dir, lo_face);
            boundary_box.setBig(dir, lo_face);                

            if (boundary_box.ok()) {

              const int nscbc_type = nscbc_lo[dir];

              ParallelFor( boundary_box, [=] AMREX_GPU_DEVICE( int i, int j, int k) noexcept
                {
                  const IntVect iv_face( AMREX_D_DECL(i,j,k));
                  const IntVect iv_inner = iv_face;

                  nscbc::add_lodi_rhs_to_cons< PROB::ProbClosures>(
                      iv_face,iv_inner, dir, +1,dxinv,cls_d, qbc,q,
                          rhs_bc, second_order, nscbc_type, *nscbc_parm);
                });
            }
          }

          // ==============================================================
          // High physical boundary
          // ==============================================================
          if (nscbc_hi[dir] > 0 && face_valid_box.smallEnd(dir) <= hi_face && face_valid_box.bigEnd(dir)   >= hi_face) {

            Box boundary_box = face_valid_box;
            boundary_box.setSmall(dir, hi_face);
            boundary_box.setBig(dir, hi_face);
                
            if (boundary_box.ok()) {
              const int nscbc_type = nscbc_hi[dir];

              ParallelFor( boundary_box, [=] AMREX_GPU_DEVICE( int i, int j, int k) noexcept
                {
                  const IntVect iv_face(AMREX_D_DECL(i,j,k));
                  const IntVect iv_inner = iv_face - edir;

                  nscbc::add_lodi_rhs_to_cons< PROB::ProbClosures>(
                      iv_face,iv_inner, dir, -1,dxinv,cls_d, qbc,q,
                          rhs_bc, second_order, nscbc_type, *nscbc_parm);
                });
            }
          }
        }
    }
  }
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::overwrite_nscbc_inviscid_flux(MFIter const& mfi, Array4<const Real> const& cell_prims,
    std::array<FArrayBox*, AMREX_SPACEDIM> const& fluxes)
  {
    //if (!use_nscbc) return;

    const PROB::ProbClosures* cls_d = CNS::d_prob_closures;

    const Box domain = Geom().Domain();

    const Box cell_box = mfi.tilebox();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (!nscbc_qbc[dir]) continue;

        auto const& qbc = nscbc_qbc[dir]->const_array(mfi);

        auto const& flux = fluxes[dir]->array();

        const Box face_box = amrex::surroundingNodes(cell_box, dir);

        const int lo_face = domain.smallEnd(dir);

        const int hi_face = domain.bigEnd(dir) + 1;

        const IntVect edir = IntVect::TheDimensionVector(dir);

        // --------------------------------------------------------------
        // Low NSCBC face
        // --------------------------------------------------------------
        if (nscbc_lo[dir] > 0 && face_box.smallEnd(dir) <= lo_face && face_box.bigEnd(dir)   >= lo_face) {

          Box boundary_box = face_box;

          boundary_box.setSmall(dir, lo_face);
          boundary_box.setBig(dir, lo_face);

          ParallelFor( boundary_box, [=] AMREX_GPU_DEVICE( int i, int j, int k) noexcept
            {
              const IntVect iv_face( AMREX_D_DECL(i,j,k));
              const IntVect iv_inner = iv_face;

              nscbc::hllc_primitive_flux< PROB::ProbClosures>( iv_face, iv_inner, dir, +1, qbc,cell_prims, flux,*cls_d);
            });
        }
        // --------------------------------------------------------------
        // High NSCBC face
        // --------------------------------------------------------------
        if (nscbc_hi[dir] > 0 && face_box.smallEnd(dir) <= hi_face && face_box.bigEnd(dir)   >= hi_face) {

          Box boundary_box =face_box;

          boundary_box.setSmall(dir, hi_face);
          boundary_box.setBig(dir, hi_face);

          ParallelFor( boundary_box,[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
              const IntVect iv_face(AMREX_D_DECL(i,j,k));
              const IntVect iv_inner = iv_face - edir;

              nscbc::hllc_primitive_flux< PROB::ProbClosures>( iv_face, iv_inner, dir, -1, qbc, cell_prims, flux, *cls_d);
            });
        }
    }
  }

