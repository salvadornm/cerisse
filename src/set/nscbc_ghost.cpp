#include <CNS.h>
#include <prob.h>
#include <nscbc.h>

using namespace amrex;

//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void decode_nscbc_owner(int owner,int& dir, int& side_sign) noexcept
{
    // owner: // 1,3,5 -> low side. // 2,4,6 -> high side
    dir = (owner - 1) / 2;
    side_sign = (owner % 2 == 1) ? +1 : -1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int number_of_external_directions(amrex::IntVect const& iv, amrex::IntVect const& dom_lo, amrex::IntVect const& dom_hi) noexcept
{
    int count = 0;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (iv[dir] < dom_lo[dir] || iv[dir] > dom_hi[dir]) {
            ++count;
        }
    }

    return count;
}

} // namespace

//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::compute_nscbc_ghostcell_rhs( amrex::MultiFab& stage_state)
{
    const auto* cls_d = d_prob_closures;
    const auto& cls_h = *h_prob_closures;
    
    // NEW
    const auto nslo = CNS::nscbc_lo;
    const auto nshi = CNS::nscbc_hi;
    const auto parm_lo = CNS::nscbc_parm_lo;
    const auto parm_hi = CNS::nscbc_parm_hi;

    const int nprim  = cls_h.NPRIM;
    const int ng     = cls_h.NGHOST;

    /*
     * Fill ordinary inter-FAB ghosts and then install the current persistent
     * NSCBC state in every local physical-boundary slab.
     */
    fill_nscbc_ghost_cells(stage_state);

    MultiFab q( stage_state.boxArray(), stage_state.DistributionMap(), nprim, ng, MFInfo().SetArena(The_Async_Arena()), Factory());
    q.setVal(Real(0.0)); // set to 0

    for (MFIter mfi(stage_state, false); mfi.isValid(); ++mfi) {
        cls_h.cons2prims(mfi,stage_state.array(mfi),q.array(mfi));
    }

    q.FillBoundary(Geom().periodicity());
    clear_nscbc_ghost_rhs();

    const auto dxinv = Geom().InvCellSizeArray();

    const Box domain = Geom().Domain();
    const auto dom_lo = domain.smallEnd();
    const auto dom_hi = domain.bigEnd();

    for (MFIter mfi(q, false); mfi.isValid(); ++mfi) {

        auto const& qp   = q.const_array(mfi);
        auto const& rhs  = nscbc_shell.rhs->array(mfi);
        auto const& own  = nscbc_shell.owner->const_array(mfi);

        const Box bx = mfi.fabbox();

        const IntVect fab_lo = bx.smallEnd();
        const IntVect fab_hi = bx.bigEnd();

        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k) noexcept
            {
                const IntVect iv(AMREX_D_DECL(i,j,k));
                const int face = own(iv,0);

                if (face == 0) {return;}

                int dir;
                int side_sign;

                decode_nscbc_owner(face, dir, side_sign);

                const int nscbc_type = (side_sign > 0) ? nslo[dir] : nshi[dir];

                const auto& parm = (side_sign > 0) ? parm_lo[dir] : parm_hi[dir];

                const int external_dirs = number_of_external_directions( iv,dom_lo,dom_hi);
                
                //const bool face_interior = (external_dirs == 1);

                bool transverse_stencil_available = (external_dirs == 1);

                // check if transverse terms possible
                for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
                    if (tdir == dir) {continue;}

                    if (iv[tdir] <= fab_lo[tdir] || iv[tdir] >= fab_hi[tdir]) {
                    transverse_stencil_available = false;
                    }
                }

                // transverse_stencil_available = false; //temp                
                // temp: only activates transverse terms in interior points
                // palce holder for corners and edges

                nscbc::add_lodi_ghost_rhs_to_cons <PROB::ProbClosures>(
                        iv,dir,side_sign,dxinv,
                        cls_d,qp,rhs,nscbc_order,nscbc_type, parm, transverse_stencil_available); //NEW
            });
    }
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Build the ownership mask for the persistent NSCBC ghost shell.
//
// owner = 0 : not an NSCBC physical ghost cell
// owner = 1 : x-low
// owner = 2 : x-high
// owner = 3 : y-low
// owner = 4 : y-high
// owner = 5 : z-low
// owner = 6 : z-high
//
// At edges and corners, a cell can lie outside the domain in more than one
// direction. Until the coupled Lodato edge/corner treatment is implemented,
// ownership is assigned using the deterministic priority
//
//     x > y > z.
//
// Thus, the first active external direction owns the cell.
//------------------------------------------------------------------------------
void CNS::build_nscbc_shell_owner()
{
    using namespace amrex;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( nscbc_shell.owner != nullptr, "CNS::build_nscbc_shell_owner(): owner iMultiFab is not defined");

    iMultiFab& owner_mf = *nscbc_shell.owner;

    owner_mf.setVal(0);

    const Box domain = Geom().Domain();

    const auto dom_lo = domain.smallEnd();
    const auto dom_hi = domain.bigEnd();

    const auto nslo = CNS::nscbc_lo;
    const auto nshi = CNS::nscbc_hi;

    for (MFIter mfi(owner_mf, false); mfi.isValid(); ++mfi) {

        auto const& owner = owner_mf.array(mfi);

        // Includes the valid region and all allocated ghost cells.
        const Box bx = mfi.fabbox();

        ParallelFor(
            bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const IntVect iv(AMREX_D_DECL(i, j, k));

                int shell_owner = 0;

                // Deterministic priority:
                // x-low, x-high, y-low, y-high, z-low, z-high.
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

                    if (iv[dir] < dom_lo[dir] && nslo[dir] > 0) {
                        shell_owner = 2*dir + 1;
                        break;
                    }

                    if (iv[dir] > dom_hi[dir] && nshi[dir] > 0) {
                        shell_owner = 2*dir + 2;
                        break;
                    }
                }

                owner(iv, 0) = shell_owner;
            });
    }

    Gpu::streamSynchronize();
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int number_of_external_directions(amrex::IntVect const& iv, amrex::IntVect const& dom_lo, amrex::IntVect const& dom_hi)
{
    int count = 0;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (iv[dir] < dom_lo[dir] || iv[dir] > dom_hi[dir]) {
            ++count;
        }
    }

    return count;
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::fill_nscbc_ghost_cells2(MultiFab& state) const //old version
{
    AMREX_ALWAYS_ASSERT(nscbc_shell.defined());

    state.FillBoundary(Geom().periodicity());

    MultiFab const& shell = *nscbc_shell.state;
    iMultiFab const& owner = *nscbc_shell.owner;

    const int ncons = PROB::ProbClosures::NCONS;

    for (MFIter mfi(state, false); mfi.isValid(); ++mfi) {

        auto const& S   = state.array(mfi);
        auto const& G   = shell.const_array(mfi);
        auto const& own = owner.const_array(mfi);

        const Box bx = mfi.fabbox();

        ParallelFor(bx, ncons, [=] AMREX_GPU_DEVICE( int i, int j, int k, int n) noexcept
            {
                if (own(i,j,k,0) != 0) { S(i,j,k,n) = G(i,j,k,n);}
            });
    }
}
//------------------------------------------------------------------------------
// Overwrite NSCBC physical ghost cells with the persistent ghost-shell state.
//
// Ordinary inter-box and periodic ghost cells are synchronized first.
// Only slabs outside active NSCBC physical boundaries are then copied.
//
// Edge and corner cells may belong to more than one slab, but every copy reads
// the same cell-centred value from nscbc_shell.state, so repeated writes are
// harmless.
//------------------------------------------------------------------------------
void CNS::fill_nscbc_ghost_cells(amrex::MultiFab& state) const
{
    using namespace amrex;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nscbc_shell.state != nullptr, "CNS::fill_nscbc_ghost_cells(): shell state is not defined");

    MultiFab const& shell = *nscbc_shell.state;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(state.boxArray() == shell.boxArray(),"CNS::fill_nscbc_ghost_cells(): incompatible BoxArray");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(state.DistributionMap() == shell.DistributionMap(),"CNS::fill_nscbc_ghost_cells(): incompatible DistributionMap");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(state.nComp() == shell.nComp(),"CNS::fill_nscbc_ghost_cells(): incompatible component count");

    const int ng = amrex::min(state.nGrow(), shell.nGrow());

    const int ncomp = state.nComp();

    const Box domain = Geom().Domain();
    const auto dom_lo = domain.smallEnd();
    const auto dom_hi = domain.bigEnd();

    /*
     * Synchronize internal box boundaries and periodic boundaries.
     * NSCBC physical ghosts are overwritten below.
     */
    state.FillBoundary(Geom().periodicity());

    for (MFIter mfi(state, false); mfi.isValid(); ++mfi) {

        auto const& dst = state.array(mfi);
        auto const& src = shell.const_array(mfi);

        const Box valid = mfi.validbox();

        /*
         * Restrict work to ghost cells allocated by this FAB.
         * The intersection is useful when a factory supplies a FAB box
         * different from grow(valid,ng).
         */
        const Box available = (amrex::grow(valid, ng) & mfi.fabbox());

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

            // ==============================================================
            // Low physical boundary
            // ==============================================================
            if (nscbc_lo[dir] > 0 && valid.smallEnd(dir) == dom_lo[dir])
            {
                Box ghost_box = available;

                ghost_box.setSmall(dir, dom_lo[dir] - ng);
                ghost_box.setBig  (dir, dom_lo[dir] - 1);

                if (ghost_box.ok()) {
                    ParallelFor(ghost_box, ncomp, [=] AMREX_GPU_DEVICE( int i, int j, int k, int n) noexcept
                        {
                            dst(i,j,k,n) = src(i,j,k,n);
                        });
                }
            }
            // ==============================================================
            // High physical boundary
            // ==============================================================
            if (nscbc_hi[dir] > 0 && valid.bigEnd(dir) == dom_hi[dir])
            {
                Box ghost_box = available;

                ghost_box.setSmall(dir, dom_hi[dir] + 1);
                ghost_box.setBig  (dir, dom_hi[dir] + ng);

                if (ghost_box.ok()) {
                    ParallelFor(ghost_box,ncomp,[=] AMREX_GPU_DEVICE(int i,int j,int k,int n) noexcept
                        {
                            dst(i,j,k,n) = src(i,j,k,n);
                        });
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
// Set the NSCBC ghost-shell RHS to zero.
// Only the shell cells are later overwritten by compute_nscbc_ghostcell_rhs().  and
// Clearing the complete MultiFab is cheap and prevents stale RHS values
//------------------------------------------------------------------------------
void CNS::clear_nscbc_ghost_rhs()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( nscbc_shell.rhs != nullptr,
        "CNS::clear_nscbc_ghost_rhs(): ghost-shell RHS is not defined");

    nscbc_shell.rhs->setVal(amrex::Real(0.0));
}
//------------------------------------------------------------------------------
// Copy one NSCBC ghost-shell MultiFab into another.
//
// The source and destination must have identical layouts.
//------------------------------------------------------------------------------
void CNS::copy_nscbc_shell(amrex::MultiFab& dst, amrex::MultiFab const& src)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst.boxArray() == src.boxArray(),
        "copy_nscbc_shell: incompatible BoxArray");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst.DistributionMap() == src.DistributionMap(),
        "copy_nscbc_shell: incompatible DistributionMap");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst.nComp() == src.nComp(),
        "copy_nscbc_shell: incompatible number of components");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst.nGrow() == src.nGrow(),
        "copy_nscbc_shell: incompatible number of ghost cells");

    amrex::MultiFab::Copy( dst, src, 0, 0,src.nComp(),src.nGrow());   
}


