#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

#include "CNS.H"
#include "prob.H"

using namespace amrex;

struct CnsFillExtDir
{
  ProbParm const* lprobparm;

  /**
   * \param iv    index of the ghost cell
   * \param[out] dest data array
   * \param dcomp starting component of the data requested
   * \param numcomp number of components needed
   * \param geom  geometry data
   * \param time  time
   * \param bcr   boundary conditions holder
   */
  AMREX_GPU_DEVICE void operator()(const IntVect& iv, Array4<Real> const& dest,
                                   const int dcomp, const int numcomp,
                                   GeometryData const& geom, const Real time,
                                   const BCRec* bcr, const int bcomp,
                                   const int /*orig_comp*/) const
  {
    // Get BC data
    const BCRec& bc = bcr[bcomp];

    // Get geom data
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const Real* prob_lo = geom.ProbLo();
    const Real* dx = geom.CellSize();
    const Real x[AMREX_SPACEDIM] = {
      AMREX_D_DECL(prob_lo[0] + static_cast<Real>(iv[0] + 0.5) * dx[0],
                   prob_lo[1] + static_cast<Real>(iv[1] + 0.5) * dx[1],
                   prob_lo[2] + static_cast<Real>(iv[2] + 0.5) * dx[2])};

    Real s_int[LEN_STATE] = {0.0};
    Real s_refl[LEN_STATE] = {0.0};
    Real s_ext[LEN_STATE] = {0.0};

    // Fill internal state data then run bcnormal to fill external state data
    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
      if // lo
        ((bc.lo(idir) == BCType::ext_dir) && (iv[idir] < domlo[idir])) {
        //  Ghost | Ghost | Ghost || Real | Real | Real
        //  ^cell to be filled       ^firs         ^refl
        IntVect refl(iv);
        refl[idir] = 2 * domlo[idir] - iv[idir] -
                     1; // reflection image of the ghost cell across boundary
        IntVect firs(iv);
        firs[idir] = domlo[idir];              // interior first cell
        for (int nc = 0; nc < numcomp; ++nc) { // Fill internal state data
          s_int[dcomp + nc] = dest(firs, dcomp + nc);
          s_refl[dcomp + nc] = dest(refl, dcomp + nc);
        }

        bcnormal(x, s_int, s_refl, s_ext, idir, 1, time, geom,
                 *lprobparm); // Call bcnormal from prob.H

        for (int nc = 0; nc < numcomp; ++nc) {
          dest(iv, dcomp + nc) =
            s_ext[dcomp + nc]; // Only take the wanted components
        }
      } else if // hi
        ((bc.hi(idir) == BCType::ext_dir) && (iv[idir] > domhi[idir])) {
        IntVect refl(iv);
        refl[idir] = 2 * domhi[idir] - iv[idir] + 1;
        IntVect firs(iv);
        firs[idir] = domhi[idir];
        for (int nc = 0; nc < numcomp; ++nc) {
          s_int[dcomp + nc] = dest(firs, dcomp + nc);
          s_refl[dcomp + nc] = dest(refl, dcomp + nc);
        }

        bcnormal(x, s_int, s_refl, s_ext, idir, -1, time, geom, *lprobparm);

        for (int nc = 0; nc < numcomp; ++nc) {
          dest(iv, dcomp + nc) = s_ext[dcomp + nc];
        }
      }
    }
  }
};

struct CnsReactFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(const amrex::IntVect& /*iv*/,
                  amrex::Array4<amrex::Real> const& /*dest*/, const int /*dcomp*/,
                  const int /*numcomp*/, amrex::GeometryData const& /*geom*/,
                  const amrex::Real /*time*/, const amrex::BCRec* /*bcr*/,
                  const int /*bcomp*/, const int /*orig_comp*/) const
  {
    // do something for external Dirichlet (BCType::ext_dir)
  }
};

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the descriptor set up in
// CNS::variableSetUp.

void cns_bcfill(Box const& bx, FArrayBox& data, const int dcomp, const int numcomp,
                Geometry const& geom, const Real time, const Vector<BCRec>& bcr,
                const int bcomp, const int scomp)
{
  ProbParm const* lprobparm = CNS::d_prob_parm;

  GpuBndryFuncFab<CnsFillExtDir> gpu_hyp_bndry_func(CnsFillExtDir{lprobparm});
  gpu_hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

void cns_react_bcfill(Box const& bx, FArrayBox& data, const int dcomp,
                      const int numcomp, Geometry const& geom, const Real time,
                      const Vector<BCRec>& bcr, const int bcomp, const int scomp)
{
  GpuBndryFuncFab<CnsReactFillExtDir> gpu_react_bndry_func(CnsReactFillExtDir{});
  gpu_react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
