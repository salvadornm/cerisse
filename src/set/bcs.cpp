
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <CNS.h>
#include <prob.h>
using namespace amrex;

// This is called per boundary point
struct CnsFillExtDir {
  // create pointers to device (for gpu) parms
  PROB::ProbParm* lprobparm = CNS::d_prob_parm;
  PROB::ProbClosures* lclosures = CNS::d_prob_closures;

  AMREX_GPU_DEVICE
  void operator()(const IntVect& iv, Array4<Real> const& dest, const int dcomp,
                  const int numcomp, GeometryData const& geom, const Real time,
                  const BCRec* bcr, const int bcomp,
                  const int /*orig_comp*/) const {
    // Get BC data
    const BCRec& bc = bcr[bcomp];

    // Get geom data
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const Real* prob_lo = geom.ProbLo();
    const Real* prob_hi = geom.ProbHi();
    const Real* dx = geom.CellSize();
    const Real x[AMREX_SPACEDIM] = {
        AMREX_D_DECL(prob_lo[0] + static_cast<Real>(iv[0] + 0.5) * dx[0],
                     prob_lo[1] + static_cast<Real>(iv[1] + 0.5) * dx[1],
                     prob_lo[2] + static_cast<Real>(iv[2] + 0.5) * dx[2])};

    Real s_int[PROB::ProbClosures::NCONS] = {0.0};
    Real s_refl[PROB::ProbClosures::NCONS] = {0.0};
    Real s_ext[PROB::ProbClosures::NCONS] = {0.0};

    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
      if ((bc.lo(idir) == BCType::ext_dir) && (iv[idir] < domlo[idir])) {
        //  Ghost | Ghost | Ghost || Real | Real | Real
        //  ^cell to be filled (ext)  ^firs         ^refl
        IntVect refl(iv);
        refl[idir] = 2 * domlo[idir] - iv[idir] -
                     1;  // reflection image of the ghost cell across boundary
        IntVect firs(iv);
        firs[idir] = domlo[idir];               // interior first cell
        for (int nc = 0; nc < numcomp; ++nc) {  // Fill internal state data
          s_int[dcomp + nc] = dest(firs, dcomp + nc);
          s_refl[dcomp + nc] = dest(refl, dcomp + nc);
        }

        Real di = Real(firs[idir] + 0.5) * dx[idir] - prob_lo[idir];
        Real de = prob_lo[idir] - x[idir];
        Real dratio = de / di;  // wall-ghost/wall-first internal distance ratio
        // For uniform grid, first ghost point dratio=1, second ghost point
        // dratio=3, 5...

        bcnormal(x, dratio, s_int, s_refl, s_ext, idir, 1, time, geom,
                 *lclosures, *lprobparm);  // Call bcnormal from prob.H

        for (int nc = 0; nc < numcomp; ++nc) {
          dest(iv, dcomp + nc) =
              s_ext[dcomp + nc];  // Only take the wanted components
        }
      } else if  // hi
          ((bc.hi(idir) == BCType::ext_dir) && (iv[idir] > domhi[idir])) {
        IntVect refl(iv);
        refl[idir] = 2 * domhi[idir] - iv[idir] + 1;
        IntVect firs(iv);
        firs[idir] = domhi[idir];
        for (int nc = 0; nc < numcomp; ++nc) {
          s_int[dcomp + nc] = dest(firs, dcomp + nc);
          s_refl[dcomp + nc] = dest(refl, dcomp + nc);
        }

        Real di = prob_hi[idir] - Real(firs[idir] + 0.5) * dx[idir];
        Real de = x[idir] - prob_hi[idir];
        Real dratio = de / di;  // wall-ghost/wall-first internal distance ratio

        bcnormal(x, dratio, s_int, s_refl, s_ext, idir, -1, time, geom,
                 *lclosures, *lprobparm);

        for (int nc = 0; nc < numcomp; ++nc) {
          dest(iv, dcomp + nc) = s_ext[dcomp + nc];
        }
      }
    }
  }
};

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the descriptor set up
// in CNS::variableSetUp. This is called once per fab (if fab has domain
// boundary and NCONS=NVAR boundary to set)
void cns_bcfill(Box const& bx, FArrayBox& data, const int dcomp,
                const int numcomp, Geometry const& geom, const Real time,
                const Vector<BCRec>& bcr, const int bcomp, const int scomp) {
  GpuBndryFuncFab<CnsFillExtDir> gpu_bndry_func(CnsFillExtDir{});
  gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
  Gpu::streamSynchronize();

  // TODO : pass 0,1,Nghost internal points to bcnormal. The current approach
}