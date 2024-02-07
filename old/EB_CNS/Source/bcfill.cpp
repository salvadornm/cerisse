#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

#include "CNS.H"
#include "prob.H"

using namespace amrex;

struct CnsFillExtDir
{
  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;

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
          s_ext[dcomp + nc] = dest(iv, dcomp + nc); // this is needed for turbinflow
        }

        bcnormal(x, s_int, s_refl, s_ext, idir, 1, time, geom, *lprobparm
#ifdef USE_PMFDATA
                 ,
                 lpmfdata
#endif
        ); // Call bcnormal from prob.H

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
          s_ext[dcomp + nc] = dest(iv, dcomp + nc); // this is needed for turbinflow
        }

        bcnormal(x, s_int, s_refl, s_ext, idir, -1, time, geom, *lprobparm
#ifdef USE_PMFDATA
                 ,
                 lpmfdata
#endif
        );

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
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata = nullptr;
#ifdef USE_PMFDATA
  lpmfdata = CNS::pmf_data.getDeviceData();
#endif

  // Fill turb_inflow if requested (from PeleLM)
  if (CNS::turb_inflow.is_initialized() && scomp < AMREX_SPACEDIM) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(), dir) & bx);
      if (bcr[1].lo()[dir] == BCType::ext_dir && bndryBoxLO.ok()) {
        // Create box with ghost cells and set them to zero
        amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
        int Grow = 4; // Being conservative
        for (int n = 0; n < AMREX_SPACEDIM; n++) { growVect[n] = Grow; }
        growVect[dir] = 0;
        amrex::Box modDom = geom.Domain();
        modDom.grow(growVect);
        auto bndryBoxLO_ghost = amrex::Box(amrex::adjCellLo(modDom, dir, Grow) & bx);
        for (int nf = 0; nf <= NUM_FIELD; ++nf) {
          data.setVal<amrex::RunOn::Host>(0.0, bndryBoxLO_ghost, nf * NVAR + UMX,
                                          AMREX_SPACEDIM);
          CNS::turb_inflow.add_turb(bndryBoxLO_ghost, data, nf * NVAR + UMX, geom,
                                    time, dir, amrex::Orientation::low);
        }
      }

      auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(), dir) & bx);
      if (bcr[1].hi()[dir] == BCType::ext_dir && bndryBoxHI.ok()) {
        // Create box with ghost cells and set them to zero
        amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
        int Grow = 4;
        for (int n = 0; n < AMREX_SPACEDIM; n++) { growVect[n] = Grow; }
        growVect[dir] = 0;
        amrex::Box modDom = geom.Domain();
        modDom.grow(growVect);
        auto bndryBoxHI_ghost = amrex::Box(amrex::adjCellHi(modDom, dir, Grow) & bx);
        for (int nf = 0; nf <= NUM_FIELD; ++nf) {
          data.setVal<amrex::RunOn::Host>(0.0, bndryBoxHI_ghost, nf * NVAR + UMX,
                                          AMREX_SPACEDIM);
          CNS::turb_inflow.add_turb(bndryBoxHI_ghost, data, nf * NVAR + UMX, geom,
                                    time, dir, amrex::Orientation::high);
        }
      }
    }
  }

  GpuBndryFuncFab<CnsFillExtDir> gpu_hyp_bndry_func(
    CnsFillExtDir{lprobparm, lpmfdata});
  gpu_hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);

  // NSCBC (done after bcnormal so that we can use the ghost cells as target values)
  // Ref: Motheau et al. AIAA J. Vol. 55, No. 10: pp. 3399-3408, 2017.
  if (CNS::do_nscbc) {
    const int* lo_bc = CNS::phys_bc.lo();
    const int* hi_bc = CNS::phys_bc.hi();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      // Lo side
      auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(), dir) & bx);
      if ((/*lo_bc[dir] == 1 ||*/ lo_bc[dir] == 2) && // inflow or outflow
          bndryBoxLO.ok()) {
        for (int nf = 0; nf <= NUM_FIELD; ++nf) {
          CNS::apply_nscbc(bx, data.array(), nf * NVAR, geom, dir, 1, lo_bc[dir]);
        }
      }

      // Hi side
      auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(), dir) & bx);
      if ((/*hi_bc[dir] == 1 ||*/ hi_bc[dir] == 2) && // inflow or outflow
          bndryBoxHI.ok()) {
        for (int nf = 0; nf <= NUM_FIELD; ++nf) {
          CNS::apply_nscbc(bx, data.array(), nf * NVAR, geom, dir, -1, hi_bc[dir]);
        }
      }
    }
  }
}

void cns_react_bcfill(Box const& bx, FArrayBox& data, const int dcomp,
                      const int numcomp, Geometry const& geom, const Real time,
                      const Vector<BCRec>& bcr, const int bcomp, const int scomp)
{
  GpuBndryFuncFab<CnsReactFillExtDir> gpu_react_bndry_func(CnsReactFillExtDir{});
  gpu_react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
