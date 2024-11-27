
#ifndef Rusanov_H_
#define Rusanov_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

template <typename cls_t>
class rusanov_t {
 public:

  AMREX_GPU_HOST_DEVICE
  rusanov_t() {}

  AMREX_GPU_HOST_DEVICE
  ~rusanov_t() {}

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void inline eflux_ibm(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& flx,
                    const Array4<Real>& rhs,
                    const cls_t* cls, const Array4<bool>& ibMarkers) {
#else
  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& flx,
                    const Array4<Real>& rhs,
                    const cls_t* cls) {
#endif

    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls->NGHOST);
    const Box& bxgnodal = mfi.grownnodaltilebox(
        -1, 0);  // extent is 0,N_cell+1 in all directions -- -1 means for all
                 // directions. amrex::surroundingNodes(bx) does the same
    // copy cons
    FArrayBox cell_fluxesf(bxg, cls_t::NCONS, The_Async_Arena()); // cell center fluxes; flx is interface fluxes
    FArrayBox consf(bxg, cls_t::NCONS, The_Async_Arena());
    FArrayBox lambda_maxf(bxg, 1, The_Async_Arena());

    Array4<Real> cell_fluxes = cell_fluxesf.array();
    Array4<Real> cons   = consf.array();
    Array4<Real> lambda_max = lambda_maxf.array();


    // copy conservative variables from rhs to cons and clear rhs
    ParallelFor(bxg, cls_t::NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      cons(i,j,k,n) = rhs(i,j,k,n);
      rhs(i, j, k, n)=0.0;
      });

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
      ParallelFor(bxg,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    GpuArray<Real, cls_t::NWAVES> eigenvalues = cls->cons2eigenvals(i, j, k, cons, vdir);

                    lambda_max(i,j,k,0) = -1e40;
                    for (int n = 0; n < cls_t::NWAVES; n++) {
                      lambda_max(i, j, k, 0) = std::max(lambda_max(i, j, k, 0), std::abs(eigenvalues[n]));
                    }
                    cls->prims2fluxes(i, j, k, prims, cell_fluxes, vdir);
                  });

      // compute interface fluxes at i-1/2, j-1/2, k-1/2
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k, vdir, cons, cell_fluxes, lambda_max, flx, cls);
                  });

      // add flux derivative to rhs = -(fi+1 - fi)/dx = (fi - fi+1)/dx
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& cell_flux, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    int il= i-vdir[0]; int jl= j-vdir[1]; int kl= k-vdir[2];

    Real lambda = std::max(lambda_max(i,j,k,0),lambda_max(il,jl,kl,0));
    AMREX_ASSERT_WITH_MESSAGE(lambda >= 0.0_rt, "lambda is negative");
    for (int n = 0; n < cls_t::NCONS; n++) {
      flx(i, j, k, n) = 0.5_rt*((cell_flux(i,j,k,n) + cell_flux(il,jl,kl,n)) 
                      - lambda*(cons(i,j,k,n) - cons(il,jl,kl,n)));
      };
    }
  };
#endif