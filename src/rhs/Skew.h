
#ifndef Skew_H_
#define Skew_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

template <bool isAD, bool isIB, int order, typename cls_t>
class skew_t {
 public:

  AMREX_GPU_HOST_DEVICE
  skew_t() {
    // initialize coefficients skew symmetric TODO
    coeffs(0, 0) = 1.0;
    coeffs(0, 1) = 0.0;
    coeffs(0, 2) = 0.0;
    coeffs(1, 0) = 4.0 / 3;
    coeffs(1, 1) = -2.0 / 12;
    coeffs(1, 2) = 0.0;
    coeffs(2, 0) = 6.0 / 4;
    coeffs(2, 1) = -6.0 / 20;
    coeffs(2, 2) = 2.0 / 60;
  }

  int order_skew = order;


  AMREX_GPU_HOST_DEVICE
  ~skew_t() {}

  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& flx,
                    const Array4<Real>& rhs,
                    const cls_t* cls) {
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

      int Qdir =  cls.QRHO + dir

      // compute interface fluxes at i-1/2, j-1/2, k-1/2
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });

      // add dissipative fluxes            

      // add flux derivative to rhs = -(fi+1 - fi)/dx = (fi - fi+1)/dx
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });
    }


    if constexpr (isIB) {
      eflux_ibm();
    }

    if constexpr (isAD) {
      art_dissipation_flux();
    }


  }

  // compute flux in each direction at i-1/2 //  
  // f = U*V (where U= (rho, rho ux, rho uy, rho uz rho e) and V = (ux,uy,uz)    //
  // f = 1/2 ( U + Ui-1/2) * 1/2 *(V + Vi-1/2)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    int il= i-vdir[0]; int jl= j-vdir[1]; int kl= k-vdir[2];


     
    // skew 2
    for (int n = 0; n < cls_t::NCONS; n++) {
      flx(i, j, k, n) = Real(0.25)*( cons(il,jl,kl,n) + cons(i,j,k,n) )*(prims(il,jl,kl,Qdir) + prims(i,j,k,Qdir) )
      };

    // how do I do 2,4,6  

    }

  // immersed boundary ghost point
  void inline eflux_ibm() { amrex::Print() << "IBM eflux (skew)" << std::endl; }

  // artificial dissipation flux
  void inline art_dissipation_flux() {
    amrex::Print() << "AD eflux (skew)" << std::endl;
  }


  };




#endif