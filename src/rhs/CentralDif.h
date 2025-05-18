#ifndef CentralDif_H_
#define CentralDif_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

#include "diff_ops.H"
template <bool isAD, bool isIB, int order, typename cls_t>
class centraldif_t {
  public:

  AMREX_GPU_HOST_DEVICE
  centraldif_t() {
    // initialize coefficients for flux interpolation based on order
    calc_CDcoeffs<order>(INTcoef,CDcoef);
  }

  AMREX_GPU_HOST_DEVICE
  ~centraldif_t() {}

  // vars accessed by functions 
  int order_sch=order;  
  int halfsten = order / 2;

  typedef Array1D<Real, 0, order> arrayNumCoef;
  arrayNumCoef CDcoef,INTcoef;

  
  //////////////////////////////////////////////////////
  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& cons, const cls_t* cls) {
                      
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls->NGHOST);
    const Box& bxgnodal = mfi.grownnodaltilebox(
        -1, 0);  // extent is 0,N_cell+1 in all directions -- -1 means for all
                 // directions. amrex::surroundingNodes(bx) does the same

    // ---------------------------------------------------------------------  //
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      auto const& flx = flxt[dir]->array(); 

      // compute interface fluxes at i-1/2, j-1/2, k-1/2
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k,dir, vdir, cons, prims, flx, cls);
                  });
    }
  }  

  // compute flux in each direction at f[i-1/2]   stored in i,j,k
  // central formulation following f[i-1/2] = 1/2 (fi + fi-1)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int dir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    Real flux_l[cls_t::NCONS];
    
    // prepare flux arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   
        
    for (int l = 0; l < order; l++) {  

      IntVect iv(AMREX_D_DECL(il, jl, kl));

      // evaluate flux from primitive
      cls->prims2flux(iv,dir,prims,flux_l);

      // compute flux
      for (int n=0;n< cls_t::NCONS;n++){
        flx(i, j, k,n) += flux_l[n]*INTcoef(l);
      }  

      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];

    }

  }
  ////////////////////////////////////////////////////////////////////////////////////////


  };




#endif