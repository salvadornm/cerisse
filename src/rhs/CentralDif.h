#ifndef CentralDif_H_
#define CentralDif_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

template <bool isAD, bool isIB, int order, typename cls_t>
class centraldif_t {
  public:

  AMREX_GPU_HOST_DEVICE
  centraldif_t() {
    
    // initialize coefficients for symmetrical interpolation based on order
    switch (order)
    {
    case 2:
      coefdif(0,0) =  Real(0.5);
      coefdif(1,0) =  Real(0.5);
      break;
    case 4:
      coefdif(0,0) = -Real(1.0/12.0);
      coefdif(1,0) =  Real(7.0/12.0);
      coefdif(2,0) =  Real(7.0/12.0);
      coefdif(3,0) = -Real(1.0/12.0);                     
      break; 
    case 6:             
      coefdif(0,0) = Real(1.0/60.0);
      coefdif(1,0) = -Real(8.0/60.0);
      coefdif(2,0) = Real(37.0/60.0);
      coefdif(3,0) = Real(37.0/60.0);        
      coefdif(4,0) = -Real(8.0/60.0);
      coefdif(5,0) = Real(1.0/60.0);
      break;
    default: //second order
      coefdif(0,0) =  Real(0.5);
      coefdif(1,0) =  Real(0.5);
      break;
    }
  
  }

  AMREX_GPU_HOST_DEVICE
  ~centraldif_t() {}

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
    FArrayBox cell_fluxesf(bxg, cls_t::NCONS, The_Async_Arena());
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

    // ---------------------------------------------------------------------  //
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      int Qdir =  cls_t::QRHO + dir + 1; 

      // compute interface fluxes at i-1/2, j-1/2, k-1/2
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });

      // add dissipative fluxes            

      // add flux derivative to rhs, i.e.  rhs[n] += (fi[n] - fi+1[n])/dx
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });
    }

  }  

  // compute flux in each direction at i-1/2 //  
  // central formulation following f = 1/2 (fi + fi+1)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    Real rho,eint,kin,rhou,UN;
    Real flux[order][cls_t::NCONS];
     
    // prepare flux arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   
        
    for (int l = 0; l < order; l++) {  

      rho  = prims(il,jl,kl,cls_t::QRHO); UN   = prims(il,jl,kl,Qdir);
      rhou = rho*UN;
      flux[l][cls_t::URHO] = rhou;
      flux[l][cls_t::UMX]  = rhou*prims(il,jl,kl,cls_t::QU);
      flux[l][cls_t::UMY]  = rhou*prims(il,jl,kl,cls_t::QV);
      flux[l][cls_t::UMZ]  = rhou*prims(il,jl,kl,cls_t::QW);
      flux[l][Qdir-1]     +=  prims(il,jl,kl,cls_t::QPRES);
      kin  = prims(il,jl,kl,cls_t::QU)*prims(il,jl,kl,cls_t::QU);
      kin += prims(il,jl,kl,cls_t::QV)*prims(il,jl,kl,cls_t::QV);
      kin += prims(il,jl,kl,cls_t::QW)*prims(il,jl,kl,cls_t::QW);
      eint= prims(il, jl, kl, cls_t::QEINT) + 0.5_rt*kin; 
      flux[l][cls_t::UET]  = rhou*eint + UN*prims(il,jl,kl,cls_t::QPRES);     

      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
    }

    // compute fluxes
    for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
      flx(i, j, k, nvar) =  0.0_rt;     
      for (int l = 0; l < order; l++) { 
        flx(i, j, k, nvar) +=  flux[l][nvar]*coefdif(l,0);     
      }
    }  

    // adding high freq damping?

    // for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
    //   flx(i, j, k, nvar) += 
    // }  

    
  }

  // immersed boundary ghost point  NOT READY YET
  void inline eflux_ibm() { amrex::Print() << "IBM eflux (central dif)" << std::endl; }

  // artificial dissipation flux NOT READY YET
  void inline art_dissipation_flux() {
    amrex::Print() << "AD eflux (central dif)" << std::endl;
  }

  typedef Array2D<Real, 0, order, 0, order> arrCoeff_t;
  arrCoeff_t coefdif;

  int halfsten = order / 2;

  };




#endif