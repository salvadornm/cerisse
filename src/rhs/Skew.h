#ifndef Skew_H_
#define Skew_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

template <bool isAD, bool isIB, int order, typename cls_t>
class skew_t {
  public:

  AMREX_GPU_HOST_DEVICE
  skew_t() {
    
    // initialize coefficients for skew-symmetric
    switch (order)
    {
    case 2:
      coefskew(0,0) =  Real(0.25);
      coefskew(1,0) =  Real(0.25);
      coefskew(0,1) =  Real(0.25);
      coefskew(1,1) =  Real(0.25);
      break;
    case 4:
      coefskew(1,1) =  Real(1.0/3.0) - Real(1.0/24.0); 
      coefskew(2,1) =  Real(1.0/3.0);       
      coefskew(1,2) =  Real(1.0/3.0); 
      coefskew(2,2) =  Real(1.0/3.0) - Real(1.0/24.0);     
      coefskew(0,0) =  - Real(1.0/24.0);
      coefskew(0,1) =  - Real(1.0/24.0);
      coefskew(1,3) =  - Real(1.0/24.0);
      coefskew(2,0) =  - Real(1.0/24.0);                     
      coefskew(3,1) =  - Real(1.0/24.0);
      coefskew(3,2) =  - Real(1.0/24.0);
      // FV correction (following Ducros)
      coefskew(1,1) +=  Real(1.0/6.0) - Real(1.0/12.0);  
      coefskew(2,2) +=  Real(1.0/6.0) - Real(1.0/12.0);  
      coefskew(1,2) -= Real(1.0/12.0);  
      coefskew(2,1) -=  Real(1.0/6.0) - Real(1.0/12.0);  
      break; 
    case 6:    //TODO         
      // coefskew(0,0) = -Real(1.0/12.0);
      // coefskew(1,0) =  Real(7.0/12.0);
      // coefskew(2,0) =  Real(7.0/12.0);
      // coefskew(3,0) = -Real(1.0/12.0);    
      break;
    default: //second order
      coefskew(0,0) =  Real(0.25);
      coefskew(1,0) =  Real(0.25);
      coefskew(0,1) =  Real(0.25);
      coefskew(1,1) =  Real(0.25);
      break;
    }
  
  }

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

      // add flux derivative to rhs, i.e.  rhs[n] + = (fi[n] - fi+1[n])/dx
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });

    // if constexpr (isIB) {
    //   eflux_ibm();
    // }

    // if constexpr (isAD) {
    //   art_dissipation_flux();
    // }

    }
  }

  // compute flux in each direction at i-1/2 //  
  // skew-symmetric formulation following f = U*V
  // U vector of conservative  vars (rho, rho ux, rho uy, rho uz rho e) //
  // V velocity vector               (ux,uy,uz)    //
  // fi = 1/2 ( U + Ui-1) * 1/2 *(V + Vi-1)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    Real rho,eint,kin,fluxskew[cls_t::NCONS];
    Real V[order];
    Real U[order][cls_t::NCONS];
    
    // prepare arrays
    int il= i-order*vdir[0]; int jl= j-order*vdir[1]; int kl= k-order*vdir[2];   
    for (int l = 0; l < order; l++) {  
      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
      rho =prims(il,jl,kl,cls_t::QRHO);
      V[l] = prims(il,jl,kl,Qdir); 
      U[l][cls_t::QRHO] = rho;
      U[l][cls_t::QU] = rho*prims(il,jl,kl,cls_t::QU);
      U[l][cls_t::QV] = rho*prims(il,jl,kl,cls_t::QV);
      U[l][cls_t::QW] = rho*prims(il,jl,kl,cls_t::QW);      
      kin  = prims(il,jl,kl,cls_t::QU)*prims(il,jl,kl,cls_t::QU);
      kin += prims(il,jl,kl,cls_t::QV)*prims(il,jl,kl,cls_t::QV);
      kin += prims(il,jl,kl,cls_t::QW)*prims(il,jl,kl,cls_t::QW);
      eint= prims(il, jl, kl, cls_t::QEINT) + Real(0.5)*kin;
      U[l][cls_t::NCONS-1] = rho*eint  + prims(il,jl,kl,cls_t::QPRES);
    }

    // compute fluxes
    for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
      fluxskew[nvar] =  Real(0.0);     
      for (int l = 0; l < order; l++) { 
        for (int m = 0; m < order; m++) { 
          fluxskew[nvar]+= coefskew(l,m)*U[l][nvar]*V[m]; // +coef*U*V
          }
      }
    }  

    // pressure flux (linear interpolation pressure)
    const Real P  = (prims(il,jl,kl,cls_t::QPRES) + prims(i,j,k,cls_t::QPRES))*Real(0.5);
    fluxskew[Qdir] += P;
  
    // fluxskew (primitive notartion)-->flux (conervative notation)  //
    flx(i, j, k, cls_t::URHO) = fluxskew[cls_t::QRHO];
    flx(i, j, k, cls_t::UMX)  = fluxskew[cls_t::QU];
    flx(i, j, k, cls_t::UMX)  = fluxskew[cls_t::QV];
    flx(i, j, k, cls_t::UMY)  = fluxskew[cls_t::QW];
    flx(i, j, k, cls_t::UET)  = fluxskew[cls_t::NCONS-1];
        
    }

  // immersed boundary ghost point  NOT READY YET
  void inline eflux_ibm() { amrex::Print() << "IBM eflux (skew)" << std::endl; }

  // artificial dissipation flux
  void inline art_dissipation_flux() {
    amrex::Print() << "AD eflux (skew)" << std::endl;
  }

  typedef Array2D<Real, 0, order, 0, order> arrCoeff_t;
  arrCoeff_t coefskew;

  };




#endif