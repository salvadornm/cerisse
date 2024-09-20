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
    for (int l = 0; l < order; l++) {       
      for (int m = 0; m < order; m++) { 
        coefskew(l,m) = Real(0.0); coefP(l,m)= Real(0.0);
      }
    }      

    switch (order)
    {
    case 2:
      coefskew(0,0) =  Real(0.25);
      coefskew(1,0) =  Real(0.25);
      coefskew(0,1) =  Real(0.25);
      coefskew(1,1) =  Real(0.25);
      //
      coefP(0,0)    =  Real(0.5);
      coefP(1,0)    =  Real(0.5);
      break;
    case 4:
      coefskew(0,0) =  - Real(1.0/24.0);
      coefskew(0,1) =  Real(0.0);
      coefskew(0,2) =  - Real(1.0/24.0);
      coefskew(0,3) =   Real(0.0);      

      coefskew(1,0) =  Real(0.0);
      coefskew(1,1) =  Real(1.0/3.0) - Real(1.0/24.0);  
      coefskew(1,2) =  Real(1.0/3.0);
      coefskew(1,3) =  - Real(1.0/24.0);
     
      coefskew(2,0) =  - Real(1.0/24.0);      
      coefskew(2,1) =  Real(1.0/3.0);        
      coefskew(2,2) =  Real(1.0/3.0) - Real(1.0/24.0); 
      coefskew(2,3) =  - Real(0.0);
     
      coefskew(3,0) =  - Real(0.0);
      coefskew(3,1) =  - Real(1.0/24.0);
      coefskew(3,2) =  - Real(0.0);
      coefskew(3,3) =  - Real(1.0/24.0);
      
      // FV correction (following Ducros)
      coefskew(1,1) += Real(1.0/12.0);  
      coefskew(1,2) -= Real(1.0/12.0);  
      coefskew(2,1) -= Real(1.0/12.0);  
      coefskew(2,2) += Real(1.0/12.0);  

      // pressure interpolation
      coefP(0,0) = -Real(1.0/12.0);
      coefP(1,0) =  Real(7.0/12.0);
      coefP(2,0) =  Real(7.0/12.0);
      coefP(3,0) = -Real(1.0/12.0);   
      break; 
    case 6:                   
      coefskew(0,0) = Real(1.0/120.0);
      coefskew(0,3) = Real(1.0/120.0);

      coefskew(1,1) = -Real(8.0/120.0);
      coefskew(1,3) = -Real(9.0/120.0);    
      coefskew(1,4) = Real(1.0/120.0);    
      
      coefskew(2,3) = Real(45.0/120.0);
      coefskew(2,4) = -Real(9.0/120.0);    
      coefskew(2,5) = Real(1.0/120.0);    

      coefskew(3,0) = Real(1.0/120.0);
      coefskew(3,1) = -Real(9.0/120.0);    
      coefskew(3,2) = Real(45.0/120.0);    
      coefskew(3,3) = Real(37.0/120.0);    
      
      coefskew(4,1) = Real(1.0/120.0);
      coefskew(4,2) = -Real(9.0/120.0);    
      coefskew(4,4) = Real(8.0/120.0);    

      coefskew(5,2) = Real(1.0/120.0);
      coefskew(5,5) = Real(7.0/120.0);

      // no FV correction
        
      // pressure interpolation 
      coefP(0,0) = Real(1.0/60.0);
      coefP(1,0) = -Real(8.0/60.0);
      coefP(2,0) = Real(37.0/60.0);
      coefP(3,0) = Real(37.0/60.0);        
      coefP(4,0) = -Real(8.0/60.0);
      coefP(5,0) = Real(1.0/60.0);
      break;
    default: //second order
      coefskew(0,0) =  Real(0.25);
      coefskew(1,0) =  Real(0.25);
      coefskew(0,1) =  Real(0.25);
      coefskew(1,1) =  Real(0.25);
      coefP(0,0)    =  Real(0.5);
      coefP(1,0)    =  Real(0.5);
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
    
    Real rho,eint,kin;
    Real V[order],P[order];
    Real U[order][cls_t::NCONS];


    // for (int l = 0; l < order; l++) { 
    //     for (int m = 0; m < order; m++) { 
    //       std::cout << l << " , " << m <<  std::endl;
    //       std::cout  << " coefskew " << coefskew(l,m) << std::endl;
    //     }
    //     std::cout << " coefP " << coefP(l,0) << std::endl;
    // }  
      
    //  exit(0); 


    // prepare arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   

    for (int l = 0; l < order; l++) {  
      rho =prims(il,jl,kl,cls_t::QRHO);
      V[l] = prims(il,jl,kl,Qdir); 
      U[l][cls_t::URHO] = rho;
      U[l][cls_t::UMX] = rho*prims(il,jl,kl,cls_t::QU);
      U[l][cls_t::UMY] = rho*prims(il,jl,kl,cls_t::QV);
      U[l][cls_t::UMZ] = rho*prims(il,jl,kl,cls_t::QW);      
      kin  = prims(il,jl,kl,cls_t::QU)*prims(il,jl,kl,cls_t::QU);
      kin += prims(il,jl,kl,cls_t::QV)*prims(il,jl,kl,cls_t::QV);
      kin += prims(il,jl,kl,cls_t::QW)*prims(il,jl,kl,cls_t::QW);
      eint = prims(il, jl, kl, cls_t::QEINT) + Real(0.5)*kin;
      P[l] = prims(il,jl,kl,cls_t::QPRES);
      U[l][cls_t::UET] = rho*eint  + P[l];

      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
    }


    // compute fluxes
    for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
      flx(i, j, k, nvar) =  0.0_rt;     
      for (int l = 0; l < order; l++) { 
        for (int m = 0; m < order; m++) { 
          flx(i, j, k, nvar) += coefskew(l,m)*U[l][nvar]*V[m];  
        } 
      }
    }  

    // pressure flux (symmetric interpolation)    
    for (int l = 0; l < order; l++) { 
      flx(i, j, k,Qdir-1) +=  P[l]*coefP(l,0);     
    }
 
  }

  // immersed boundary ghost point  NOT READY YET
  void inline eflux_ibm() { amrex::Print() << "IBM eflux (skew)" << std::endl; }

  // artificial dissipation flux
  void inline art_dissipation_flux() {
    amrex::Print() << "AD eflux (skew)" << std::endl;
  }

  typedef Array2D<Real, 0, order, 0, order> arrCoeff_t;
  arrCoeff_t coefskew,coefP;


  int halfsten = order / 2;


  };




#endif