#ifndef Skew_H_
#define Skew_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

//-------------------
// discontinuity sensor function
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real disconSensor(Real pp, Real pl,
                                                      Real pr) {
  Real pjst = pr + 2.0_rt * pp + pl;
  Real ptvd = std::abs(pr - pp) + std::abs(pp - pl);
  return std::abs(2.0_rt * (pr - 2.0_rt * pp + pl) /
                  (pjst + ptvd + Real(1.0e-40)));
  }


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


    // compute lambda (only if AD)
    if constexpr (isAD) {
    ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->ComputeLambda(i, j, k, prims, lambda_max,cls);
                  });                    
    } 
    // ---------------------------------------------------------------------  //
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      int Qdir =  cls_t::QRHO + dir + 1; 

      // compute interface fluxes at  flx[i] => f[i-1/2]
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });

      // dissipative fluxes 
      if constexpr (isAD) {      
        ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->fluxdissip_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });       
      }

      // add flux derivative to rhs, i.e.  rhs + = (flx[i] - flx[i+1])/dx
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });


    // if constexpr (isIB) {
    //   eflux_ibm();
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

    // prepare arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   

    for (int l = 0; l < order; l++) {  
      V[l] = prims(il,jl,kl,Qdir); 
      U[l][cls_t::URHO] = prims(il,jl,kl,cls_t::QRHO);
      U[l][cls_t::UMX]  = cons(il,jl,kl,cls_t::UMX);
      U[l][cls_t::UMY]  = cons(il,jl,kl,cls_t::UMY);
      U[l][cls_t::UMZ]  = cons(il,jl,kl,cls_t::UMZ);      
      P[l] = prims(il,jl,kl,cls_t::QPRES);
      U[l][cls_t::UET] = cons(il,jl,kl,cls_t::UET)  + P[l];

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

  
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void fluxdissip_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda, const Array4<Real>& flx,
    const cls_t* cls) const {
         
    int ir = i+vdir[0];   int jr = j+vdir[1];   int kr = k+vdir[2];   
    int il = i-vdir[0];   int jl = j-vdir[1];   int kl = k-vdir[2];   
    int ill= il-vdir[0]; int jll = jl-vdir[1]; int kll= kl-vdir[2];   

    // calculate sensor based on P   
    Real sen_num= Real(0.0),sen_denom=Real(1.0e-16);
    // do it clever in a loop NVARSEN
    int nv = cls_t::QPRES;
    
    // sensor
    Real p0 =  prims(ill,jll,kll,nv);
    Real p1 =  prims(il,jl,kl,nv);
    Real p2 =  prims(i,j,k,nv);
    Real p3 =  prims(ir,jr,kr,nv);    
    Real sen  = std::max(disconSensor(p0,p1,p2), disconSensor(p1,p2,p3) );
    sen_num += sen*sen;sen_denom +=sen;
    sen = sen_num/sen_denom;

    // spectral radius Jacobian matrix (u + c)
    Real rr = std::max(lambda(il, jl, kl, Qdir), lambda(i, j, k, Qdir));    
    rr = std::abs( prims(i,j,k,Qdir) ) + prims(i, j, k, cls_t::QC);           
    Real eps2 = Cshock*rr*sen;
    Real eps4 = std::max(0.0, Cdamp*rr - eps2);
        
    for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {   
      // shock capturing
      // 2nd order     
      // flx(i, j, k, nvar) -=  eps2*(cons(i, j, k, nvar) - cons(il, jl, kl, nvar));

      // 4th order
      flx(i, j, k, nvar) -=  eps2*(-cons(ir, jr, kr, nvar) + 7.0* cons(i, j, k, nvar)
                               -5.0*cons(il, jl, kl, nvar) - cons(ill, jll, kll, nvar)  );


      // high freq damping   
      flx(i, j, k, nvar) += eps4*(cons(ir, jr, kr, nvar) - 3.0 * cons(i, j, k, nvar) +
        3.0 * cons(il, jl, kl, nvar) - cons(ill, jll, kll, nvar));
        
    }

  } 
  
  // compute lambda = u + c
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ComputeLambda(
    int i, int j, int k, const auto& prims, const auto& lambda,
    const cls_t* cls) const {
 
    Real gamma = prims(i, j, k, cls_t::QG); 
    // speed of sound
    Real cs    = prims(i, j, k, cls_t::QC); 
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      lambda(i, j, k, dir) = std::abs( prims(i,j,k,dir)) + cs;    
    }

  }    

  // immersed boundary ghost point  NOT READY YET
  void inline eflux_ibm() { amrex::Print() << "IBM eflux (skew)" << std::endl; }

  // coefficents
  typedef Array2D<Real, 0, order, 0, order> arrCoeff_t;
  arrCoeff_t coefskew,coefP;

  int halfsten = order / 2;

  // parameters for damping and shock capturing   
  const Real Cshock = 0.1,Cdamp = 0.00016; //0.016

  };


#endif