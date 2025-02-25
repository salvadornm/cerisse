#ifndef Skew_H_
#define Skew_H_

#include <AMReX_FArrayBox.H>


//-------------------
// discontinuity sensor function
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real disconSensor(Real pp, Real pl,
                                                      Real pr) {
  Real pjst = pr + 2.0_rt * pp + pl;
  Real ptvd = std::abs(pr - pp) + std::abs(pp - pl);
  return std::abs(2.0_rt * (pr - 2.0_rt * pp + pl) /
                  (pjst + ptvd + Real(1.0e-40)));
  }

// int imask
template <typename param, typename cls_t>
class skew_t {
  public:

  AMREX_GPU_HOST_DEVICE
  skew_t() {
    
    // init coefficients for skew-symmetric
    for (int l = 0; l < order; l++) {       
      for (int m = 0; m < order; m++) { 
        coefskew(l,m) = Real(0.0); 
      }
    }
    // init interpolation, damping and shock coeffcients      
    for (int l = 0; l < order; l++) {       
      coefdamp(l)  = Real(0.0);
      coefP(l)     = Real(0.0);
      coefshock(l) = Real(0.0);
    }    
    
    switch (order)
    {
    case 2:
      coefskew(0,0) =  Real(0.25);
      coefskew(1,0) =  Real(0.25);
      coefskew(0,1) =  Real(0.25);
      coefskew(1,1) =  Real(0.25);
      // pressure interpolation
      coefP(0)      =  Real(0.5);
      coefP(1)      =  Real(0.5);
      // damping coefficients
      coefdamp(0)   = Real(0.0);
      coefdamp(1)   = Real(0.0);
      // shock
      coefshock(0) =  -Real(1.0);
      coefshock(1) =  +Real(1.0);            
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
      coefP(0) = -Real(1.0/12.0);
      coefP(1) =  Real(7.0/12.0);
      coefP(2) =  Real(7.0/12.0);
      coefP(3) = -Real(1.0/12.0);   

      // damping coefficients
      coefdamp(0) = -Real(1.0);
      coefdamp(1) = +Real(3.0);
      coefdamp(2) = -Real(3.0);
      coefdamp(3) = +Real(1.0);

      // shock coefficients
      coefshock(0) = -Real(1.0/6.0);
      coefshock(1) = -Real(5.0/6.0);
      coefshock(2) = +Real(7.0/6.0);
      coefshock(3) = -Real(1.0/6.0);

      break; 
    case 6:                   
      coefskew(0,0) = Real(1.0/120.0);
      coefskew(0,3) = Real(1.0/120.0);

      coefskew(1,1) = -Real(8.0/120.0);
      coefskew(1,3) = -Real(9.0/120.0);    
      coefskew(1,4) = Real(1.0/120.0);    
      
      coefskew(2,2) = Real(37.0/120.0);
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
      coefP(0) = Real(1.0/60.0);
      coefP(1) = -Real(8.0/60.0);
      coefP(2) = Real(37.0/60.0);
      coefP(3) = Real(37.0/60.0);        
      coefP(4) = -Real(8.0/60.0);
      coefP(5) = Real(1.0/60.0);

      // damping coefficients
      coefdamp(0) = +Real(1.0/12.0);
      coefdamp(1) = -Real(17.0/12.0);
      coefdamp(2) = +Real(46.0/12.0);
      coefdamp(3) = -Real(46.0/12.0);
      coefdamp(4) = +Real(17.0/12.0);
      coefdamp(5) = -Real(1.0/12.0);

      // shock coefficients
      coefshock(0) =  Real(0.0);      
      coefshock(1) = -Real(1.0/6.0);
      coefshock(2) = -Real(5.0/6.0);
      coefshock(3) = +Real(7.0/6.0);
      coefshock(4) = -Real(1.0/6.0);
      coefshock(5) =  Real(0.0);                             
      break;

    default: 

      amrex::Abort("Skew-symmetric order available 2/4/6 MUST specify one of those ");      

      break;
    }
    
    // sensor variables (denisty and pressure by default)
    NSEN[0] = cls_t::QRHO;
    NSEN[1] = cls_t::QPRES;

    // no masking
    for (int l = 0; l < AMREX_SPACEDIM; l++) {  
      mask_sen(l,1) = -1;     
      mask_sen(l,2) = 9999999;     
    }  
    
  }

  AMREX_GPU_HOST_DEVICE
  ~skew_t() {}

  // parameters for damping and shock capturing   
  static const int order = param::order;
  // default values Cshock = 0.1 Cdamp = 0.016; 
  Real Cshock = param::C2skew;
  Real Cdamp  = param::C4skew; 


#if (AMREX_USE_GPIBM || CNS_USE_EB )  
 void inline eflux_ibm(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& rhs, const cls_t* cls,const Array4<bool>& ibMarkers) {

#else
  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& rhs, const cls_t* cls) {
#endif


    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(cls->NGHOST);
    const Box& bxgnodal = mfi.grownnodaltilebox(
        -1, 0);  // extent is 0,N_cell+1 in all directions -- -1 means for all
                 // directions. amrex::surroundingNodes(bx) does the same
    
    FArrayBox consf(bxg, cls_t::NCONS, The_Async_Arena());
    FArrayBox lambda_maxf(bxg, 1, The_Async_Arena());

    // create lambda(0) and cons array
    Array4<Real> cons   = consf.array();
    Array4<Real> lambda_max = lambda_maxf.array();

    // copy conservative variables from rhs to cons and clear rhs
    ParallelFor(bxg, cls_t::NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      cons(i,j,k,n) = rhs(i,j,k,n);
      rhs(i, j, k, n)=0.0;              
      });

    int imask = 3; // reduce order 
    // get global index
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

    // masking BC here function of global geometry  (only non-periodic dirs)       
    for (int l = 0; l < AMREX_SPACEDIM; l++) {  
      if (geom.isPeriodic(l)==0) {
        mask_sen(l,1) = domlo[l] + imask;     
        mask_sen(l,2) = domhi[l] - imask;     
      }
    }  

    // ---------------------------------------------------------------------  //
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      int Qdir =  cls_t::QRHO + dir + 1; 

      auto const& flx = flxt[dir]->array(); 
  

#if (AMREX_USE_GPIBM || CNS_USE_EB )  
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir_ibm(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls,ibMarkers);
                  });
#else    
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });                  
#endif

      
      // dissipative fluxes 
      if constexpr (param::dissipation) {

#if (AMREX_USE_GPIBM || CNS_USE_EB )  

        ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->fluxdissip_dir_ibm(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls, ibMarkers);
                  });             
#else
        ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->fluxdissip_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });  
#endif
        

      }
      // add flux derivative to rhs, i.e.  rhs + = (flx[i] - flx[i+1])/dx (THIS WILL GO)
      // ParallelFor(bx, cls_t::NCONS,
      //             [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      //               rhs(i, j, k, n) +=
      //                   dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
      //             });

    }

  

  }

  // ............................................................. 
  // compute flux in each direction at i-1/2 //  
  // skew-symmetric formulation following f = U*V
  // U vector of conservative  vars (rho, rho ux, rho uy, rho uz rho e) //
  // V velocity vector               (ux,uy,uz)    //
  // fi = 1/2 ( U + Ui-1) * 1/2 *(V + Vi-1)  (example of second order)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls) const {
    
    Real V[order],P[order];
    Real U[order][cls_t::NCONS];

    // prepare arrays
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];   

    for (int l = 0; l < order; l++) {  
      V[l] = prims(il,jl,kl,Qdir);       
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {U[l][nvar] = cons(il,jl,kl,nvar);}
      P[l] = prims(il,jl,kl,cls_t::QPRES);      
      U[l][cls_t::UET] += P[l];
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
      flx(i, j, k,Qdir-1) +=  P[l]*coefP(l);     
    }
 
  }  
  // .............................................................
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir_ibm(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls, const Array4<bool>& marker) const {
        
    Real V[order],P[order];
    Real U[order][cls_t::NCONS];

    int il= i-vdir[0]; int jl= j-vdir[1]; int kl= k-vdir[2];
    const bool close_to_wall = marker(i,j,k,1) || marker(il,jl,kl,1);           // //  flux close to a GP (IBM) or a cut-cell (EB)
    const bool intersolid_flx = marker(i,j,k,0) &&  marker(il,jl,kl,0);    // inter-flux

    if (intersolid_flx) return;  // flux =0  inside solid 
      

#ifdef CNS_USE_EB   
    // in EBM wall flux will be computed afterwards
    const bool next_to_wall = marker(i,j,k,0) || marker(il,jl,kl,0);
    if (next_to_wall) return;  
#endif


    // reduce to second order scheme close to wall
    if (close_to_wall)
    {   
      //il= i-vdir[0]; jl= j-vdir[1]; kl= k-vdir[2];
      for (int l = 0; l < 2; l++) {  
        V[l] = prims(il,jl,kl,Qdir); 
        for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {U[l][nvar] = cons(il,jl,kl,nvar);}
        P[l] = prims(il,jl,kl,cls_t::QPRES);
        U[l][cls_t::UET] += P[l];
        il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
      }
      // fluxes
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {
        flx(i, j, k, nvar) = Real(0.25)*(U[0][nvar] + U[1][nvar])*(V[0]+V[1]);  
      }
      flx(i, j, k,Qdir-1) += Real(0.5)*(P[0] + P[1]);
    }
    else
    {
      il= i-halfsten*vdir[0]; jl= j-halfsten*vdir[1]; kl= k-halfsten*vdir[2];

      for (int l = 0; l < order; l++) {  
        V[l] = prims(il,jl,kl,Qdir); 
        for (int nvar = 0; nvar < cls_t::NCONS; nvar++) {U[l][nvar] = cons(il,jl,kl,nvar);}
        P[l] = prims(il,jl,kl,cls_t::QPRES);
        U[l][cls_t::UET] += P[l];
        il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
      } 
      // fluxes      
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++)  {        
        flx(i, j, k, nvar) =  0.0_rt;             
        for (int l = 0; l < order; l++) { 
          for (int m = 0; m < order; m++) { 
            flx(i, j, k, nvar) += coefskew(l,m)*U[l][nvar]*V[m];  
          } 
        }        
      }  
      // pressure flux (symmetric interpolation)    
      for (int l = 0; l < order; l++) { 
        flx(i, j, k,Qdir-1) +=  P[l]*coefP(l);     
      }    
    } 
             
  }
  // .............................................................
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void fluxdissip_dir_ibm(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda, const Array4<Real>& flx,
    const cls_t* cls,const Array4<bool>& marker) const {

    int il= i-vdir[0]; int jl= j-vdir[1]; int kl= k-vdir[2];
    const bool close_to_wall  = marker(i,j,k,1) || marker(il,jl,kl,1);          //  flux close to a GP (IBM) or a cut-cell (EB)
    const bool intersolid_flx = marker(i,j,k,0) &&  marker(il,jl,kl,0);  

    if (intersolid_flx) return;  // flux =0  inside solid

#ifdef CNS_USE_EB   
    // in EBM wall flux will be computed afterwards
    const bool next_to_wall = marker(i,j,k,0) || marker(il,jl,kl,0);
    if (next_to_wall) return;  
#endif

    int i0[3],i1[3],i2[3],i3[3];
 
    const int idir = Qdir -1;

    // calculate sensor    
    Real p0,p1,p2,p3;
    Real sen_num= Real(0.0),sen_denom=Real(1.0e-16);
    // loop over sensor variables
    int nv = 0;
    Real sen = Real(0.0);
    if (close_to_wall)
    {          

      const bool wall_right = marker(i,j,k,0) || marker(i+vdir[0],j+vdir[1],k+vdir[2],0); // wall i(IB) or i+1 (EB)

      if (wall_right)       // solid towards the right use i-2,i-1,i(GP/CUT)
      {
        i1[0] = i - 2*vdir[0]; i1[1] = j - 2*vdir[1];i1[2] = k - 2*vdir[2];
        for (int l=0;l<3;l++) {i2[l]=i1[l]+vdir[l];i3[l]=i2[l]+vdir[l];}
      } 
      else                  // solid towards the left, use i-1(GP/CUT),i,i+1
      {
        i1[0] = i - vdir[0]; i1[1] = j - vdir[1];i1[2] = k - vdir[2];
        for (int l=0;l<3;l++) {i2[l]=i1[l]+vdir[l];i3[l]=i2[l]+vdir[l];}
      }

      for (int l=0;l<NVARSEN;l++)
      { 
        nv = NSEN[l];       
        p1 =  prims(i1[0],i1[1],i1[2],nv);
        p2 =  prims(i2[0],i2[1],i2[2],nv);
        p3 =  prims(i3[0],i2[1],i2[2],nv);    
        sen  = disconSensor(p1,p2,p3);
        sen_num += sen*sen;sen_denom +=sen;
        sen = sen_num/sen_denom;
      }
    } 
    else
    {
      i0[0] = i - 2*vdir[0]; i0[1] = j - 2*vdir[1];i0[2] = k - 2*vdir[2];
      for (int l=0;l<3;l++) {
        i1[l]=i0[l]+vdir[l];i2[l]=i1[l]+vdir[l];i3[l]=i2[l]+vdir[l];
      }
      for (int l=0;l<NVARSEN;l++)
      { 
        nv = NSEN[l];       
        p0 =  prims(i0[0],i0[1],i0[2],nv);
        p1 =  prims(i1[0],i1[1],i1[2],nv);
        p2 =  prims(i2[0],i2[1],i2[2],nv);
        p3 =  prims(i3[0],i3[1],i3[2],nv);        
        sen  = std::max(disconSensor(p0,p1,p2), disconSensor(p1,p2,p3) );
        sen_num += sen*sen;sen_denom +=sen;
        sen = sen_num/sen_denom;
      }
    }  

    // reduce order close to BC by making sensor  = 1   
    // sen = (i < mask_sen(idir,1)) ? 1.0 : sen;  
    // sen = (i > mask_sen(idir,2)) ? 1.0 : sen;  
    
    // spectral radius Jacobian matrix (u + c)  
    //Real rr = std::max(lambda(il, jl, kl, 0), lambda(i, j, k, 0));          
    Real rr = std::abs( prims(i,j,k,Qdir) ) + prims(i, j, k, cls_t::QC);              
    Real eps2 = Cshock*rr*sen;
    Real eps4 = std::max(0.0, Cdamp*rr - eps2);

    // shock capturing and damping 
    if (close_to_wall)
    {
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++) { 
        flx(i, j, k, nvar) -=  eps2*( cons(i,j,k,nvar) -  cons(il,jl,kl,nvar)) ;         
      } 
    } 
    else
    {
      int ii= i-halfsten*vdir[0]; int jj= j-halfsten*vdir[1]; int kk= k-halfsten*vdir[2];                    
      for (int l = 0; l < order; l++) { 
        for (int nvar = 0; nvar < cls_t::NCONS; nvar++) { 
          flx(i, j, k, nvar) -=  eps2*coefshock(l)*cons(ii,jj,kk,nvar);
          flx(i, j, k, nvar) +=   eps4*coefdamp(l)*cons(ii,jj,kk,nvar);     
        }
        ii +=  vdir[0];jj +=  vdir[1];kk +=  vdir[2];      
      }
    }   
  }

  // .............................................................
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void fluxdissip_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda, const Array4<Real>& flx,
    const cls_t* cls) const {
         
    int ir = i+vdir[0];   int jr = j+vdir[1];   int kr = k+vdir[2];   
    int il = i-vdir[0];   int jl = j-vdir[1];   int kl = k-vdir[2];   
    int ill= il-vdir[0]; int jll = jl-vdir[1]; int kll= kl-vdir[2];   
 
    const int idir = Qdir -1;

    // calculate sensor    
    Real p0,p1,p2,p3;
    Real sen_num= Real(0.0),sen_denom=Real(1.0e-16);
    // loop over sensor variables
    int nv = 0;
    Real sen = Real(0.0);
    for (int l=0;l<NVARSEN;l++)
    { 
      nv = NSEN[l];       
      p0 =  prims(ill,jll,kll,nv);
      p1 =  prims(il,jl,kl,nv);
      p2 =  prims(i,j,k,nv);
      p3 =  prims(ir,jr,kr,nv);    
      sen  = std::max(disconSensor(p0,p1,p2), disconSensor(p1,p2,p3) );
      sen_num += sen*sen;sen_denom +=sen;
      sen = sen_num/sen_denom;
    }

    // reduce order close to BC by making sensor  = 1   
    // sen = (i < mask_sen(idir,1)) ? 1.0 : sen;  
    // sen = (i > mask_sen(idir,2)) ? 1.0 : sen;  
    
    // spectral radius Jacobian matrix (u + c)
    //Real rr = std::max(lambda(il, jl, kl, 0), lambda(i, j, k, 0));    
    Real rr = std::abs( prims(i,j,k,Qdir) ) + prims(i, j, k, cls_t::QC);           
    Real eps2 = Cshock*rr*sen;
    Real eps4 = std::max(0.0, Cdamp*rr - eps2);

    // shock capturing and damping 
    int ii= i-halfsten*vdir[0]; int jj= j-halfsten*vdir[1]; int kk= k-halfsten*vdir[2];                    
    for (int l = 0; l < order; l++) { 
      for (int nvar = 0; nvar < cls_t::NCONS; nvar++) { 
        flx(i, j, k, nvar) -=  eps2*coefshock(l)*cons(ii,jj,kk,nvar);
        flx(i, j, k, nvar) +=   eps4*coefdamp(l)*cons(ii,jj,kk,nvar);     
      }
      ii +=  vdir[0];jj +=  vdir[1];kk +=  vdir[2];      
    }   
  }
  // .............................................................  

  
  // coefficents
  //static const int ordermax = 6;
  typedef Array2D<Real, 0, order, 0, order> arrCoeff_t;
  arrCoeff_t coefskew;
  typedef Array1D<Real, 0, order> arrayNumCoef;
  arrayNumCoef coefdamp,coefshock,coefP;

  int halfsten = order / 2;

  // sensor variables
  const int NVARSEN=2;
  int NSEN[2];
  // masking sensor  
  typedef Array2D<int, 0, AMREX_SPACEDIM, 0, 2> arrIntCoeff_t;
  arrIntCoeff_t mask_sen;

  //-----------------------------------------------------------------------------------
  };


#endif