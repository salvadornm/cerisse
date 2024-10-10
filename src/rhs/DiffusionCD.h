#ifndef DiffusionCD_H_
#define DiffusionCD_H_

#include <AMReX_FArrayBox.H>

/////////////////////////////////////////////////////////////////////////////
// This template adds diffusion of heat only due to temperatute gradients
/////////////////////////////////////////////////////////////////////////////

template <typename param, typename cls_t>
class diffusiveheat_t {
 public:
  
  AMREX_GPU_HOST_DEVICE
  diffusiveheat_t() {

  switch (order_sch)
    {
    case 2:
      // interpolation
      INTcoef(0) =  Real(0.5);
      INTcoef(1) =  Real(0.5);
      // derivatives coefficients
      CDcoef(0)   = -Real(1.0/2.0);
      CDcoef(1)   = Real(1.0/2.0);      
      break;
    case 4:
      // interpolation
      INTcoef(0) = -Real(1.0/12.0);
      INTcoef(1) =  Real(7.0/12.0);
      INTcoef(2) =  Real(7.0/12.0);
      INTcoef(3) = -Real(1.0/12.0);   
      // derivatives coefficients
      CDcoef(0)   =   Real(1.0/12);
      CDcoef(1)   = - Real(8.0/12);      
      CDcoef(2)   = + Real(8.0/12);      
      CDcoef(3)   = - Real(1.0/12);      
      break;
    case 6: 
      // interpolation 
      INTcoef(0) = Real(1.0/60.0);
      INTcoef(1) = -Real(8.0/60.0);
      INTcoef(2) = Real(37.0/60.0);
      INTcoef(3) = Real(37.0/60.0);        
      INTcoef(4) = -Real(8.0/60.0);
      INTcoef(5) = Real(1.0/60.0);
      // derivatives coefficients
      CDcoef(0)   =   Real(1.0/60);
      CDcoef(1)   = - Real(9.0/60);      
      CDcoef(2)   = + Real(45.0/60);      
      CDcoef(3)   = - Real(45.0/60); 
      CDcoef(4)   = + Real(9.0/60); 
      CDcoef(5)   = - Real(1.0/60);                  
      break;
    default: // crash
    
      amrex::Abort("Diffusive order available 2/4/6 MUST specify one of those ");      

      break;
    }
    
    // the derivatives are in the faces, so coefficients are doubled
    for(int l=0;l<param::order;l++) {CDcoef(l)=2.0_rt*CDcoef(l);}
      
  }

  AMREX_GPU_HOST_DEVICE
  ~diffusiveheat_t() {}

  // vars accessed by functions 
  int order_sch=param::order;
  typedef Array1D<Real, 0, param::order> arrayNumCoef;
  arrayNumCoef CDcoef,INTcoef;
  int halfsten = param::order / 2;


#ifdef AMREX_USE_GPIBM   
  void inline dflux_ibm(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, const Array4<Real>& flx,
            const Array4<Real>& rhs, const cls_t* cls,const Array4<bool>& ibMarkers) {
#else
  void inline dflux(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, const Array4<Real>& flx,
            const Array4<Real>& rhs, const cls_t* cls) {
#endif

    // mesh sizes
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

    // grid
    const Box& bx = mfi.tilebox();    
    //const Box& bxg1 = amrex::grow(bx, 1);
    const Box& bxg = mfi.growntilebox(cls->NGHOST);     // to handle high-order 
    const Box& bxgnodal = mfi.grownnodaltilebox(-1, 0); // to handle fluxes

    // allocate 2 arrays for transport properties
    FArrayBox coeffs(bxg, 2, The_Async_Arena());
    const auto& mu_arr  = coeffs.array(0);
    const auto& lam_arr = coeffs.array(1);

    // calculate all transport properties and store in array (up to ghost points)
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          mu_arr(i,j,k)  = cls->visc(prims(i,j,k,cls_t::QT));
          lam_arr(i,j,k) = cls->cond(prims(i,j,k,cls_t::QT));          
        });
    
    // loop over directions -----------------------------------------------
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      // computr diffusion fluxes
#ifdef AMREX_USE_GPIBM  
      amrex::Abort(" IBM +diffusion only No ready yet")
      amrex::ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
                    this->flux_diffdir_ibm(i, j, k,dir,vdir, flx, prims, mu_arr, lam_arr,
                          dxinv, cls,ibMarkers);
                  });
#else    
      amrex::ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
                    this->flux_diffdir(i, j, k,dir,vdir, flx, prims, mu_arr, lam_arr,
                          dxinv, cls);
                  });                  
#endif
  
      // add flux derivative to rhs, i.e.  rhs + = (flx[i] - flx[i+1])/dx
      amrex::ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });

    }  

    // end loop  ------------------------------------------------------
  }
  // @brief subroutine to compute heat flux lamf*dT/dx depend direction  
  // ............................................................. 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_diffdir(
      int i, int j, int k, int dir,const GpuArray<int, 3>& vdir,
      amrex::Array4<amrex::Real> const& fx,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      const cls_t* cls) const {

    using amrex::Real;    

    Real lam = 0.0_rt; Real dTdr= 0.0_rt;
    int il= i-halfsten*vdir[0]; int jl= j-halfsten*vdir[1]; int kl= k-halfsten*vdir[2];     
    for (int l = 0; l < order_sch; l++) { 
      lam  += INTcoef(l)*lam_arr(il,jl,kl);     
      dTdr += CDcoef(l)*prims(il,jl,kl,cls_t::QT);
      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];     
    }  
    fx(i, j, k, cls_t::UET) +=   -lam* dTdr*dxinv[dir];
  }
};




#endif
