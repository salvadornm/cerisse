#ifndef DiffusionCD_H_
#define DiffusionCD_H_

#include <AMReX_FArrayBox.H>

/////////////////////////////////////////////////////////////////////////////
// This template adds diffusion of heat only due to temperatute gradients
/////////////////////////////////////////////////////////////////////////////

template <typename param, typename cls_t>
class diffusiveheat_t {
 public:
  
  //AMREX_GPU_HOST_DEVICE
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
      CDcoef(4)   = - Real(9.0/60); 
      CDcoef(5)   = - Real(1.0/60);                  
      break;
    default: // crash
    
      amrex::Abort("Diffusive order available 2/4/6 MUST specify one of those ");      

      break;
    }
  
  }

  // AMREX_GPU_HOST_DEVICE
  ~diffusiveheat_t() {}

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
    const Box& bxg1 = amrex::grow(bx, 1);

    // allocate arrays
    FArrayBox coeffs(bxg1, 3 + NUM_SPECIES, The_Async_Arena());
    const auto& mu_arr  = coeffs.array(0);
    const auto& lam_arr = coeffs.array(1);

    // calculate all transport properties and store in array
    amrex::ParallelFor(
        bxg1, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          mu_arr(i,j,k)  = cls->visc(prims(i,j,k,cls_t::QT));
          lam_arr(i,j,k) = cls->cond(prims(i,j,k,cls_t::QT));
        });

    // SNM
    // amrex::Print() << " Diffusion Simple order " << order_sch << std::endl;  
    // int ii=0,jj=0,kk=0;
    // Real dt = Real(1e-4);
    // Real Cp = cls->cp; //cls_t::cp
    // Real lam = lam_arr(ii,jj,kk); Real rho = prims(ii,jj,kk,cls_t::QRHO);

    // Real DIF = dxinv[0]*dxinv[0]*dt*lam/(Cp*rho);

    // amrex::Print() << " DIF " << DIF << std::endl;

    // exit(0);

    //
 


  // _>>>>>>
  //  // loop over directions
  //   for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
  //     GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};



    // compute diffusion fluxes in x-direction
    int cdir = 0;
    const Box& xflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        xflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
          this->cns_diff_x(i, j, k, flx, prims, mu_arr, lam_arr,
                           dxinv, *cls);
        });


    amrex::ParallelFor(
        bx, cls_t::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) +=
              dxinv[cdir] * (flx(i, j, k, n) - flx(i + 1, j, k, n));

        });



#if AMREX_SPACEDIM >= 2
    // y-direction
    cdir = 1;
    const Box& yflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        yflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
          this->cns_diff_y(i, j, k, flx, prims, mu_arr,lam_arr,
                           dxinv, *cls);
        });
    amrex::ParallelFor(
        bx, cls_t::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) +=
              dxinv[cdir] * (flx(i, j, k, n) - flx(i, j + 1, k, n));
        });
#endif

#if AMREX_SPACEDIM == 3
    // z-direction
    cdir = 2;
    const Box& zflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        zflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
          this->cns_diff_z(i, j, k, flx, prims, mu_arr,  lam_arr,
                           dxinv, *cls);
        });
    amrex::ParallelFor(
        bx, cls_t::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) +=
              dxinv[cdir] * (flx(i, j, k, n) - flx(i, j, k + 1, n));
        });
#endif
  }

  
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_x(
      int i, int j, int k, amrex::Array4<amrex::Real> const& fx,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      cls_t const& cls) const {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};

  
    Real dTdx = (prims(i,j,k,cls.QT) - prims(i-1,j,k,cls.QT))*dxinv[0];
    Real lam  = 0.5_rt*(lam_arr(i,j,k)+lam_arr(i-1,j,k));
    fx(i, j, k, cls.UET) =   -lam* dTdx;

  }

#if AMREX_SPACEDIM >= 2
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_y(
      int i, int j, int k, amrex::Array4<amrex::Real> const& fy,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      cls_t const& cls) const noexcept {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    
    Real dTdy = (prims(i,j,k,cls.QT) - prims(i,j-1,k,cls.QT))*dxinv[1];
    Real lam  = 0.5_rt*(lam_arr(i,j,k)+lam_arr(i,j-1,k));
    fy(i, j, k, cls.UET) +=   -lam * dTdy;

  }
#endif
#if AMREX_SPACEDIM == 3
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_z(
      int i, int j, int k, amrex::Array4<amrex::Real> const& fz,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      cls_t const& cls) const noexcept {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    constexpr int order = order_sch;

    Real dTdz = normal_diff<order>(iv, 2, cls.QT, prims, dxinv);
    Real lam  = interp<order>(iv, 2, 0, lam_arr) ;
    fz(i, j, k, cls.UET) += - lam * dTdz;
  }
#endif
  

  // private differencing function ------------------------------------------
  /// @brief central difference of given order 
  template <int order>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    // x 2
    Real sum = 0.0; amrex::IntVect iv0 = iv-halfsten*ivn;
    for (int l = 0; l < order; l++) { 
      sum +=  CDcoef(l)*q(iv0+l*ivn,comp);     
    }
    return(2.0_rt*sum)* dxinv[idir];
  }
  
  /// @brief general interpolation with variable order using polynomial
  template <int order>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real interp(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q) noexcept{
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    Real sum = 0.0; amrex::IntVect iv0= iv-halfsten*ivn;
    for (int l = 0; l < order; l++) { 
      sum +=  INTcoef(l)*q(iv0+l*ivn,comp);     
    }
    return(sum);
  }

  /// @brief 2nd order central difference
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff2(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    return (q(iv, comp) - q(iv - ivn, comp)) * dxinv[idir];
  }

  /// @brief 2nd interpolation 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real interp2(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q) noexcept {
    const auto ivn = amrex::IntVect::TheDimensionVector(idir);
    return  (q(iv, comp) + q(iv - ivn, comp));
  }
//  ------------------------------------------



};




#endif
