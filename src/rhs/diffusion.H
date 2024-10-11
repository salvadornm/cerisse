#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

#include "diff_ops.H"


template <typename param, typename cls_t>
class viscous_t {

  public:
  
  AMREX_GPU_HOST_DEVICE
  diffusiveheat_t() {}

  AMREX_GPU_HOST_DEVICE
  ~diffusiveheat_t() {}

  // vars accessed by functions 
  int order_sch=param::order;

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
    const Box& bxg = mfi.growntilebox(cls->NGHOST);     // to handle high-order 
    const Box& bxgnodal = mfi.grownnodaltilebox(-1, 0); // to handle fluxes

    // allocate arrays for transport properties
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
    
#else    
      amrex::ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
                    this->cns_diff(i, j, k,dir, prims,flx,coeffs, dxinv, cls);
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

  
  /**
  * @brief Compute diffusion fluxes.
  *
  * @param i,j,k  x, y, z index.
  * @param dir    direction, 0:x, 1:y, 2:z.
  * @param q      primitive variables.
  * @param[out] flx  output diffusion fluxes.
  * @param coeffs transport coefficients.
  * @param dxinv  1/dx
  * @param cls_t  ProbClosures 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
  cns_diff(const int i, const int j, const int k, const int d1,
        amrex::Array4<const amrex::Real> const& q,
        amrex::Array4<amrex::Real> const& flx,
        amrex::Array4<const amrex::Real> const& coeffs,
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
        const cls_t* cls) noexcept
    {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(d1);

    const int d2 = d1 == 0 ? 1 : 0;
    const int d3 = d1 == 2 ? 1 : 2;
    AMREX_D_TERM(const int QU1 = cls_t::QU + d1;, const int QU2 = cls_t::QU + d2;
               , const int QU3 = cls_t::QU + d3;)
    AMREX_D_TERM(const int UM1 = cls_t::UMX + d1;, const int UM2 = cls_t::UMX + d2;
               , const int UM3 = cls_t::UMX + d3;)

    constexpr int order = 4;

    // Aij = dA_i/dx_j
    const Real dTdn = normal_diff<order>(iv, d1, cls_t::QT, q, dxinv);
    const Real u11  = normal_diff<order>(iv, d1, QU1, q, dxinv);
#if (AMREX_SPACEDIM >= 2)
    const Real u21  = normal_diff<order>(iv, d1, QU2, q, dxinv);
    const Real u12  = tangent_diff<order>(iv, d1, d2, QU1, q, dxinv);
    const Real u22  = tangent_diff<order>(iv, d1, d2, QU2, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real u31  = normal_diff<order>(iv, d1, QU3, q, dxinv);
    const Real u13  = tangent_diff<order>(iv, d1, d3, QU1, q, dxinv);
    const Real u33  = tangent_diff<order>(iv, d1, d3, QU3, q, dxinv);
#endif
  const Real divu   = AMREX_D_TERM(u11, +u22, +u33);

  // CHECK  ....  need CMU,CXI,CLAM
  const Real muf    = interp<order>(iv, d1, CMU, coeffs);
  const Real xif    = interp<order>(iv, d1, CXI, coeffs);
  //

  AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - (2.0 / 3.0) * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

  AMREX_D_TERM(flx(iv, UM1) -= tau11;, flx(iv, UM2) -= tau12;, flx(iv, UM3) -= tau13;)
  flx(iv, cls_t::UET) -= 0.5 * (AMREX_D_TERM((q(iv, QU1) + q(ivm, QU1)) * tau11,
                                            +(q(iv, QU2) + q(ivm, QU2)) * tau12,
                                            +(q(iv, QU3) + q(ivm, QU3)) * tau13) +
                           (coeffs(iv, CLAM) + coeffs(ivm, CLAM)) * dTdn);

  // Species transport
  cns_diff_species(iv, d1, q, coeffs, dxinv, flx);
}

  

//---------------------------------------------
#endif