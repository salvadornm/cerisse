#ifndef VISCOUS_H_
#define VISCOUS_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

#include <Constants.h>
#include <TransPele.h>

#include "diff_ops.H"
template <typename param, typename cls_t>
class viscous_t {

  public:
  
  AMREX_GPU_HOST_DEVICE
  viscous_t() {
    calc_CDcoeffs<param::order>(INTcoef,CDcoef);
  }

  AMREX_GPU_HOST_DEVICE
  ~viscous_t() {}

  // vars accessed by functions 
  int order_sch=param::order;
  typedef Array1D<Real, 0, param::order> arrayNumCoef;
  arrayNumCoef CDcoef,INTcoef;
  int halfsten = param::order / 2;

#if (AMREX_USE_GPIBM || CNS_USE_EB )  
  void inline dflux_ibm(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,            
            const Array4<Real>& rhs, const cls_t* cls,const Array4<bool>& ibMarkers) {
#else
  void inline dflux(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt, 
            const Array4<Real>& rhs, const cls_t* cls) {
#endif

    // mesh sizes
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

    // grid
   // const Box& bx = mfi.tilebox();        
    const Box& bxg = mfi.growntilebox(cls->NGHOST);     // to handle high-order 
    const Box& bxgnodal = mfi.grownnodaltilebox(-1, 0); // to handle fluxes    

    // allocate arrays for transport properties  
    FArrayBox coeffs(bxg, cls_t::NCOEF, The_Async_Arena());
    const int CMU    = cls_t::CMU;
    const int CLAM   = cls_t::CLAM;
    const int CXI    = cls_t::CXI;
    const int CRHOD  = cls_t::CRHOD;
        
    const auto& mu_arr   = coeffs.array(CMU);     // dynamic viscosity
    const auto& lam_arr  = coeffs.array(CLAM);    // thermal conductivity 
    const auto& xi_arr   = coeffs.array(CXI);     // bulk viscosity
    const auto& rhoD_arr = coeffs.array(CRHOD);   // species diffusivity (times rho)

    // pointer to array of transport coefficients    
    const amrex::Array4<const amrex::Real>& coeftrans = coeffs.array();
    
    // calculate all transport properties and store in array (up to ghost points)
#ifdef USE_PELEPHYSICS

    FArrayBox qfab(bxg, cls_t::NPRIM, The_Async_Arena());     // prep space q-arrays
    auto const& q = qfab.array();
    // fill it with prims data
    amrex::ParallelFor( bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
      for (int n=0;n<cls_t::NPRIM; n++){ q(i,j,k,n) = prims(i,j,k,n);}   
    });  
    auto const& q_y   = qfab.const_array(cls_t::QFS);  // array: species mass fraction
    auto const& q_T   = qfab.const_array(cls_t::QT);   // array: temperature 
    auto const& q_rho = qfab.const_array(cls_t::QRHO); // array: density
     
    BL_PROFILE("PelePhysics::get_transport_coeffs()");
    Array4<Real> chi; // dummy Soret effect coef
    
    // temp snm
    //pele::physics::transport::TransportParams< pele::physics::PhysicsType::transport_type> trans_parms;
    trans_parms.allocate(); 

    auto const* ltransparm = trans_parms.device_trans_parm();
    
    amrex::launch(bxg, [=] AMREX_GPU_DEVICE(Box const& tbx) {

            auto trans = pele::physics::PhysicsType::transport();                      
            trans.get_transport_coeffs(tbx, q_y, q_T, q_rho, 
                rhoD_arr, chi, mu_arr,xi_arr, lam_arr, ltransparm);
          });

    // change units
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
        mu_arr(i,j,k)  = mu_arr(i,j,k) *visc_cgs2si;
        lam_arr(i,j,k) = lam_arr(i,j,k)*cond_cgs2si;    
        for (int n=0;n<NUM_SPECIES; n++){        
          rhoD_arr(i,j,k,n) = rhoD_arr(i,j,k,n)*rhodiff_cgs2si;
        }   
        xi_arr(i,j,k)  = xi_arr(i,j,k)*visc_cgs2si;
        });

#else
    amrex::ParallelFor(
        bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
        mu_arr(i,j,k)  = cls->visc(prims(i,j,k,cls_t::QT));
        lam_arr(i,j,k) = cls->cond(prims(i,j,k,cls_t::QT));       
        xi_arr(i,j,k)  = 0.0;
        });
#endif     

    // TEST properties    
    // amrex::ParallelFor( bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {              
    //   printf(" i = %d j=%d k= %d \n");
    //   const Real rho= prims(i,j,k,cls_t::QRHO);
    //   std::cout << " T " << prims(i,j,k,cls_t::QT) << std::endl;
    //   std::cout << " rho " << rho << std::endl;      
    //   for (int n=0;n<NUM_SPECIES; n++){
    //     std::cout << " spec= " << n << " Y " << prims(i,j,k,cls_t::QFS +n) << std::endl;
    //   }
    //   std::cout << " visc " << mu_arr(i,j,k) << std::endl;
    //   std::cout << " cond " << lam_arr(i,j,k) << std::endl;
    //   std::cout << " xi   " << xi_arr(i,j,k) << std::endl;
    //   for (int n=0;n<NUM_SPECIES; n++){        
    //     printf(" spec= %d rhoD= %f D= %f \n",n,rhoD_arr(i,j,k,n),rhoD_arr(i,j,k,n)/rho );
    //   }      
    // });    
    //



    // if (LES)   Pseudo-code for LES
    // {
    //  cls->compute_sgsterms(prims(UX-UZ)) 
    //  muarr+= cls->musgs()
    //  lamarr+=  ls->lamsgs()
    // }
    
    // loop over directions -----------------------------------------------
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
     // GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
      auto const& flx = flxt[dir]->array(); 
  
      // compute diffusion fluxes
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
      amrex::ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                                       
                    this->cns_diff_ibm(i, j, k,dir, prims,flx,coeftrans, dxinv, cls,ibMarkers);
                  });                      
#else    
      amrex::ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {                   
                    this->cns_diff(i, j, k,dir, prims,flx,coeftrans, dxinv, cls);
                  });                  
#endif
        
    }  
    // end loop  ------------------------------------------------------
  } 

 
  // ----------------------------------------------------------------------------------------------  
  /**
  * @brief Compute diffusion fluxes (viscosity + heat + diffusion).
  *
  * @param i,j,k  x, y, z index.cls_t::CLAM
  * @param d1    direction, 0:x, 1:y, 2:z (dir)
  * @param q      primitive variables.
  * @param[out] flx  output diffusion fluxes.
  * @param coeffs transport coefficients.
  * @param dxinv  1/dx
  * @param cls_t  ProbClosures 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff(
      const int i, const int j, const int k, const int d1,
      amrex::Array4<const amrex::Real> const& q,
      amrex::Array4<amrex::Real> const& flx,
      amrex::Array4<const amrex::Real> const& coeffs,
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
      const cls_t* cls) const {
    
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(d1);

    const int d2 = d1 == 0 ? 1 : 0;
    const int d3 = d1 == 2 ? 1 : 2;
    AMREX_D_TERM(const int QU1 = cls_t::QU + d1;,  const int QU2 = cls_t::QU + d2;
               , const int QU3 = cls_t::QU + d3;)
    AMREX_D_TERM(const int UM1 = cls_t::UMX + d1;, const int UM2 = cls_t::UMX + d2;
               , const int UM3 = cls_t::UMX + d3;)

    // Aij = dA_i/dx_j
    const Real dTdn = normal_diff<param::order>(iv, d1, cls_t::QT, q, dxinv);
    const Real u11  = normal_diff<param::order>(iv, d1, QU1, q, dxinv);
#if (AMREX_SPACEDIM >= 2)
    const Real u21  = normal_diff<param::order>(iv, d1, QU2, q, dxinv);
    const Real u12  = tangent_diff<param::order>(iv, d1, d2, QU1, q, dxinv);
    const Real u22  = tangent_diff<param::order>(iv, d1, d2, QU2, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real u31  = normal_diff<param::order>(iv, d1, QU3, q, dxinv);
    const Real u13  = tangent_diff<param::order>(iv, d1, d3, QU1, q, dxinv);
    const Real u33  = tangent_diff<param::order>(iv, d1, d3, QU3, q, dxinv);
#endif
    const Real divu   = AMREX_D_TERM(u11, +u22, +u33);
    
    const Real muf    = interp<param::order>(iv, d1, cls_t::CMU, coeffs);
    const Real xif    = interp<param::order>(iv, d1, cls_t::CXI, coeffs);
    const Real lamf   = interp<param::order>(iv, d1, cls_t::CLAM, coeffs);
    
    AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - (2.0 / 3.0) * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

    // momentum
    AMREX_D_TERM(flx(iv, UM1) -= tau11;, flx(iv, UM2) -= tau12;, flx(iv, UM3) -= tau13;)
   
    // energy
  
    flx(iv, cls_t::UET) -= 0.5 * (AMREX_D_TERM((q(iv, QU1) + q(ivm, QU1)) * tau11,
                                              +(q(iv, QU2) + q(ivm, QU2)) * tau12,
                                              +(q(iv, QU3) + q(ivm, QU3)) * tau13))
                                              + lamf* dTdn;
    //flx(iv, cls_t::UET) -= lamf* dTdn; // temp snm
    

#if NUM_SPECIES > 1    
    // --------------------------------------------------------------------------
    // diffusion species  
    Real ymass[order_sch][NUM_SPECIES],xmole[order_sch][NUM_SPECIES];
    Real hi[order_sch][NUM_SPECIES];
    Real yaux[NUM_SPECIES],xaux[NUM_SPECIES],haux[NUM_SPECIES];

    auto thermo = typename cls_t::multispecies_pele_gas_t();
    
    amrex::IntVect ivp(iv -halfsten*amrex::IntVect::TheDimensionVector(d1));
    for (int l = 0; l < order_sch; l++) {
      ivp +=  amrex::IntVect::TheDimensionVector(d1);
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        yaux[n] = q(ivp,  cls_t::QFS + n); 
      }
      // calculate xmol and specific enthalpies
      thermo.Y2X(yaux, xaux);    
      thermo.RTY2Hi(q(ivp,  cls_t::QRHO), q(ivp,   cls_t::QT), yaux, haux);

      // store in temp array 
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        xmole[l][n] = xaux[n];
        ymass[l][n] = yaux[n]; 
        hi[l][n]    = haux[n];
      }
     
    }
    
    // const Real dpdx  = normal_diff<param::order>(iv, d1, cls_t::QPRES, q, dxinv); 
    // const Real pface = interp<param::order>(iv, d1, cls_t::QPRES, q);
    // const Real dlnp = dpdx/pface;
    Real dlnp = 0.0;  

     
    Real Vc = 0.0;    
    Real Yf[NUM_SPECIES],hf[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; n++) {
      Real Xface = 0.0, Yface = 0.0, hface = 0.0, dXdx = 0.0;
      for (int l = 0; l < param::order; l++) {
        Xface += xmole[l][n]*INTcoef(l);
        Yface += ymass[l][n]*INTcoef(l);
        hface += hi[l][n]*INTcoef(l);
        dXdx  += xmole[l][n]*CDcoef(l);
      }      
      dXdx  /= dxinv[d1];
      Yf[n] = Yface; hf[n] = hface;

      const Real rhoD_f = interp<param::order>(iv, d1, cls_t::CRHOD + n, coeffs); 

      const Real Vd = -rhoD_f * (dXdx + (Xface - Yface) * dlnp);
      Vc += Vd;
      flx(iv, cls_t::UFS + n) += Vd; 
      flx(iv, cls_t::UET)     += Vd * hface;
     }
    // Add correction velocity to fluxes so sum(Vd) = 0
    for (int n = 0; n < NUM_SPECIES; ++n) {       
      flx(iv, cls_t::UFS + n)-= Yf[n] * Vc;
      flx(iv, cls_t::UET)    -= Yf[n] * hf[n] * Vc;
    }

    // --------------------------------------------------------------------------
#endif                                           


      
      
       
  
    }
  // ---------------------------------------------------------------------------------------------  
  /**
  * @brief Compute diffusion fluxes in IB/EB.
  *
  * @param i,j,k  x, y, z index.cls_t::CLAM
  * @param dir    direction, 0:x, 1:y, 2:z.
  * @param q      primitive variables.
  * @param[out] flx  output diffusion fluxes.
  * @param coeffs transport coefficients.
  * @param dxinv  1/dx
  * @param cls_t  ProbClosures 
  * @param marker  geometry markers  (sld and cutcells) 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_ibm(
      const int i, const int j, const int k, const int d1,
      amrex::Array4<const amrex::Real> const& q,
      amrex::Array4<amrex::Real> const& flx,
      amrex::Array4<const amrex::Real> const& coeffs,
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
      const cls_t* cls, const Array4<bool>& marker) const {
    
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(d1);

    const int d2 = d1 == 0 ? 1 : 0;
    const int d3 = d1 == 2 ? 1 : 2;
    AMREX_D_TERM(const int QU1 = cls_t::QU + d1;, const int QU2 = cls_t::QU + d2;
               , const int QU3 = cls_t::QU + d3;)
    AMREX_D_TERM(const int UM1 = cls_t::UMX + d1;, const int UM2 = cls_t::UMX + d2;
               , const int UM3 = cls_t::UMX + d3;)

    // ivm  is iv -1 
    const bool close_to_wall  = marker(iv,1) || marker(ivm,1);      //  flux close to a GP (IBM) or a cut-cell (EB)
    const bool intersolid_flx = marker(iv,0) &&  marker(ivm,0);    // inter-flux

    if (intersolid_flx) return;  // flux =0  inside solid   

    Real u11,dTdn,u21,u12,u22,u31,u13,u33,muf,xif,lamf;
    
    // reduce to second order close to walls
    if (close_to_wall) {        
      constexpr int order = 2;
      dTdn = normal_diff<order>(iv, d1, cls_t::QT, q, dxinv);
      u11  = normal_diff<order>(iv, d1, QU1, q, dxinv);
#if (AMREX_SPACEDIM >= 2)
      u21  = normal_diff<order>(iv, d1, QU2, q, dxinv);
      u12  = tangent_diff<order>(iv, d1, d2, QU1, q, dxinv);
      u22  = tangent_diff<order>(iv, d1, d2, QU2, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
      Real u31  = normal_diff<order>(iv, d1, QU3, q, dxinv);
      Real u13  = tangent_diff<order>(iv, d1, d3, QU1, q, dxinv);
      Real u33  = tangent_diff<order>(iv, d1, d3, QU3, q, dxinv);
#endif  
      // properties
      muf  = interp<order>(iv, d1, cls_t::CMU, coeffs);
      xif  = interp<order>(iv, d1, cls_t::CXI, coeffs);
      lamf = interp<order>(iv, d1, cls_t::CLAM, coeffs);
    }
    else
    {
      dTdn = normal_diff<param::order>(iv, d1, cls_t::QT, q, dxinv);
      u11  = normal_diff<param::order>(iv, d1, QU1, q, dxinv);
#if (AMREX_SPACEDIM >= 2)
      u21  = normal_diff<param::order>(iv, d1, QU2, q, dxinv);
      u12  = tangent_diff<param::order>(iv, d1, d2, QU1, q, dxinv);
      u22  = tangent_diff<param::order>(iv, d1, d2, QU2, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
      u31  = normal_diff<param::order>(iv, d1, QU3, q, dxinv);
      u13  = tangent_diff<param::order>(iv, d1, d3, QU1, q, dxinv);
      u33  = tangent_diff<param::order>(iv, d1, d3, QU3, q, dxinv);
#endif  
      // properties
      muf  = interp<param::order>(iv, d1, cls_t::CMU, coeffs);
      xif  = interp<param::order>(iv, d1, cls_t::CXI, coeffs);
      lamf = interp<param::order>(iv, d1, cls_t::CLAM, coeffs);
    }

    const Real divu   = AMREX_D_TERM(u11, +u22, +u33);
    
    AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - (2.0 / 3.0) * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

    // momentum
    AMREX_D_TERM(flx(iv, UM1) -= tau11;, flx(iv, UM2) -= tau12;, flx(iv, UM3) -= tau13;)
   
    // energy
    flx(iv, cls_t::UET) -= 0.5 * (AMREX_D_TERM((q(iv, QU1) + q(ivm, QU1)) * tau11,
                                              +(q(iv, QU2) + q(ivm, QU2)) * tau12,
                                              +(q(iv, QU3) + q(ivm, QU3)) * tau13)) +
                                              + lamf* dTdn;

  }
  // ---------------------------------------------------------------------------------------------


  }; 

//---------------------------------------------
#endif