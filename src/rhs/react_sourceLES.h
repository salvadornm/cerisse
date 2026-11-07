#ifndef REACT_SOURCELES_H
#define REACT_SOURCELES_H

#include <PelePhysics.H>
#include <ReactorBase.H>
#include <Constants.h>
#include <CNSconstants.h>

// use ::  .., reactor_sourceLES<user_source_t, ProbClosures> >;
template <typename source_t, typename cls_t >
class reactor_sourceLES_t {
 public:
  bool m_initialized = false;
  std::unique_ptr<pele::physics::reactions::ReactorBase> m_reactor;

  // factor to multiply reaction (for transients)
  amrex::Real reaction_relax=1.0;
  // max temperature for reaction (to avoid instability in early stages of ATF)
  amrex::Real max_react_temp = 3500.0; 

  // reactor types (hardcoded) and specific options
  inline static constexpr int therm_reactor_type = 1; // 1: U  2:H
  inline static constexpr bool reactor_constant_pressure = false;
  inline static constexpr bool pass_source_term   = false; 
  inline static constexpr bool check_problem_cell = false; 
  /// 

  reactor_sourceLES_t() {
    std::string reactor_type;
    {
      amrex::ParmParse pp("cns");
      pp.get("reactor_type", reactor_type);
      //pp.get("reaction_relax", reaction_relax);
      if (!pp.query("reaction_relax", reaction_relax)) {
        amrex::Print() << " using no relaxation in chem source term  \n ";   
      }
      if (!pp.query("max_react_temp", max_react_temp)) {
        amrex::Print() << " using default max_react_temp = " << max_react_temp << " K \n ";   
      }


    }
    m_reactor = pele::physics::reactions::ReactorBase::create(reactor_type);
    
    if (!m_reactor) {
      amrex::Abort("reactor_sourceLES (): Unknown reactor type " + reactor_type);
    }
    m_reactor->init(therm_reactor_type, 1); // create reactor solver
    m_initialized = true;
  };

  ~reactor_sourceLES_t() {
    if (m_initialized) m_reactor->close();
  }

  /**
   * @brief Calculate chemical reaction source term, adding to the
   * right-hand-side (rhs) array, as well as calling the pass source term (defiend in prob)
   *
   * @tparam cls_t The problem closure class typename.
   * @param mfi    The MFIter object representing the current grid patch.
   * @param prims  The input primitive variables array.
   * @param rhs    The output array where the updated right-hand-side will be
   * stored.
   * @param cls    The problem closure object (for indicies).
   * @param dt     The time step size. (react() requires it to be non-const)
   */

#if (AMREX_USE_GPIBM || CNS_USE_EB )     
  void inline src(const Geometry& geomdata, const amrex::MFIter& mfi,
                  const amrex::Array4<const amrex::Real>& prims,
                  const amrex::Array4<amrex::Real>& rhs, const cls_t* cls_d,
                  amrex::Real dt, amrex::Real real_time, const Array4<uint8_t>& marker) {
#else
  void inline src(const Geometry& geomdata, const amrex::MFIter& mfi,
                  const amrex::Array4<const amrex::Real>& prims,
                  const amrex::Array4<amrex::Real>& rhs, const cls_t* cls_d,
                  amrex::Real dt, amrex::Real real_time) {
#endif


    if (!m_initialized) amrex::Abort("reactor_t not initialised");

    BL_PROFILE("reactor_sourceLES_t::src()");

    // put here because this is a .h file
    using amrex::Array4;
    using amrex::Box;
    using amrex::FArrayBox;
    using amrex::IArrayBox;
    using amrex::Real;

    const Box bx = mfi.tilebox();
    // mesh sizes
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geomdata.InvCellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dx = geomdata.CellSizeArray(); 

       

    // TODO: stochastic fields indexing

    ///////////////////// Prepare for react /////////////////////
    FArrayBox tempf(bx, 2 * NUM_SPECIES + 4, The_Async_Arena());
    auto const& rY = tempf.array(0);
    auto const& rEi = tempf.array(NUM_SPECIES);
    auto const& T = tempf.array(NUM_SPECIES + 1);
    auto const& rYsrc  = tempf.array(NUM_SPECIES + 2);
    auto const& rEisrc = tempf.array(2 * NUM_SPECIES + 2);
    auto const& fc     = tempf.array(2 * NUM_SPECIES + 3);  // number of RHS eval (not used)
    IArrayBox maskf(bx, 1, The_Async_Arena());
    maskf.setVal<RunOn::Gpu>(1);

    IArrayBox maskf_noreact(bx, 1, The_Async_Arena());
    amrex::Gpu::DeviceScalar<int> skip_react(0);
    int* skip_react_ptr = skip_react.dataPtr();
    auto const& mask_noreact = maskf_noreact.array();


    auto const& mask = maskf.array();  // 1: do reaction, -1: skip reaction

    const Real o_dt = 1.0 / dt;

    
    // parameters pased by prob.h 
    constexpr bool do_react       = source_t::do_reactions;  
    constexpr bool mask_closewall = source_t::mask_cells_boundary;  

    // parameters pased by input
    const Real maxT = max_react_temp;
    const Real relax = reaction_relax;

    // models defined by source_t  (ATF, PaSR, LES) if not defined, default to false
    constexpr bool useATF = []{
    if constexpr (requires { source_t::use_ATF; })
        return source_t::use_ATF;
    else
        return false;
    }();

    constexpr bool usePaSR = []{
    if constexpr (requires { source_t::use_PaSR; })
        return source_t::use_PaSR;
    else
        return false;
    }();

    constexpr bool useLES = []{
    if constexpr (requires { source_t::use_LES; })
        return source_t::use_LES;
    else
        return false;
    }();
    //

    // use cls_d to compute sensor and filter width for ATF and LES models (if needed)
    const auto& cls = *cls_d; 

    // prepare input -------------------------------------------------------------------------
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

      // [rY, rEi, T, rYsrc, rEisrc] convert to CGS  to enter PelePhysics

      Real rho = prims(i, j, k, cls_t::QRHO);
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        rY(i, j, k, ns)    =   rho* prims(i, j, k, cls_t::QFS + ns) * rho_si2cgs;        
      }      
      rEi(i, j, k) = rho * prims(i, j, k, cls_t::QEINT) * rhoenergy_si2cgs;
      T(i, j, k) = prims(i, j, k, cls_t::QT);
      // Enthalpy (if reactor_type 2)
      if (therm_reactor_type==2){
        rEi(i,j,k)  +=  prims(i, j, k, cls_t::QPRES)*pres_si2cgs;   // rEi stores rhoHi
      }          
      
      // communicate source terms to solver
      if (pass_source_term){
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          rYsrc(i, j, k, ns) = rhs(i, j, k, cls_t::UFS + ns) * rho_si2cgs;
        }        
        Real mx = rho * prims(i, j, k, cls_t::QU);
        Real my = rho * prims(i, j, k, cls_t::QV);
        Real mz = rho * prims(i, j, k, cls_t::QW);
        Real rke = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          rho += rhs(i, j, k, cls_t::UFS + ns) * dt;
        }
        mx += rhs(i, j, k, cls_t::UMX) * dt;
        my += rhs(i, j, k, cls_t::UMY) * dt;
        mz += rhs(i, j, k, cls_t::UMZ) * dt;
        Real rke_new = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
        rEisrc(i, j, k) = (rhs(i, j, k, cls_t::UET) - (rke_new - rke) * o_dt ) * rhoenergy_si2cgs;
      }      
      else {
        rEisrc(i, j, k) = 0.0;
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          rYsrc(i, j, k, ns) = 0.0;       
        } 
      }     
      // --------  MASK -------------
      // fill mask  with max/min temperature and solid boundaries (if EB or GPIBM)    
      mask(i, j, k) = (T(i, j, k) > CNSConstants::min_react_temp) ? 1 : -1;
      if (T(i,j,k) > maxT) mask(i,j,k)  = -1;
#if (AMREX_USE_GPIBM || CNS_USE_EB )        
      mask(i, j, k) = marker(i, j, k, 0) ? -1 : mask(i, j, k);
      //remove cells close to solid from chemistry (input by prob)
      if (mask_closewall) { if (marker(i, j, k, 1)) mask(i,j,k) = -1; }
#endif
           
    });

    /////////////////////////// React ///////////////////////////
    Real current_time = 0.0;

    if (do_react) {
    // Not necessary to start a stream here, however pelePhysics function only takes a stream.     
#ifdef AMREX_USE_GPU
    m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time,
                     amrex::Gpu::gpuStream());
#else
    m_reactor->react(bx, rY, rYsrc, T, rEi, rEisrc, fc, mask, dt, current_time);
#endif
    amrex::Gpu::Device::streamSynchronize();  // Important
      } // end do react

    // Convert to SI units
    ParallelFor(bx, [rY,rEi] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        for (int n = 0; n < NUM_SPECIES; n++) {
          rY(i,j,k,n) *= rho_cgs2si;
        }
        rEi(i,j,k) *= rhoenergy_cgs2si;
      });
    /////////////////////////////////////////////////////////////
   
    // pre-multiplying factor array
    FArrayBox turbcombf(bx, 1, The_Async_Arena());
    auto const& wfactor = turbcombf.array(0);
    amrex::ParallelFor(
      bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {      
          wfactor(i,j,k) = relax; // initialize wfactor to 1 (no modification )
    });      
    // LES
    Real Delta;
    if constexpr(useLES) {
      Delta = cls.calc_delta(dx);
    } else {
      Delta = dx[0];
    }
    /// ATF
    if constexpr(useATF) {                      
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {      
          const auto& cls = *cls_d;
          const Real omega = cls.flame_sensor(i,j,k,prims);  // compute sensor for ATF
          const Real F     = cls.thickening(omega); 
          // ATF models: (1) Classic Ducros (2) Rathore transformation (default)
          if (source_t::ATF_model == 1) {
            // compute wrinkling  only inside flame
            Real E = 1.0; 
            if (F > 1.001) { // avoid computing u_sgs when F ~ 1 (outside flame)
              const Real usgs  = cls.usgs_cell(i,j,k,prims, dxinv, Delta); // compute u_sgs for model 1
              E     = cls.efficiency(usgs,Delta);                                        
            } 
            wfactor(i,j,k) *= E/F;    // modify reaction source term by ATF factor and efficiency             
          }
          else {                    
            wfactor(i,j,k) *= 1.0/F;  // modify reaction source term by ATF factor             
          }
      }); 
    }     
    /// PASR
    if constexpr(usePaSR) {                      
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {      
          const auto& cls = *cls_d;
          const Real tau_chem = cls.tau_chem(i,j,k,rY,prims,prims(i, j, k, cls_t::QRHO),dt);
          const Real tau_mix  = cls.tau_sgs(i,j,k,prims,dxinv,Delta);
          const Real Dasgs = tau_mix/tau_chem; // Damkohler number based on sgs mixing time scale
          wfactor(i,j,k)  *= 1.0/(1.0 + Dasgs); // modify reaction source term by PaSR efficiency
        }); 
    }

    //////////////////////// Unpack dat + Update RHS ////////////////////////
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const auto& cls = *cls_d;

      // only update rhs if mask != -1 (valid cell for reaction)
      if (mask(i, j, k) != -1) {               
        
        const Real rr = wfactor(i,j,k); // reaction rate multiplyier (relaxation, comb models)

        const Real rho = prims(i, j, k, cls_t::QRHO);          
        // Option 1: constant pressure ----------- (h,P constant)
        if (reactor_constant_pressure) {
          Real Yk[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; ++n) { Yk[n] = rY(i, j, k, n)/rho;}            
          // enthalpy
          const Real h = prims(i, j, k, cls_t::QEINT) + prims(i,j,k,cls_t::QPRES)/rho;
          // recalculate Temperature
	        cls.RHY2T(rho, h, Yk, T(i,j,k));
          // recalculate rho based on new T, Y and P (unchanged during reaction)
          Real rhonew;
          cls.PYT2R(prims(i,j,k,cls_t::QPRES),Yk,T(i,j,k),rhonew);
          // calculate density change drho/dt
          const Real drhodt = (rhonew - rho)* o_dt;
          // mass rhs:   rho dY/dt + Y drho/dt
          for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            rhs(i, j, k, cls_t::UFS + ns) +=  rho*(Yk[ns]- prims(i, j, k, cls_t::QFS + ns))* o_dt*rr;
            rhs(i, j, k, cls_t::UFS + ns) +=  prims(i, j, k, cls_t::QFS + ns)* drhodt*rr;                         
          }
          // energy rhs:  h drho/dt
          rhs(i, j, k, cls_t::UET) +=  h*drhodt;
        }
        else {
        // Option 2: constant volume ----------- (rho, e constant)          
          for (int ns = 0; ns < NUM_SPECIES; ++ns) {                          
            Real Wchem = (rY(i, j, k, ns) - rho * prims(i, j, k, cls_t::QFS + ns)) * o_dt;
            if (pass_source_term){
              rhs(i, j, k, cls_t::UFS + ns) = Wchem*rr;
            }
            else {
              rhs(i, j, k, cls_t::UFS + ns) += Wchem*rr;
            }
          }
        }    
      } //end mask      
    });    
    ///

    // clear memory
    tempf.clear();
    turbcombf.clear();
     
    // call user source term (passed as argument)
    //  - assume source_t is a user_source_t is lightweight (no persistent state, just logic),
    source_t{}.rsrc(geomdata,mfi, prims, rhs, cls_d, dt, real_time);

  }
};

#endif
