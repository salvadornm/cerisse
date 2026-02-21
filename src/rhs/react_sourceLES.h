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
  Real reaction_relax=1.0;

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
      pp.get("reaction_relax", reaction_relax);

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
      
    // TODO: do not work in fine covered box ??
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

    
    // parameters pased by prob.h (only in LES)
    constexpr bool do_react       = source_t::do_reactions;  
    constexpr bool mask_closewall = source_t::mask_cells_boundary;  

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
        rEisrc(i, j, k) = (rhs(i, j, k, cls_t::UET) - (rke_new - rke) / dt) * rhoenergy_si2cgs;
      }      
      else {
        rEisrc(i, j, k) = 0.0;
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          rYsrc(i, j, k, ns) = 0.0;       
        } 
      }     
      // fill mask      
      mask(i, j, k) = (T(i, j, k) > CNSConstants::min_react_temp) ? 1 : -1;

      // mask solid boundaries
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

    /// Compute LES properties
    // if (LES)
    // {
    //   // do stuff compute taus sgs, Efficiency ...
    // }


    // relaxation factor
    const Real rr = reaction_relax;


    // ATF options    (by default no ATF: o_F=1)
    constexpr amrex::Real Fthick = (source_t::ATF ? source_t::thickfactor : amrex::Real(1.0));
    const amrex::Real o_F = rr / Fthick;


    //////////////////////// Unpack dat + Update RHS ////////////////////////
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const auto& cls = *cls_d;

      // only update rhs if mask != -1 (valid cell for reaction)
      if (mask(i, j, k) != -1) {
        
        bool problem_cell = false;      

        // check for problematic cell after reaction
        if (check_problem_cell) {
          for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            problem_cell |=
              (rY(i, j, k, ns) < -1e-5 || rY(i, j, k, ns) > 1.0 + 1e-5 ||
               std::isnan(rY(i, j, k, ns)));
          }
          problem_cell =  problem_cell || (T(i, j, k) <= 0.0);        
        }

        // Monitor problem cell, do not add reaction source         
        if (!problem_cell)
        {
          // Update species source terms ----------------------
          const Real rho = prims(i, j, k, cls_t::QRHO);

          // Option 1: constant pressure ----------- (h,P constant)
          if (reactor_constant_pressure) 
          {
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
              rhs(i, j, k, cls_t::UFS + ns) +=  rho*(Yk[ns]- prims(i, j, k, cls_t::QFS + ns))* o_dt;
              rhs(i, j, k, cls_t::UFS + ns) +=  prims(i, j, k, cls_t::QFS + ns)* drhodt;                         
            }
            // energy rhs:  h drho/dt
            rhs(i, j, k, cls_t::UET) +=  h*drhodt;
          }
          else
          // Option 2: constant volume ----------- (rho, e constant)
          {
            for (int ns = 0; ns < NUM_SPECIES; ++ns) {                          
              Real Wchem = (rY(i, j, k, ns) - rho * prims(i, j, k, cls_t::QFS + ns)) * o_dt;
              if (pass_source_term){
                rhs(i, j, k, cls_t::UFS + ns) = Wchem*o_F;
              }
              else {
                rhs(i, j, k, cls_t::UFS + ns) += Wchem*o_F;
              }
            }
          }
      

        }
      }
    });
    ///



    //


    // clear memory
    tempf.clear();

     
     // if LES multiply by something


    // call user source term (passed as argument)
    //  - assume source_t is a user_source_t is lightweight (no persistent state, just logic),
    source_t{}.rsrc(geomdata,mfi, prims, rhs, cls_d, dt, real_time);

  }
};

#endif
