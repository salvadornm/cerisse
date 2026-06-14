#ifndef CNS_H_
#define CNS_H_

#include <AMReX_AmrLevel.H>
#if CNS_USE_EB
#include <AMReX_EBFluxRegister.H>
using FluxReg = amrex::EBFluxRegister;
#else
#include <AMReX_YAFluxRegister.H>
#include <AMReX_Math.H>
using FluxReg = amrex::YAFluxRegister;
#endif
#include <prob.h>
#include <CNSconstants.h>

#include <Utilities.h>


// using namespace amrex;

class CNS : public amrex::AmrLevel {
 public:
  // Init --------------------------------------------------------------------
  CNS();
  CNS(amrex::Amr& papa, int lev, const amrex::Geometry& level_geom,
      const amrex::BoxArray& bl, const amrex::DistributionMapping& dm,
      amrex::Real time);
  ~CNS();

  CNS(const CNS& rhs) = delete;
  CNS& operator=(const CNS& rhs) = delete;

  // Read parameters
  static void read_params();

  // Define data descriptors.
  static void variableSetUp();

  // Cleanup data descriptors at end of run.
  static void variableCleanUp();

  // Initialize data on this level from another CNS (during regrid).
  void init(amrex::AmrLevel& old) override;

  // Initialize data on this level after regridding if old level did not
  // previously exist
  void init() override;

  // Initialize grid data at problem start-up.
  virtual void initData() override;

  // Do work after init().
  virtual void post_init(amrex::Real stop_time) override;
  // -------------------------------------------------------------------------

  // Time-stepping -----------------------------------------------------------
  void compute_rhs(amrex::MultiFab& S, amrex::Real dt,
                   FluxReg* fr_as_crse, FluxReg* fr_as_fine);

#if NUM_SPECIES > 1                   
  void clip_species_state(amrex::MultiFab& S);                   
#endif  

  // void computeTemp(amrex::MultiFab& State, int ng);

  GpuArray<Real,AMREX_SPACEDIM> maxEigen();

  // Compute initial time step.
  amrex::Real initialTimeStep();

  void computeInitialDt(int finest_level, int sub_cycle,
                        amrex::Vector<int>& n_cycle,
                        const amrex::Vector<amrex::IntVect>& ref_ratio,
                        amrex::Vector<amrex::Real>& dt_level,
                        amrex::Real stop_time) override;

  void computeNewDt(int finest_level, int sub_cycle,
                    amrex::Vector<int>& n_cycle,
                    const amrex::Vector<amrex::IntVect>& ref_ratio,
                    amrex::Vector<amrex::Real>& dt_min,
                    amrex::Vector<amrex::Real>& dt_level, amrex::Real stop_time,
                    int post_regrid_flag) override;

  // Advance grids at this level in time.
  Real advance(amrex::Real time, amrex::Real dt, int iteration,
               int ncycle) override;

  // Do work after timestep().
  virtual void post_timestep(int iteration) override;

  virtual void postCoarseTimeStep(Real time) override;

  virtual void post_restart() override;

  void set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint) override;

  // -------------------------------------------------------------------------

  // Gridding ----------------------------------------------------------------
  virtual void post_regrid(int lbase, int new_finest) override;

#ifdef AMREX_USE_GPIBM
  void rebuildIBM();
#endif

  // Error estimation for regridding.
  // virtual void errorEst (int lev, TagBoxArray& tags, Real time, int ngrow);
  virtual void errorEst(amrex::TagBoxArray& tb, int clearval, int tagval,
                        amrex::Real time, int n_error_buf = 0,
                        int ngrow = 0) override;

  // init
  CNS& getLevel(int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

  enum StateDataType { State_Type = 0, Stats_Type, Cost_Type };

  void buildMetrics();

  static AMREX_FORCE_INLINE void rz_sanity_check(amrex::Geometry const& geom)
  {
    // RZ axis sanity check: for RZ the axis must be at r=0 and r-direction
    // cannot be periodic.
    if (geom.IsRZ()) {
#if (AMREX_SPACEDIM != 2)
      amrex::Abort("RZ requires AMREX_SPACEDIM=2 (axisymmetric r-z)");
#endif
      if (geom.isPeriodic(0)) {
        amrex::Abort(
            "RZ requires geometry.is_periodic[0]=0 (non-periodic r-direction)");
      }

      const amrex::Real rlo = geom.ProbLo(0);
      if (amrex::Math::abs(rlo) > amrex::Real(1.e-14)) {
        amrex::Abort("RZ requires geometry.prob_lo[0]=0 (axis at r=0)");
      }
    }
  }

  int okToContinue() override;

  void avgDown();

  void printTotal() const;

  virtual void writePlotFile(const std::string& dir, std::ostream& os,
                             VisMF::How how = VisMF::NFiles) override;

  virtual void writePlotFilePost(const std::string& dir,
                                 std::ostream& os) override;

#if AMREX_USE_GPIBM
  virtual void writeSurfFile();
#endif

  // diagnostics
  static bool record_probe;
  void setupTimeProbe();
  void recordTimeProbe();
  static int time_probe_lev;
  static int time_probe_int;
  static amrex::Vector<std::string> time_probe_names;
  static amrex::Vector<std::string> time_probe_reductions;
  static amrex::Vector<amrex::Box> time_probe_boxes;

  // Parameters
  static int num_state_data_types;
  FluxReg flux_reg;
  static int do_reflux;

  static bool verbose;
  // static amrex::IntVect hydro_tile_size;
  static amrex::Real cfl;

  static int refine_max_dengrad_lev;
  static amrex::Real refine_dengrad;

  // Statistics
  static amrex::Real time_stats;
  static amrex::Real time_stat_level[10];
  static bool compute_stats, record_stats;
  static int INDEX_THERM;
  void setupStats();
  void computeStats();
  void debugStats(int is, int js, int ks);

  static amrex::Real gravity;

  static amrex::Real dt_constant;
  static bool dt_dynamic;
  static int nstep_screen_output;
  static int dist_linear;
  static int order_rk;
  static int stages_rk;

  // When true, the end-of-step IBM abort check uses "approaching the
  // smallr / ei_min clipping floors" as the failure criterion instead of
  // "already non-positive". Catches silent clipping in cons2prims that
  // would otherwise mask a numerical breakdown.
  static bool strict_positivity;

  // Pass 2 (flood-fill interior solid + zero momentum) always runs for FSI.
  // For static geometry, it is opt-in: some complex-geometry / multi-body
  // cases are actually *destabilised* by the flood-fill because the averaged
  // neighbour values spread post-shock states into the body, which the next
  // step's WENO stencil then reads back and oscillates on. Enable via
  // cns.pass2_static = 1 only when you explicitly want that behaviour.
  static bool pass2_static;

  // Utility-variables
  static bool use_utility;
  static Utility utilidades;

  // IBM-specific keywords
  static bool ib_move;
  static bool plot_surf;
  static int  surf_int;
  static std::string surf_filename;

  // EB-specific keywords
  static amrex::Real eb_weight;
  static bool eb_redistribution;
  static std::string eb_redistribution_type;

  // NSCBC-specific keywords  
  static bool use_nscbc;
  static int nscbc_order;
  static amrex::GpuArray<int, AMREX_SPACEDIM> nscbc_lo;
  static amrex::GpuArray<int, AMREX_SPACEDIM> nscbc_hi;
  struct NSCBCParm {
    amrex::Real Lchar   = 1.0;
    amrex::Real Mmax    = 0.1;
    amrex::Real Ptarget = 101325.0;
    amrex::Real sigma   = 0.28;

    amrex::Real utarget = 0.0;
    amrex::Real vtarget = 0.0;
    amrex::Real wtarget = 0.0;
    amrex::Real Ttarget = 300.0;
    amrex::Real eta     = 1.0;
  };

  static NSCBCParm h_nscbc_parm;
  static NSCBCParm* d_nscbc_parm;

  // LES-variables
  static bool use_LES;

 public:
  PROB::ProbRHS prob_rhs{};   // per-level RHS object (Euler + diffusive + source functors)
  static PROB::ProbClosures* h_prob_closures;
  static PROB::ProbClosures* d_prob_closures;
  static PROB::ProbParm* h_prob_parm;       // host-resident objects used on CPU and as sources for copies
  static PROB::ProbParm* d_prob_parm;       // device-resident objects used on GPU (copied from host at initialization)
  static BCRec* h_phys_bc;
  static BCRec* d_phys_bc;
};

void cns_bcfill(amrex::Box const& bx, amrex::FArrayBox& data, const int dcomp,
                const int numcomp, amrex::Geometry const& geom,
                const amrex::Real time, const amrex::Vector<amrex::BCRec>& bcr,
                const int bcomp, const int scomp);

// declare main IB class instance
#ifdef AMREX_USE_GPIBM  
namespace IBM{
  inline PROB::ProbIB ib;
}
#endif

// declare main EB class instance
#ifdef CNS_USE_EB
namespace EBM{
  inline PROB::ProbEB eb;
}
#endif


#endif