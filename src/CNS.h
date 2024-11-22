#ifndef CNS_H_
#define CNS_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <prob.h>

using namespace amrex;

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
                   amrex::FluxRegister* fr_as_crse,
                   amrex::FluxRegister* fr_as_fine);

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

  // -------------------------------------------------------------------------

  // Gridding ----------------------------------------------------------------
  virtual void post_regrid(int lbase, int new_finest) override;

  // Error estimation for regridding.
  // virtual void errorEst (int lev, TagBoxArray& tags, Real time, int ngrow);
  virtual void errorEst(amrex::TagBoxArray& tb, int clearval, int tagval,
                        amrex::Real time, int n_error_buf = 0,
                        int ngrow = 0) override;

  // init
  CNS& getLevel(int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

  // enum StateVariable {
  //     Density = 0, Xmom, Ymom, Zmom, Etot
  // };


  enum StateDataType { State_Type = 0, Stats_Type, Cost_Type };

  void buildMetrics();

  void avgDown();

  void printTotal() const;

  virtual void writePlotFile(const std::string& dir, std::ostream& os,
                             VisMF::How how = VisMF::NFiles) override;

  virtual void writePlotFilePost(const std::string& dir,
                                 std::ostream& os) override;

  // diagnostics
  static bool record_probe;
  void setupTimeProbe();
  void recordTimeProbe();
  static int time_probe_lev;
  static int time_probe_int;
  static amrex::Vector<std::string> time_probe_names;
  static amrex::Vector<amrex::Box> time_probe_boxes;


  // Parameters
  static int num_state_data_types;
  std::unique_ptr<amrex::FluxRegister> flux_reg;
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
  static bool ib_move;
  static bool plot_surf;

 public:
  /*static inline*/ PROB::ProbRHS prob_rhs;
  static PROB::ProbClosures* h_prob_closures;
  static PROB::ProbClosures* d_prob_closures;
  static PROB::ProbParm* h_prob_parm;
  static PROB::ProbParm* d_prob_parm;
  static BCRec* h_phys_bc;
  static BCRec* d_phys_bc;
};

void cns_bcfill(amrex::Box const& bx, amrex::FArrayBox& data, const int dcomp,
                const int numcomp, amrex::Geometry const& geom,
                const amrex::Real time, const amrex::Vector<amrex::BCRec>& bcr,
                const int bcomp, const int scomp);

#ifdef AMREX_USE_GPIBM
  // declare main IB class instance
namespace IBM{
  inline PROB::ProbIB ib;
}
#endif

#ifdef CNS_USE_EB
  // declare main EB class instance
namespace EBM{
  inline PROB::ProbEB eb;
}
#endif


#endif