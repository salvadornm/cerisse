#ifndef CNS_H
#define CNS_H

#include <AMReX_AmrLevel.H>

#if CNS_USE_EB
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFluxRegister.H>
#else
#include <AMReX_YAFluxRegister.H>
#endif
#include <AMReX_ErrorList.H>
#include <PMFData.H>
#include <PelePhysics.H>
#include <ReactorBase.H>
#include <turbinflow.H>

#include "LES.H"
#include "constants.H"
#include "index_macros.H"
#include "prob_parm.H"
#include "random.H"

using namespace amrex;

class CNS : public amrex::AmrLevel
{
public:
  ////////////////////////////////////////////////////////////////////////////
  //    Constructors and destructors                                        //
  ////////////////////////////////////////////////////////////////////////////

  CNS() {}; // Default constructor, builds invalid object

  CNS(amrex::Amr& papa, int lev, const amrex::Geometry& level_geom,
      const amrex::BoxArray& bl, const amrex::DistributionMapping& dm,
      amrex::Real time);

  virtual ~CNS();

  CNS(const CNS& rhs) = delete;

  CNS& operator=(const CNS& rhs) = delete;

  ////////////////////////////////////////////////////////////////////////////
  //    AmrLevel virtual functions                                          //
  ////////////////////////////////////////////////////////////////////////////

  // Advance grids at this level in time
  virtual amrex::Real advance(amrex::Real time, amrex::Real dt, int iteration,
                              int ncycle) override;

  // Write checkpoint
  virtual void checkPoint(const std::string& dir, std::ostream& os,
                          amrex::VisMF::How how = amrex::VisMF::NFiles,
                          bool dump_old = true) override;

  virtual void computeInitialDt(int finest_level, int sub_cycle,
                                amrex::Vector<int>& n_cycle,
                                const amrex::Vector<amrex::IntVect>& ref_ratio,
                                amrex::Vector<amrex::Real>& dt_level,
                                amrex::Real stop_time) override;

  virtual void computeNewDt(int finest_level, int sub_cycle,
                            amrex::Vector<int>& n_cycle,
                            const amrex::Vector<amrex::IntVect>& ref_ratio,
                            amrex::Vector<amrex::Real>& dt_min,
                            amrex::Vector<amrex::Real>& dt_level,
                            amrex::Real stop_time, int post_regrid_flag) override;

  // Error estimation for regridding
  virtual void errorEst(amrex::TagBoxArray& tb, int clearval, int tagval,
                        amrex::Real time, int n_error_buf = 0,
                        int ngrow = 0) override;

  // Initialize data on this level after regridding if old level did not
  // previously exist
  virtual void init() override;

  // Initialize data on this level from another CNS (during regrid)
  virtual void init(amrex::AmrLevel& old) override;

  // Initialize grid data at problem start-up
  virtual void initData() override;

  virtual int okToContinue() override;

  // After a full time step
  virtual void postCoarseTimeStep(amrex::Real time) override;

  // Do work after init()
  virtual void post_init(amrex::Real stop_time) override;

  virtual void post_regrid(int lbase, int new_finest) override;

  virtual void post_restart() override;

  // Do work after timestep()
  virtual void post_timestep(int iteration) override;

  // Restart from a checkpoint file
  virtual void restart(amrex::Amr& papa, std::istream& is,
                       bool bReadSpecial = false) override;

  virtual void setPlotVariables() override;

  // Provide state data type to AmrLevel for load balancing. -1 means none.
  virtual int WorkEstType() override { return do_load_balance ? Cost_Type : -1; }

  // Write a plotfile to specified directory
  virtual void writePlotFile(const std::string& dir, std::ostream& os,
                             amrex::VisMF::How how) override;

  virtual std::string thePlotFileType() const override { return {"HyperCLaw-V1.1"}; }

  ////////////////////////////////////////////////////////////////////////////
  //    Public functions                                                    //
  ////////////////////////////////////////////////////////////////////////////

  // Define data descriptors
  static void variableSetUp();

  // Cleanup data descriptors at end of run
  static void variableCleanUp();

  void compute_dSdt_box(const amrex::Box& bx, amrex::Array4<amrex::Real const>& Sfab,
                        amrex::Array4<amrex::Real>& dSdtfab,
                        const std::array<FArrayBox*, AMREX_SPACEDIM>& flux,
                        amrex::Array4<const amrex::Real>& shock_sensor);

#if CNS_USE_EB
  void compute_dSdt_box_eb(const amrex::Box& bx,
                           amrex::Array4<amrex::Real const> const& Sfab,
                           amrex::Array4<amrex::Real> const& dSdtfab,
                           const std::array<FArrayBox*, AMREX_SPACEDIM>& flux,
                           amrex::Array4<amrex::EBCellFlag const> const& flag,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                        amrex::Array4<amrex::Real const> const& apy,
                                        amrex::Array4<amrex::Real const> const& apz),
                           AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                        amrex::Array4<amrex::Real const> const& fcy,
                                        amrex::Array4<amrex::Real const> const& fcz),
                           amrex::Array4<amrex::Real const> const& bcent,
                           int as_crse,
                           amrex::Array4<amrex::Real> const& drho_as_crse,
                           amrex::Array4<int const> const& rrflag_as_crse,
                           int as_fine, amrex::Array4<amrex::Real> const& dm_as_fine,
                           amrex::Array4<int const> const& lev_mask, amrex::Real dt,
                           amrex::Array4<const amrex::Real>& shock_sensor);
#endif

  void computeAvg(amrex::MultiFab& State);

  void compute_pdf_model(amrex::MultiFab& S, amrex::Real dt, int iteration);

  // Estimate timestep size based on CFL number
  amrex::Real estTimeStep();

  void enforce_consistent_state();

  void enforce_consistent_state(amrex::MultiFab& State);

  static void apply_nscbc(const amrex::Box& bx,
                          const amrex::Array4<amrex::Real>& data, int start_comp,
                          const amrex::Geometry& geom, int dir, int lo_hi,
                          int bc_type);

  // These derive functions need to access variables in the class
  static void cns_dertranscoef(const amrex::Box& bx, amrex::FArrayBox& derfab,
                               int dcomp, int ncomp, const amrex::FArrayBox& datafab,
                               const amrex::Geometry& geomdata, amrex::Real time,
                               const int* bcrec, int level);

  // static void cns_derextsrc (const amrex::Box& bx,
  //                            amrex::FArrayBox& derfab,
  //                            int dcomp,
  //                            int ncomp,
  //                            const amrex::FArrayBox& datafab,
  //                            const amrex::Geometry& geomdata,
  //                            amrex::Real time,
  //                            const int* bcrec,
  //                            int level);

  ////////////////////////////////////////////////////////////////////////////
  //    Public data                                                         //
  ////////////////////////////////////////////////////////////////////////////

  static ProbParm* h_prob_parm;
  static ProbParm* d_prob_parm;
#ifdef USE_PROB_PARM_HOST
  static ProbParmHost* prob_parm_host;
#endif

  static pele::physics::transport::TransportParams<
    pele::physics::PhysicsType::transport_type>
    trans_parms;
  static pele::physics::turbinflow::TurbInflow turb_inflow;
  static pele::physics::PMF::PmfData pmf_data;

  // Eddy viscosity type LES

  static bool do_les;
  static bool do_pasr;
  static amrex::Real Cs;
  static amrex::Real C_I;
  static amrex::Real Pr_T;
  static amrex::Real Sc_T;
  static std::unique_ptr<LESModel> les_model;
  static amrex::Real Cm;

  // NSCBC

  static bool do_nscbc;
  static amrex::Real nscbc_relax_p;
  static amrex::Real nscbc_relax_u;
  static amrex::Real nscbc_relax_T;
  static amrex::Real ambient_p;
  static amrex::BCRec phys_bc; // this has to be public for cns_bcfill to access

  // Reconstruction

  static amrex::Real teno_cutoff;

#ifndef AMREX_USE_CUDA
protected:
#endif
  ////////////////////////////////////////////////////////////////////////////
  //    Protected functions                                                 //
  ////////////////////////////////////////////////////////////////////////////

  // Average down for all the state types
  void avgDown();

  // Build EB info
  void buildMetrics();

  void compute_dSdt(const amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt,
#if CNS_USE_EB
                    amrex::EBFluxRegister* fr_as_crse,
                    amrex::EBFluxRegister* fr_as_fine,
#else
                    amrex::YAFluxRegister* fr_as_crse,
                    amrex::YAFluxRegister* fr_as_fine,
#endif
                    bool write_to_flux_register = true);

  void fill_ext_src(int i, int j, int k, amrex::Real time,
                    amrex::GeometryData const& geomdata,
                    amrex::Array4<const amrex::Real> const& state,
                    amrex::Array4<amrex::Real> const& ext_src,
                    ProbParm const& prob_parm);

#ifdef USE_FULL_PROB_POST_TIMESTEP
  void full_prob_post_timestep(int iteration);
#endif

#ifdef USE_PROB_POST_COARSETIMESTEP
  void prob_post_coarsetimestep(amrex::Real time);
#endif

  CNS& getLevel(int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

  void printTotalandCheckNan() const;

  void checkRuntimeMessages();

  // void writeTimeStat (); // TODO

  void react_state(amrex::Real time, amrex::Real dt, bool init_react = false);

  static void read_params();

  static void read_errtags();

  void set_typical_values_chem();

#if CNS_USE_EB
  const amrex::MultiFab& volFrac() const { return *volfrac; }
#endif

  // void fillFineMask();

  ////////////////////////////////////////////////////////////////////////////
  //    Private data                                                        //
  ////////////////////////////////////////////////////////////////////////////

  enum StateDataType {
    State_Type = 0,
    Reactions_Type,
    Cost_Type,
    num_state_data_types
  };

  //
  // Advance related

#if CNS_USE_EB
  amrex::EBFluxRegister flux_reg;
#else
  amrex::YAFluxRegister flux_reg;
#endif

  amrex::MultiFab Sborder; // A state array with ghost zones. It is put here inside
                           // the class so alloc overhead only happen once.
  amrex::MultiFab shock_sensor_mf; // for hybrid scheme switch
  // amrex::iMultiFab ifine_mask; // fine covered cells plus 2 buffers because they
  // are used for interpolation in fine level ghost cells

  //
  // EB related

#if CNS_USE_EB
  const amrex::MultiFab* volfrac;
  const amrex::MultiCutFab* bndrycent;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;
  amrex::iMultiFab level_mask;
#endif

  //
  // PelePhysics related

  static amrex::Vector<std::string> spec_names;
  static std::string chem_integrator;
  std::unique_ptr<pele::physics::reactions::ReactorBase> reactor;

  //
  // PDF related

  UniqueRand WienerProcess;

  //
  // Control parameters

  static amrex::Real cfl;
  static amrex::Real fixed_dt;
  static amrex::Real dt_cutoff;
  static amrex::Real dt_max_change;
  static int rk_order;
  static int amr_interp_order;
  static int verbose;
  static bool do_reflux;
  static bool do_load_balance;
  static bool signalStopJob;
  static int check_message_int;

  static bool do_hydro;
  static bool do_visc;
  static bool do_react;
  static bool do_ext_src;

  static amrex::Real min_react_temp;
  static amrex::Real clip_temp;
  static bool use_typical_vals_chem;
  static int reset_typical_vals_int;
  static bool update_heat_release;
  static int rk_reaction_iter;

  // Reconstruction

  static bool recon_char_var;
  static int char_sys;
  static amrex::Real plm_theta;
  static int recon_scheme;
  static bool use_hybrid_scheme;

  // LES-PDF stochastic fields

  static bool do_restart_fields;
  static bool do_psgs;
  static bool do_pd_model;
  static bool do_vpdf;
  static bool do_spdf;
  static bool do_species_langevin;

  // EB

#if CNS_USE_EB
  static int eb_recon_mode;
  static amrex::Real eb_weight;
  static bool eb_no_slip;
  static std::string redistribution_type;
#endif
  static bool eb_isothermal;
  static bool eb_wall_model;
  static amrex::Real eb_wall_temp;

  // Buffer layer

  static amrex::RealBox buffer_box;

  // AMR controls

  static bool refine_cutcells;
  static int refine_cutcells_max_lev;
  static amrex::Vector<amrex::RealBox> refine_boxes;
  static amrex::Vector<int> refine_boxes_max_lev;
  static amrex::RealBox* dp_refine_boxes;
  static amrex::Vector<amrex::Real> refine_dengrad;
  static int refine_dengrad_max_lev;
  static amrex::Vector<amrex::Real> refine_velgrad;
  static int refine_velgrad_max_lev;
  static amrex::Vector<amrex::Real> refine_presgrad;
  static int refine_presgrad_max_lev;
  static amrex::Vector<amrex::Real> refine_magvort;
  static int refine_magvort_max_lev;
  static amrex::Vector<amrex::Real> refine_tke;
  static int refine_tke_max_lev;
  static amrex::Vector<amrex::AMRErrorTag> errtags;
};

void cns_bcfill(amrex::Box const& bx, amrex::FArrayBox& data, const int dcomp,
                const int numcomp, amrex::Geometry const& geom,
                const amrex::Real time, const amrex::Vector<amrex::BCRec>& bcr,
                const int bcomp, const int scomp);

void cns_null_bcfill(amrex::Box const& bx, amrex::FArrayBox& data, const int dcomp,
                     const int numcomp, amrex::Geometry const& geom,
                     const amrex::Real time, const amrex::Vector<amrex::BCRec>& bcr,
                     const int bcomp, const int scomp);

#endif