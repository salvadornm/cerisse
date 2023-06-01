#ifndef CNS_H
#define CNS_H

#include <AMReX_AmrLevel.H>

#if CNS_USE_EB
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFluxRegister.H>
#else
#include <AMReX_YAFluxRegister.H>
#endif
// #include <AMReX_ErrorList.H>

#include "PelePhysics.H"
#include "ReactorBase.H"

#include "prob_parm.H"
#include "index_macros.H"
#include "parm.H"
#include "random.H"

using namespace amrex;

class CNS : public amrex::AmrLevel
{
public:
  ////////////////////////////////////////////////////////////////////////////
  //    Constructors and destructors                                        //
  ////////////////////////////////////////////////////////////////////////////

  CNS () {}; // Default constructor, builds invalid object

  CNS (amrex::Amr&            papa,
       int                    lev,
       const amrex::Geometry& level_geom,
       const amrex::BoxArray& bl,
       const amrex::DistributionMapping& dm,
       amrex::Real            time);
  
  virtual ~CNS ();

  CNS (const CNS& rhs) = delete;

  CNS& operator= (const CNS& rhs) = delete;

  ////////////////////////////////////////////////////////////////////////////
  //    AmrLevel virtual functions                                          //
  ////////////////////////////////////////////////////////////////////////////
  
  // Advance grids at this level in time
  virtual amrex::Real advance (amrex::Real time,
                               amrex::Real dt,
                               int iteration,
                               int ncycle) override;

  // Write checkpoint
  virtual void checkPoint (const std::string& dir,
                           std::ostream&      os,
                           amrex::VisMF::How  how = amrex::VisMF::NFiles,
                           bool               dump_old = true) override;

  virtual void computeInitialDt (int                 finest_level,
                                 int                 sub_cycle,
                                 amrex::Vector<int>& n_cycle,
                                 const amrex::Vector<amrex::IntVect>& ref_ratio,
                                 amrex::Vector<amrex::Real>& dt_level,
                                 amrex::Real         stop_time) override;

  virtual void computeNewDt (int                 finest_level,
                             int                 sub_cycle,
                             amrex::Vector<int>& n_cycle,
                             const amrex::Vector<amrex::IntVect>& ref_ratio,
                             amrex::Vector<amrex::Real>& dt_min,
                             amrex::Vector<amrex::Real>& dt_level,
                             amrex::Real         stop_time,
                             int                 post_regrid_flag) override;

  // Error estimation for regridding
  virtual void errorEst (amrex::TagBoxArray& tb,
                          int                 clearval,
                          int                 tagval,
                          amrex::Real         time,
                          int                 n_error_buf = 0,
                          int                 ngrow = 0) override;

  // Initialize data on this level after regridding if old level did not 
  // previously exist
  virtual void init () override;

  // Initialize data on this level from another CNS (during regrid)
  virtual void init (amrex::AmrLevel& old) override;

  // Initialize grid data at problem start-up
  virtual void initData () override;

  virtual int okToContinue () override;

  // After a full time step
  virtual void postCoarseTimeStep (amrex::Real time) override;

  // Do work after init()
  virtual void post_init (amrex::Real stop_time) override;

  virtual void post_regrid (int lbase, int new_finest) override;

  virtual void post_restart () override;

  // Do work after timestep()
  virtual void post_timestep (int iteration) override;

  // Restart from a checkpoint file
  virtual void restart (amrex::Amr&   papa,
                        std::istream& is,
                        bool          bReadSpecial = false) override;

  virtual void setPlotVariables() override;

  virtual int WorkEstType () override { return Cost_Type; }

  // Write a plotfile to specified directory
  virtual void writePlotFile (const std::string& dir,
                              std::ostream&      os,
                              amrex::VisMF::How  how) override;

  virtual std::string thePlotFileType () const override {
    return {"HyperCLaw-V1.1"};
  }

  ////////////////////////////////////////////////////////////////////////////
  //    Public functions                                                    //
  ////////////////////////////////////////////////////////////////////////////
  
  // Define data descriptors
  static void variableSetUp ();

  // Cleanup data descriptors at end of run
  static void variableCleanUp ();

  void compute_dSdt_box (const amrex::Box& bx,
                         amrex::Array4<amrex::Real const>& Sfab,
                         amrex::Array4<amrex::Real      >& dSdtfab,
                         const std::array<FArrayBox*, AMREX_SPACEDIM>& flux);

#if CNS_USE_EB
  void compute_dSdt_box_eb (const amrex::Box& bx,
                            amrex::Array4<amrex::Real       const> const& Sfab,
                            amrex::Array4<amrex::Real            > const& dSdtfab,
                            const std::array<FArrayBox*, AMREX_SPACEDIM>& flux,
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            amrex::Array4<amrex::Real       const> const& vfrac,
                            AMREX_D_DECL(
                              amrex::Array4<amrex::Real     const> const& apx,
                              amrex::Array4<amrex::Real     const> const& apy,
                              amrex::Array4<amrex::Real     const> const& apz),
                            AMREX_D_DECL(
                              amrex::Array4<amrex::Real     const> const& fcx,
                              amrex::Array4<amrex::Real     const> const& fcy,
                              amrex::Array4<amrex::Real     const> const& fcz),
                            amrex::Array4<amrex::Real       const> const& bcent,
                            int                                           as_crse,
                            amrex::Array4<amrex::Real            > const& drho_as_crse,
                            amrex::Array4<int               const> const& rrflag_as_crse,
                            int                                           as_fine,
                            amrex::Array4<amrex::Real            > const& dm_as_fine,
                            amrex::Array4<int               const> const& lev_mask,
                            amrex::Real                                   dt);

  void cns_flux_redistribute (const amrex::Box& bx,
                              amrex::Array4<amrex::Real            > const& dqdt,
                              amrex::Array4<amrex::Real       const> const& divc,
                              amrex::Array4<amrex::Real       const> const& q,
                              amrex::Array4<amrex::Real       const> const& vfrac,
                              amrex::Array4<amrex::EBCellFlag const> const& flag,
                              int as_crse,
                              amrex::Array4<amrex::Real            > const& drho_as_crse,
                              amrex::Array4<int               const> const& rrflag_as_crse,
                              int as_fine,
                              amrex::Array4<amrex::Real            > const& dm_as_fine,
                              amrex::Array4<int               const> const& lev_mask,
                              amrex::Real dt);
#endif

  void computeAvg (amrex::MultiFab& State);

  void compute_pdf_model (amrex::MultiFab& S, amrex::Real dt, int iteration);
  
  // Estimate timestep size based on CFL number
  amrex::Real estTimeStep ();

  void enforce_consistent_state ();

  // These derive functions need to access variables in the class
  static void cns_dertranscoef (const amrex::Box& bx, 
                                amrex::FArrayBox& derfab, 
                                int dcomp, 
                                int ncomp,
                                const amrex::FArrayBox& datafab, 
                                const amrex::Geometry& geomdata,
                                amrex::Real time, 
                                const int* bcrec, 
                                int level);
  
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
  static pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type> trans_parms;

protected:

  ////////////////////////////////////////////////////////////////////////////
  //    Protected functions                                                 //
  ////////////////////////////////////////////////////////////////////////////
  
  // Average down for all the state types
  void avgDown ();

  // Build EB info
  void buildMetrics ();

  void compute_dSdt (const amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt,
#if CNS_USE_EB
                      amrex::EBFluxRegister* fr_as_crse, amrex::EBFluxRegister* fr_as_fine,
#else
                      amrex::YAFluxRegister* fr_as_crse, amrex::YAFluxRegister* fr_as_fine,
#endif
                      bool write_to_flux_register);  

  void fill_ext_src (int i, int j, int k, amrex::Real time,
                      amrex::GeometryData const& geomdata, 
                      amrex::Array4<const amrex::Real> const& state, 
                      amrex::Array4<amrex::Real> const& ext_src, 
                      Parm const& parm, ProbParm const& prob_parm);

  CNS& getLevel (int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

  void printTotalandCheckNan () const;

  // void writeTimeStat (); // TODO

  void react_state (amrex::Real time, amrex::Real dt, bool init_react = false);

  static void read_params ();

  void set_typical_values_chem ();

#if CNS_USE_EB
  const amrex::MultiFab& volFrac () const { return *volfrac; }
#endif

  ////////////////////////////////////////////////////////////////////////////
  //    Private data                                                        //
  ////////////////////////////////////////////////////////////////////////////

  enum StateDataType { State_Type = 0, Reactions_Type, Cost_Type };
  static int num_state_data_types;

  static amrex::BCRec phys_bc;

  //
  // EB related

#if CNS_USE_EB
  const amrex::MultiFab* volfrac;
  const amrex::MultiCutFab* bndrycent;
  std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> areafrac;
  std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> facecent;
  amrex::iMultiFab level_mask;
#endif

#if CNS_USE_EB
  amrex::EBFluxRegister flux_reg;
#else
  amrex::YAFluxRegister flux_reg;
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
  static amrex::Real dt_cutoff;
  static int verbose;
  static int do_reflux;
  static bool signalStopJob;

  static bool do_visc;
  static bool do_react;
  static bool do_ext_src;
  
  static amrex::Real min_react_temp;
  static amrex::Real clip_temp;
  static bool use_typical_vals_chem;
  static int reset_typical_vals_int;
  static bool update_heat_release;
  static int rk_reaction_iter;

  static amrex::Real plm_theta;
  static int recon_scheme;

  static bool do_restart_fields;
  static bool do_psgs;
  static bool do_vpdf;
  static bool do_spdf;

#if CNS_USE_EB
  static int eb_weights_type;
  static int do_reredistribution;
  static bool eb_no_slip;
  static bool eb_isothermal;
  static amrex::Real eb_wall_temp;
#endif
  
  static int refine_cutcells;  
  static int refine_cutcells_max_lev;
  static amrex::Vector<amrex::RealBox> refine_boxes;
  static amrex::Vector<int> refine_boxes_max_lev;
  static amrex::RealBox* dp_refine_boxes;
  static amrex::Vector<amrex::Real> refine_dengrad;
  static int refine_dengrad_max_lev;
  static amrex::Vector<amrex::Real> refine_velgrad;
  static int refine_velgrad_max_lev;
  static amrex::Vector<amrex::Real> refine_magvort;
  static int refine_magvort_max_lev;
  static amrex::Vector<amrex::Real> refine_tke;
  static int refine_tke_max_lev;
  // static amrex::Vector<amrex::AMRErrorTag> errtagger;

  static Parm* h_parm;
  static Parm* d_parm;
};

#endif