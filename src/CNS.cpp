#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>
#include <CNS.h>
#include <CNS_K.h>
#include <prob.h>

#include <nscbc.h>

// 
#include "mandebug.h"


#ifdef CNS_USE_FSI
#include <fsi/Kinematics.h>
#include <fsi/RigidBodyProperties.h>
#endif

#ifdef USE_PELEPHYSICS
#include "TransPele.h"

pele::physics::PeleParams<
  pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type
  >> trans_parms;
#endif


using namespace amrex;

bool CNS::verbose = true;
bool CNS::record_probe = false;
bool CNS::dt_dynamic = false;
bool CNS::ib_move = false;
bool CNS::plot_surf = false;

int CNS::surf_int = 10000000;
std::string CNS::surf_filename = "surfplot";

// utilities
bool CNS::use_utility = false; 
Utility CNS::utilidades;

Real CNS::eb_weight = 0.0;
bool CNS::eb_redistribution = false;
std::string CNS::eb_redistribution_type = "NoRedist";

int CNS::nstep_screen_output = 10;
int CNS::order_rk = 2;
int CNS::stages_rk = 2;
bool CNS::strict_positivity = false;
bool CNS::pass2_static = false;
int CNS::do_reflux = 0; // default reflux is off
int CNS::refine_max_dengrad_lev = -1;
Real CNS::cfl = 0.0_rt;
Real CNS::dt_constant = 0.0_rt;
Real CNS::refine_dengrad = 1.0e10;
int  CNS::INDEX_THERM = 0;
bool CNS::compute_stats = false;
bool CNS::record_stats = false;
Real CNS::time_stats = 0.0;
Real CNS::time_stat_level[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

PROB::ProbClosures *CNS::h_prob_closures = nullptr;
PROB::ProbClosures *CNS::d_prob_closures = nullptr;
PROB::ProbParm *CNS::h_prob_parm = nullptr;
PROB::ProbParm *CNS::d_prob_parm = nullptr;
BCRec *CNS::h_phys_bc = nullptr;
BCRec *CNS::d_phys_bc = nullptr;

bool CNS::use_nscbc = false;
int CNS::nscbc_order = 2;
amrex::GpuArray<int, AMREX_SPACEDIM> CNS::nscbc_lo = {AMREX_D_DECL(0,0,0)};
amrex::GpuArray<int, AMREX_SPACEDIM> CNS::nscbc_hi = {AMREX_D_DECL(0,0,0)};

CNS::NSCBCParm CNS::h_nscbc_parm{};
CNS::NSCBCParm* CNS::d_nscbc_parm = nullptr;

int CNS::second_order_flux_method = CNS::NONE;

// needed for CNSBld - derived from LevelBld (abstract class, pure virtual
// functions must be implemented)

CNS::CNS() {}

CNS::CNS(Amr &papa, int lev, const Geometry &level_geom, const BoxArray &bl,
         const DistributionMapping &dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time) {
  if (do_reflux && level > 0) {
    flux_reg.define(bl, papa.boxArray(level - 1), dm,
                    papa.DistributionMap(level - 1), level_geom,
                    papa.Geom(level - 1), papa.refRatio(level - 1), level,
                    PROB::ProbClosures::NCONS);
  }

#ifdef AMREX_USE_GPIBM
  IBM::ib.build_mf(grids, dmap, level);
  IBM::ib.computeMarkers(level);
#endif

#ifdef CNS_USE_EB
  EBM::eb.build_mf(grids, dmap, level);
#endif

  buildMetrics();


  // build NSCBC ghost state if needed
  // if (use_nscbc)
  // {

  //   amrex::Print() << "CNS::CNS() building NSCBC ghost state MultiFab " << std::endl;

    
  // }


  rz_sanity_check(Geom());
};

CNS::~CNS() {}
// -----------------------------------------------------------------------------

// ------------------------------------------------------------------------------------//
// Read cerisse-specific parameters from input file, starting by cns.
//
void CNS::read_params() {

  ParmParse pp("cns");

  pp.query("nstep_screen_output", nstep_screen_output);
  pp.query("verbose", verbose);

  pp.query("record_probe", record_probe);
  pp.query("record_stats", record_stats);

  Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    h_phys_bc->setLo(i, lo_bc[i]);
    h_phys_bc->setHi(i, hi_bc[i]);
  }

  // Read NSBC type 0 (N/A) 1:Inflow 2:Outflow (constant P)
  Vector<int> nslo(AMREX_SPACEDIM, 0); Vector<int> nshi(AMREX_SPACEDIM, 0);

  bool has_nslo = pp.queryarr("nscbc_lo", nslo, 0, AMREX_SPACEDIM);
  bool has_nshi = pp.queryarr("nscbc_hi", nshi, 0, AMREX_SPACEDIM);

  use_nscbc = false;

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    nscbc_lo[d] = nslo[d];
    nscbc_hi[d] = nshi[d];

    if (nscbc_lo[d] != 0 || nscbc_hi[d] != 0) {
      use_nscbc = true;
    }
  }


  if (use_nscbc) {
    amrex::Print() << "Using NSCBC/LODI flags: lo="
                 << AMREX_D_TERM(nscbc_lo[0], << " " << nscbc_lo[1], << " " << nscbc_lo[2])
                 << " hi="
                 << AMREX_D_TERM(nscbc_hi[0], << " " << nscbc_hi[1], << " " << nscbc_hi[2])
                 << "\n";
   
    // read NSBC parameters and store in host_nsbc_parm
    pp.query("nscbc_Lchar",   h_nscbc_parm.Lchar);
    pp.query("nscbc_Mmax",    h_nscbc_parm.Mmax);
    pp.query("nscbc_Ptarget", h_nscbc_parm.Ptarget);
    pp.query("nscbc_sigma",   h_nscbc_parm.sigma);
    pp.query("nscbc_utarget", h_nscbc_parm.utarget);
    pp.query("nscbc_vtarget", h_nscbc_parm.vtarget);
    pp.query("nscbc_wtarget", h_nscbc_parm.wtarget);
    pp.query("nscbc_Ttarget", h_nscbc_parm.Ttarget);
    pp.query("nscbc_eta",     h_nscbc_parm.eta);
    pp.query("nscbc_use_transverse", h_nscbc_parm.use_transverse);
    pp.query("nscbc_beta_transverse",h_nscbc_parm.beta_transverse);
    pp.query("nscbc_order", nscbc_order);
    

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(d_nscbc_parm, &h_nscbc_parm, sizeof(NSCBCParm));
#endif

    // check parameters validity
    const NSCBCParm* nscbc_parm =d_nscbc_parm;
    nscbc::check_nscbc<PROB::ProbClosures>(nslo,nshi,*nscbc_parm);
  } 
  
  // second order flux close to BC 
  pp.query("second_order_flux_method", second_order_flux_method);

  if (second_order_flux_method != KEEP && second_order_flux_method != SKEW && second_order_flux_method != NONE) {
    amrex::Abort( "cns.second_order_flux_method must be 0 (NONE), 1 (KEEP) or 2 (SKEW)");
  }
  //


  pp.query("do_reflux", do_reflux);

  if (!pp.query("order_rk", order_rk)) {
    amrex::Abort(
        "Need to specify SSPRK scheme order of accuracy, order_rk={-2, 1, 2, "
        "3}");
  }

  if (!pp.query("stages_rk", stages_rk)) {
    amrex::Abort("Need to specify SSPRK number of stages, stages_rk");
  } else {
    if (order_rk == 1 && stages_rk != 1) {
      amrex::Abort("Forward Euler number of stages must be 1");
    }

    if (order_rk == 2 && stages_rk < order_rk) {
      amrex::Abort(
          "SSPRK2 number of stages must equal or greater than order of "
          "accuracy");
    }
    if (order_rk == 3 && !(stages_rk == 4 || stages_rk == 3)) {
      amrex::Abort("SSPRK3 number of stages must equal 3 or 4");
    }
  }

  pp.query("strict_positivity", strict_positivity);
  if (strict_positivity) {
    amrex::Print() << "  cns.strict_positivity = 1 (abort if state approaches smallr/ei_min floors)\n";
  }

  pp.query("pass2_static", pass2_static);
  if (pass2_static) {
    amrex::Print() << "  cns.pass2_static = 1 (Pass 2 flood-fill runs for static geometry too)\n";
  }

  //  Utilities options ----------------------------------------------------
  // specific keywords for Utilities
  pp.query("use_utility",use_utility);
  if (use_utility)
  {
    amrex::Print() << " Using Utilities " << std::endl;      
    ParmParse pp_util("util");

    // PMF
    bool use_PMF=false;
    pp_util.query("use_PMF",use_PMF); 
    if (use_PMF){
#ifdef USE_PELEPHYSICS      
      amrex::Print() << " Reading PMF from file.. " << std::endl;      
      CNS::utilidades.initPMF();
#else      
      amrex::Abort("using PMF files need PelePhysics");
#endif
    }  

    // Read from file
    bool use_turb_file  = false;
    pp_util.query("use_turb_file",use_turb_file); 
    if (use_turb_file){
      std::string turbfilename;
      pp_util.query("turb_file",turbfilename); 
      amrex::Print() << " Reading turbulence from file: " << turbfilename << std::endl; 
      CNS::utilidades.initTurbulenceFile(turbfilename);
    }

  }

#ifdef AMREX_USE_GPIBM
  // IBM-specific input keywords (ib.* namespace)
  ParmParse ppib("ib");
  if (!ppib.query("move", ib_move)) {
    amrex::Abort("ib.move not specified (0=false, 1=true)");
  }
  if (!ppib.query("plot_surf", plot_surf)) {
    amrex::Abort("ib.plot_surf not specified (0=false, 1=true)");
  }
  if (plot_surf) {
    ppib.query("surf_int", surf_int);
    ppib.query("surf_file", surf_filename);
  }
#endif
  
#if CNS_USE_EB 
  // specific keywords for EB boundaries
  ParmParse ppeb2("eb2");  
  ppeb2.query("eb_weight",eb_weight); 
  ppeb2.query("redistribution_type", eb_redistribution_type);
  if (eb_redistribution_type != "StateRedist" && eb_redistribution_type != "FluxRedist" &&
      eb_redistribution_type != "NoRedist"    && eb_redistribution_type != "NewRedist") {
    amrex::Abort( " input file: redistribution_type must be StateRedist/FluxRedist/NewRedist/NoRedist");
  }
  if (eb_redistribution_type != "NoRedist") eb_redistribution =true;
  // This communicates to the class (not very elegant)
  EBM::eb.eb_weight = eb_weight;
  EBM::eb.redistribution_type = eb_redistribution_type; 

#endif


#ifdef USE_PELEPHYSICS
  // One-time transport parameter initialization (host->device)
  static bool trans_inited = false;
  if (!trans_inited) {
    trans_parms.initialize();
    trans_inited = true;
  }
#endif


#if AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(d_prob_closures, h_prob_closures,
                          sizeof(PROB::ProbClosures));
  amrex::Gpu::htod_memcpy(d_prob_parm, h_prob_parm, sizeof(PROB::ProbParm));
  amrex::Gpu::htod_memcpy(d_phys_bc, h_phys_bc, sizeof(BCRec));
#endif
}

// ------------------------------------------------------------------------------------//
void CNS::init(AmrLevel &old) {
  auto &oldlev = dynamic_cast<CNS &>(old);

  // amrex::Print( ) << " oo CNS::init (AMR recast) -----  " << std::endl;  

  Real dt_new = parent->dtLevel(level);
  Real cur_time = oldlev.state[State_Type].curTime();
  Real prev_time = oldlev.state[State_Type].prevTime();
  Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  MultiFab &S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0,PROB::ProbClosures::NCONS);

  if (use_nscbc) {
    define_nscbc_ghost_shell();
    initialise_nscbc_ghost_shell(S_new);
  }

  if (compute_stats){
    MultiFab &Sstat_new = get_new_data(Stats_Type);
    FillPatch(old, Sstat_new, 0, cur_time, Stats_Type, 0,PROB::ProbClosures::NSTAT);
  }

}
// ------------------------------------------------------------------------------------//
void CNS::init() {

  // amrex::Print( ) << " oo CNS::init -----  " << std::endl;  

  Real dt = parent->dtLevel(level);
  Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  Real prev_time = getLevel(level - 1).state[State_Type].prevTime();
  Real dt_old = (cur_time - prev_time) /
                static_cast<Real>(parent->MaxRefRatio(level - 1));
  setTimeLevel(cur_time, dt_old, dt);

  MultiFab &S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0,PROB::ProbClosures::NCONS);
};

// ------------------------------------------------------------------------------------//
void CNS::initData() {
  BL_PROFILE("CNS::initData()");

  const auto geomdata = geom.data();
  MultiFab &S_new = get_new_data(State_Type);
  auto const &sma = S_new.arrays();

  PROB::ProbClosures const *lclosures = d_prob_closures;
  PROB::ProbParm const *lprobparm = d_prob_parm;

  // Initialise problem by calling user-given prob.h
#if USE_UTILITY
  amrex::ParallelFor(
    S_new, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      prob_initdata(i, j, k, sma[box_no], geomdata, *lclosures, *lprobparm, use_utility ? &CNS::utilidades : nullptr);
  });
#else
  amrex::ParallelFor(
      S_new, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        prob_initdata(i, j, k, sma[box_no], geomdata, *lclosures, *lprobparm);        
  });
#endif
            
  // Initialise stats 
  if (compute_stats) {
    setupStats();
  }


  // Initialise coefficients   
   prob_rhs.init_coeffs();

  // Initialise NSBC face values 
  if (use_nscbc) {
    define_nscbc_ghost_shell();
    initialise_nscbc_ghost_shell(S_new);
  }


}

// ------------------------------------------------------------------------------------//
void CNS::buildMetrics() {
  if (verbose) {
    const Real *dx = geom.CellSize();
    amrex::Print() << "Mesh size (dx,dy,dz) = ";
    amrex::Print() << AMREX_D_TERM(dx[0], << "  " << dx[1], << "  " << dx[2]) << "  \n";
  }  
}
//------------------------------------------------------------------------------------//
void CNS::post_init(Real /*stop_time*/) {

  //amrex::Print() << " CNS::post_init level= "  << level << std::endl;

  if (level > 0) {
    return;
  };

  for (int k = parent->finestLevel() - 1; k >= 0; --k) {
    getLevel(k).avgDown();
  }

  if (verbose) {
    printTotal();
  }
  
  // Set up diagnostics
  if (record_probe) {
    setupTimeProbe();
  }
  
#if CNS_USE_EB
  EBM::eb.check_geometry(level);
#endif
  

}
// -----------------------------------------------------------------------------

// Time-stepping ---------------------------------------------------------------
void CNS::computeInitialDt(int finest_level, int sub_cycle,
                           Vector<int> &n_cycle,  // no. of subcycling steps
                           const Vector<IntVect> &ref_ratio,
                           Vector<Real> &dt_level, Real stop_time) {
  BL_PROFILE("CNS::computeInitialDt()");
  //amrex::Print() << " oo CNS::computeInitialDt  " << std::endl;

  Real dt0 = std::numeric_limits<Real>::max();
  Vector<GpuArray<Real,AMREX_SPACEDIM>> eigenvals_level;
  Vector<Real> CFL_level;
  eigenvals_level.resize(finest_level + 1);
  CFL_level.resize(finest_level + 1);

  // Compute max eigenvalues in all directions on each level
  for (int i = 0; i <= finest_level; i++) {
    eigenvals_level[i] = getLevel(i).maxEigen();
  }

  if (dt_dynamic) {
    // Dynamic dt
    for (int i = 0; i <= finest_level; i++) {
      const GpuArray<Real, AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
#if (AMREX_SPACEDIM == 1)
      dt_level[i] = cfl * dx[0] / eigenvals_level[i][0];
#elif (AMREX_SPACEDIM == 2)
      dt_level[i] =
          cfl * amrex::min( dx[0] / eigenvals_level[i][0],
                            dx[1] / eigenvals_level[i][1] );                                        
#else
      dt_level[i] =
          cfl * amrex::min(AMREX_D_DECL(dx[0] / eigenvals_level[i][0],
                                        dx[1] / eigenvals_level[i][1],
                                        dx[2] / eigenvals_level[i][2]));
#endif
    }
    // Find min dt across all levels
    int nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt0 = std::min(dt0, nfactor * dt_level[i]);
    }
  } 
  else {
    // If constant dt
    dt0 = dt_constant;
  }
  // Set dt for all levels
  int nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    const GpuArray<Real,AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
    nfactor *= n_cycle[i];
    dt_level[i] = dt0 / nfactor;
#if (AMREX_SPACEDIM == 1)
    CFL_level[i] = dt_level[i] * eigenvals_level[i][0] / dx[0];
#elif (AMREX_SPACEDIM == 2)
     CFL_level[i] =
        dt_level[i] * amrex::max( eigenvals_level[i][0] / dx[0],
                                  eigenvals_level[i][1] / dx[1]);
#else
    CFL_level[i] =
        dt_level[i] * amrex::max(AMREX_D_DECL(eigenvals_level[i][0] / dx[0],
                                              eigenvals_level[i][1] / dx[1],
                                              eigenvals_level[i][2] / dx[2]));
#endif                                          
  }
  // Print
  if (ParallelDescriptor::IOProcessor()) {
    for (int i = 0; i <= finest_level; i++) {
      printf("[computeInitialDt] Level %d, Max CFL= %f   Max eigenvalues =( ",i,CFL_level[i]);
      for (int dim=0; dim < AMREX_SPACEDIM; dim++)
      { printf(" %f ",eigenvals_level[i][dim]);}
      printf(" ) \n");
      
    }
  }
}

// ------------------------------------------------------------------------------------//
// Called at the end of a coarse grid timecycle or after regrid, to compute the
// dt (time step) for all levels, for the next step.
// Output dt_level
// ------------------------------------------------------------------------------------//
void CNS::computeNewDt(int finest_level, int sub_cycle, Vector<int> &n_cycle,
                       const Vector<IntVect> &ref_ratio, Vector<Real> &dt_min,
                       Vector<Real> &dt_level, Real stop_time,
                       int post_regrid_flag) {
  BL_PROFILE("CNS::computeNewDt()");

  //amrex::Print() << " oo CNS::computeNewDt " << std::endl;

  Real dt0 = std::numeric_limits<Real>::max();
  Vector<GpuArray<Real,AMREX_SPACEDIM>> eigenvals_level;
  Vector<Real> CFL_level;
  eigenvals_level.resize(finest_level + 1);
  CFL_level.resize(finest_level + 1);

  // Compute max eigenvalues in all directions on each level
  for (int i = 0; i <= finest_level; i++) {
    eigenvals_level[i] = getLevel(i).maxEigen();
  }

  if (dt_dynamic) {
    // Estimate timestep across all points, levels, and procs
    for (int i = 0; i <= finest_level; i++) {
      const GpuArray<Real,AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
#if (AMREX_SPACEDIM == 1)
      dt_min[i] = cfl * dx[0] / eigenvals_level[i][0];
#elif (AMREX_SPACEDIM ==2) 
      dt_min[i] = cfl * amrex::min( dx[0] / eigenvals_level[i][0],
                                    dx[1] / eigenvals_level[i][1] );   
#else     
      dt_min[i] = cfl * amrex::min(AMREX_D_DECL(dx[0] / eigenvals_level[i][0],
                                                dx[1] / eigenvals_level[i][1],
                                                dx[2] / eigenvals_level[i][2]));
#endif                                          
    };

    // Limit dt
    if (post_regrid_flag == 1) {
      // Limit dt's by pre-regrid dt
      for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = std::min(dt_min[i], dt_level[i]);
      }
    } else {
      // Limit dt's by change_max * old dt
      static Real change_max = 1.1; //////// PARAMETER /////////////
      for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
      }
    }
    // Find the minimum over all levels
    int nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt0 = std::min(dt0, nfactor * dt_min[i]);
    }
    // Limit dt0 by the value of stop_time.
    const Real eps = 0.001_rt * dt0; ////////// PARAMETER /////////////
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0_rt) {
      if ((cur_time + dt0) > (stop_time - eps)) {
        dt0 = stop_time - cur_time;
      }
    }
    // Set dt at all levels
    nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt_level[i] = dt0 / nfactor;
    }
  } else {
    // If constant dt
    dt0 = dt_constant;
  }
  // Set dt for all levels
  int nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    const GpuArray<Real,AMREX_SPACEDIM> dx = parent->Geom(i).CellSizeArray();
    nfactor *= n_cycle[i];
    dt_level[i] = dt0 / nfactor;
#if (AMREX_SPACEDIM == 1)
    CFL_level[i] = dt_level[i] * eigenvals_level[i][0] / dx[0];
#elif (AMREX_SPACEDIM ==2) 
    CFL_level[i] =
        dt_level[i] * amrex::max(eigenvals_level[i][0] / dx[0],
                                 eigenvals_level[i][1] / dx[1]);    
#else    
    CFL_level[i] =
        dt_level[i] * amrex::max(AMREX_D_DECL(eigenvals_level[i][0] / dx[0],
                                              eigenvals_level[i][1] / dx[1],
                                              eigenvals_level[i][2] / dx[2]));
#endif                                          
  }
}

// Returns maximum eigenvalue in each direction
 [[nodiscard]] GpuArray<Real,AMREX_SPACEDIM> CNS::maxEigen() {
  BL_PROFILE("CNS::maxEigen()");

  PROB::ProbClosures const *d_cls = d_prob_closures;

  // Get multifabs
  MultiFab& consmf = get_new_data(State_Type);

  GpuArray<Real,AMREX_SPACEDIM> h_max_eigenvals;

  // Use ReduceOps for proper GPU-parallel reduction (replaces the original
  // AsyncArray approach which had a data race on the device max-reduction).
#if (AMREX_SPACEDIM == 1)
  ReduceOps<ReduceOpMax> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (MFIter mfi(consmf, false); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    const Array4<Real>& cons = consmf.array(mfi);
    reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
      GpuArray<int, 3> vdir0 = {1, 0, 0};
      auto temp0 = d_cls->cons2eigenvals(i, j, k, cons, vdir0);
      Real maxe0 = Real(0.0);
      for (int iw = 0; iw < PROB::ProbClosures::NWAVES; iw++)
        maxe0 = amrex::max(maxe0, std::abs(temp0[iw]));
      return {maxe0};
    });
  }
  auto hv = reduce_data.value(reduce_op);
  h_max_eigenvals[0] = amrex::get<0>(hv);

#elif (AMREX_SPACEDIM == 2)
  ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_GPIBM
  // IBM-aware CFL: skip solid cells when computing max eigenvalue.
  // Solid cells don't participate in time integration (their RHS is zeroed
  // in compute_rhs), so their density/pressure values — which come from
  // extrapolation, not from the flow — should not constrain the timestep.
  // Without this, a solid cell with extrapolated ρ near zero produces a
  // huge sound speed → dt → 0 → simulation stalls.
  auto& ib_mf = *IBM::ib.bmf_a[level];
#endif

  for (MFIter mfi(consmf, false); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    const Array4<Real>& cons = consmf.array(mfi);
#ifdef AMREX_USE_GPIBM
    const auto& ibm = ib_mf.const_array(mfi);
#endif
    reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
#ifdef AMREX_USE_GPIBM
      // Skip solid cells — they don't evolve and shouldn't constrain CFL
      if (ibm(i,j,k,0) != 0) return {Real(0.0), Real(0.0)};
#endif
      GpuArray<int, 3> vdir0 = {1, 0, 0};
      GpuArray<int, 3> vdir1 = {0, 1, 0};
      auto temp0 = d_cls->cons2eigenvals(i, j, k, cons, vdir0);
      auto temp1 = d_cls->cons2eigenvals(i, j, k, cons, vdir1);
      Real maxe0 = Real(0.0), maxe1 = Real(0.0);
      for (int iw = 0; iw < PROB::ProbClosures::NWAVES; iw++) {
        maxe0 = amrex::max(maxe0, std::abs(temp0[iw]));
        maxe1 = amrex::max(maxe1, std::abs(temp1[iw]));
      }
      return {maxe0, maxe1};
    });
  }
  auto hv = reduce_data.value(reduce_op);
  h_max_eigenvals[0] = amrex::get<0>(hv);
  h_max_eigenvals[1] = amrex::get<1>(hv);

#else // 3D
  ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_GPIBM
  auto& ib_mf = *IBM::ib.bmf_a[level];
#endif

  for (MFIter mfi(consmf, false); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    const Array4<Real>& cons = consmf.array(mfi);
#ifdef AMREX_USE_GPIBM
    const auto& ibm = ib_mf.const_array(mfi);
#endif
    reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
#ifdef AMREX_USE_GPIBM
      if (ibm(i,j,k,0) != 0) return {Real(0.0), Real(0.0), Real(0.0)};
#endif
      GpuArray<int, 3> vdir0 = {1, 0, 0};
      GpuArray<int, 3> vdir1 = {0, 1, 0};
      GpuArray<int, 3> vdir2 = {0, 0, 1};
      auto temp0 = d_cls->cons2eigenvals(i, j, k, cons, vdir0);
      auto temp1 = d_cls->cons2eigenvals(i, j, k, cons, vdir1);
      auto temp2 = d_cls->cons2eigenvals(i, j, k, cons, vdir2);
      Real maxe0 = Real(0.0), maxe1 = Real(0.0), maxe2 = Real(0.0);
      for (int iw = 0; iw < PROB::ProbClosures::NWAVES; iw++) {
        maxe0 = amrex::max(maxe0, std::abs(temp0[iw]));
        maxe1 = amrex::max(maxe1, std::abs(temp1[iw]));
        maxe2 = amrex::max(maxe2, std::abs(temp2[iw]));
      }
      return {maxe0, maxe1, maxe2};
    });
  }
  auto hv = reduce_data.value(reduce_op);
  h_max_eigenvals[0] = amrex::get<0>(hv);
  h_max_eigenvals[1] = amrex::get<1>(hv);
  h_max_eigenvals[2] = amrex::get<2>(hv);
#endif

  // Communicate across processors
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    ParallelDescriptor::ReduceRealMax(h_max_eigenvals[idir]);
  }

  return h_max_eigenvals;
}
// ------------------------------------------------------------------------------------//
void CNS::post_timestep(int /* iteration */) {
  BL_PROFILE("post_timestep");
  
  //amrex::Print() << " oo CNS::post_timestep "  << std::endl;

  if (do_reflux && level < parent->finestLevel()) {
    CNS& fine_level = getLevel(level + 1);
    MultiFab& S_crse = get_new_data(State_Type);
    const int ncomp = PROB::ProbClosures::NCONS;
#if CNS_USE_EB
    MultiFab& S_fine = fine_level.get_new_data(State_Type);
    const MultiFab& volfrac_crse = *EBM::eb.volmf_a[level];
    const MultiFab& volfrac_fine = *EBM::eb.volmf_a[level + 1];
    fine_level.flux_reg.Reflux(S_crse, volfrac_crse, S_fine, volfrac_fine, 0, 0, ncomp);
#else
    fine_level.flux_reg.Reflux(S_crse, 0, 0, ncomp);
#endif
  }

  if (level < parent->finestLevel()) {
    avgDown();
    getLevel(level + 1).resetFillPatcher();
  }

  // Record time statistics
  if (record_probe) {
    recordTimeProbe();    
  }

  // Record statistics
  if (record_stats) {
    time_stat_level[level] += parent->dtLevel(level);
    computeStats();
  }
}

// ------------------------------------------------------------------------------------//
void CNS::postCoarseTimeStep(Real time) {

 // amrex::Print() << " oo CNS::postCoarseTimeStep " << time <<  std::endl;

#ifdef AMREX_USE_GPIBM
  // Surface output at user-specified step interval
  const int istep = parent->levelSteps(0);
  if (plot_surf && (istep % surf_int == 0)) {
      for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
          dynamic_cast<CNS&>(parent->getLevel(lev)).writeSurfFile();
      }
  }

#ifdef CNS_USE_FSI
  // Calculate and print FSI loads and properties
  {
    auto& ib = IBM::ib;
    if (ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "\n=== FSI Loads (Step " << istep
                       << ", Time " << time << ") ===\n";
    }
    for (int i = 0; i < ib.ngeom; ++i) {
        auto props = FSI::RigidBodyProperties::readOrCompute(ib.geom_a[i], i);
        auto loads = FSI::Kinematics::computeLoads(i, props.xcenter);
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Geometry " << i << ":\n"
                           << "  Mass: " << props.mass << "\n"
                           << "  Center of Mass: " << props.xcenter << "\n"
                           << "  Inertia Tensor:\n";
            for (int r = 0; r < 3; ++r) {
                amrex::Print() << "    [ " << props.inertia[r][0] << ", "
                               << props.inertia[r][1] << ", "
                               << props.inertia[r][2] << " ]\n";
            }
            amrex::Print() << "  Fluid Force: " << loads.force << "\n"
                           << "  Fluid Moment (about CM): " << loads.moment << "\n"
                           << "----------------------------------------\n";
        }
    }
  }
#endif  // CNS_USE_FSI
#endif  // AMREX_USE_GPIBM

  if (verbose && ((this->nStep() % nstep_screen_output) == 0)) {
    printTotal();
  }
}
// -----------------------------------------------------------------------------

// Gridding -------------------------------------------------------------------
// Called for each level from 0,1...nlevs-1

void CNS::post_regrid(int lbase, int new_finest) {

#ifdef AMREX_USE_GPIBM
  rebuildIBM();

  // ==========================================================================
  // FSI post-regrid state cleanup.
  //
  // After regrid, AMReX::FillPatch has populated the new-grid conservative
  // state by interpolating from the old grids. If the old grids had any
  // stale or extreme values in solid cells (which don't evolve during RK),
  // those values can contaminate freshly-created fine cells near the solid
  // boundary. This causes downstream WENO reconstruction to explode.
  //
  // Two passes:
  //   [A] Flood-fill solid cells from valid fluid/ghost neighbors so they
  //       carry bounded, physically plausible data.
  //   [B] Zero momentum in interior solid cells (same rationale as
  //       end-of-step pass in advance.cpp).
  //
  // Static geometry doesn't need this: solid cells never transition and
  // their data stays consistent with the (unchanging) flow around them.
  // ==========================================================================
  if (ib_move) {
    MultiFab& S = get_new_data(State_Type);
    auto& ib_mf = *IBM::ib.bmf_a[level];
    const int ncons = d_prob_closures->NCONS;
    constexpr int MAX_FLOOD_ITER = 32;

    for (MFIter mfi(S, false); mfi.isValid(); ++mfi) {
      const Box& bx  = mfi.tilebox();
      const Box& bxg = mfi.growntilebox(d_prob_closures->NGHOST);
      auto const& state = S.array(mfi);
      auto const& mk    = ib_mf.const_array(mfi);

      // Tag array: 0 = fluid or ghost point (valid)
      //            1 = interior solid (needs fixing)
      //            2 = solid, already fixed this iteration
      BaseFab<int> tagfab(bxg, 1, The_Managed_Arena());
      auto const& tag = tagfab.array();

      ParallelFor(bxg, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tag(i,j,k) = (mk(i,j,k,0) != 0 && mk(i,j,k,1) == 0) ? 1 : 0;
      });

      // [A] Iterative flood fill — propagate valid data inward one ring
      //     per iteration until all solid cells are reached (or budget runs out).
      for (int iter = 0; iter < MAX_FLOOD_ITER; ++iter) {
        Gpu::DeviceScalar<int> d_nfixed(0);
        int* p_nfixed = d_nfixed.dataPtr();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          if (tag(i,j,k) != 1) return;

          Real sum[PROB::ProbClosures::NCONS] = {};
          int count = 0;
          for (int dj = -1; dj <= 1; ++dj) {
            for (int di = -1; di <= 1; ++di) {
#if (AMREX_SPACEDIM == 3)
              for (int dk = -1; dk <= 1; ++dk) {
#else
              { int dk = 0;
#endif
                if (di == 0 && dj == 0 && dk == 0) continue;
                const int ii = i+di, jj = j+dj, kk = k+dk;
                if (!bxg.contains(IntVect(AMREX_D_DECL(ii,jj,kk)))) continue;
                if (tag(ii,jj,kk) == 0 || tag(ii,jj,kk) == 2) {
                  for (int n = 0; n < ncons; ++n)
                    sum[n] += state(ii,jj,kk,n);
                  count++;
                }
              }
            }
          }
          if (count > 0) {
            const Real inv = Real(1.0) / count;
            for (int n = 0; n < ncons; ++n)
              state(i,j,k,n) = sum[n] * inv;
            tag(i,j,k) = 2;
            Gpu::Atomic::Add(p_nfixed, 1);
          }
        });

        Gpu::streamSynchronize();
        if (d_nfixed.dataValue() == 0) break;
      }

      // [B] Zero momentum in interior solid cells to prevent spurious
      //     velocity amplification from acoustic-phase averaging.
      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (mk(i,j,k,0) == 0) return;  // fluid
        if (mk(i,j,k,1) != 0) return;  // ghost point

        using PC = PROB::ProbClosures;
        const Real rho = state(i,j,k, PC::URHO);
        if (rho <= Real(0)) return;

        const Real mx = state(i,j,k, PC::UMX);
        const Real my = state(i,j,k, PC::UMY);
#if (AMREX_SPACEDIM == 3)
        const Real mz = state(i,j,k, PC::UMZ);
        const Real ke = Real(0.5) * (mx*mx + my*my + mz*mz) / rho;
#else
        const Real ke = Real(0.5) * (mx*mx + my*my) / rho;
#endif
        state(i,j,k, PC::UMX) = Real(0);
        state(i,j,k, PC::UMY) = Real(0);
#if (AMREX_SPACEDIM == 3)
        state(i,j,k, PC::UMZ) = Real(0);
#endif
        state(i,j,k, PC::UET) -= ke;
      });
    }
  }  // end if (ib_move)
#endif  // AMREX_USE_GPIBM

#ifdef CNS_USE_EB
  EBM::eb.destroy_mf(level);
  EBM::eb.build_mf(grids, dmap, level);

  // update volfrac and relevant EB data
  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  EBM::eb.volmf_a[level]  = &(ebfactory.getVolFrac()); 
  EBM::eb.normmcf_a[level] = &(ebfactory.getBndryNormal());
  EBM::eb.areamcf_a[level] = ebfactory.getAreaFrac();  
  EBM::eb.ebflags_a[level] = &(ebfactory.getMultiEBCellFlagFab());
  EBM::eb.bcareamcf_a[level] = &(ebfactory.getBndryArea());
  EBM::eb.bndrycent_a[level] = &(ebfactory.getBndryCent());
  EBM::eb.volcent_a[level] = &(ebfactory.getCentroid());
  

  // Level mask for redistribution (stored as object not pointer)
  EBM::eb.level_mask_a[level].clear();
  EBM::eb.level_mask_a[level].define(grids, dmap, 1, 3);
  EBM::eb.level_mask_a[level].BuildMask(
        geom.Domain(), geom.periodicity(), CNSConstants::level_mask_covered,
        CNSConstants::level_mask_notcovered, CNSConstants::level_mask_physbnd,
        CNSConstants::level_mask_interior);

  // EBM::eb.facecent  = ebfactory.getFaceCent();

  // Calculate markers  
  EBM::eb.computeMarkers(level);
  
  EBM::eb.check_geometry(level);


#endif

}

void CNS::errorEst(TagBoxArray &tags, int /*clearval*/, int /*tagval*/,
                   Real time, int /*n_error_buf*/, int /*ngrow*/) {

 // amrex::Print() << " oo CNS::errorEst " << std::endl;

  // MF without ghost points filled (why?)
  MultiFab sdata(get_new_data(State_Type).boxArray(),
                 get_new_data(State_Type).DistributionMap(), PROB::ProbClosures::NCONS, PROB::ProbClosures::NGHOST,
                 MFInfo(), Factory());

  // filling ghost points (copied from PeleC)
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this, sdata, PROB::ProbClosures::NGHOST, cur_time, State_Type, 0, PROB::ProbClosures::NCONS);
  const auto geomdata = geom.data();

  
  // fill ghost points of stats (SNM needed? )
  // if (compute_stats) {
  //   MultiFab sdata_stat(get_new_data(Stats_Type).boxArray(),
  //                get_new_data(Stats_Type).DistributionMap(), PROB::ProbClosures::NSTAT, PROB::ProbClosures::NGHOST,
  //                MFInfo(), Factory());
  //   FillPatch(*this, sdata_stat, PROB::ProbClosures::NGHOST, cur_time, Stats_Type, 0, PROB::ProbClosures::NSTAT);
  
  // }

#ifdef AMREX_USE_GPIBM
  auto &ibdata = (*IBM::ib.bmf_a[level]);
#elif defined CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif
  for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    auto const &tagfab = tags.array(mfi);
    auto const &sdatafab = sdata.array(mfi);
#ifdef AMREX_USE_GPIBM
    auto const &ibfab = ibdata.array(mfi); // was const_array for user_tagging???
#elif defined CNS_USE_EB
    auto const& flag = flags.const_array(mfi);
#endif
    int lev = level;
    int nt_lev = nStep();
    PROB::ProbParm const *lprobparm = d_prob_parm;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // the user_tagging function is defined in prob.h
#ifdef AMREX_USE_GPIBM
      user_tagging(i, j, k, nt_lev, tagfab, sdatafab, ibfab, geomdata,*lprobparm, lev);
#elif defined CNS_USE_EB
      user_tagging(i, j, k, nt_lev, tagfab, sdatafab, flag, geomdata, *lprobparm, lev);
#else
      user_tagging(i, j, k, nt_lev, tagfab, sdatafab, geomdata ,*lprobparm, lev);
#endif
    });
  }
}

void CNS::post_restart() {

// recreate markers
amrex::Print() << " recreate markers " << std::endl;


#ifdef AMREX_USE_GPIBM
  rebuildIBM();
#endif

#ifdef CNS_USE_EB
  EBM::eb.destroy_mf(level);
  EBM::eb.build_mf(grids, dmap, level);

  // update volfrac and relevant EB data
  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

  EBM::eb.volmf_a[level]  = &(ebfactory.getVolFrac()); 
  EBM::eb.normmcf_a[level] = &(ebfactory.getBndryNormal());
  EBM::eb.areamcf_a[level] = ebfactory.getAreaFrac();  
  EBM::eb.ebflags_a[level] = &(ebfactory.getMultiEBCellFlagFab());
  EBM::eb.bcareamcf_a[level] = &(ebfactory.getBndryArea());
  EBM::eb.bndrycent_a[level] = &(ebfactory.getBndryCent());
  EBM::eb.volcent_a[level] = &(ebfactory.getCentroid());
  

  // Level mask for redistribution (stored as object not pointer)
  EBM::eb.level_mask_a[level].clear();
  EBM::eb.level_mask_a[level].define(grids, dmap, 1, 3);
  EBM::eb.level_mask_a[level].BuildMask(
        geom.Domain(), geom.periodicity(), CNSConstants::level_mask_covered,
        CNSConstants::level_mask_notcovered, CNSConstants::level_mask_physbnd,
        CNSConstants::level_mask_interior);

  // EBM::eb.facecent  = ebfactory.getFaceCent();

  // Calculate markers  
  EBM::eb.computeMarkers(level);

  EBM::eb.check_geometry(level);

#endif

  // WARNING !! UBC not read from restart, it will be rebuild    
  if (use_nscbc) {
    MultiFab& S_new = get_new_data(State_Type);
    define_nscbc_ghost_shell();
    initialise_nscbc_ghost_shell(S_new);
  }

  // Initialise stats arrays to zero when restarting from a checkpoint
  // that was written without statistics (NSTAT was 0 in the old build).
  if (compute_stats) {
    if (!state[Stats_Type].hasNewData()) {
      const Real cur_time = state[State_Type].curTime();
      const Real dt_old   = cur_time - state[State_Type].prevTime();
      state[Stats_Type].define(geom.Domain(), grids, dmap,
                               desc_lst[Stats_Type], cur_time, dt_old,
                               Factory());
    }
    setupStats();
    time_stat_level[level] = 0.0;
  }

  // Set up diagnostics after restart
  if (record_probe) {
    setupTimeProbe();
  }




}

void CNS::set_state_in_checkpoint(Vector<int>& state_in_checkpoint) {
  // This is only called when the checkpoint has fewer state types than
  // the current code.  Mark Stats_Type as absent so AMReX skips reading it.
  if (compute_stats) {
    state_in_checkpoint[Stats_Type] = 0;
  }
}

int CNS::okToContinue()
{
  if (level > 0) { return 1; }

  int test = 1;
  MultiFab &S = get_new_data(State_Type);
  if (S.contains_nan(0, S.nComp())) {
    test = 0;
  }

  return test;
}

// 
void CNS::avgDown() {
  BL_PROFILE("CNS::avgDown()");

  if (level == parent->finestLevel()) return;

  auto &fine_lev = getLevel(level + 1);

  MultiFab &S_crse = get_new_data(State_Type);
  MultiFab &S_fine = fine_lev.get_new_data(State_Type);

  // Standard avgDown first.
  amrex::average_down(S_fine, S_crse, fine_lev.geom, geom, 0, S_fine.nComp(),
                      parent->refRatio(level));
  Gpu::streamSynchronize();  // ensure GPU average_down is complete before CPU read

#ifdef AMREX_USE_GPIBM
  // IBM-aware avgDown correction.
  // Standard average_down has already run above. Now correct coarse cells
  // that overlap MIXED fine regions (some solid + some fluid fine sub-cells)
  // by re-averaging using only fluid fine cells. Pure-fluid coarse cells
  // keep the standard average; pure-solid coarse cells keep their own value.
  if (IBM::ib.bmf_a[level + 1] != nullptr)
  {
    auto& fine_ibmf = *IBM::ib.bmf_a[level + 1];
    const int nc = S_fine.nComp();
    const IntVect rr = parent->refRatio(level);

    // Build coarsened box array on fine's DistributionMap
    BoxArray cba = S_fine.boxArray();
    cba.coarsen(rr);

    // Use ncons+1 components: first ncons are data, last is validity flag
    const int ncp1 = nc + 1;
    MultiFab S_corr(cba, S_fine.DistributionMap(), ncp1, 0);
    S_corr.setVal(Real(0.0));

    // Populate S_corr with fluid-only averages for MIXED cells
    // Note: use GPU ParallelFor because in CUDA builds, MultiFab data
    // lives in device memory and can't be accessed via raw CPU loops.
    for (MFIter fmfi(S_fine, false); fmfi.isValid(); ++fmfi) {
      const Box fbx = fmfi.tilebox();
      const Box cbx = amrex::coarsen(fbx, rr);
      auto const fine = S_fine.const_array(fmfi);
      auto const fmk  = fine_ibmf.const_array(fmfi);
      auto const corr = S_corr.array(fmfi);
      const int ncomp = nc;
      const IntVect ratio = rr;

      amrex::ParallelFor(cbx,
      [=] AMREX_GPU_DEVICE (int ci, int cj, int ck) noexcept
      {
        int n_fluid = 0, n_solid = 0;
        Real sum[PROB::ProbClosures::NCONS] = {};
#if (AMREX_SPACEDIM == 2)
        for (int fj = cj*ratio[1]; fj < (cj+1)*ratio[1]; ++fj)
        for (int fi = ci*ratio[0]; fi < (ci+1)*ratio[0]; ++fi) {
          if (fmk(fi,fj,0,0) == 0) {
            for (int n = 0; n < ncomp; ++n) sum[n] += fine(fi,fj,0,n);
            ++n_fluid;
          } else ++n_solid;
        }
#else
        for (int fk = ck*ratio[2]; fk < (ck+1)*ratio[2]; ++fk)
        for (int fj = cj*ratio[1]; fj < (cj+1)*ratio[1]; ++fj)
        for (int fi = ci*ratio[0]; fi < (ci+1)*ratio[0]; ++fi) {
          if (fmk(fi,fj,fk,0) == 0) {
            for (int n = 0; n < ncomp; ++n) sum[n] += fine(fi,fj,fk,n);
            ++n_fluid;
          } else ++n_solid;
        }
#endif
        // Only overwrite MIXED cells
        if (n_fluid > 0 && n_solid > 0) {
          Real inv = Real(1.0) / n_fluid;
          for (int n = 0; n < ncomp; ++n) corr(ci,cj,ck,n) = sum[n] * inv;
          corr(ci,cj,ck,ncomp) = Real(1.0);  // mark valid
        }
      });
    }

    // ParallelCopy S_corr to a MultiFab aligned with S_crse
    MultiFab S_corr_aligned(S_crse.boxArray(), S_crse.DistributionMap(), ncp1, 0);
    S_corr_aligned.setVal(Real(0.0));
    S_corr_aligned.ParallelCopy(S_corr, 0, 0, ncp1);

    // Apply corrections: where validity flag == 1, overwrite S_crse
    for (MFIter mfi(S_crse, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& crse = S_crse.array(mfi);
      auto const& corr = S_corr_aligned.const_array(mfi);
      const int ncomp = nc;

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (corr(i,j,k,ncomp) > Real(0.5)) {
          for (int n = 0; n < ncomp; ++n)
            crse(i,j,k,n) = corr(i,j,k,n);
        }
      });
    }
  }
#endif

  if (compute_stats) {
    MultiFab &Sstat_crse = get_new_data(Stats_Type);
    MultiFab &Sstat_fine = fine_lev.get_new_data(Stats_Type);
    amrex::average_down(Sstat_fine, Sstat_crse, fine_lev.geom, geom, 0, Sstat_fine.nComp(),
                      parent->refRatio(level));
  }

}

void CNS::printTotal() const {
  // Get conservatives multifab
  const MultiFab& consmf = get_new_data(State_Type);

  // Volume-weighted integral of conserved variables (works for Cartesian and RZ)
  MultiFab volume(consmf.boxArray(), consmf.DistributionMap(), 1, 0);
  geom.GetVolume(volume, consmf.boxArray(), consmf.DistributionMap(), 0);

  std::array<Real, PROB::ProbClosures::NCONS> tot{};
  for (int comp = 0; comp < PROB::ProbClosures::NCONS; ++comp) {
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    auto const& a = consmf.const_arrays();
    auto const& v = volume.const_arrays();
    reduce_op.eval(consmf, IntVect(0), reduce_data,
                   [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept -> Real {
                     return a[box_no](i, j, k, comp) * v[box_no](i, j, k);
                   });
    Gpu::streamSynchronize();
    auto const& hv = reduce_data.value(reduce_op);
    tot[comp] = amrex::get<0>(hv);
  }

  // Communicate across processors
  ParallelDescriptor::ReduceRealSum(tot.data(), PROB::ProbClosures::NCONS,
                                    ParallelDescriptor::IOProcessorNumber());
  // Print
  Vector<std::string> names= PROB::ProbClosures::get_cons_vars_names();
  for (int comp = 0; comp < PROB::ProbClosures::NCONS; ++comp) {
    amrex::Print().SetPrecision(17) << "   Total " << names[comp] << " = " << tot[comp] << "\n";
  }
}

void CNS::variableCleanUp() {
  delete h_prob_closures;
  delete h_phys_bc;

#ifdef AMREX_USE_GPU
  The_Arena()->free(d_prob_closures);
  The_Arena()->free(d_phys_bc);
#endif
  desc_lst.clear();
  derive_lst.clear();
}

// Plotting
//------------------------------------------------------------------------------
void CNS::writePlotFile(const std::string &dir, std::ostream &os,
                        VisMF::How how) {
  int i, n;
  //
  // The list of indices of State to write to plotfile.
  // first component of pair is state_type,
  // second component of pair is component # within the state_type
  //
  std::vector<std::pair<int, int>> plot_var_map;
  for (int typ = 0; typ < desc_lst.size(); typ++) {
    for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
      if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
          desc_lst[typ].getType() == IndexType::TheCellType()) {
        plot_var_map.push_back(std::pair<int, int>(typ, comp));
      }
    }
  }

  int num_derive = 0;
  std::vector<std::string> derive_names;
  const std::list<DeriveRec> &dlist = derive_lst.dlist();
  for (auto const &d : dlist) {
    if (parent->isDerivePlotVar(d.name())) {
      derive_names.push_back(d.name());
      num_derive += d.numDerive();
    }
  }

  int n_data_items = plot_var_map.size() + num_derive;

//----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
  n_data_items += 2;
#endif
#ifdef CNS_USE_EB
  n_data_items += 1;
#endif
  //------------------------------------------------------------------------------

  // get the time from the first State_Type
  // if the State_Type is ::Interval, this will get t^{n+1/2} instead of t^n
  Real cur_time = state[0].curTime();

  if (level == 0 && ParallelDescriptor::IOProcessor()) {
    //
    // The first thing we write out is the plotfile type.
    //
    os << thePlotFileType() << '\n';

    if (n_data_items == 0)
      amrex::Error("Must specify at least one valid data item to plot");

    os << n_data_items << '\n';

    //
    // Names of variables
    //
    for (i = 0; i < static_cast<int>(plot_var_map.size()); i++) {
      int typ = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      os << desc_lst[typ].name(comp) << '\n';
    }

    // derived
    for (auto const &dname : derive_names) {
      const DeriveRec *rec = derive_lst.get(dname);
      for (i = 0; i < rec->numDerive(); ++i) {
        os << rec->variableName(i) << '\n';
      }
    }

    //----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
    os << "sld\n";
    os << "ghs\n";
#endif
#ifdef CNS_USE_EB
    os << "vfrac\n";
#endif

    //------------------------------------------------------------------------------

    os << AMREX_SPACEDIM << '\n';
    os << parent->cumTime() << '\n';
    int f_lev = parent->finestLevel();
    os << f_lev << '\n';
    for (i = 0; i < AMREX_SPACEDIM; i++) os << Geom().ProbLo(i) << ' ';
    os << '\n';
    for (i = 0; i < AMREX_SPACEDIM; i++) os << Geom().ProbHi(i) << ' ';
    os << '\n';
    for (i = 0; i < f_lev; i++) os << parent->refRatio(i)[0] << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) os << parent->Geom(i).Domain() << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) os << parent->levelSteps(i) << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) {
      for (int k = 0; k < AMREX_SPACEDIM; k++)
        os << parent->Geom(i).CellSize()[k] << ' ';
      os << '\n';
    }
    os << (int)Geom().Coord() << '\n';
    os << "0\n";  // Write bndry data.
  }
  // Build the directory to hold the MultiFab at this level.
  // The name is relative to the directory containing the Header file.
  //
  static const std::string BaseName = "/Cell";
  char buf[64];
  snprintf(buf, sizeof buf, "Level_%d", level);
  std::string sLevel = buf;
  //
  // Now for the full pathname of that directory.
  //
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.size() - 1] != '/') {
    FullPath += '/';
  }
  FullPath += sLevel;
  //
  // Only the I/O processor makes the directory if it doesn't already exist.
  //
  if (!levelDirectoryCreated) {
    if (ParallelDescriptor::IOProcessor()) {
      if (!amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier();
  }

  if (ParallelDescriptor::IOProcessor()) {
    os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
    os << parent->levelSteps(level) << '\n';

    for (i = 0; i < grids.size(); ++i) {
      RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
      for (n = 0; n < AMREX_SPACEDIM; n++)
        os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
    }
    //
    // The full relative pathname of the MultiFabs at this level.
    // The name is relative to the Header file containing this name.
    // It's the name that gets written into the Header.
    //
    if (n_data_items > 0) {
      std::string PathNameInHeader = sLevel;
      PathNameInHeader += BaseName;
      os << PathNameInHeader << '\n';
    }
  }
  //
  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  int cnt = 0;
  const int nGrow = 0;
  MultiFab plotMF(grids, dmap, n_data_items, nGrow, MFInfo(), Factory());
  MultiFab *this_dat = 0;
  //
  // Cull data from state variables -- use no ghost cells.
  //
  for (i = 0; i < static_cast<int>(plot_var_map.size()); i++) {
    int typ = plot_var_map[i].first;
    int comp = plot_var_map[i].second;
    this_dat = &state[typ].newData();
    MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
    cnt++;
  }

  // derived
  if (derive_names.size() > 0) {
    for (auto const &dname : derive_names) {
      derive(dname, cur_time, plotMF, cnt);
      cnt += derive_lst.get(dname)->numDerive();
    }  //exit(1);
  }

  //------------------------------------------------------------------------------
  // additional plotting ...
  //------------------------------------------------------------------------------
  
#ifdef AMREX_USE_GPIBM
  plotMF.setVal(0.0_rt, cnt, 2, nGrow);
  IBM::ib.bmf_a[level]->copytoRealMF(plotMF, 0, cnt);
  cnt+=2;
#endif

#ifdef CNS_USE_EB
  plotMF.setVal(0.0_rt, cnt, 0, nGrow); 
 // EBM::eb.bmf_a[level]->copytoRealMF(plotMF, 0, cnt);  // boolean 
  const MultiFab *vfrac = EBM::eb.volmf_a[level];
  MultiFab::Copy(plotMF, *vfrac, 0, cnt, 1, 0);
  cnt++;
#endif

  //------------------------------------------------------------------------------

  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  if (AsyncOut::UseAsyncOut()) {
    VisMF::AsyncWrite(plotMF, TheFullPath);
  } else {
    VisMF::Write(plotMF, TheFullPath, how, true);
  }

  levelDirectoryCreated = false;  // ---- now that the plotfile is finished
}

// This is called once per level on write timestep.
// surf_int must be the same as plot_int
void CNS::writePlotFilePost(const std::string &dir, std::ostream &os) {

#if AMREX_USE_GPIBM
  // writeSurfFile();
#endif

}

// this subroutine is called from the main loop 
// should be called per level
#if AMREX_USE_GPIBM

void CNS::rebuildIBM() {
  IBM::ib.destroy_mf(level);
  IBM::ib.build_mf(grids, dmap, level);
  IBM::ib.computeMarkers(level);
  IBM::ib.initialiseGPs(level);
  // Surface indices depend on all levels having valid bmf_a, so we can
  // only rebuild them once the entire regrid cascade is complete (i.e.,
  // when the finest level calls rebuildIBM).
  if (level == parent->finestLevel()) {
     for (int lev = parent->finestLevel(); lev >= 0; --lev) {
        IBM::ib.computeSurfIndices(lev);
     }
  }
}
//------------------------------------------------------------------------------
//
//
//------------------------------------------------------------------------------
void CNS::writeSurfFile() {
      
  // calculate and  write surface data  
  int istep = parent->levelSteps(0);

  if (plot_surf && (istep % surf_int == 0))  {
     
    MultiFab& Sdata = get_new_data(State_Type); 

    int ncons = CNS::d_prob_closures->NCONS;
    int nghost= CNS::d_prob_closures->NGHOST;

    Real time = parent->cumTime();

    if (this->level == parent->maxLevel()) {
      Print() << "Computing surface properties ";
      Print() << " at time= " << time << " and step= " << istep << std::endl;
    }
    
    FillPatch(*this, Sdata, nghost, time, State_Type, 0, ncons);

    // Convert conservative to primitive variables for surface interpolation
    int nprim = PROB::ProbClosures::NPRIM;
    MultiFab prims_mf(Sdata.boxArray(), Sdata.DistributionMap(),
                      nprim, nghost, MFInfo().SetArena(The_Async_Arena()));
    for (MFIter mfi(Sdata, false); mfi.isValid(); ++mfi) {
      CNS::h_prob_closures->cons2prims(mfi, Sdata.array(mfi), prims_mf.array(mfi));
    }

    const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
    const PROB::ProbClosures* cls_h = CNS::h_prob_closures; 

    IBM::ib.computeSURFs(prims_mf,cls_d,this->level); // computed at each level. From low to high.

    // Only gather and write on the finest level to ensure all levels are processed
    if (this->level == parent->finestLevel()){
      // collect data to rank 0
      IBM::ib.gatherSurfData(); 

      if (amrex::ParallelDescriptor::IOProcessor()){
        IBM::ib.plotSURF(time, istep, surf_filename); 
      } 
    }

  }
}
#endif

//------------------------------------------------------------------------------
//
//
//------------------------------------------------------------------------------
void CNS::define_nscbc_ghost_shell()
{
    const int ncons  = PROB::ProbClosures::NCONS;
    const int nghost = PROB::ProbClosures::NGHOST;

    nscbc_shell.state = std::make_unique<MultiFab>(grids, dmap, ncons, nghost, MFInfo(), Factory());

    nscbc_shell.rhs = std::make_unique<MultiFab>(grids, dmap, ncons, nghost, MFInfo(), Factory());

    nscbc_shell.owner = std::make_unique<iMultiFab>(grids, dmap, 1, nghost);

    nscbc_shell.state->setVal(Real(0.0));
    nscbc_shell.rhs->setVal(Real(0.0));
    nscbc_shell.owner->setVal(0);

    build_nscbc_shell_owner();
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::initialise_nscbc_ghost_shell(
    MultiFab const& cell_state)
{
    MultiFab& shell = *nscbc_shell.state;

    const Box domain = Geom().Domain();
    const int ng     = shell.nGrow();
    const int ncons  = PROB::ProbClosures::NCONS;

    /*
     * Construct a synchronized source with enough ghost cells.
     * This is essential at tangential FAB interfaces on a physical boundary.
     */
    MultiFab source(cell_state.boxArray(),cell_state.DistributionMap(),ncons,ng,MFInfo(),Factory());

    source.setVal(Real(0.0));

    MultiFab::Copy(source,cell_state,0, 0,ncons,0);             // copy valid cells

    source.FillBoundary(Geom().periodicity());

    shell.setVal(Real(0.0));

    for (MFIter mfi(shell, false); mfi.isValid(); ++mfi) {

        auto const& U   = source.const_array(mfi);
        auto const& G   = shell.array(mfi);
        auto const& own = nscbc_shell.owner->const_array(mfi);

        const Box bx = mfi.fabbox();

        ParallelFor( bx, ncons, [=] AMREX_GPU_DEVICE( int i, int j, int k, int n) noexcept
            {
                const IntVect iv(AMREX_D_DECL(i,j,k));
                if (own(iv,0) == 0) {return;}

                IntVect src = iv;

                /*
                 * Extrapolate only in physical-boundary directions.
                 * Tangential inter-FAB coordinates remain unchanged and
                 * are supplied by source.FillBoundary().
                 */
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    src[dir] = amrex::max( domain.smallEnd(dir), amrex::min(domain.bigEnd(dir), src[dir]));
                }

                G(iv,n) = U(src,n);
            });
    }
    Gpu::streamSynchronize();
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::saxpy_nscbc_shell(
    MultiFab& dst,
    Real a,
    MultiFab const& src)
{
    const int ncons = PROB::ProbClosures::NCONS;

    for (MFIter mfi(dst, false); mfi.isValid(); ++mfi) {

        auto const& d   = dst.array(mfi);
        auto const& s   = src.const_array(mfi);
        auto const& own =
            nscbc_shell.owner->const_array(mfi);

        const Box bx = mfi.fabbox();

        ParallelFor(
            bx, ncons,
            [=] AMREX_GPU_DEVICE(
                int i, int j, int k, int n) noexcept
            {
                if (own(i,j,k,0) != 0) {
                    d(i,j,k,n) += a*s(i,j,k,n);
                }
            });
    }
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
void CNS::lincomb_nscbc_shell(
    MultiFab& dst,
    Real a,
    MultiFab const& A,
    Real b,
    MultiFab const& B)
{
    const int ncons = PROB::ProbClosures::NCONS;

    for (MFIter mfi(dst, false); mfi.isValid(); ++mfi) {

        auto const& d   = dst.array(mfi);
        auto const& x   = A.const_array(mfi);
        auto const& y   = B.const_array(mfi);
        auto const& own =
            nscbc_shell.owner->const_array(mfi);

        const Box bx = mfi.fabbox();

        ParallelFor(
            bx, ncons,
            [=] AMREX_GPU_DEVICE(
                int i, int j, int k, int n) noexcept
            {
                if (own(i,j,k,0) != 0) {
                    d(i,j,k,n) =
                        a*x(i,j,k,n)
                      + b*y(i,j,k,n);
                }
            });
    }
}

