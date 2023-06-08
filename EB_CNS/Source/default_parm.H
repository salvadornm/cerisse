// #ifndef CNS_DEFAULT_PARM_H //No need this as this file will be included only once, in CNS.cpp
// #define CNS_DEFAULT_PARM_H

BCRec    CNS::phys_bc; // holds BC info

int      CNS::verbose   = 2;
Real     CNS::cfl       = 0.3;
Real     CNS::dt_cutoff = 5.e-20;

int             CNS::do_reflux = 1;

int             CNS::refine_cutcells = 1;
int             CNS::refine_cutcells_max_lev = 10;
Vector<RealBox> CNS::refine_boxes;
Vector<int>     CNS::refine_boxes_max_lev;
RealBox*        CNS::dp_refine_boxes;
int             CNS::refine_dengrad_max_lev = -1;
Vector<Real>    CNS::refine_dengrad;
int             CNS::refine_velgrad_max_lev = -1;
Vector<Real>    CNS::refine_velgrad;
int             CNS::refine_magvort_max_lev = -1;
Vector<Real>    CNS::refine_magvort;
int             CNS::refine_tke_max_lev = -1;
Vector<Real>    CNS::refine_tke;
// TaggingParm* CNS::tag_parm;
// Vector<AMRErrorTag> CNS::errtagger;

bool     CNS::do_visc        = true;  // diffusion is on by default
bool     CNS::do_ext_src     = false; // call external source term in prob.cpp

bool     CNS::do_react       = false; // reaction is off by default
std::string CNS::chem_integrator = "Reactor_Null";
pele::physics::transport::TransportParams<
    pele::physics::PhysicsType::transport_type> CNS::trans_parms;
amrex::Vector<std::string> CNS::spec_names;
bool     CNS::use_typical_vals_chem = false; // tell chem_integrator typical value of Temp
int      CNS::reset_typical_vals_int = 10;   // interval to reset the typical value
int      CNS::rk_reaction_iter = 0;          // iterate to tightly couple chemistry
Real     CNS::min_react_temp = 300.0;        // turn off reactor below this temperature (only work with CVODE)
Real     CNS::clip_temp = std::numeric_limits<Real>::epsilon();
bool     CNS::update_heat_release = true;    // write to heat release in I_R? switch it off for some performance gain (approx 3-5%)

int      CNS::char_sys       = 1;   // Thermodynamic system used in characteristic reconstruction. 0:speed of sound, 1:gamma
int      CNS::recon_scheme   = 5;   // 1: basic Godunov; 2: MUSCL; 3: WENO-Z3; 4: WENO-JS5; 5: WENO-Z5; 6: TENO5
Real     CNS::plm_theta      = 2.0; // [1,2] 1: minmod; 2: van Leer's MC (higher sharper)

bool     CNS::do_restart_fields = false;
bool     CNS::do_psgs = false;
bool     CNS::do_vpdf = false;
bool     CNS::do_spdf = false;

#if CNS_USE_EB
int   CNS::eb_weights_type     = 3; // [0,1,2,3] 0: weights are all 1; 1: eint; 2: cell mass; 3: volfrac
int   CNS::do_reredistribution = 1;
bool  CNS::eb_no_slip          = true;
bool  CNS::eb_isothermal       = false;
Real  CNS::eb_wall_temp; // EB wall temperature
#endif

bool  CNS::signalStopJob = false;

// #endif