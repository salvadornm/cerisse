// #ifndef CNS_DEFAULT_PARM_H //No need this as this file will be included only once,
// in CNS.cpp #define CNS_DEFAULT_PARM_H

BCRec CNS::phys_bc; // holds BC info

int CNS::verbose = 2;
Real CNS::cfl = 0.3;
Real CNS::fixed_dt = -1;
Real CNS::dt_cutoff = 5.e-20;
int CNS::rk_order = 2;
int CNS::amr_interp_order = 2;

bool CNS::do_reflux = true;
bool CNS::do_load_balance = true;

int CNS::refine_cutcells = 1;
int CNS::refine_cutcells_max_lev = 10;
Vector<RealBox> CNS::refine_boxes;
Vector<int> CNS::refine_boxes_max_lev;
RealBox* CNS::dp_refine_boxes;
int CNS::refine_dengrad_max_lev = -1;
Vector<Real> CNS::refine_dengrad;
int CNS::refine_velgrad_max_lev = -1;
Vector<Real> CNS::refine_velgrad;
int CNS::refine_presgrad_max_lev = -1;
Vector<Real> CNS::refine_presgrad;
int CNS::refine_magvort_max_lev = -1;
Vector<Real> CNS::refine_magvort;
int CNS::refine_tke_max_lev = -1;
Vector<Real> CNS::refine_tke;
// TaggingParm* CNS::tag_parm;
// Vector<AMRErrorTag> CNS::errtagger;

RealBox CNS::buffer_box;

bool CNS::do_hydro = true;   // hydrodynamic fluxes
bool CNS::do_visc = true;    // diffusion is on by default
bool CNS::do_ext_src = true; // call external source term in prob.cpp

bool CNS::do_react = false; // reaction is off by default
std::string CNS::chem_integrator = "Reactor_Null";
pele::physics::transport::TransportParams<pele::physics::PhysicsType::transport_type>
  CNS::trans_parms;
amrex::Vector<std::string> CNS::spec_names;
pele::physics::turbinflow::TurbInflow CNS::turb_inflow;
pele::physics::PMF::PmfData CNS::pmf_data;

bool CNS::do_nscbc = false;
amrex::Real CNS::nscbc_relax_p = 0.25; // see https://github.com/AMReX-Combustion/PeleC/blob/48f1b3d3720771351a4475051c63da5ae6eb01e9/Docs/sphinx/BoundaryConditions.rst
amrex::Real CNS::nscbc_relax_u = 0.2;
amrex::Real CNS::nscbc_relax_T = 0.2;
amrex::Real CNS::ambient_p = 1013250.0; // for outflow

bool CNS::use_typical_vals_chem = false; // tell chem_integrator typical value of Temp
int CNS::reset_typical_vals_int = 10; // interval to reset the typical value
int CNS::rk_reaction_iter = 0;        // iterate to tightly couple chemistry
Real CNS::min_react_temp =
  300.0; // turn off reactor below this temperature (only work with CVODE)
Real CNS::clip_temp = std::numeric_limits<Real>::epsilon();
bool CNS::update_heat_release = true; // write to heat release in I_R? switch it off
                                      // for some performance gain (approx 3-5%)

bool CNS::recon_char_var = true; // Reconstruct the thermodynamic variable (sos or
                                 // gamma) to cell interfaces
int CNS::char_sys = 0;            // Thermodynamic system used in characteristic
                                  // reconstruction. 0:speed of sound, 1:gamma
int CNS::recon_scheme = 5;        // 1: basic Godunov; 2: MUSCL; 3: WENO-Z3;
                                  // 4: WENO-JS5; 5: WENO-Z5; 6: TENO5
Real CNS::plm_theta = 2.0; // [1,2] 1: minmod; 2: van Leer's MC (higher sharper)
bool CNS::use_hybrid_scheme = false; // use central scheme in smooth region
Real CNS::teno_cutoff = 1e-4; // sharp cutoff for smoothness indicator in TENO5, set
                              // to < 0 triggers adaptively calculating cutoff

bool CNS::do_restart_fields = false;
bool CNS::do_psgs = false;
bool CNS::do_pd_model = false;
bool CNS::do_vpdf = false;
bool CNS::do_spdf = false;
bool CNS::do_species_langevin = false;

bool CNS::do_les = false;        // Eddy viscosity type subgrid-scale model
bool CNS::do_pasr = false;       // Partially stirred reactor
std::string CNS::les_model_name; // Avaliable models are "Smagorinsky" or "WALE"
std::unique_ptr<LESModel> CNS::les_model;
Real CNS::Cs = 0.1;     // Smagorinsky constant
Real CNS::C_I = 0.0066; // Yoshizawa constant
Real CNS::Pr_T = 0.7;   // Turbulent Prandtl number = mu_T/(cp*kappa_T)
Real CNS::Sc_T = 0.7;   // Turbulent Schmidt number = mu_T/(rho*D_T)
Real CNS::Cm = 0.1;     // PaSR mixing timescale coefficient

#if CNS_USE_EB
int CNS::eb_recon_mode = 0; //0: fill covered stencil with first cell values; 1: switch to 2nd order scheme
amrex::Real CNS::eb_weight = 0.5;
bool CNS::eb_no_slip = true;
bool CNS::eb_isothermal = false;
Real CNS::eb_wall_temp; // EB wall temperature
std::string CNS::redistribution_type =
  "StateRedist"; // StateRedist / FluxRedist / NoRedist
#endif

bool CNS::signalStopJob = false;

// #endif