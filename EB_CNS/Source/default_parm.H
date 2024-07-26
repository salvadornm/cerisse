// This file specifics default values for static parameters in the CNS class.
// It is included only once, in CNS.cpp.

// Control parameters

int  CNS::verbose       = 2;      // Verbosity (0: silent, 1: basic, 2: dianostic, 3: warning/error)
Real CNS::cfl           = 0.3;    // Courant number (must be < 1)
Real CNS::fixed_dt      = -1;     // Fixed time step size (if < 0, dt is computed from the CFL condition)
Real CNS::dt_cutoff     = 5.e-20; // Minimum allowed time step size
Real CNS::dt_max_change = 1.2;    // dt^n+1 <= dt_max_change * dt^n
int  CNS::rk_order      = 2;      // Runge-Kutta time integration order (note: order 1 is unstable for WENO, and order > 2 does not support reaction and SF)
int  CNS::rk_reaction_iter = 0;   // How many times to iterate to tightly couple chemistry
int  CNS::amr_interp_order = 2;   // Interpolation order for AMR coarse-fine boundary
bool CNS::do_reflux       = true; // Correct flux at coarse-fine interfaces?
bool CNS::do_load_balance = true; // Load balance?
bool CNS::signalStopJob  = false; // Signal to stop the job?
int  CNS::check_message_int = 100; // Check runtime message interval

// Physics parameters

bool CNS::do_hydro   = true;  // Compute hydrodynamic/hyprobolic terms?
bool CNS::do_visc    = true;  // Compute diffusiive terms?
bool CNS::do_ext_src = true;  // Call external source term in prob.cpp?
bool CNS::do_react   = false; // Compute reaction source term?
Real CNS::clip_temp  = std::numeric_limits<Real>::epsilon(); // Clip minimum temperature [K]
Real CNS::min_react_temp = 300.0;     // Turn off reactor below this temperature (can save some work with CVODE)
bool CNS::update_heat_release = true; // Write to heat release in I_R? Switch off for some performance gain (approx 3-5%)
BCRec CNS::phys_bc; // Holds BC info

// EB controls

#if CNS_USE_EB
int  CNS::eb_recon_mode = 0;     // 0: fill covered stencil with first cell values; 1: switch to 2nd order scheme close to EB
Real CNS::eb_weight     = 0.5;   // Weight for FluxRedist
bool CNS::eb_no_slip    = true;  // No-slip EB walls?
bool CNS::eb_isothermal = false; // Isothermal EB walls?
Real CNS::eb_wall_temp;          // EB wall temperature (if eb_isothermal = true)
std::string CNS::redistribution_type = "StateRedist"; // StateRedist / FluxRedist / NoRedist
#endif

// Numerics controls

int  CNS::recon_scheme      = 5;     // 1: Godunov, 2: MUSCL, 3: WENO-Z3, 4: WENO-JS5, 5: WENO-Z5, 6: TENO5
int  CNS::char_sys          = 0;     // What to use in characteristic transformation (0: c, 1: gamma)
bool CNS::recon_char_var    = true;  // Reconstruct the thermodynamic variable to cell interfaces? A must for high-order
Real CNS::plm_theta         = 2.0;   // [1,2] 1: minmod; 2: van Leer's MC (higher sharper)
bool CNS::use_hybrid_scheme = false; // Use central scheme in smooth region?
Real CNS::teno_cutoff       = 1e-4;  // TENO's cutoff value. < 0 triggers adaptively calculating cutoff
bool CNS::use_typical_vals_chem = false; // tell chem_integrator typical value of Temp
int  CNS::reset_typical_vals_int = 10;   // interval to reset the typical value

RealBox CNS::buffer_box; // Increase viscosity in this box (usually near outflow) to dampen pressure waves

bool CNS::do_nscbc = false;
Real CNS::nscbc_relax_p = 0.25;  // see https://github.com/AMReX-Combustion/PeleC/blob/48f1b3d3720771351a4475051c63da5ae6eb01e9/Docs/sphinx/BoundaryConditions.rst
Real CNS::nscbc_relax_u = 0.2;
Real CNS::nscbc_relax_T = 0.2;
Real CNS::ambient_p = 1013250.0; // for outflow

// LES

bool CNS::do_les  = false;       // Enable eddy viscosity type SGS model?
bool CNS::do_pasr = false;       // Enable partially stirred reactor model?
std::string CNS::les_model_name; // Avaliable models are "Smagorinsky" or "WALE"
std::unique_ptr<LESModel> CNS::les_model;
Real CNS::Cs   = 0.1;    // Smagorinsky constant
Real CNS::C_I  = 0.0066; // Yoshizawa constant (not used)
Real CNS::Pr_T = 0.7;    // Turbulent Prandtl number = mu_T/(cp*kappa_T)
Real CNS::Sc_T = 0.7;    // Turbulent Schmidt number = mu_T/(rho*D_T)
Real CNS::Cm   = 0.1;    // PaSR mixing timescale coefficient

// Stochastic fields

bool CNS::do_restart_fields = false; // Restart stochastic fields from mean?
bool CNS::do_vpdf           = false; // Enable simplified Langevin model?
bool CNS::do_spdf           = false; // Enable IEM model for species and energy?
bool CNS::do_psgs           = false; // Enable pressure gradient correction?
bool CNS::do_pd_model       = false; // Enable pressure-dilatation model?
bool CNS::do_species_langevin = false; // For testing only

// Tagging refinement criteria (refine_<var> and refine_<var>_max_lev)

#if CNS_USE_EB
bool CNS::refine_cutcells = true; // Refine near EB cut cells?
int  CNS::refine_cutcells_max_lev = 10;
#endif
Vector<RealBox> CNS::refine_boxes; // Refine in user-specified boxes
Vector<int>     CNS::refine_boxes_max_lev;
RealBox*        CNS::dp_refine_boxes;
Vector<Real> CNS::refine_dengrad;  // Refine based on density gradient
int          CNS::refine_dengrad_max_lev = -1;
Vector<Real> CNS::refine_velgrad;  // Refine based on velocity gradient
int          CNS::refine_velgrad_max_lev = -1;
Vector<Real> CNS::refine_presgrad; // Refine based on pressure gradient
int          CNS::refine_presgrad_max_lev = -1;
Vector<Real> CNS::refine_magvort;  // Refine based on vorticity magnitude 
int          CNS::refine_magvort_max_lev = -1;
Vector<Real> CNS::refine_tke;     // Refine based on turbulent kinetic energy (only for SF)
int          CNS::refine_tke_max_lev = -1;

Vector<AMRErrorTag> CNS::errtags; // AMReX built in error tagging tools

// PelePhysics related

std::string CNS::chem_integrator = "Reactor_Null"; // see PelePhysics/Reactions
pele::physics::transport::TransportParams<
  pele::physics::PhysicsType::transport_type> CNS::trans_parms; // Transport parameters
amrex::Vector<std::string> CNS::spec_names; // Species names
pele::physics::turbinflow::TurbInflow CNS::turb_inflow; // Turbulent inflow
pele::physics::PMF::PmfData CNS::pmf_data; // Premixed flame data holder