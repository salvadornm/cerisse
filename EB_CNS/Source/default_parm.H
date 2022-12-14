// #ifndef CNS_DEFAULT_PARM_H //No need this as this file will be included only once, in CNS.cpp
// #define CNS_DEFAULT_PARM_H

BCRec    CNS::phys_bc;

int      CNS::verbose = 2;
IntVect  CNS::hydro_tile_size {AMREX_D_DECL(1024,16,16)};
Real     CNS::cfl       = 0.3;
Real     CNS::dt_cutoff = 5.e-20;
int      CNS::do_reflux = 1;
int      CNS::refine_cutcells        = 1;
int      CNS::refine_max_dengrad_lev = -1;
Real     CNS::refine_dengrad         = 1.e10;
Vector<RealBox> CNS::refine_boxes;
RealBox* CNS::dp_refine_boxes;

bool     CNS::do_visc        = true;  // diffusion is on by default
// bool     CNS::use_const_visc = false; // diffusion does not use constant viscosity by default
bool     CNS::do_ext_src     = true;
bool     CNS::do_react       = false; // reaction is off by default
std::string CNS::chem_integrator = "Reactor_Null";
pele::physics::transport::TransportParams<
  pele::physics::PhysicsType::transport_type> CNS::trans_parms;
bool     CNS::use_typical_vals_chem = true; // tell chem_integrator typical value of Temp
int      CNS::reset_typical_vals_int = 10;  // interval to reset the typical value
int      CNS::rk_reaction_iter = 0;

bool     CNS::do_restart_fields = false;

int      CNS::recon_scheme   = 4;     // 1: basic Godunov; 2: MUSCL; 3: WENO-JS5; 4: WENO-Z5
Real     CNS::plm_theta      = 2.0;   // [1,2] 1: minmod; 2: van Leer's MC (higher sharper)

bool     CNS::do_vpdf = true;
bool     CNS::do_spdf = false;

int      CNS::eb_weights_type = 0;   // [0,1,2,3] 0: weights are all 1; 1: eint; 2: cell mass; 3: volfrac
int      CNS::do_reredistribution = 1;

bool     CNS::signalStopJob = false;

// #endif