
#include <CNS.H>
#include <CNS_K.H>
#include <tagging.H>
#include <parm.H>
#include <prob.H>

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>

#include <climits>

using namespace amrex;

// Default values
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
bool     CNS::use_const_visc = false; // diffusion does not use constant viscosity by default

bool     CNS::do_react       = false; // reaction is off by default
std::string CNS::chem_integrator = "Reactor_Null";
pele::physics::transport::TransportParams<
  pele::physics::PhysicsType::transport_type> CNS::trans_parms;
amrex::Vector<std::string> CNS::spec_names;
bool     CNS::use_typical_vals_chem = true; // tell chem_integrator typical value of Temp
int      CNS::reset_typical_vals_int = 10;  // interval to reset the typical value
int      CNS::rk_reaction_iter = 0;

bool     CNS::do_restart_fields = false;

int      CNS::plm_iorder     = 2;     // [1,2] 1: slopes are zero'd,  2: second order slopes
Real     CNS::plm_theta      = 2.0;   // [1,2] 1: minmod; 2: van Leer's MC
// Real     CNS::gravity        = 0.0;

int      CNS::eb_weights_type = 0;   // [0,1,2,3] 0: weights are all 1
int      CNS::do_reredistribution = 1;

bool     CNS::signalStopJob = false;

CNS::CNS () {}

CNS::CNS (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time) : AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    if (do_reflux && level > 0) {
        flux_reg.define(bl, papa.boxArray(level-1),
                        dm, papa.DistributionMap(level-1),
                        level_geom, papa.Geom(level-1),
                        papa.refRatio(level-1), level, LEN_STATE);
    }

    buildMetrics();

    // Initialize reaction
    get_new_data(Reactions_Type).setVal(0.0);
    if (do_react) {
        if (chem_integrator == "ReactorNull") {
            Print() << "WARNING: turning on reactions while using ReactorNull. "
                       "Make sure this is intended." << std::endl;
        }

        reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
        reactor->init(1, 1);
    }
}

CNS::~CNS () 
{
    if (do_react) {
        reactor->close();
    }
}

void
CNS::init (AmrLevel& old)
{
    BL_PROFILE("CNS::init(old)");

    auto& oldlev = dynamic_cast<CNS&>(old);

    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev.state[State_Type].curTime();
    Real prev_time = oldlev.state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old, S_new, 0, cur_time, State_Type, 0, LEN_STATE);

    amrex::MultiFab& React_new = get_new_data(Reactions_Type);
    if (do_react) {
        FillPatch(
        old, React_new, 0, cur_time, Reactions_Type, 0, React_new.nComp());
    } else {
        React_new.setVal(0);
    }

    MultiFab& C_new = get_new_data(Cost_Type);
    FillPatch(old, C_new, 0, cur_time, Cost_Type, 0, 1);
}

void
CNS::init ()
{
    // This version inits the data on a new level that did not
    // exist before regridding.
    BL_PROFILE("CNS::init()");

    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
    setTimeLevel(cur_time,dt_old,dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, LEN_STATE);

    MultiFab& C_new = get_new_data(Cost_Type);
    FillCoarsePatch(C_new, 0, cur_time, Cost_Type, 0, 1);
}

void
CNS::initData ()
{
    BL_PROFILE("CNS::initData()");

    const auto geomdata = geom.data();
    MultiFab& S_new = get_new_data(State_Type);
    S_new.setVal(0.0);

    Parm const* lparm = d_parm;
    ProbParm const* lprobparm = d_prob_parm; // <T> const* = pointer to constant <T>; const <T>* == <T> const*

    if (verbose != 0) {
        Print() << "Initializing the data at level " << level << std::endl;
    }

    auto const& sarrs = S_new.arrays();
    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
        cns_initdata(i, j, k, sarrs[box_no], geomdata, *lparm, *lprobparm);
        // Verify that the sum of (rho Y)_i = rho at every cell
        // cns_check_initial_species(i, j, k, sarrs[box_no]);
    });
    amrex::Gpu::synchronize();

    // Compute the initial temperature (will override what was set in initdata)
    computeTemp(S_new, 0);

    // set_body_state(S_new);
    // amrex::Real cur_time = state[State_Type].curTime();
    // const amrex::StateDescriptor* desc = state[State_Type].descriptor();
    // const auto& bcs = desc->getBCs();
    // InitialRedistribution(cur_time, bcs, S_new);

    // get_new_data(Reactions_Type).setVal(0.0);

    get_new_data(Cost_Type).setVal(1.0);
}

void
CNS::computeInitialDt (int                    finest_level,
                       int                    /*sub_cycle*/,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/,
                       Vector<Real>&          dt_level,
                       Real                   stop_time)
{
    // Grids have been constructed, compute dt for all levels.
    if (level > 0) return;

    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        dt_level[i] = getLevel(i).estTimeStep();
        n_factor *= n_cycle[i];
        dt_0 = std::min<amrex::Real>(dt_0, n_factor*dt_level[i]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
        dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::computeNewDt (int                    finest_level,
                   int                    /*sub_cycle*/,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& /*ref_ratio*/,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag)
{
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    if (level > 0) return;

    for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1) {
        // Limit dt's by pre-regrid dt
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = std::min(dt_min[i], dt_level[i]);
        }
    } else {
        // Limit dt's by change_max * old dt (TODO: Do we really need this?)
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++) {
            if ((verbose != 0) && ParallelDescriptor::IOProcessor()) {
                if (dt_min[i] > change_max * dt_level[i]) {
                    Print() << " ... Compute new dt : limiting dt at level "
                            << i << '\n';
                    Print() << "     new dt computed: " << dt_min[i] << '\n';
                    Print() << "     but limiting to: "
                            << change_max * dt_level[i] << " = " << change_max
                            << " * " << dt_level[i] << '\n';
                }
            }

            dt_min[i] = std::min(dt_min[i], change_max*dt_level[i]);
        }
    }

    // Find the minimum over all levels
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor*dt_min[i]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::post_regrid (int /*lbase*/, int /*new_finest*/)
{
    if (do_react && use_typical_vals_chem) {
        set_typical_values_chem();
    }
}

void
CNS::post_timestep (int /*iteration*/)
{
    BL_PROFILE("CNS::post_timestep()");

    if (do_reflux && level < parent->finestLevel()) {
        CNS& fine_level = getLevel(level+1);
        MultiFab& S_crse = get_new_data(State_Type);
        MultiFab& S_fine = fine_level.get_new_data(State_Type);
        fine_level.flux_reg.Reflux(S_crse, *volfrac, S_fine, *fine_level.volfrac);

        if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Reflux() at level " << level
                           << " : time = " << state[State_Type].curTime() << std::endl;
        }
    }

    if (level < parent->finestLevel()) {
        avgDown();
    }

    // Re-compute temperature and reset typical_values_chem
    amrex::MultiFab& S_new = get_new_data(State_Type);
    computeTemp(S_new, 0);
    if ( do_react && use_typical_vals_chem 
      && parent->levelSteps(0) % reset_typical_vals_int == 0) {
        set_typical_values_chem();
    }
}

void
CNS::postCoarseTimeStep (Real /*time*/)
{
    // This only computes sum on level 0
    // if (verbose >= 2) { 
    printTotal(); //must do because this checks for nan as well
    // }
}

void
CNS::printTotal () const
{
    // Print total and check for nan
    BL_PROFILE("CNS::printTotal()");

    const MultiFab& S_new = get_new_data(State_Type);
    MultiFab mf(grids, dmap, 1, 0);
    std::array<Real,6> tot;
    for (int comp = 0; comp < 6; ++comp) {
        MultiFab::Copy(mf, S_new, comp, 0, 1, 0);
        MultiFab::Multiply(mf, *volfrac, 0, 0, 1, 0);
        tot[comp] = mf.sum(0,true) * geom.ProbSize();
    }
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealSum(tot.data(), 5, ParallelDescriptor::IOProcessorNumber());
            amrex::Print().SetPrecision(17) << "\n[CNS] Total mass       = " << tot[0] << "\n"
                                            <<   "      Total x-momentum = " << tot[1] << "\n"
                                            <<   "      Total y-momentum = " << tot[2] << "\n"
                                            <<   "      Total z-momentum = " << tot[3] << "\n"
                                            <<   "      Total energy     = " << tot[4] << "\n"
                                            <<   "      Total int energy = " << tot[5] << "\n";
#ifdef BL_LAZY
        });
#endif

    // Nan detector
    if (std::isnan(tot[4])) signalStopJob = true;
    amrex::ParallelDescriptor::ReduceBoolOr(signalStopJob);
}

void
CNS::post_init (Real /*stop_time*/)
{
    if (level > 0) return;
    for (int k = parent->finestLevel()-1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    if (do_react) {        
        if (use_typical_vals_chem) {
            set_typical_values_chem();
        }

        react_source(parent->cumTime(), parent->dtLevel(level), true);
    }

    if (verbose >= 2) {
        printTotal();
    }
}

void
CNS::post_restart ()
{
    // Initialize reactor
    if (do_react) {
        if (chem_integrator == "ReactorNull") {
            Print() << "WARNING: turning on reactions while using ReactorNull. "
                       "Make sure this is intended." << std::endl;
        }

        reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
        reactor->init(1, 1);

        if (use_typical_vals_chem) {
            set_typical_values_chem();
        }
    }
    
    MultiFab& S_new = get_new_data(State_Type);

    // Populate fields (when restarting from a different number of fields)
    if (do_restart_fields) {
        Print() << "Resetting field state data ..." << std::endl;
        
        // Move aux variables
        MultiFab::Copy(S_new, S_new, NVAR, UFA, NUM_AUX, 0);

        // Copy mean to fields
        for (int nf = 1; nf < NUM_FIELD+1; ++nf) {
            MultiFab::Copy(S_new, S_new, 0, nf*NVAR, NVAR, 0);
        }

        if (do_react) {
            Print() << "Resetting field reaction data ..." << std::endl;

            MultiFab& I_R = get_new_data(Reactions_Type);

            // Copy mean to fields
            for (int nf = 1; nf < NUM_FIELD+1; ++nf) {
                MultiFab::Copy(I_R, I_R, 0, nf*NREACT, NREACT, 0);
            }
        }
    }

    Parm const* lparm = d_parm;
    ProbParm const* lprobparm = d_prob_parm;
    const auto geomdata = geom.data();
    auto const& sarrs = S_new.arrays();
    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
        // Modify restarted data and/or add turbulence
        cns_prob_post_restart(i, j, k, sarrs[box_no], geomdata, *lparm, *lprobparm);
    });
}

int
CNS::okToContinue () 
{
    if (level > 0) {
        return 1;
    }

    int test = 1;
    if (signalStopJob) {
        test = 0;

        amrex::Print() << "Signalling a stop of the run due to signalStopJob = true."
                       << std::endl;
    } else if (parent->dtLevel(0) < dt_cutoff) {
        test = 0;

        amrex::Print() << "Signalling a stop of the run because dt < dt_cutoff."
                       << std::endl;
    }

    return test;
}

void
CNS::errorEst (TagBoxArray& tags, int /*clearval*/, int /*tagval*/, 
               Real /*time*/, int /*n_error_buf = 0*/, int /*ngrow*/)
{
    BL_PROFILE("CNS::errorEst()");

    if (refine_cutcells) {
        const MultiFab& S_new = get_new_data(State_Type);
        TagCutCells(tags, S_new);
    }

    if (!refine_boxes.empty())
    {
        const int n_refine_boxes = refine_boxes.size();
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        auto boxes = dp_refine_boxes;

        auto const& tagma = tags.arrays();
        ParallelFor(tags,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
            RealVect pos {AMREX_D_DECL((i+0.5)*dx[0]+problo[0],
                                       (j+0.5)*dx[1]+problo[1],
                                       (k+0.5)*dx[2]+problo[2])};
            for (int irb = 0; irb < n_refine_boxes; ++irb) {
                if (boxes[irb].contains(pos)) {
                    tagma[box_no](i,j,k) = TagBox::SET;
                }
            }
        });
        Gpu::streamSynchronize();
    }

    if (level < refine_max_dengrad_lev) {
        const MultiFab& S_new = get_new_data(State_Type);
        const Real cur_time = state[State_Type].curTime();
        MultiFab rho(S_new.boxArray(), S_new.DistributionMap(), 1, 1);
        FillPatch(*this, rho, rho.nGrow(), cur_time, State_Type, URHO, 1, 0);

        const char   tagval = TagBox::SET;
        // const char clearval = TagBox::CLEAR;
        const Real dengrad_threshold = refine_dengrad;

        auto const& tagma = tags.arrays();
        auto const& rhoma = rho.const_arrays();
        ParallelFor(rho,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
            cns_tag_graderror(i, j, k, tagma[box_no], rhoma[box_no], dengrad_threshold, tagval);
        });
        Gpu::streamSynchronize();
    }
}

void
CNS::read_params ()
{
    ParmParse pp("cns");

    pp.query("v", verbose);

    Vector<int> tilesize(AMREX_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }

    pp.query("cfl", cfl);
    pp.query("dt_cutoff", dt_cutoff);

    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    // pp.query("do_reflux", do_reflux); // How can you not do reflux?!

    pp.query("do_visc", do_visc);

    // if (do_visc) {
    //     pp.query("use_const_visc",use_const_visc);
    //     if (use_const_visc) {
    //         pp.get("const_visc_mu",h_parm->const_visc_mu);
    //         pp.get("const_visc_ki",h_parm->const_visc_ki);
    //         pp.get("const_lambda" ,h_parm->const_lambda);
    //     }
    // } else {
    //    use_const_visc = true;
    //    h_parm->const_visc_mu = 0.0;
    //    h_parm->const_visc_ki = 0.0;
    //    h_parm->const_lambda  = 0.0;
    // }

    pp.query("refine_cutcells", refine_cutcells);

    pp.query("refine_max_dengrad_lev", refine_max_dengrad_lev);
    pp.query("refine_dengrad", refine_dengrad);

    int irefbox = 0;
    Vector<Real> refboxlo, refboxhi;
    while (pp.queryarr(("refine_box_lo_"+std::to_string(irefbox)).c_str(), refboxlo))
    {
        pp.getarr(("refine_box_hi_"+std::to_string(irefbox)).c_str(), refboxhi);
        refine_boxes.emplace_back(refboxlo.data(), refboxhi.data());
        ++irefbox;
    }
    if (!refine_boxes.empty()) {
#ifdef AMREX_USE_GPU
        dp_refine_boxes = (RealBox*)The_Arena()->alloc(sizeof(RealBox)*refine_boxes.size());
        Gpu::htod_memcpy_async(dp_refine_boxes, refine_boxes.data(), sizeof(RealBox)*refine_boxes.size());
#else
        dp_refine_boxes = refine_boxes.data();
#endif
    }

    pp.query("do_react", do_react);
    pp.query("chem_integrator", chem_integrator);
    pp.query("rk_reaction_iter", rk_reaction_iter);
    pp.query("use_typical_vals_chem", use_typical_vals_chem);
    pp.query("reset_typical_vals_int", reset_typical_vals_int);

    pp.query("do_restart_fields", do_restart_fields);

    // pp.query("gravity", gravity);

    // pp.query("eos_gamma", h_parm->eos_gamma);
    // pp.query("eos_mu"   , h_parm->eos_mu);
    // pp.query("Pr"       , h_parm->Pr);
    // pp.query("C_S"      , h_parm->C_S);
    // pp.query("T_S"      , h_parm->T_S);

    // h_parm->Initialize();
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_parm, h_parm+1, d_parm);

    // eb_weights_type:
    //   0 -- weights = 1
    //   1 -- use_total_energy_as_eb_weights-
    //   2 -- use_mass_as_eb_weights
    //   3 -- use_volfrac_as_eb_weights
    pp.query("eb_weights_type", eb_weights_type);
    if (eb_weights_type < 0 || eb_weights_type > 3) {
        amrex::Abort("CNS: eb_weights_type must be 0, 1, 2 or 3");
    }

    pp.query("do_reredistribution", do_reredistribution);
    if (do_reredistribution != 0 && do_reredistribution != 1) {
        amrex::Abort("CNS: do_reredistibution must be 0 or 1");
    }

    amrex::Gpu::streamSynchronize();
}

void
CNS::avgDown ()
{
    BL_PROFILE("CNS::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level+1);
          MultiFab& S_crse =          get_new_data(State_Type);
    const MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
    geom.GetVolume(volume);

    amrex::EB_average_down(S_fine, S_crse, volume, fine_lev.volFrac(),
                           0, S_fine.nComp(), fine_ratio);
    const int nghost = 0;
    computeTemp(S_crse, nghost);

    if (do_react) {
              MultiFab& R_crse =          get_new_data(Reactions_Type);
        const MultiFab& R_fine = fine_lev.get_new_data(Reactions_Type);
        amrex::EB_average_down(R_fine, R_crse, volume, fine_lev.volFrac(),
                               0, S_fine.nComp(), fine_ratio);
    }
}

void
CNS::buildMetrics ()
{
    BL_PROFILE("CNS::buildMetrics()");

    // make sure dx == dy == dz
//     const Real* dx = geom.CellSize();
// #if (AMREX_SPACEDIM == 2)
//     if (std::abs(dx[0]-dx[1]) > 1.e-12*dx[0])
//         amrex::Abort("CNS: must have dx == dy\n");
// #else
//     if (std::abs(dx[0]-dx[1]) > 1.e-12*dx[0] || std::abs(dx[0]-dx[2]) > 1.e-12*dx[0])
//         amrex::Abort("CNS: must have dx == dy == dz\n");
// #endif

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());

    volfrac = &(ebfactory.getVolFrac());
    bndrycent = &(ebfactory.getBndryCent());
    areafrac = ebfactory.getAreaFrac();
    facecent = ebfactory.getFaceCent();

    Parm const* l_parm = d_parm;

    level_mask.clear();
    level_mask.define(grids,dmap,1,3);
    level_mask.BuildMask(geom.Domain(), geom.periodicity(),
                         l_parm->level_mask_covered,
                         l_parm->level_mask_notcovered,
                         l_parm->level_mask_physbnd,
                         l_parm->level_mask_interior);
}

Real
CNS::estTimeStep ()
{
    BL_PROFILE("CNS::estTimeStep()");

    const auto dx = geom.CellSizeArray();
    const MultiFab& S = get_new_data(State_Type);
    Parm const* lparm = d_parm;

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    Real estdt = std::numeric_limits<Real>::max();

    // Reduce min operation
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& flag = flags[mfi];
        auto const& s_arr = S.array(mfi);
        if (flag.getType(bx) != FabType::covered)
        {
          reduce_op.eval(bx, reduce_data, [=]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return cns_estdt(i,j,k,s_arr,dx,*lparm);
            });
        }
    } // mfi

    ReduceTuple host_tuple = reduce_data.value();
    estdt = amrex::min(estdt,amrex::get<0>(host_tuple));

    estdt *= cfl;
    ParallelDescriptor::ReduceRealMin(estdt);
    return estdt;
}

// Real
// CNS::initialTimeStep ()
// {
//     return estTimeStep();
// }

void
CNS::computeTemp (MultiFab& State, int ng)
{
    BL_PROFILE("CNS::computeTemp()");

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(State.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    // This will reset Eint and compute Temperature
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(State, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        auto s_arr = State.array(mfi);

        if (flags[mfi].getType(bx) != FabType::covered) {
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cns_compute_temperature(i, j, k, s_arr);
            });
        }
    }
    // auto const& sarrs = State.arrays();
    // auto const& flagarrs = flags.const_arrays();
    // const amrex::IntVect ngs(ng);
    // amrex::ParallelFor(
    //     State, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
    //         if (!flagarrs[nbx](i, j, k).isCovered()) {
    //             cns_compute_temperature(i, j, k, sarrs[nbx]);
    //         }
    //     });
    amrex::Gpu::synchronize();
}
