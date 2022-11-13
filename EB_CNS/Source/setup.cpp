#include "CNS.H"
#include "derive.H"
#include "prob.H"

using namespace amrex;

int CNS::num_state_data_types = 0;
Parm* CNS::h_parm = nullptr;
Parm* CNS::d_parm = nullptr;
ProbParm* CNS::h_prob_parm = nullptr;
ProbParm* CNS::d_prob_parm = nullptr;

using BndryFunc = StateDescriptor::BndryFunc;

// Components are:
//  Interior,        Inflow,          Outflow,          Symmetry,             SlipWall,             NoSlipWall,           UserBC
static int scalar_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::ext_dir
};

static int norm_vel_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_odd,  BCType::reflect_odd,  BCType::reflect_odd,  BCType::ext_dir
};

static int tang_vel_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd,  BCType::ext_dir
};

static int react_src_bc[] = {
    BCType::int_dir, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even
};

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        bc.setLo(dir, scalar_bc[lo_bc[dir]]);
        bc.setHi(dir, scalar_bc[hi_bc[dir]]);
    }
}

static
void
set_x_vel_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    AMREX_D_TERM(
        bc.setLo(0, norm_vel_bc[lo_bc[0]]);  bc.setHi(0, norm_vel_bc[hi_bc[0]]); ,
        bc.setLo(1, tang_vel_bc[lo_bc[1]]);  bc.setHi(1, tang_vel_bc[hi_bc[1]]); ,
        bc.setLo(2, tang_vel_bc[lo_bc[2]]);  bc.setHi(2, tang_vel_bc[hi_bc[2]]); 
    )
}

static
void
set_y_vel_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    AMREX_D_TERM(
        bc.setLo(0, tang_vel_bc[lo_bc[0]]);  bc.setHi(0, tang_vel_bc[hi_bc[0]]); ,
        bc.setLo(1, norm_vel_bc[lo_bc[1]]);  bc.setHi(1, norm_vel_bc[hi_bc[1]]); ,
        bc.setLo(2, tang_vel_bc[lo_bc[2]]);  bc.setHi(2, tang_vel_bc[hi_bc[2]]); 
    )
}

static
void
set_z_vel_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    AMREX_D_TERM(
        bc.setLo(0,tang_vel_bc[lo_bc[0]]);  bc.setHi(0,tang_vel_bc[hi_bc[0]]); ,
        bc.setLo(1,tang_vel_bc[lo_bc[1]]);  bc.setHi(1,tang_vel_bc[hi_bc[1]]); ,
        bc.setLo(2,norm_vel_bc[lo_bc[2]]);  bc.setHi(2,norm_vel_bc[hi_bc[2]]);
    )
}

static 
void
set_react_src_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        bc.setLo(dir, react_src_bc[lo_bc[dir]]);
        bc.setHi(dir, react_src_bc[hi_bc[dir]]);
    }
}

void
CNS::variableSetUp ()
{
    h_parm = new Parm{};          //This lives on host
    h_prob_parm = new ProbParm{}; //This lives on host
    d_parm = (Parm*)The_Arena()->alloc(sizeof(Parm)); //static_cast, This lives on device
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm)); //This lives on device
    trans_parms.allocate(); //PelePhysics trans_parms

    read_params();

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
    Interpolater* interp = &eb_mf_cell_cons_interp;
#else
    Interpolater* interp = &cell_cons_interp; //or pc_interp, lincc_interp
#endif

    // Setup State_Type
    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, NUM_GROW, LEN_STATE,
                           interp, state_data_extrap, store_in_checkpoint);

    amrex::Vector<amrex::BCRec> bcs(LEN_STATE);
    amrex::Vector<std::string>  name(LEN_STATE);
    amrex::BCRec bc;

    pele::physics::eos::speciesNames<
        pele::physics::PhysicsType::eos_type>(spec_names);
    
    if (amrex::ParallelDescriptor::IOProcessor()) {
        // Print species names
        Print() << NUM_SPECIES << " Species: ";
        for (int i = 0; i < NUM_SPECIES; i++) {
            Print() << spec_names[i] << ' ' << ' ';
        }
        Print() << std::endl;

        // Print number of fields and aux
        Print() << NUM_FIELD << " Fields" << std::endl
                << NUM_AUX   << " Auxiliary Variables" << std::endl 
                << LEN_STATE << LEN_REACT << LEN_PRIM << LEN_COEF
                << std::endl;
    }

    int cnt = 0; // variable counter
    // Mean field
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density"; cnt++;
    set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";    cnt++;
    set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";    cnt++;
    set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";    cnt++;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";   cnt++;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";   cnt++;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";    cnt++;
    for (int i = 0; i < NUM_SPECIES; ++i) {
        set_scalar_bc(bc, phys_bc); bcs[cnt] = bc; name[cnt] = "rho_" + spec_names[i];
        cnt++;
    }

    // For each field
    for (int nf = 0; nf < NUM_FIELD; ++nf) {
        set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density_Field" + std::to_string(nf); cnt++;
        set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom_Field"    + std::to_string(nf); cnt++;
        set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom_Field"    + std::to_string(nf); cnt++;
        set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom_Field"    + std::to_string(nf); cnt++;
        set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E_Field"   + std::to_string(nf); cnt++;
        set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e_Field"   + std::to_string(nf); cnt++;
        set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp_Field"    + std::to_string(nf); cnt++;
        for (int i = 0; i < NUM_SPECIES; ++i) {
            set_scalar_bc(bc, phys_bc); bcs[cnt] = bc; name[cnt] = "rho_" + spec_names[i] + "_Field" + std::to_string(nf);
            cnt++;
        }
    }

    // Get AUX names
    for (int i = 0; i < NUM_AUX; ++i) {
        set_scalar_bc(bc, phys_bc); bcs[cnt] = bc; name[cnt] = "aux_" + std::to_string(i);
        cnt++;
    }

    StateDescriptor::BndryFunc bndryfunc(cns_bcfill); // <--- need for all fields
    bndryfunc.setRunOnGPU(true);

    desc_lst.setComponent(State_Type, 0, name, bcs, bndryfunc);

    // Setup React_Type
    store_in_checkpoint = do_react;
    desc_lst.addDescriptor(Reactions_Type, amrex::IndexType::TheCellType(),
                           amrex::StateDescriptor::Point, 0, LEN_REACT, 
                           interp, state_data_extrap, store_in_checkpoint);

    amrex::Vector<amrex::BCRec> react_bcs(LEN_REACT);
    amrex::Vector<std::string>  react_name(LEN_REACT);

    cnt = 0;
    // Mean field
    for (int i = 0; i < NUM_SPECIES; ++i) {
        set_react_src_bc(bc, phys_bc); 
        react_bcs[cnt] = bc; 
        react_name[cnt] = "rho_omega_" + spec_names[i];
        cnt++;
    }
    set_react_src_bc(bc, phys_bc); react_bcs[cnt] = bc; react_name[cnt] = "rhoe_dot";    cnt++;
    set_react_src_bc(bc, phys_bc); react_bcs[cnt] = bc; react_name[cnt] = "heatRelease"; cnt++;

    // For each field
    for (int nf = 0; nf < NUM_FIELD; ++nf) {
        for (int i = 0; i < NUM_SPECIES; ++i) {            
            set_react_src_bc(bc, phys_bc); 
            react_bcs[cnt] = bc; 
            react_name[cnt] = "rho_omega_" + spec_names[i] + "_Field" + std::to_string(nf);
            cnt++; 
        }
        set_react_src_bc(bc, phys_bc); react_bcs[cnt] = bc; react_name[cnt] = "rhoe_dot_Field"    + std::to_string(nf); cnt++; 
        set_react_src_bc(bc, phys_bc); react_bcs[cnt] = bc; react_name[cnt] = "heatRelease_Field" + std::to_string(nf); cnt++;
    }

    amrex::StateDescriptor::BndryFunc bndryfunc2(cns_react_bcfill); // <--- need for all fields
    bndryfunc2.setRunOnGPU(true);

    desc_lst.setComponent(Reactions_Type, 0, react_name, react_bcs, bndryfunc2);

    // Setup Cost_Type
    desc_lst.addDescriptor(Cost_Type, IndexType::TheCellType(), StateDescriptor::Point,
                           0, 1, &pc_interp);
    desc_lst.setComponent(Cost_Type, 0, "Cost", bc, bndryfunc);

    num_state_data_types = desc_lst.size();
    Print() << desc_lst.size() << " Data Types:" << std::endl;
    for (int typ = 0; typ < desc_lst.size(); typ++) {
        Print() << typ << " - ";
        const StateDescriptor& desc = desc_lst[typ];
        for (int n = 0; n < desc.nComp(); n++) {
            Print() << desc.name(n) << " ";
        }
        Print() << "(" << desc.nComp() << ")" << std::endl;
    }

    StateDescriptor::setBndryFuncThreadSafety(true);

    // DEFINE DERIVED QUANTITIES (Derive from MEAN field)
    // Pressure
    derive_lst.add("pressure", IndexType::TheCellType(), 1,
                   cns_derpres, amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("pressure", desc_lst, State_Type, URHO, NVAR);

    // Velocities
    derive_lst.add("x_velocity", IndexType::TheCellType(), 1,
                   cns_dervel, amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, UMX, 1);

#if (AMREX_SPACEDIM >= 2)
    derive_lst.add("y_velocity", IndexType::TheCellType(), 1,
                   cns_dervel, amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, UMY, 1);
#endif

#if (AMREX_SPACEDIM == 3)
    derive_lst.add("z_velocity", IndexType::TheCellType(), 1,
                   cns_dervel, amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, URHO, 1);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, UMZ, 1);
#endif

    // Mach number
    derive_lst.add("MachNumber", amrex::IndexType::TheCellType(), 1, 
                   cns_dermachnumber, amrex::DeriveRec::TheSameBox);
    derive_lst.addComponent("MachNumber", desc_lst, State_Type, URHO, NVAR);

    // Vorticity
    derive_lst.add("magvort", amrex::IndexType::TheCellType(), 1, 
                   cns_dermagvort, amrex::DeriveRec::GrowBoxByOne);
    derive_lst.addComponent("magvort", desc_lst, State_Type, URHO, AMREX_SPACEDIM + 1);

    // Numerical schlieren
    derive_lst.add("divu", amrex::IndexType::TheCellType(), 1, 
                   cns_derdivu, amrex::DeriveRec::GrowBoxByOne);
    derive_lst.addComponent("divu", desc_lst, State_Type, URHO, AMREX_SPACEDIM + 1);

    derive_lst.add("divrho", amrex::IndexType::TheCellType(), 1, 
                   cns_derdivrho, amrex::DeriveRec::GrowBoxByOne);
    derive_lst.addComponent("divrho", desc_lst, State_Type, URHO, 1);
}

void
CNS::variableCleanUp ()
{
    delete h_parm;
    delete h_prob_parm;
    The_Arena()->free(d_parm);
    The_Arena()->free(d_prob_parm);
    desc_lst.clear();
    derive_lst.clear();

#ifdef AMREX_USE_GPU
    The_Arena()->free(dp_refine_boxes);
#endif
}
