#include "prob.H"

using namespace amrex;

extern "C" {
void 
amrex_probinit (const int* /*init*/,
                const int* /*name*/,
                const int* /*namelen*/,
                const amrex_real* /*problo*/,
                const amrex_real* /*probhi*/)
{
    std::string fuel;
    amrex::ParmParse pp("prob");
    pp.get("fuel", fuel);

    // amrex::Real T_l = 300;
    // amrex::Real p_l = 1013250; //1 atm
    // amrex::Real air_N2_O2_ratio = 3.76;
    // amrex::Real stoich_O2;

    //stoich_O2 = 0.25*gas.n_atoms('H2','H')
    //stoich_O2 = gas.n_atoms('CH4','C') + 0.25*gas.n_atoms('CH4','H')

    //u = 233.6766568266508 for H2
    //u = 38.103032815126364 for CH4

    // if (fuel == "CH4") {
    //     CNS::h_prob_parm->rho_l = 1.1220692645119495143e-3;
    //     CNS::h_prob_parm->rho_r = CNS::h_prob_parm->rho_l; //0.15893692830758096276e-3;
    //     CNS::h_prob_parm->u_l = 0.0; //38.103032815126364;
    //     CNS::h_prob_parm->u_r = -(284.4938475139773 - 38.08304243848445); //269.00131049297897;
    //     // CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l*(-344913.69742286915e7);
    //     // CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r*(-890656.7590928635e7);
    //     CNS::h_prob_parm->massfrac_l[CH4_ID] = 0.05518666598235154;
    //     CNS::h_prob_parm->massfrac_l[O2_ID] = 0.2201412376866278;
    //     CNS::h_prob_parm->massfrac_l[N2_ID] = 0.7246720963310207;
    //     CNS::h_prob_parm->massfrac_r[H2O_ID] = /*0.12394;*/0.1205124213818434;
    //     CNS::h_prob_parm->massfrac_r[CO2_ID] = /*0.15139;*/0.1356685616068702;
    //     CNS::h_prob_parm->massfrac_r[CO_ID] = 9.764257515157837e-03;
    //     CNS::h_prob_parm->massfrac_r[N2_ID] = /*0.72467;*/0.7246541035377703;
    //     CNS::h_prob_parm->massfrac_r[H2_ID] = 2.842357047863630e-04;
    //     auto eos = pele::physics::PhysicsType::eos();
    //     amrex::Real massfrac[NUM_SPECIES];
    //     for (int n = 0; n < NUM_SPECIES; ++n) massfrac[n] = CNS::h_prob_parm->massfrac_l[n];
    //     eos.TY2E(300.0, massfrac, CNS::h_prob_parm->rhoe_l);
    //     // for (int n = 0; n < NUM_SPECIES; ++n) massfrac[n] = CNS::h_prob_parm->massfrac_r[n];
    //     // eos.TY2E(2077.3674796042046, massfrac, CNS::h_prob_parm->rhoe_r);
    //     CNS::h_prob_parm->rhoe_l *= CNS::h_prob_parm->rho_l;
    //     CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rhoe_l;
    //     // CNS::h_prob_parm->rhoe_r *= CNS::h_prob_parm->rho_r;
    //     amrex::Print() << CNS::h_prob_parm->rhoe_l << CNS::h_prob_parm->rhoe_r << std::endl;
    // // } else if (fuel == "H2") {
    //     CNS::h_prob_parm->rho_l = 0.8494721394837322e-3;
    //     CNS::h_prob_parm->rho_r = 0.158936928307581e-3;
    //     CNS::h_prob_parm->u_l = 233.6766568266508;
    //     CNS::h_prob_parm->u_r = 1477.9485731992894;
    //     CNS::h_prob_parm->rhoe_l = CNS::h_prob_parm->rho_l*(-116671.84351781593e7);
    //     CNS::h_prob_parm->rhoe_r = CNS::h_prob_parm->rho_r*(-748306.9422316741e7);
    //     CNS::h_prob_parm->massfrac_l[H2_ID] = 0.028522;
    //     CNS::h_prob_parm->massfrac_l[O2_ID] = 0.22635;
    //     CNS::h_prob_parm->massfrac_l[N2_ID] = 0.74512;
    //     CNS::h_prob_parm->massfrac_r[H2O_ID] = 0.25488;
    //     CNS::h_prob_parm->massfrac_r[N2_ID] = 0.74512;
    // } else {
    //     amrex::Abort("Unknown fuel type. Choose between CH4 and H2");
    // }

    amrex::Real p = 1013250.0;
    amrex::Real T = 300.0;

    // CNS::h_prob_parm->u_l = 38.08304243848445; 
    // CNS::h_prob_parm->u_r = 284.4938475139773;
    CNS::h_prob_parm->u_l = 40.0;
    CNS::h_prob_parm->u_r = 40.0; //38.08304243848445;

    CNS::h_prob_parm->massfrac_l[CH4_ID] = 0.05518666598235154;
    CNS::h_prob_parm->massfrac_l[O2_ID] = 0.2201412376866278;
    CNS::h_prob_parm->massfrac_l[N2_ID] = 0.7246720963310207;

    // CNS::h_prob_parm->massfrac_r[H2O_ID] = 0.1205124213818434;
    // CNS::h_prob_parm->massfrac_r[CO2_ID] = 0.1356685616068702;
    // CNS::h_prob_parm->massfrac_r[CO_ID] = 9.764257515157837e-03;
    // CNS::h_prob_parm->massfrac_r[N2_ID] = 0.7246541035377703;
    // CNS::h_prob_parm->massfrac_r[H2_ID] = 2.842357047863630e-04;
    // CNS::h_prob_parm->massfrac_r[O2_ID] = 6.769637978540457e-03;
    CNS::h_prob_parm->massfrac_r[CH4_ID] = 0.05518666598235154;
    CNS::h_prob_parm->massfrac_r[O2_ID] = 0.2201412376866278;
    CNS::h_prob_parm->massfrac_r[N2_ID] = 0.7246720963310207;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2RE(p, CNS::h_prob_parm->massfrac_l.begin(), T, CNS::h_prob_parm->rho_l, CNS::h_prob_parm->rhoe_l);
    CNS::h_prob_parm->rhoe_l *= CNS::h_prob_parm->rho_l;
    eos.PYT2RE(p, CNS::h_prob_parm->massfrac_r.begin(), 1800.0, CNS::h_prob_parm->rho_r, CNS::h_prob_parm->rhoe_r);
    CNS::h_prob_parm->rhoe_r *= CNS::h_prob_parm->rho_r;

    amrex::Print() << CNS::h_prob_parm->rhoe_l << CNS::h_prob_parm->rhoe_r << std::endl;
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm+1, CNS::d_prob_parm);
// #ifdef AMREX_USE_GPU
//     // Cannot use Gpu::copy because ProbParm is not trivailly copyable.
//     Gpu::htod_memcpy_async(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
// #else
//     std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
// #endif
}
}

/**
 * \brief Fill external source term.
 * 
 * \warning This function will overwrite dsdt fab. Always use += instead of =.
 *  
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param time      time.
 * @param geomdata  domain geometry data.
 * @param state     state data.
 * @param ext_src   external source term.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
void
CNS::fill_ext_src (int i, int j, int k, 
                   amrex::Real /*time*/,
                   amrex::GeometryData const& /*geomdata*/, 
                   amrex::Array4<const amrex::Real> const& state, 
                   amrex::Array4<amrex::Real> const& ext_src, 
                   Parm const& /*parm*/,
                   ProbParm const& /*prob_parm*/)
{
    // Example
    // const Real g = 9.81 * 100;
    // ext_src(i, j, k, UMZ)   += g * state(i, j, k, URHO);
    // ext_src(i, j, k, UEDEN) += g * state(i, j, k, UMZ);
}