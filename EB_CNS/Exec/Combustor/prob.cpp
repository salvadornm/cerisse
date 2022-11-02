#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/,
                    const int* /*name*/,
                    const int* /*namelen*/,
                    const amrex_real* /*problo*/,
                    const amrex_real* /*probhi*/)
{
    ParmParse pp("prob");
    
    pp.query("inflow_T"   , CNS::h_prob_parm->inflow_T);
    pp.query("inflow_p"   , CNS::h_prob_parm->inflow_p);
    pp.query("inflow_mach", CNS::h_prob_parm->inflow_mach);
    pp.query("interior_T" , CNS::h_prob_parm->interior_T);
    pp.query("interior_P" , CNS::h_prob_parm->interior_p);

    CNS::h_prob_parm->massfrac[N2_ID] = 0.7;
    CNS::h_prob_parm->massfrac[O2_ID] = 0.3;

#ifdef AMREX_USE_GPU
    // Cannot use Gpu::copy because ProbParm is not trivailly copyable.
    Gpu::htod_memcpy_async(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#endif

    Gpu::HostVector<Real> inflow_state(NUM_STATE);

    auto eos = pele::physics::PhysicsType::eos();
    Real rho, eint, cs;
    eos.PYT2RE(CNS::h_prob_parm->inflow_p, CNS::h_prob_parm->massfrac.begin(), CNS::h_prob_parm->inflow_T,   
               rho, eint);
    eos.RTY2Cs(rho, CNS::h_prob_parm->inflow_T, CNS::h_prob_parm->massfrac.begin(), 
               cs);
    Real v = CNS::h_prob_parm->inflow_mach * cs;

    inflow_state[URHO] = rho;
    inflow_state[UMX]  = 0.0;
    inflow_state[UMY]  = 0.0;
    inflow_state[UMZ]  = rho * v;
    inflow_state[UEDEN] = rho * eint + 0.5 * rho * v * v;
    inflow_state[UEINT] = rho * eint;
    inflow_state[UTEMP] = CNS::h_prob_parm->inflow_T;
    for (int i = 0; i < NUM_SPECIES; ++i) {
        inflow_state[UFS + i] = rho * CNS::h_prob_parm->massfrac[i];
    }

    Gpu::copyAsync(Gpu::hostToDevice, inflow_state.data(),
                   inflow_state.data() + NUM_STATE,
                   CNS::h_prob_parm->inflow_state);
    Gpu::streamSynchronize();
}
}
