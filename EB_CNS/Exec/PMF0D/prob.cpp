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
    amrex::Real p = 1.5e7;
    amrex::Real T = 1500.0;
    // amrex::ParmParse pp("prob");
    // pp.query("p", p);
    // pp.query("T", T);

    CNS::h_prob_parm->massfrac[H2_ID] = 0.022949;
    CNS::h_prob_parm->massfrac[O2_ID] = 0.22765;
    CNS::h_prob_parm->massfrac[N2_ID] = 0.7494;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2RE(p, CNS::h_prob_parm->massfrac.begin(), T, 
               CNS::h_prob_parm->rho, CNS::h_prob_parm->eint);
    
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm+1, CNS::d_prob_parm);
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
                   amrex::Real time,
                   amrex::GeometryData const& /*geomdata*/, 
                   amrex::Array4<const amrex::Real> const& state, 
                   amrex::Array4<amrex::Real> const& ext_src, 
                   Parm const& /*parm*/,
                   ProbParm const& /*prob_parm*/)
{
}