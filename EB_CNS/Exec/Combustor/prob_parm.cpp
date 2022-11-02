#include "CNS.H"
#include "index_macros.H"
#include "prob_parm.H"

#include <AMReX_Arena.H>

ProbParm::ProbParm ()
{
    inflow_state = (amrex::Real*)The_Arena()->alloc(sizeof(Real)*NUM_STATE);
}

ProbParm::~ProbParm ()
{
    The_Arena()->free(inflow_state);
}
