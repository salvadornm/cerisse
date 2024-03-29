#include "prob_parm.H"

#include <AMReX_Arena.H>

#include "CNS.H"
#include "index_macros.H"

ProbParm::ProbParm()
{
  inflow_state = (amrex::Real*)The_Arena()->alloc(sizeof(Real) * NVAR);
}

ProbParm::~ProbParm() { The_Arena()->free(inflow_state); }
