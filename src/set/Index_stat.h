#ifndef INDEXSTAT_H_
#define INDEXSTAT_H_

#include <AMReX_GpuContainers.H>

#include <Index.h>

struct indicies_stat_t : indicies_t {

// Statistic space

static constexpr int record_velocity   = 1;
static constexpr int record_PTrho      = 0;

// velocity
#if (AMREX_SPACEDIM == 1)
static constexpr int NSTAT_VEL=2;              
#elif (AMREX_SPACEDIM == 2)
static constexpr int NSTAT_VEL=5;              
#elif  (AMREX_SPACEDIM == 3)
static constexpr int NSTAT_VEL=9;              
#else
static constexpr int NSTAT_VEL=0;              
#endif

// P,T,rho
static constexpr int NSTAT_THERM=6;              

static constexpr int NSTAT = NSTAT_VEL*record_velocity + NSTAT_THERM*record_PTrho;              


};

#endif