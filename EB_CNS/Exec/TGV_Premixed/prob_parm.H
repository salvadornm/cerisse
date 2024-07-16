#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  int mode = 0; // 0: TGV, 1: HIT
  int inres; // HIT input resolution

  amrex::Real T0 = 800.0; // [K]
  amrex::Real p0 = 2.0e5; // [Ba]

  amrex::Real L;
  amrex::Real rho0;
  amrex::Real v0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {0.0};

  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_kinput = nullptr;
  amrex::Real ksgs_avg = 0.0;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_uinput;
  amrex::Vector<amrex::Real> h_vinput;
  amrex::Vector<amrex::Real> h_winput;
  amrex::Vector<amrex::Real> h_kinput;
  amrex::Gpu::DeviceVector<amrex::Real> uinput;
  amrex::Gpu::DeviceVector<amrex::Real> vinput;
  amrex::Gpu::DeviceVector<amrex::Real> winput;
  amrex::Gpu::DeviceVector<amrex::Real> kinput;
  ProbParmHost() : uinput(0), vinput(0), winput(0), kinput(0) {}
};

#endif
