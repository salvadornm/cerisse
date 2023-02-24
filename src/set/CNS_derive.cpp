#include "CNS_derive.H"
#include "CNS.H"
#include "CNS_parm.H"

using namespace amrex;

void derpres (const Box& bx, FArrayBox& pfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& rhoefab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const rhoe = rhoefab.array();
    auto       p    = pfab.array();
    Parm const* parm = CNS::d_parm;
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // p(i,j,k,dcomp) = (parm->eos_gamma-1.);//*rhoe(i,j,k);
        p(i,j,k,dcomp) = (1.4-1.)*rhoe(i,j,k);
    });
}

void dervelx(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,     
                const FArrayBox& datfab, const Geometry& /*geomdata*/, Real /*time*/, const int* /*bcrec*/, const int /*level*/)
{
    auto const dat = datfab.const_array();
    auto vel = derfab.array();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      vel(i, j, k, dcomp) = 0.0;//dat(i, j, k, 0);/// dat(i, j, k, 0);
    });
}

void dervely(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,     
                const FArrayBox& datfab, const Geometry& /*geomdata*/, Real /*time*/, const int* /*bcrec*/, const int /*level*/)
{
    auto const dat = datfab.const_array();
    auto vel = derfab.array();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      vel(i, j, k, dcomp) = 0.0;//dat(i, j, k, 0);/// dat(i, j, k, 0);
    });
}

void dervelz(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,     
                const FArrayBox& datfab, const Geometry& /*geomdata*/, Real /*time*/, const int* /*bcrec*/, const int /*level*/)
{
    auto const dat = datfab.const_array();
    auto vel = derfab.array();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      vel(i, j, k, dcomp) = 0.0;//dat(i, j, k, 0);/// dat(i, j, k, 0);
    });
}
