#include "CNS.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

struct CnsFillExtDir
{
    Real* inflow_state = nullptr;

    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real /*time*/,
                     const BCRec* bcr, const int bcomp,
                     const int /*orig_comp*/) const
    {
        // Print() << "(" << dcomp << "," << numcomp << ") ";

        const Box& domain_box = geom.Domain();

        const BCRec& bc = bcr[bcomp];
        // const int* bc = bcr[bcomp]->data();

        // AMREX_D_TERM(int i = iv[0];,
        //              int j = iv[1];,
        //              int k = iv[2];);

        // The combustor problem is hard-wired for inflow at low-z right now
        // Print() << i << j << k << dcomp << std::endl;
        if ((bc.lo(2) == BCType::ext_dir) && (iv[2] < domain_box.smallEnd(2)))
        // if ((bc[2] == amrex::BCType::ext_dir) && (iv[2] < domain_box.smallEnd(2)))
        {
            for (int nf = 0; nf <= NUM_FIELD; ++nf) {
                for (int nc = 0; nc < NVAR; ++nc) {
                    // Print() << nf << " " << NVAR << " " << dcomp + nc << std::endl;
                    if (nf*NVAR + nc < numcomp)
                        dest(iv, dcomp + nf*NVAR + nc) = inflow_state[dcomp + nc];
                }
            }
        }
    }
};

struct CnsReactFillExtDir
{
    AMREX_GPU_DEVICE
    void operator()(const amrex::IntVect& /*iv*/,
                    amrex::Array4<amrex::Real> const& /*dest*/,
                    const int /*dcomp*/,
                    const int /*numcomp*/,
                    amrex::GeometryData const& /*geom*/,
                    const amrex::Real /*time*/,
                    const amrex::BCRec* /*bcr*/,
                    const int /*bcomp*/,
                    const int /*orig_comp*/) const
    {
        // do something for external Dirichlet (BCType::ext_dir)
    }
};

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the descriptor set up in CNS::variableSetUp.

void cns_bcfill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
    GpuBndryFuncFab<CnsFillExtDir> gpu_bndry_func(CnsFillExtDir{CNS::h_prob_parm->inflow_state});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}

void cns_react_bcfill(Box const& bx, FArrayBox& data,
                      const int dcomp, const int numcomp,
                      Geometry const& geom, const Real time,
                      const Vector<BCRec>& bcr, const int bcomp,
                      const int scomp)
{
    GpuBndryFuncFab<CnsReactFillExtDir> gpu_react_bndry_func(CnsReactFillExtDir{});
    gpu_react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}