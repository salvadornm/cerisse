#include "derive.H"
#include "CNS.H"
#include "parm.H"

using namespace amrex;

void cns_derpres(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const dat = datafab.array();
    auto pfab = derfab.array(dcomp);
    // Parm const* parm = CNS::d_parm;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // p(i,j,k,dcomp) = (parm->eos_gamma-1.)*rhoe(i,j,k); // change to pelephysics
        const amrex::Real rho = dat(i, j, k, URHO);
        const amrex::Real rhoinv = 1.0 / rho;
        amrex::Real massfrac[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; ++n) {
            massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
        }
        auto eos = pele::physics::PhysicsType::eos();
        eos.RTY2P(rho, dat(i, j, k, UTEMP), massfrac, pfab(i, j, k));
    });
}

void cns_dervel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const FArrayBox& datfab, const Geometry& /*geomdata*/,
                Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const dat = datfab.array();
    auto vel = derfab.array();
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vel(i,j,k,dcomp) = dat(i,j,k,1) / dat(i,j,k,0);
    });
}

void cns_dermagvort(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datfab.const_array();
  auto vort = derfab.array(dcomp);

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3);
  amrex::Elixir local_eli = local.elixir();
  auto larr = local.array();

  const auto& flag_fab = amrex::getEBCellFlagFab(datfab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;

  // Convert momentum to velocity
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoinv = 1.0 / dat(i, j, k, URHO);
    AMREX_D_TERM(larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoinv; ,
                 larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoinv; ,
                 larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoinv;)
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0); , 
               const amrex::Real dy = geomdata.CellSize(1); ,
               const amrex::Real dz = geomdata.CellSize(2););

  // Calculate vorticity
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im; int ip; , 
                 int jm; int jp; , 
                 int km; int kp;)

    // if fab is all regular -> call regular idx and weights
    // otherwise
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip); ,
                 get_idx(j, 1, all_regular, flags(i, j, k), jm, jp); ,
                 get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip); ,
                 const amrex::Real wj = get_weight(jm, jp); ,
                 const amrex::Real wk = get_weight(km, kp);)

    AMREX_D_TERM(
        vort(i, j, k) = 0.0 * dx;
        ,
        const amrex::Real vx = wi * (larr(ip, j, k, 1) - larr(im, j, k, 1)) / dx;
        const amrex::Real uy = wj * (larr(i, jp, k, 0) - larr(i, jm, k, 0)) / dy;
        const amrex::Real v3 = vx - uy;
        ,
        const amrex::Real wx = wi * (larr(ip, j, k, 2) - larr(im, j, k, 2)) / dx;
        const amrex::Real wy = wj * (larr(i, jp, k, 2) - larr(i, jm, k, 2)) / dy;
        const amrex::Real uz = wk * (larr(i, j, kp, 0) - larr(i, j, km, 0)) / dz;
        const amrex::Real vz = wk * (larr(i, j, kp, 1) - larr(i, j, km, 1)) / dz;
        const amrex::Real v1 = wy - vz; 
        const amrex::Real v2 = uz - wx;);
    vort(i, j, k) = sqrt(AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void cns_derdivu(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datfab.const_array();
  auto divu = derfab.array(dcomp);

  const auto& flag_fab = amrex::getEBCellFlagFab(datfab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0); , 
               const amrex::Real dy = geomdata.CellSize(1); ,
               const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im; int ip;, 
                 int jm; int jp;, 
                 int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip); ,
                 get_idx(j, 1, all_regular, flags(i, j, k), jm, jp); ,
                 get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip); ,
                 const amrex::Real wj = get_weight(jm, jp); ,
                 const amrex::Real wk = get_weight(km, kp);)

    AMREX_D_TERM(
      const amrex::Real uhi = dat(ip, j, k, UMX) / dat(ip, j, k, URHO);
      const amrex::Real ulo = dat(im, j, k, UMX) / dat(im, j, k, URHO);
      ,
      const amrex::Real vhi = dat(i, jp, k, UMY) / dat(i, jp, k, URHO);
      const amrex::Real vlo = dat(i, jm, k, UMY) / dat(i, jm, k, URHO);
      , 
      const amrex::Real whi = dat(i, j, kp, UMZ) / dat(i, j, kp, URHO);
      const amrex::Real wlo = dat(i, j, km, UMZ) / dat(i, j, km, URHO););
    divu(i, j, k) = AMREX_D_TERM(
      wi * (uhi - ulo) / dx, +wj * (vhi - vlo) / dy, +wk * (whi - wlo) / dz);
  });
}

void cns_derdivrho(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datfab.const_array();
  auto divrho = derfab.array(dcomp);
    
  const auto& flag_fab = amrex::getEBCellFlagFab(datfab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0); , 
               const amrex::Real dy = geomdata.CellSize(1); ,
               const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im; int ip;, 
                 int jm; int jp;, 
                 int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip); ,
                 get_idx(j, 1, all_regular, flags(i, j, k), jm, jp); ,
                 get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip); ,
                 const amrex::Real wj = get_weight(jm, jp); ,
                 const amrex::Real wk = get_weight(km, kp);)

    AMREX_D_TERM(
      const amrex::Real xhi = dat(ip, j, k, URHO);
      const amrex::Real xlo = dat(im, j, k, URHO);
      ,
      const amrex::Real yhi = dat(i, jp, k, URHO);
      const amrex::Real ylo = dat(i, jm, k, URHO);
      , 
      const amrex::Real zhi = dat(i, j, kp, URHO);
      const amrex::Real zlo = dat(i, j, km, URHO););
    divrho(i, j, k) = AMREX_D_TERM(
      wi * (xhi - xlo) / dx, +wj * (yhi - ylo) / dy, +wk * (zhi - zlo) / dz);
  });
}

void cns_dermachnumber(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& /*geomdata*/,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datfab.const_array();
  auto mach = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    const amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real cs;
    eos.RTY2Cs(rho, T, massfrac, cs);
    const amrex::Real momxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX);
    const amrex::Real momysq = dat(i, j, k, UMY) * dat(i, j, k, UMY);
    const amrex::Real momzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ);
    mach(i, j, k) = sqrt(momxsq + momysq + momzsq) / dat(i, j, k, URHO) / cs;
  });
}