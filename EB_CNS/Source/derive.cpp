#include "derive.H"
#include "CNS.H"
#include "parm.H"

using namespace amrex;

void cns_dertemp (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datafab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.array();
  auto Tfab = derfab.array(dcomp);
  
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(amrex::Real vx = dat(i, j, k, UMX) * rhoinv; ,
                 amrex::Real vy = dat(i, j, k, UMY) * rhoinv; ,
                 amrex::Real vz = dat(i, j, k, UMZ) * rhoinv;);
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vx*vx, +vy*vy, +vz*vz));
    
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, massfrac, Tfab(i, j, k));
  });
}

void cns_derpres (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datafab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.array();
  auto pfab = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(amrex::Real vx = dat(i, j, k, UMX) * rhoinv; ,
                 amrex::Real vy = dat(i, j, k, UMY) * rhoinv; ,
                 amrex::Real vz = dat(i, j, k, UMZ) * rhoinv;);
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vx*vx, +vy*vy, +vz*vz));
    amrex::Real T;
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2P(rho, T, massfrac, pfab(i, j, k));
  });
}

void cns_dervel (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.array();
  auto vel = derfab.array();
  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    vel(i,j,k,dcomp) = dat(i,j,k,1) / dat(i,j,k,0);
  });
}

void cns_dermagvort (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datafab, const Geometry& geomdata,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto vort = derfab.array(dcomp);

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3);
  amrex::Elixir local_eli = local.elixir();
  auto larr = local.array();

  // const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  // const auto& typ = flag_fab.getType(bx);
  // if (typ == amrex::FabType::covered) {
  //   derfab.setVal<amrex::RunOn::Device>(0.0, bx);
  //   return;
  // }
  // const auto& flags = flag_fab.const_array();
  // const bool all_regular = typ == amrex::FabType::regular;

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
    // AMREX_D_TERM(int im; int ip; , 
    //              int jm; int jp; , 
    //              int km; int kp;)
    // if fab is all regular -> call regular idx and weights
    // otherwise
    // AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip); ,
    //              get_idx(j, 1, all_regular, flags(i, j, k), jm, jp); ,
    //              get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    // AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip); ,
    //              const amrex::Real wj = get_weight(jm, jp); ,
    //              const amrex::Real wk = get_weight(km, kp);)
    AMREX_D_TERM(int im = i-1; int ip = i+1; amrex::Real wi = 0.5; , 
                 int jm = j-1; int jp = j+1; amrex::Real wj = 0.5; , 
                 int km = k-1; int kp = k+1; amrex::Real wk = 0.5;)

    AMREX_D_TERM(
      vort(i, j, k) = 0.0;
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

void cns_derdivu (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datafab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto divu = derfab.array(dcomp);

  // const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  // const auto& typ = flag_fab.getType(bx);
  // if (typ == amrex::FabType::covered) {
  //   derfab.setVal<amrex::RunOn::Device>(0.0, bx);
  //   return;
  // }
  // const auto& flags = flag_fab.const_array();
  // const bool all_regular = typ == amrex::FabType::regular;

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0); , 
               const amrex::Real dy = geomdata.CellSize(1); ,
               const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // AMREX_D_TERM(int im; int ip;, 
    //              int jm; int jp;, 
    //              int km; int kp;)
    // AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip); ,
    //              get_idx(j, 1, all_regular, flags(i, j, k), jm, jp); ,
    //              get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    // AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip); ,
    //              const amrex::Real wj = get_weight(jm, jp); ,
    //              const amrex::Real wk = get_weight(km, kp);)
    AMREX_D_TERM(int im = i-1; int ip = i+1; amrex::Real wi = 0.5; , 
                 int jm = j-1; int jp = j+1; amrex::Real wj = 0.5; , 
                 int km = k-1; int kp = k+1; amrex::Real wk = 0.5;)

    AMREX_D_TERM(const amrex::Real uhi = dat(ip, j, k, UMX) / dat(ip, j, k, URHO);
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

void cns_derdivrho (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto divrho = derfab.array(dcomp);
  
  // const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  // const auto& typ = flag_fab.getType(bx);
  // if (typ == amrex::FabType::covered) {
  //   derfab.setVal<amrex::RunOn::Device>(0.0, bx);
  //   return;
  // }
  // const auto& flags = flag_fab.const_array();
  // const bool all_regular = typ == amrex::FabType::regular;

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0); , 
               const amrex::Real dy = geomdata.CellSize(1); ,
               const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // AMREX_D_TERM(int im; int ip;, 
    //              int jm; int jp;, 
    //              int km; int kp;)
    // AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip); ,
    //              get_idx(j, 1, all_regular, flags(i, j, k), jm, jp); ,
    //              get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    // AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip); ,
    //              const amrex::Real wj = get_weight(jm, jp); ,
    //              const amrex::Real wk = get_weight(km, kp);)
    AMREX_D_TERM(int im = i-1; int ip = i+1; amrex::Real wi = 0.5; , 
                 int jm = j-1; int jp = j+1; amrex::Real wj = 0.5; , 
                 int km = k-1; int kp = k+1; amrex::Real wk = 0.5;)

    AMREX_D_TERM(const amrex::Real xhi = dat(ip, j, k, URHO);
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

void cns_dermachnumber (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                        const FArrayBox& datafab, const Geometry& /*geomdata*/,
                        Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto mach = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv; ,
                 const amrex::Real vysq = dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv; ,
                 const amrex::Real vzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real cs, T;
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2Cs(rho, T, massfrac, cs);

    mach(i, j, k) = sqrt(AMREX_D_TERM(vxsq, +vysq, +vzsq)) / cs;
  });
}

// void cns_derextsrc (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
//                     const FArrayBox& datafab, const Geometry& geomdata,
//                     Real time, const int* /*bcrec*/, int /*level*/)
// {
//   auto const dat = datafab.array();
//   auto Sfab = derfab.array(dcomp);
  
//   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
//     CNS::fill_ext_src(i, j, k, time, geomdata, dat, Sfab, parm, prob_parm);
//   });
// }

void cns_dercp (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const FArrayBox& datafab, const Geometry& /*geomdata*/,
                Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto cp_arr = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv; ,
                 const amrex::Real vysq = dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv; ,
                 const amrex::Real vzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T, cp;
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2Cp(rho, T, massfrac, cp);

    cp_arr(i, j, k) = cp;
  });
}

void cns_dercv (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const FArrayBox& datafab, const Geometry& /*geomdata*/,
                Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto cv_arr = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv; ,
                 const amrex::Real vysq = dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv; ,
                 const amrex::Real vzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T, cv;
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2Cv(rho, T, massfrac, cv);

    cv_arr(i, j, k) = cv;
  });
}

void CNS::cns_dertranscoef (
  const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
  const amrex::FArrayBox& datafab, const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/) 
{
  auto const dat = datafab.const_array();
  auto Diag_arr = derfab.array(dcomp);
  auto mu_arr = derfab.array(dcomp+NUM_SPECIES);
  auto xi_arr = derfab.array(dcomp+NUM_SPECIES+1);
  auto lam_arr = derfab.array(dcomp+NUM_SPECIES+2);
  auto const* ltransparm = trans_parms.device_trans_parm();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq = dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv; ,
                 const amrex::Real vysq = dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv; ,
                 const amrex::Real vzsq = dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5* (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }    
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T;
    eos.REY2T(rho, ei, massfrac, T);

    amrex::Real mu, xi, lam;
    amrex::Real Ddiag[NUM_SPECIES];
    bool get_xi = true, get_mu = true, get_lam = true, get_Ddiag = true;
    auto trans = pele::physics::PhysicsType::transport();
    trans.transport(get_xi, get_mu, get_lam, get_Ddiag, T, rho, massfrac, 
                    Ddiag, mu, xi, lam, ltransparm);

    mu_arr(i, j, k) = mu;
    xi_arr(i, j, k) = xi;
    lam_arr(i, j, k) = lam;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Diag_arr(i, j, k, n) = Ddiag[n];
    }
  });
}

void cns_dermassfrac (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datafab, const Geometry& /*geomdata*/,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto massfrac = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac(i, j, k, n) = dat(i, j, k, UFS + n) * rhoinv;
    }
  });
}

void cns_dermolefrac (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datafab, const Geometry& /*geomdata*/,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto molefrac = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    amrex::Real mass[NUM_SPECIES], mole[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      mass[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.Y2X(mass, mole);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      molefrac(i, j, k, n) = mole[n];
    }
  });
}

void cns_dervaru (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datafab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto rs = derfab.array(dcomp); //reynolds stresses

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    unsigned int ctr = 0; // counter
    for (int dir1 = 0; dir1 < AMREX_SPACEDIM; ++dir1) {
    for (int dir2 = 0; dir2 < AMREX_SPACEDIM; ++dir2) {
      if (dir1 >= dir2) {
        rs(i, j, k, ctr) = 0.0;
        amrex::Real meanu = dat(i, j, k, UMX+dir1) / dat(i, j, k, URHO);
        amrex::Real meanv = dat(i, j, k, UMX+dir2) / dat(i, j, k, URHO);
        for (int nf = 1; nf < NUM_FIELD; ++nf) {
          amrex::Real u = dat(i, j, k, nf*NVAR + UMX+dir1) / dat(i, j, k, nf*NVAR + URHO);
          amrex::Real v = dat(i, j, k, nf*NVAR + UMX+dir2) / dat(i, j, k, nf*NVAR + URHO);
          rs(i, j, k, ctr) += (u - meanu) * (v - meanv);
        }
        rs(i, j, k, ctr) /= (amrex::Real)(NUM_FIELD - 1); //over N-1 instead of N is more accurate

        ctr++;
      }
    }
    }
  });
}

void cns_dertke (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto tke = derfab.array(dcomp); //trace of reynolds stresses

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    tke(i, j, k) = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Real meanu = dat(i, j, k, UMX+dir) / dat(i, j, k, URHO);
      amrex::Real u;
      for (int nf = 1; nf < NUM_FIELD; ++nf) {
        u = dat(i, j, k, nf*NVAR + UMX+dir) / dat(i, j, k, nf*NVAR + URHO);
        tke(i, j, k) += (u - meanu) * (u - meanu);
      }
    }
    tke(i, j, k) /= amrex::Real(NUM_FIELD - 1); //over N-1 instead of N is more accurate    
    tke(i, j, k) *= 0.5;
  });
}

void cns_dervary (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datafab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto vary = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = 0; n < NUM_SPECIES; ++n) { 
      vary(i, j, k, n) = 0.0;       
      for (int nf = 1; nf < NUM_FIELD; ++nf) {
        amrex::Real y = dat(i, j, k, nf*NVAR + UFS+n) / dat(i, j, k, nf*NVAR + URHO);
        vary(i, j, k, n) += y * y;
      }
      amrex::Real meany = dat(i, j, k, UFS+n) / dat(i, j, k, URHO);
      vary(i, j, k, n) = vary(i, j, k, n) / NUM_FIELD - meany * meany;
    }
  });
}

void cns_dervelgrad (const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datafab, const Geometry& geomdata, 
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto gradu = derfab.array(dcomp);

  AMREX_D_TERM(const Real dx = geomdata.CellSize(0); , 
               const Real dy = geomdata.CellSize(1); ,
               const Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(int im = i-1; int ip = i+1; , 
                 int jm = j-1; int jp = j+1; , 
                 int km = k-1; int kp = k+1;)

    AMREX_D_TERM(      Real uhi = dat(ip, j, k, UMX) / dat(ip, j, k, URHO);
                       Real ulo = dat(im, j, k, UMX) / dat(im, j, k, URHO);
                 const Real dudx = 0.5 * (uhi - ulo) / dx;
                 ,
                            uhi = dat(i, jp, k, UMX) / dat(i, jp, k, URHO);
                            ulo = dat(i, jm, k, UMX) / dat(i, jm, k, URHO);
                       Real vhi = dat(i, jp, k, UMY) / dat(i, jp, k, URHO);
                       Real vlo = dat(i, jm, k, UMY) / dat(i, jm, k, URHO);
                 const Real dudy = 0.5 * (uhi - ulo) / dy;
                 const Real dvdy = 0.5 * (vhi - vlo) / dy;
                            vhi = dat(ip, j, k, UMY) / dat(ip, j, k, URHO);
                            vlo = dat(im, j, k, UMY) / dat(im, j, k, URHO);
                 const Real dvdx = 0.5 * (vhi - vlo) / dx;
                 , 
                            uhi = dat(i, j, kp, UMX) / dat(i, j, kp, URHO);
                            ulo = dat(i, j, km, UMX) / dat(i, j, km, URHO);
                            vhi = dat(i, j, kp, UMY) / dat(i, j, kp, URHO);
                            vlo = dat(i, j, km, UMY) / dat(i, j, km, URHO);
                       Real whi = dat(i, j, kp, UMZ) / dat(i, j, kp, URHO);
                       Real wlo = dat(i, j, km, UMZ) / dat(i, j, km, URHO);
                 const Real dudz = 0.5 * (uhi - ulo) / dz;
                 const Real dvdz = 0.5 * (vhi - vlo) / dz;
                 const Real dwdz = 0.5 * (whi - wlo) / dz;
                            whi = dat(ip, j, k, UMZ) / dat(ip, j, k, URHO);
                            wlo = dat(im, j, k, UMZ) / dat(im, j, k, URHO);
                 const Real dwdx = 0.5 * (whi - wlo) / dx;
                            whi = dat(i, jp, k, UMZ) / dat(i, jp, k, URHO);
                            wlo = dat(i, jm, k, UMZ) / dat(i, jm, k, URHO);
                 const Real dwdy = 0.5 * (whi - wlo) / dy;
                 );


    gradu(i, j, k) = sqrt(AMREX_D_TERM(dudx*dudx, 
                                      +dudy*dudy +dvdx*dvdx+dvdy*dvdy, 
                                      +dudz*dudz +dvdz*dvdz +dwdx*dwdx+dwdy*dwdy+dwdz*dwdz));
  });
}