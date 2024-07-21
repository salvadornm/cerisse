#include "derive.H"

#include "CNS.H"

using namespace amrex;

void cns_derpres(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto pfab = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(amrex::Real vx = dat(i, j, k, UMX) * rhoinv;
                 , amrex::Real vy = dat(i, j, k, UMY) * rhoinv;
                 , amrex::Real vz = dat(i, j, k, UMZ) * rhoinv;);
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv -
                     0.5 * (AMREX_D_TERM(vx * vx, +vy * vy, +vz * vz));
    amrex::Real T = dat(i, j, k, UTEMP);
    auto eos = pele::physics::PhysicsType::eos();
    // eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2P(rho, T, massfrac, pfab(i, j, k));
  });
}

void cns_dereint(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto eifab = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(amrex::Real vx = dat(i, j, k, UMX) * rhoinv;
                 , amrex::Real vy = dat(i, j, k, UMY) * rhoinv;
                 , amrex::Real vz = dat(i, j, k, UMZ) * rhoinv;);
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv -
                     0.5 * (AMREX_D_TERM(vx * vx, +vy * vy, +vz * vz));

    eifab(i, j, k) = ei;
  });
}

void cns_dervel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const FArrayBox& datafab, const Geometry& /*geomdata*/,
                Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto vel = derfab.array();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(vel(i, j, k, dcomp) = dat(i, j, k, UMX) / dat(i, j, k, URHO);
                 , vel(i, j, k, dcomp + 1) = dat(i, j, k, UMY) / dat(i, j, k, URHO);
                 , vel(i, j, k, dcomp + 2) = dat(i, j, k, UMZ) / dat(i, j, k, URHO);)
  });
}

void cns_dermagvort(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto vort = derfab.array(dcomp);

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3);
  amrex::Elixir local_eli = local.elixir();
  auto larr = local.array();

#if CNS_USE_EB
  const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    // derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      vort(i, j, k) = 0.0;
    });
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;
#endif

  // Convert momentum to velocity
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoinv = 1.0 / dat(i, j, k, URHO);
    AMREX_D_TERM(larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoinv;
                 , larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoinv;
                 , larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoinv;)
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate vorticity
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if CNS_USE_EB
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip);
                 , const amrex::Real wj = get_weight(jm, jp);
                 , const amrex::Real wk = get_weight(km, kp);)
#else
    AMREX_D_TERM(int im = i - 1; int ip = i + 1; amrex::Real wi = 0.5;
                 , int jm = j - 1; int jp = j + 1; amrex::Real wj = 0.5;
                 , int km = k - 1; int kp = k + 1; amrex::Real wk = 0.5;)
#endif

    AMREX_D_TERM(
      vort(i, j, k) = 0.0;
      , const amrex::Real vx = wi * (larr(ip, j, k, 1) - larr(im, j, k, 1)) / dx;
      const amrex::Real uy = wj * (larr(i, jp, k, 0) - larr(i, jm, k, 0)) / dy;
      const amrex::Real v3 = vx - uy;
      , const amrex::Real wx = wi * (larr(ip, j, k, 2) - larr(im, j, k, 2)) / dx;
      const amrex::Real wy = wj * (larr(i, jp, k, 2) - larr(i, jm, k, 2)) / dy;
      const amrex::Real uz = wk * (larr(i, j, kp, 0) - larr(i, j, km, 0)) / dz;
      const amrex::Real vz = wk * (larr(i, j, kp, 1) - larr(i, j, km, 1)) / dz;
      const amrex::Real v1 = wy - vz; const amrex::Real v2 = uz - wx;);
    vort(i, j, k) = sqrt(AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

void cns_derdivu(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& geomdata, Real /*time*/,
                 const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto divu = derfab.array(dcomp);

#if CNS_USE_EB
  const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    // derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      divu(i, j, k) = 0.0;
    });
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;
#endif

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if CNS_USE_EB
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip);
                 , const amrex::Real wj = get_weight(jm, jp);
                 , const amrex::Real wk = get_weight(km, kp);)
#else
    AMREX_D_TERM(int im = i - 1; int ip = i + 1; amrex::Real wi = 0.5;
                 , int jm = j - 1; int jp = j + 1; amrex::Real wj = 0.5;
                 , int km = k - 1; int kp = k + 1; amrex::Real wk = 0.5;)
#endif

    AMREX_D_TERM(const amrex::Real uhi = dat(ip, j, k, UMX) / dat(ip, j, k, URHO);
                 const amrex::Real ulo = dat(im, j, k, UMX) / dat(im, j, k, URHO);
                 , const amrex::Real vhi = dat(i, jp, k, UMY) / dat(i, jp, k, URHO);
                 const amrex::Real vlo = dat(i, jm, k, UMY) / dat(i, jm, k, URHO);
                 , const amrex::Real whi = dat(i, j, kp, UMZ) / dat(i, j, kp, URHO);
                 const amrex::Real wlo = dat(i, j, km, UMZ) / dat(i, j, km, URHO););
    divu(i, j, k) = AMREX_D_TERM(wi * (uhi - ulo) / dx, +wj * (vhi - vlo) / dy,
                                 +wk * (whi - wlo) / dz);
  });
}

void cns_derdivrho(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datafab, const Geometry& geomdata, Real /*time*/,
                   const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto divrho = derfab.array(dcomp);

#if CNS_USE_EB
  const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    // derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      divrho(i, j, k) = 0.0;
    });
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;
#endif

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if CNS_USE_EB
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip);
                 , const amrex::Real wj = get_weight(jm, jp);
                 , const amrex::Real wk = get_weight(km, kp);)
#else
    AMREX_D_TERM(int im = i - 1; int ip = i + 1; amrex::Real wi = 0.5;
                 , int jm = j - 1; int jp = j + 1; amrex::Real wj = 0.5;
                 , int km = k - 1; int kp = k + 1; amrex::Real wk = 0.5;)
#endif

    AMREX_D_TERM(const amrex::Real xhi = dat(ip, j, k, URHO);
                 const amrex::Real xlo = dat(im, j, k, URHO);
                 , const amrex::Real yhi = dat(i, jp, k, URHO);
                 const amrex::Real ylo = dat(i, jm, k, URHO);
                 , const amrex::Real zhi = dat(i, j, kp, URHO);
                 const amrex::Real zlo = dat(i, j, km, URHO););
    divrho(i, j, k) = AMREX_D_TERM(wi * (xhi - xlo) / dx, +wj * (yhi - ylo) / dy,
                                   +wk * (zhi - zlo) / dz);
  });
}

void cns_dermachnumber(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datafab, const Geometry& /*geomdata*/,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto mach = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    const amrex::Real q2 = (AMREX_D_TERM(dat(i, j, k, UMX) * dat(i, j, k, UMX),
                                         +dat(i, j, k, UMY) * dat(i, j, k, UMY),
                                         +dat(i, j, k, UMZ) * dat(i, j, k, UMZ))) *
                           rhoinv * rhoinv;
    const amrex::Real ei = dat(i, j, k, UEDEN) * rhoinv - 0.5 * q2;
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real cs;
    auto eos = pele::physics::PhysicsType::eos();
    eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2Cs(rho, T, massfrac, cs);

    mach(i, j, k) = sqrt(q2) / cs;
  });
}

// void CNS::cns_derextsrc (
//     const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
//     const FArrayBox& datafab, const Geometry& /*geomdata*/,
//     Real time, const int* /*bcrec*/, int /*level*/)
// {
//   auto const sarr = datafab.const_array();
//   auto srcarr = derfab.array(dcomp);

//   ProbParm const* lprobparm = d_prob_parm;
//   const auto geomdata = geom.data();

//   amrex::ParallelFor(bx,
//   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
//     CNS::fill_ext_src(i, j, k, time, geomdata, sarr, srcarr, *lprobparm);
//   });
// }

void cns_dercp(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const FArrayBox& datafab, const Geometry& /*geomdata*/, Real /*time*/,
               const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto cp_arr = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq =
                   dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv;
                 , const amrex::Real vysq =
                     dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv;
                 , const amrex::Real vzsq =
                     dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei =
      dat(i, j, k, UEDEN) * rhoinv - 0.5 * (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real /*T = 150.0,*/ cp;
    // eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2Cp(rho, T, massfrac, cp);

    cp_arr(i, j, k) = cp;
  });
}

void cns_dercv(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const FArrayBox& datafab, const Geometry& /*geomdata*/, Real /*time*/,
               const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto cv_arr = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq =
                   dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv;
                 , const amrex::Real vysq =
                     dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv;
                 , const amrex::Real vzsq =
                     dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei =
      dat(i, j, k, UEDEN) * rhoinv - 0.5 * (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T = dat(i, j, k, UTEMP);
    amrex::Real cv;
    // eos.REY2T(rho, ei, massfrac, T);
    eos.RTY2Cv(rho, T, massfrac, cv);

    cv_arr(i, j, k) = cv;
  });
}

void CNS::cns_dertranscoef(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                           int /*ncomp*/, const amrex::FArrayBox& datafab,
                           const amrex::Geometry& /*geomdata*/, amrex::Real /*time*/,
                           const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto Diag_arr = derfab.array(dcomp);
  auto mu_arr = derfab.array(dcomp + NUM_SPECIES);
  auto xi_arr = derfab.array(dcomp + NUM_SPECIES + 1);
  auto lam_arr = derfab.array(dcomp + NUM_SPECIES + 2);
  auto const* ltransparm = trans_parms.device_trans_parm();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, URHO);
    const amrex::Real rhoinv = 1.0 / rho;
    AMREX_D_TERM(const amrex::Real vxsq =
                   dat(i, j, k, UMX) * dat(i, j, k, UMX) * rhoinv * rhoinv;
                 , const amrex::Real vysq =
                     dat(i, j, k, UMY) * dat(i, j, k, UMY) * rhoinv * rhoinv;
                 , const amrex::Real vzsq =
                     dat(i, j, k, UMZ) * dat(i, j, k, UMZ) * rhoinv * rhoinv;);
    const amrex::Real ei =
      dat(i, j, k, UEDEN) * rhoinv - 0.5 * (AMREX_D_TERM(vxsq, +vysq, +vzsq));
    amrex::Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, UFS + n) * rhoinv;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real T = dat(i, j, k, UTEMP);
    // eos.REY2T(rho, ei, massfrac, T);

    amrex::Real mu, xi, lam;
    amrex::Real Ddiag[NUM_SPECIES];
    const bool get_xi = true, get_mu = true, get_lam = true, get_Ddiag = true,
               get_chi = false;
    auto trans = pele::physics::PhysicsType::transport();
    trans.transport(get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, massfrac,
                    Ddiag, nullptr, mu, xi, lam, ltransparm);

    mu_arr(i, j, k) = mu;
    xi_arr(i, j, k) = xi;
    lam_arr(i, j, k) = lam;
    for (int n = 0; n < NUM_SPECIES; ++n) { Diag_arr(i, j, k, n) = Ddiag[n]; }
  });
}

void cns_dermassfrac(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
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

void cns_dermolefrac(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
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
    for (int n = 0; n < NUM_SPECIES; ++n) { molefrac(i, j, k, n) = mole[n]; }
  });
}

#if NUM_FIELD > 0 // Cannot compute subgrid-variance and tke without SF
void cns_dervaru(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto rs = derfab.array(dcomp); // reynolds stresses

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    unsigned int ctr = 0; // counter
    for (int dir1 = 0; dir1 < AMREX_SPACEDIM; ++dir1) {
      for (int dir2 = 0; dir2 < AMREX_SPACEDIM; ++dir2) {
        if (dir1 >= dir2) {
          rs(i, j, k, ctr) = 0.0;
          amrex::Real meanu = dat(i, j, k, UMX + dir1) / dat(i, j, k, URHO);
          amrex::Real meanv = dat(i, j, k, UMX + dir2) / dat(i, j, k, URHO);
          for (int nf = 1; nf < NUM_FIELD; ++nf) {
            amrex::Real u =
              dat(i, j, k, nf * NVAR + UMX + dir1) / dat(i, j, k, nf * NVAR + URHO);
            amrex::Real v =
              dat(i, j, k, nf * NVAR + UMX + dir2) / dat(i, j, k, nf * NVAR + URHO);
            rs(i, j, k, ctr) += (u - meanu) * (v - meanv);
          }
          rs(i, j, k, ctr) /=
            (amrex::Real)(NUM_FIELD - 1); // over N-1 instead of N is more accurate

          ctr++;
        }
      }
    }
  });
}

void cns_dertke(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                const FArrayBox& datafab, const Geometry& /*geomdata*/,
                Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto tke = derfab.array(dcomp); // trace of reynolds stresses

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    tke(i, j, k) = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Real meanu = dat(i, j, k, UMX + dir) / dat(i, j, k, URHO);
      amrex::Real u;
      for (int nf = 1; nf < NUM_FIELD; ++nf) {
        u = dat(i, j, k, nf * NVAR + UMX + dir) / dat(i, j, k, nf * NVAR + URHO);
        tke(i, j, k) += (u - meanu) * (u - meanu);
      }
    }
    tke(i, j, k) /=
      amrex::Real(NUM_FIELD - 1); // over N-1 instead of N is more accurate
    tke(i, j, k) *= 0.5;
  });
}

void cns_dervary(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto vary = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = 0; n < NUM_SPECIES; ++n) {
      vary(i, j, k, n) = 0.0;
      amrex::Real meany = dat(i, j, k, UFS + n) / dat(i, j, k, URHO);
      for (int nf = 1; nf < NUM_FIELD; ++nf) {
        amrex::Real y =
          dat(i, j, k, nf * NVAR + UFS + n) / dat(i, j, k, nf * NVAR + URHO);
        vary(i, j, k, n) += (y - meany) * (y - meany);
      }
      vary(i, j, k, n) /= amrex::Real(NUM_FIELD - 1);
    }
  });
}

void cns_dervarp(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto varp = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    varp(i, j, k) = 0.0;
    amrex::Real p[NUM_FIELD + 1];
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      amrex::Real rho = dat(i, j, k, nf * NVAR + URHO);
      amrex::Real rhoinv = 1.0 / rho;
      AMREX_D_TERM(amrex::Real ux = dat(i, j, k, nf * NVAR + UMX) * rhoinv;
                  , amrex::Real uy = dat(i, j, k, nf * NVAR + UMY) * rhoinv;
                  , amrex::Real uz = dat(i, j, k, nf * NVAR + UMZ) * rhoinv;);
      amrex::Real ei =
        dat(i, j, k, nf * NVAR + UEDEN) * rhoinv -
        amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));
      ;
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n)
        Y[n] = dat(i, j, k, nf * NVAR + UFS + n) * rhoinv;

      amrex::Real T = dat(i, j, k, UTEMP);
      auto eos = pele::physics::PhysicsType::eos();
      // eos.REY2T(rho, ei, Y, T);
      eos.RTY2P(rho, T, Y, p[nf]);

      if (nf > 0) varp(i, j, k) += (p[nf] - p[0]) * (p[nf] - p[0]);
    }
    varp(i, j, k) /= amrex::Real(NUM_FIELD - 1);    
  });
}

void cns_dervarT(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datafab, const Geometry& /*geomdata*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto varT = derfab.array(dcomp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    varT(i, j, k) = 0.0;
    amrex::Real T[NUM_FIELD + 1] = {0.0};
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      amrex::Real rho = dat(i, j, k, nf * NVAR + URHO);
      amrex::Real rhoinv = 1.0 / rho;
      AMREX_D_TERM(amrex::Real ux = dat(i, j, k, nf * NVAR + UMX) * rhoinv;
                  , amrex::Real uy = dat(i, j, k, nf * NVAR + UMY) * rhoinv;
                  , amrex::Real uz = dat(i, j, k, nf * NVAR + UMZ) * rhoinv;);
      amrex::Real ei =
        dat(i, j, k, nf * NVAR + UEDEN) * rhoinv -
        amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));
      ;
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n)
        Y[n] = dat(i, j, k, nf * NVAR + UFS + n) * rhoinv;

      auto eos = pele::physics::PhysicsType::eos();
      T[nf] = dat(i, j, k, nf * NVAR + UTEMP);
      // eos.REY2T(rho, ei, Y, T[nf]);

      if (nf > 0) varT(i, j, k) += (T[nf] - T[0]) * (T[nf] - T[0]);
    }
    varT(i, j, k) /= amrex::Real(NUM_FIELD - 1);    
  });
}
#endif

void cns_dervelgrad(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto gradu = derfab.array(dcomp);

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3, amrex::The_Async_Arena());
  auto larr = local.array();

#if CNS_USE_EB
  const auto& flag_fab = amrex::getEBCellFlagFab(datafab);
  const auto& typ = flag_fab.getType(bx);
  if (typ == amrex::FabType::covered) {
    // derfab.setVal<amrex::RunOn::Device>(0.0, bx);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      gradu(i, j, k) = 0.0;
    });
    return;
  }
  const auto& flags = flag_fab.const_array();
  const bool all_regular = typ == amrex::FabType::regular;
#endif

  // Convert momentum to velocity
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoinv = 1.0 / dat(i, j, k, URHO);
    AMREX_D_TERM(larr(i, j, k, 0) = dat(i, j, k, UMX) * rhoinv;
                 , larr(i, j, k, 1) = dat(i, j, k, UMY) * rhoinv;
                 , larr(i, j, k, 2) = dat(i, j, k, UMZ) * rhoinv;)
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2););

  // Calculate velocity gradients
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if CNS_USE_EB
    AMREX_D_TERM(int im; int ip;, int jm; int jp;, int km; int kp;)
    AMREX_D_TERM(get_idx(i, 0, all_regular, flags(i, j, k), im, ip);
                 , get_idx(j, 1, all_regular, flags(i, j, k), jm, jp);
                 , get_idx(k, 2, all_regular, flags(i, j, k), km, kp);)
    AMREX_D_TERM(const amrex::Real wi = get_weight(im, ip);
                 , const amrex::Real wj = get_weight(jm, jp);
                 , const amrex::Real wk = get_weight(km, kp);)
#else
    AMREX_D_TERM(int im = i - 1; int ip = i + 1; amrex::Real wi = 0.5;
                 , int jm = j - 1; int jp = j + 1; amrex::Real wj = 0.5;
                 , int km = k - 1; int kp = k + 1; amrex::Real wk = 0.5;)
#endif

    AMREX_D_TERM(
      const Real dudx = wi * (larr(ip, j, k, 0) - larr(im, j, k, 0)) / dx;
      , const Real dudy = wj * (larr(i, jp, k, 0) - larr(i, jm, k, 0)) / dy;
      const Real dvdy = wj * (larr(i, jp, k, 1) - larr(i, jm, k, 1)) / dy;
      const Real dvdx = wi * (larr(ip, j, k, 1) - larr(im, j, k, 1)) / dx;
      , const Real dudz = wk * (larr(i, j, kp, 0) - larr(i, j, km, 0)) / dz;
      const Real dvdz = wk * (larr(i, j, kp, 1) - larr(i, j, km, 1)) / dz;
      const Real dwdz = wk * (larr(i, j, kp, 2) - larr(i, j, km, 2)) / dz;
      const Real dwdx = wi * (larr(ip, j, k, 2) - larr(im, j, k, 2)) / dx;
      const Real dwdy = wj * (larr(i, jp, k, 2) - larr(i, jm, k, 2)) / dy;)

    gradu(i, j, k) = sqrt(AMREX_D_TERM(
      dudx * dudx, +dudy * dudy + dvdx * dvdx + dvdy * dvdy,
      +dudz * dudz + dvdz * dvdz + dwdx * dwdx + dwdy * dwdy + dwdz * dwdz));
  });
}

// Shock sensor 
void cns_dershocksensor(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                        const FArrayBox& datafab, const Geometry& geomdata,
                        Real time, const int* bcrec, int level)
{
  // Get divu and magvort
  FArrayBox divu_fab(bx, 1, amrex::The_Async_Arena());
  FArrayBox magvort_fab(bx, 1, amrex::The_Async_Arena());
  cns_derdivu(bx, divu_fab, 0, 1, datafab, geomdata, time, bcrec, level);
  cns_dermagvort(bx, magvort_fab, 0, 1, datafab, geomdata, time, bcrec, level);
  
  auto const divu = divu_fab.const_array();
  auto const magvort = magvort_fab.const_array();
  auto shock_sensor = derfab.array(dcomp);
  const Real* dx = geomdata.CellSize();
  auto const sarr = datafab.const_array();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // "Shock" sensor
    Real divu2 = divu(i, j, k) * divu(i, j, k);
    Real magvort2 = magvort(i, j, k) * magvort(i, j, k);    
    // shock_sensor(i, j, k) = divu2 / (divu2 + magvort2 + 1.0e-3); // Ducros
    Real magvel2 = (AMREX_D_TERM(sarr(i, j, k, UMX) * sarr(i, j, k, UMX),
                                 +sarr(i, j, k, UMY) * sarr(i, j, k, UMY),
                                 +sarr(i, j, k, UMZ) * sarr(i, j, k, UMZ))) /
                   sarr(i, j, k, URHO) / sarr(i, j, k, URHO);
    Real cellsize2 =
      std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]), 2.0 / Real(amrex::SpaceDim));
    // shock_sensor(i, j, k) = divu2 / (divu2 + magvel2 / cellsize2 + 1.0e-6); // Hendrickson, Kartha, Candler
    shock_sensor(i, j, k) = divu2 / (divu2 + magvort2 + std::sqrt(magvel2 / cellsize2)); // mixed
    // shock_sensor(i, j, k) *= divu(i, j, k) < 0.0 ? 1.0 : 0.0; // only tag compression

    // Partial density discontinuity sensor
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      AMREX_D_TERM(
        Real drYdx =
          (sarr(i + 1, j, k, UFS + ns) - sarr(i - 1, j, k, UFS + ns)) / dx[0];
        , Real drYdy =
            (sarr(i, j + 1, k, UFS + ns) - sarr(i, j - 1, k, UFS + ns)) / dx[1];
        , Real drYdz =
            (sarr(i, j, k + 1, UFS + ns) - sarr(i, j, k - 1, UFS + ns)) / dx[2];)
      Real drY2 =
        std::sqrt(AMREX_D_TERM(drYdx * drYdx, +drYdy * drYdy, +drYdz * drYdz));
      Real rY2 = sarr(i, j, k, UFS + ns) * sarr(i, j, k, UFS + ns);
      shock_sensor(i, j, k) =
        amrex::max(shock_sensor(i, j, k), drY2 / (drY2 + rY2 / cellsize2 + 1.0e-6));
    }
  });
}

void cns_derturbvisc(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datafab, const Geometry& geomdata,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
  auto const dat = datafab.const_array();
  auto mu_T = derfab.array(dcomp);
  auto xi_T = derfab.array(dcomp + 1);

  const auto dxinv = geomdata.InvCellSizeArray();

  const amrex::Box& gbx = amrex::grow(bx, 1);
  amrex::FArrayBox local(gbx, AMREX_SPACEDIM + 1,
                         amrex::The_Async_Arena()); // [rho, u, v, w]
  auto larr = local.array();

  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhoinv = 1.0 / dat(i, j, k, URHO);
    larr(i, j, k, QRHO) = dat(i, j, k, URHO);
    AMREX_D_TERM(larr(i, j, k, QU) = dat(i, j, k, UMX) * rhoinv;
                 , larr(i, j, k, QV) = dat(i, j, k, UMY) * rhoinv;
                 , larr(i, j, k, QW) = dat(i, j, k, UMZ) * rhoinv;)
  });

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if AMREX_SPACEDIM == 2
    Real delta = 1.0 / std::sqrt(dxinv[0] * dxinv[1]);
#else
    Real delta = std::pow(dxinv[0]*dxinv[1]*dxinv[2], -1./3.);
#endif

    CNS::les_model->mu_T_cc(i, j, k, larr, dxinv, delta, CNS::Cs, mu_T(i, j, k));
    CNS::les_model->xi_T_cc(i, j, k, larr, dxinv, delta, CNS::Cs, xi_T(i, j, k));
  });
}