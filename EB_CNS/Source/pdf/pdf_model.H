#ifndef CNS_PDF_MODEL_H_
#define CNS_PDF_MODEL_H_

#include "CNS.H"
#include "index_macros.H"
#include "recon.H"

/**
 * @brief Pull dp/dx of each field towards the mean <dp/dx>. This relaxation is
 * essential for stability in compressible VPDF. Ref
 * https://doi.org/10.1063/5.0039038.
 * @param bx  box.
 * @param[in,out] s  state array.
 * @param dt  timestep size.
 * @param dx  grid size.
 */
AMREX_FORCE_INLINE AMREX_GPU_DEVICE void vpdf_psgs_model(
  const amrex::Box& bx, amrex::Array4<amrex::Real>& s, const amrex::Real dt,
  const amrex::Real (*)[] /*dW[AMREX_SPACEDIM][NUM_FIELD]*/,
  const amrex::Real dx[AMREX_SPACEDIM])
{
  // Compute p
  const Box& bxg1 = amrex::grow(bx, 1);
  FArrayBox pfab(bxg1, NUM_FIELD, The_Async_Arena());
  auto const& p = pfab.array();

  amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    auto eos = pele::physics::PhysicsType::eos();

    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      amrex::Real rho = s(i, j, k, nf * NVAR + URHO);
      amrex::Real rhoinv = 1.0 / rho;
      AMREX_D_TERM(amrex::Real ux = s(i, j, k, nf * NVAR + UMX) * rhoinv;
                   , amrex::Real uy = s(i, j, k, nf * NVAR + UMY) * rhoinv;
                   , amrex::Real uz = s(i, j, k, nf * NVAR + UMZ) * rhoinv;);
      amrex::Real ei =
        s(i, j, k, nf * NVAR + UEDEN) * rhoinv -
        amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n)
        Y[n] = s(i, j, k, nf * NVAR + UFS + n) * rhoinv;

      amrex::Real T = 0, pres;
      eos.REY2T(rho, ei, Y, T);
      eos.RTY2P(rho, T, Y, pres);
      p(i, j, k, nf - 1) = pres;

      if (isnan(pres))
        std::cout << i << "," << j << "," << k << ": p is nan!!! T=" << T
                  << ", rho=" << rho << ", ei=" << ei << std::endl;
    }
  });

  // Compute (dp/dx)_sgs and (dup/dx)_sgs
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);
      amrex::Real mean_dp = 0.0;
      amrex::Real dp[NUM_FIELD];
      amrex::Real mean_dup = 0.0;
      amrex::Real dup[NUM_FIELD];

      amrex::Real um, up, pm, pp;
      for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        // cell {i-1}
        pm = p(iv - iv_dir, nf - 1);
        um =
          s(iv - iv_dir, nf * NVAR + UMX + dir) / s(iv - iv_dir, nf * NVAR + URHO);
        // cell {i+1}
        pp = p(iv + iv_dir, nf - 1);
        up =
          s(iv + iv_dir, nf * NVAR + UMX + dir) / s(iv + iv_dir, nf * NVAR + URHO);

        dp[nf - 1] = pp - pm; // p_{i+1} - p{i-1}
        mean_dp += dp[nf - 1];

        dup[nf - 1] = up * pp - um * pm; // up_{i+1} - up_{i-1}
        mean_dup += dup[nf - 1];
      }
      mean_dp /= amrex::Real(NUM_FIELD); // or use mean calculated p?
      mean_dup /= amrex::Real(NUM_FIELD);

      for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        // central difference for dpdx
        // mom_sgs = -(mean_dp - dp[nf-1]) / 0.5*dx[dir] * dt
        s(i, j, k, nf * NVAR + UMX + dir) -=
          0.5 * (mean_dp - dp[nf - 1]) / dx[dir] * dt;

        s(i, j, k, nf * NVAR + UEDEN) -=
          0.5 * (mean_dup - dup[nf - 1]) / dx[dir] * dt;
        // consistent e
        // amrex::Real mean_u = s(i, j, k, UMX + dir) / s(i, j, k, URHO);
        // s(i, j, k, nf * NVAR + UEDEN) -= mean_u * 0.5 * (mean_dp - dp[nf - 1]) /
        // dx[dir] * dt; amrex::Real u = s(i, j, k, nf * NVAR + UMX + dir) / s(i, j,
        // k, nf * NVAR + URHO); s(i, j, k, nf * NVAR + UEDEN) -= u * 0.5 * (mean_dp
        // - dp[nf - 1]) / dx[dir] * dt;
      }
    }
  });
}

/**
 * @brief Simplified Langevin Model for VPDF. Ref https://doi.org/10.1063/1.865723.
 * @param bx  box.
 * @param[in,out] s  state array.
 * @param dt  timestep size.
 * @param dW  generated white noise (@see UniqueRand).
 * @param dx  grid size.
 * @param C0,C_e  model parameters
 */
AMREX_FORCE_INLINE AMREX_GPU_DEVICE void vpdf_langevin_model(
  const amrex::Box& bx, amrex::Array4<amrex::Real>& s, const amrex::Real dt,
  const amrex::Real dW[AMREX_SPACEDIM][NUM_FIELD],
  const amrex::Real dx[AMREX_SPACEDIM], const amrex::Real C0 = 2.1,
  const amrex::Real C_e = 1.05) // 1.05? or 0.53?
{
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Compute turbulent statistics
    amrex::Real mean_u, u, tke = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      mean_u = s(i, j, k, UMX + dir) / s(i, j, k, URHO); // this is Favre-averaged
      for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        u = s(i, j, k, nf * NVAR + UMX + dir) / s(i, j, k, nf * NVAR + URHO);
        tke += (u - mean_u) * (u - mean_u);
      }
    }
    tke = 0.5 * tke / amrex::Real(NUM_FIELD - 1); // k = 0.5*(u^2 + v^2 + w^2)
    amrex::Real delta = std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]),
                                 1.0 / amrex::Real(AMREX_SPACEDIM));
    amrex::Real eps = C_e * std::pow(tke, 1.5) / delta;
    tke = amrex::max(tke, 1e-10); // prevent divide by 0

    // rhou_i = rhou_i + rho_i*Gij*(mean_u_i - u_i)*dt + rho_i*sqrt(C0*eps)*dW
    //                                                 need to /mean_rho? ^
    amrex::Real u_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        mean_u = s(i, j, k, UMX + dir) / s(i, j, k, URHO);
        u = s(i, j, k, nf * NVAR + UMX + dir) / s(i, j, k, nf * NVAR + URHO);

        // Gij = (0.5 + 0.75*C0)*eps/tke; B = sqrt(C0*eps) need to /mean_rho?
        u_sgs = (0.5 + 0.75 * C0) * eps / tke * (mean_u - u) * dt +
                sqrt(C0 * eps) * dW[dir][nf - 1];
        // u_sgs = (0.5 + 0.75 * C0) / (tke / eps + dt) * (mean_u - u) * dt +
        //         sqrt(C0 * eps) * dW[dir][nf - 1];

        s(i, j, k, nf * NVAR + UMX + dir) += s(i, j, k, nf * NVAR + URHO) * u_sgs;
      }
    }
  });
}

/**
 * @brief Interact-by-Exchange-with-Mean Model for SPDF. Ref TODO.
 * @param bx  box.
 * @param[in,out] s  state array.
 * @param dt  timestep size.
 * @param dx  grid size.
 * @param C_E,C_Y,C_e  model parameters.
 */
AMREX_FORCE_INLINE AMREX_GPU_DEVICE void spdf_iem_model(
  const amrex::Box& bx, amrex::Array4<amrex::Real>& s, const amrex::Real dt,
  const amrex::Real (*)[] /*dW[AMREX_SPACEDIM][NUM_FIELD]*/,
  const amrex::Real dx[AMREX_SPACEDIM], const amrex::Real C_E = 2.0,
  const amrex::Real C_Y = 2.0, const amrex::Real C_e = 1.05)
{
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Compute turbulent statistics
    amrex::Real mean_u, u, tke = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      mean_u = s(i, j, k, UMX + dir) / s(i, j, k, URHO); // this is Favre-averaged
      for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        u = s(i, j, k, nf * NVAR + UMX + dir) / s(i, j, k, nf * NVAR + URHO);
        tke += (u - mean_u) * (u - mean_u);
      }
    }
    tke = 0.5 * tke / amrex::Real(NUM_FIELD - 1); // k = 0.5*(u'^2 + v'^2 + w'^2)
    amrex::Real delta = std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]),
                                 1.0 / amrex::Real(AMREX_SPACEDIM));
    amrex::Real eps = C_e * std::pow(tke, 1.5) / delta;
    tke = amrex::max(tke, 1e-10); // prevent divide by 0

    // rhoQ_i = rhoQ_i + rho_i*0.5*C_Q*eps/k*(mean_y - y)*dt
    amrex::Real mean_Q, Q, Q_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        // Species (Q = Y)
        mean_Q = s(i, j, k, UFS + ns) / s(i, j, k, URHO);
        Q = s(i, j, k, nf * NVAR + UFS + ns) / s(i, j, k, nf * NVAR + URHO);
        Q_sgs = 0.5 * C_Y * eps / tke * (mean_Q - Q) * dt;
        // Q_sgs = 0.5 * C_Y  / (tke / eps + dt) * (mean_Q - Q) * dt;
        s(i, j, k, nf * NVAR + UFS + ns) += s(i, j, k, nf * NVAR + URHO) * Q_sgs;
      }

      // Energy (Q = E)
      mean_Q = s(i, j, k, UEDEN) / s(i, j, k, URHO);
      Q = s(i, j, k, nf * NVAR + UEDEN) / s(i, j, k, nf * NVAR + URHO);
      Q_sgs = 0.5 * C_E * eps / tke * (mean_Q - Q) * dt;
      // Q_sgs = 0.5 * C_E / (tke / eps + dt) * (mean_Q - Q) * dt;
      s(i, j, k, nf * NVAR + UEDEN) += s(i, j, k, nf * NVAR + URHO) * Q_sgs;
    }
  });
}

#endif