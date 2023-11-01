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
 * @warning This is moved to psgs_model.cpp
 */
AMREX_FORCE_INLINE AMREX_GPU_DEVICE void vpdf_psgs_model(
  const amrex::Box& bx, amrex::Array4<amrex::Real>& s, const amrex::Real dt,
  const amrex::Real (*)[] /*dW[AMREX_SPACEDIM][NUM_FIELD]*/,
  const amrex::Real dx[AMREX_SPACEDIM])
{
  // Compute p
  const Box& bxg = amrex::grow(bx, 1);
  FArrayBox pfab(bxg, NUM_FIELD+1, The_Async_Arena());
  auto const& p = pfab.array();

  amrex::ParallelFor(bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    auto eos = pele::physics::PhysicsType::eos();

    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
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
      p(i, j, k, nf) = pres;

      // if (isnan(pres))
      //   std::cout << i << "," << j << "," << k << ": p is nan!!! T=" << T
      //             << ", rho=" << rho << ", ei=" << ei << std::endl;
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
        pm = p(iv - iv_dir, nf);
        um =
          s(iv - iv_dir, nf * NVAR + UMX + dir) / s(iv - iv_dir, nf * NVAR + URHO);
        // cell {i+1}
        pp = p(iv + iv_dir, nf);
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
        Real p_sgs = -0.5 * (mean_dp - dp[nf - 1]) / dx[dir] * dt;
        Real pu_sgs = -0.5 * (mean_dup - dup[nf - 1]) / dx[dir] * dt;

        s(i, j, k, nf * NVAR + UMX + dir) += p_sgs;
        s(i, j, k, nf * NVAR + UEDEN) += pu_sgs;
      }
    }
  });
}

/**
 * @brief Pressure-dilatation model for VPDF.
 * @param bx  box.
 * @param[in,out] s  state array.
 * @param dt  timestep size.
 * @param dW  generated white noise (@see UniqueRand).
 * @param dx  grid size.
 * @param C_e  turbulent dissipation model parameter.
 * @param shock_threshold  threshold in shock sensor.
 */
AMREX_FORCE_INLINE AMREX_GPU_DEVICE void vpdf_pd_model(
  const amrex::Box& bx, amrex::Array4<amrex::Real>& s, const amrex::Real dt,
  const amrex::Real dW[AMREX_SPACEDIM][NUM_FIELD],
  const amrex::Real dx[AMREX_SPACEDIM],  
  const amrex::Real C_e = 1.05,
  const amrex::Real shock_threshold = 0.02)
{
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::IntVect iv(AMREX_D_DECL(i,j,k));

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
    tke = amrex::max(tke, std::numeric_limits<amrex::Real>::epsilon()); // prevent divide by 0
    amrex::Real delta = std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]),
                                 1.0 / amrex::Real(AMREX_SPACEDIM));
    amrex::Real eps = C_e * std::pow(tke, 1.5) / delta;
    
    // Pressure-dilatation (Delaure & Pope 1997)
    amrex::Real p[NUM_FIELD+1], p_var = 0.0, cs_mean = 0.0, gamma[NUM_FIELD+1], ei_mean;
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      amrex::Real rho = s(i, j, k, nf * NVAR + URHO);
      amrex::Real rhoinv = 1.0 / rho;
      AMREX_D_TERM(amrex::Real ux = s(i, j, k, nf * NVAR + UMX) * rhoinv;
                  , amrex::Real uy = s(i, j, k, nf * NVAR + UMY) * rhoinv;
                  , amrex::Real uz = s(i, j, k, nf * NVAR + UMZ) * rhoinv;);
      amrex::Real ei = s(i, j, k, nf * NVAR + UEDEN) * rhoinv 
                       - amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n)
        Y[n] = s(i, j, k, nf * NVAR + UFS + n) * rhoinv;

      amrex::Real T = 0.0;
      auto eos = pele::physics::PhysicsType::eos();
      eos.REY2T(rho, ei, Y, T);
      eos.RTY2P(rho, T, Y, p[nf]);
      eos.RTY2G(rho, T, Y, gamma[nf]);

      if (nf == 0) {
        eos.RTY2Cs(rho, T, Y, cs_mean);
        ei_mean = ei;
      }

      if (nf > 0) 
        p_var += (p[nf] - p[0]) * (p[nf] - p[0]);
    }
    p_var = p_var / amrex::Real(NUM_FIELD-1);

    // Compute shock sensor
    // bool shock_sensor[AMREX_SPACEDIM+1] = {true};
    bool shock_sensor[AMREX_SPACEDIM+1] = {false};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        amrex::Real pp[3];
        for (int i_s: {-1, 0, 1}) {
          amrex::Real rho = s(iv+i_s*iv_dir, nf*NVAR+URHO);
          amrex::Real rhoinv = 1.0 / rho;
          AMREX_D_TERM(amrex::Real ux = s(iv + i_s*iv_dir, nf*NVAR+UMX) * rhoinv;
                      , amrex::Real uy = s(iv + i_s*iv_dir, nf*NVAR+UMY) * rhoinv;
                      , amrex::Real uz = s(iv + i_s*iv_dir, nf*NVAR+UMZ) * rhoinv;);
          amrex::Real ei = s(iv + i_s*iv_dir, nf*NVAR+UEDEN) * rhoinv 
                          - amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));
          amrex::Real Y[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; ++n)
            Y[n] = s(iv + i_s*iv_dir, nf*NVAR+UFS+n) * rhoinv;

          amrex::Real T = 0.0;
          auto eos = pele::physics::PhysicsType::eos();
          eos.REY2T(rho, ei, Y, T);
          eos.RTY2P(rho, T, Y, pp[i_s+1]);
        }
        shock_sensor[dir] |= (amrex::Math::abs(pp[2]-2*pp[1]+pp[0])/(pp[2]+2*pp[1]+pp[0]) > shock_threshold);
      }
      shock_sensor[AMREX_SPACEDIM] |= shock_sensor[dir];
    }

    Real M_turb = sqrt(2.0 * tke) / cs_mean; // turbulent Mach number
    Real tau_a = M_turb * tke / eps / (1.0 + M_turb*M_turb); // acoustic timescale (Zerman 1991), 1/sqrt((M_turb/3 + 1)*54) ~ 0.13
    // Real pe2 = s(i, j, k, URHO)*s(i, j, k, URHO)*cs_mean*cs_mean*2*tke
    //           *(pow(M_turb,2)+2*pow(M_turb,4))/(1+pow(M_turb,2)+2*pow(M_turb,4)); // equlibrium pressure squared (Zerman 1991)
    Real pe2 = 0.0;
    Real Pi_d = (p_var - pe2) / (2 * s(iv, URHO) * cs_mean * cs_mean * tau_a); // ~ <p' divu'> (Zerman 1991)
    p_var = amrex::max(p_var, std::numeric_limits<amrex::Real>::epsilon());
    Pi_d = amrex::max(Pi_d, std::numeric_limits<amrex::Real>::epsilon());

    amrex::Real u_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        mean_u = s(iv, UMX + dir) / s(iv, URHO);
        u = s(iv, nf * NVAR + UMX + dir) / s(iv, nf * NVAR + URHO);

        u_sgs = 1.0 / (2.0 * tke * s(iv, URHO) / Pi_d + dt) * (mean_u - u) * dt;
        if (!shock_sensor[dir]) {u_sgs = 0.0;}

        s(iv, nf * NVAR + UMX + dir) += s(iv, nf * NVAR + URHO) * u_sgs;
        s(iv, nf * NVAR + UEDEN) += 0.5 * s(iv, nf * NVAR + URHO) * ((u + u_sgs) * (u + u_sgs) - u * u); // conservative formulation
        // s(i, j, k, nf * NVAR + UEDEN) += s(i, j, k, nf * NVAR + URHO) * u * u_sgs; // force/work formulation
      }
      
      Real dp = 1.0 / (p_var / (Pi_d * p[0]) + dt) * (p[0] - p[nf]) * dt; 
      if (!shock_sensor[AMREX_SPACEDIM]) {dp = 0.0;}
      // dp += p[0] * sqrt(pe2 / tau_a / (1.0 - 1.0/gamma_mean) / (s(i, j, k, URHO)*s(i, j, k, URHO)*cs_mean*cs_mean) / ei_mean) * dW[dir][NUM_FIELD - nf - 1]; 

      s(iv, nf * NVAR + UEDEN) += dp / (gamma[nf] - 1.0);
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
  const amrex::Real C_e = 1.05,
  const amrex::Real shock_threshold = 0.02)
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
    tke = amrex::max(tke, std::numeric_limits<amrex::Real>::epsilon()); // prevent divide by 0

    // Compute shock sensor
    amrex::IntVect iv(AMREX_D_DECL(i,j,k));
    bool shock_sensor[AMREX_SPACEDIM+1] = {false};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir){
      amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);
      for (int nf = 0; nf <= NUM_FIELD; ++nf) {
        amrex::Real pp[3];
        for (int i: {-1, 0, 1}) {
          amrex::Real rho = s(iv+i*iv_dir, nf*NVAR+URHO);
          // pp[i+1]=rho;
          amrex::Real rhoinv = 1.0 / rho;
          AMREX_D_TERM(amrex::Real ux = s(iv+i*iv_dir, nf*NVAR+UMX) * rhoinv;
                      , amrex::Real uy = s(iv+i*iv_dir, nf*NVAR+UMY) * rhoinv;
                      , amrex::Real uz = s(iv+i*iv_dir, nf*NVAR+UMZ) * rhoinv;);
          amrex::Real ei = s(iv+i*iv_dir, nf*NVAR+UEDEN) * rhoinv 
                          - amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));;
          amrex::Real Y[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; ++n)
            Y[n] = s(iv+i*iv_dir, nf*NVAR+UFS+n) * rhoinv;

          amrex::Real T = 0.0;
          auto eos = pele::physics::PhysicsType::eos();
          eos.REY2T(rho, ei, Y, T);
          eos.RTY2P(rho, T, Y, pp[i+1]);
        }
        shock_sensor[dir] |= (amrex::Math::abs(pp[2]-2*pp[1]+pp[0])/(pp[2]+2*pp[1]+pp[0]) > shock_threshold);
      }
      shock_sensor[AMREX_SPACEDIM] |= shock_sensor[dir];
    }
    
    // Real M_turb = sqrt(2.0 * tke) / cs_mean; // turbulent Mach number
    // eps *= (1 + M_turb*M_turb);    

    // rhou_i = rhou_i + rho_i*Gij*(mean_u_i - u_i)*dt + rho_i*sqrt(C0*eps)*dW
    //                                                 need to /mean_rho? ^
    amrex::Real u_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        mean_u = s(i, j, k, UMX + dir) / s(i, j, k, URHO);
        u = s(i, j, k, nf * NVAR + UMX + dir) / s(i, j, k, nf * NVAR + URHO);

        // Langevin Gij = (0.5 + 0.75*C0)*eps/tke; B = sqrt(C0*eps) need to /mean_rho?
        // u_sgs = (0.5 + 0.75 * C0) * eps / tke * (mean_u - u) * dt +
        //         sqrt(C0 * eps) * dW[dir][nf - 1];
        if (shock_sensor[dir]) {
          u_sgs = 0.0; //(0.5 + 0.75 * C0) / (tke / eps + dt) * (mean_u - u) * dt;
        } else {
          u_sgs = (0.5 + 0.75 * C0) / (tke / eps + dt) * (mean_u - u) * dt +
                  sqrt(C0 * eps) * dW[dir][nf - 1];
        }

        s(i, j, k, nf * NVAR + UMX + dir) += s(i, j, k, nf * NVAR + URHO) * u_sgs;
        s(i, j, k, nf * NVAR + UEDEN) += 0.5 * s(i, j, k, nf * NVAR + URHO) * ((u + u_sgs) * (u + u_sgs) - u * u); // conservative formulation
        // s(i, j, k, nf * NVAR + UEDEN) += s(i, j, k, nf * NVAR + URHO) * u * u_sgs; // force/work formulation
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

      // // Energy (Q = E)
      // mean_Q = s(i, j, k, UEDEN) / s(i, j, k, URHO);
      // Q = s(i, j, k, nf * NVAR + UEDEN) / s(i, j, k, nf * NVAR + URHO);
      // Q_sgs = 0.5 * C_E * eps / tke * (mean_Q - Q) * dt;
      // // Q_sgs = 0.5 * C_E / (tke / eps + dt) * (mean_Q - Q) * dt;
      // s(i, j, k, nf * NVAR + UEDEN) += s(i, j, k, nf * NVAR + URHO) * Q_sgs;
      
      // Energy 2 (Q = Ei)
      mean_Q = s(i, j, k, UEDEN) / s(i, j, k, URHO) 
        - 0.5 * (AMREX_D_TERM(s(i, j, k, UMX) * s(i, j, k, UMX), 
                            + s(i, j, k, UMY) * s(i, j, k, UMY),
                            + s(i, j, k, UMZ) * s(i, j, k, UMZ))) 
                            / s(i, j, k, URHO) / s(i, j, k, URHO);
      Q = s(i, j, k, nf * NVAR + UEDEN) / s(i, j, k, nf * NVAR + URHO)
        - 0.5 * (AMREX_D_TERM(s(i, j, k, nf * NVAR + UMX) * s(i, j, k, nf * NVAR + UMX), 
                            + s(i, j, k, nf * NVAR + UMY) * s(i, j, k, nf * NVAR + UMY),
                            + s(i, j, k, nf * NVAR + UMZ) * s(i, j, k, nf * NVAR + UMZ))) 
                            / s(i, j, k, nf * NVAR + URHO) / s(i, j, k, nf * NVAR + URHO);
      Q_sgs = 0.5 * C_E * eps / tke * (mean_Q - Q) * dt;
      // Q_sgs = 0.5 * C_E / (tke / eps + dt) * (mean_Q - Q) * dt;
      s(i, j, k, nf * NVAR + UEDEN) += s(i, j, k, nf * NVAR + URHO) * Q_sgs;
    }
  });
}

#endif