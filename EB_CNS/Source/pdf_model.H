#ifndef CNS_PDF_MODEL_H_
#define CNS_PDF_MODEL_H_

#include "CNS.H"
#include "index_macros.H"

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
vpdf_langevin_model (const amrex::Box& bx,
                     amrex::Array4<amrex::Real>& s,
                     const amrex::Real dt,
                     const amrex::Real dW[NUM_FIELD][AMREX_SPACEDIM],
                     const amrex::Real dx[AMREX_SPACEDIM])
{
  amrex::ParallelFor(bx, 
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    // Compute turbulent statistics
    amrex::Real mean_u, u, tke = 0.0;        
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO); // this is Favre-averaged
      for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
        tke += (u - mean_u)*(u - mean_u);
      }
    }
    tke = 0.5 * tke / amrex::Real(NUM_FIELD); // k = 0.5*(u^2 + v^2 + w^2)
    if (tke < 0.0) { amrex::Print() << tke; amrex::Error("TKE < 0!!!"); }
    amrex::Real deltainv = 1.0 / std::pow(AMREX_D_TERM(dx[0], * dx[1], * dx[2]), 1.0/amrex::Real(AMREX_SPACEDIM));
    amrex::Real eps = 0.53 * std::pow(tke, 1.5) * deltainv;
    eps = amrex::max(eps, 1e-10);

    // Compute pressure (store pressure of each field in p[])
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(0);
    amrex::Real mean_dp = 0.0;
    amrex::Real rho, rhoinv, ei, T, p, massfrac[NUM_SPECIES], dp[NUM_FIELD];
    auto eos = pele::physics::PhysicsType::eos();            
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      // i-1
      rho = s(iv-iv_dir, nf*NVAR + URHO);
      rhoinv = 1.0 / rho;
      AMREX_D_TERM(amrex::Real ux = s(iv-iv_dir, nf*NVAR + UMX) * rhoinv; ,
                   amrex::Real uy = s(iv-iv_dir, nf*NVAR + UMY) * rhoinv; ,
                   amrex::Real uz = s(iv-iv_dir, nf*NVAR + UMZ) * rhoinv;);
      ei = s(iv-iv_dir, nf*NVAR + UEDEN) * rhoinv - amrex::Real(0.5)*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = s(iv-iv_dir, nf*NVAR + UFS+n) * rhoinv;
      }
      eos.EY2T(ei, massfrac, T); 
      eos.RTY2P(rho, T, massfrac, p);
      if (isnan(p)) amrex::Print() << i << ": p_{i-1} is nan!!! T = " << T << 
                                                            ", ei = " << ei << 
                                                           ", rho = " << rho << " = 1/" << 1 / rho <<
                                                            ", mx = " << std::setprecision(5) << s(i-1, j, k, nf*NVAR + UMX) << 
                                                            ", ux = " << ux << " = " << s(i-1, j, k, nf*NVAR + UMX) / rho << std::endl;

      // i+1
      rho = s(iv+iv_dir, nf*NVAR + URHO);
      rhoinv = 1.0 / rho;
      AMREX_D_TERM(ux = s(iv+iv_dir, nf*NVAR + UMX)*rhoinv; ,
                    uy = s(iv+iv_dir, nf*NVAR + UMY)*rhoinv; ,
                    uz = s(iv+iv_dir, nf*NVAR + UMZ)*rhoinv;);
      ei = s(iv+iv_dir, nf*NVAR + UEDEN) * rhoinv - amrex::Real(0.5)*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));
      for (int n = 0; n < NUM_SPECIES; ++n) {
          massfrac[n] = s(iv+iv_dir, nf*NVAR + UFS+n) * rhoinv;
      }
      eos.EY2T(ei, massfrac, T);
      eos.RTY2P(rho, T, massfrac, dp[nf-1]);
      if (isnan(dp[nf-1])) amrex::Print() << i << ": p_{i+1} is nan!!! T = " << T << 
                                                                   ", ei = " << ei << 
                                                                  ", rho = " << rho <<
                                                                   ", mx = " << s(i+1, j, k, nf*NVAR + UMX) <<
                                                                   ", ux = " << ux << std::endl;
      dp[nf-1] = dp[nf-1] - p; // p_{i+1} - p{i-1}
      mean_dp += dp[nf-1];

      // // i
      // rho = s(i, j, k, nf*NVAR + URHO);
      // rhoinv = 1.0 / rho;
      // AMREX_D_TERM(amrex::Real ux = s(i, j, k, nf*NVAR + UMX)*rhoinv; ,
      //              amrex::Real uy = s(i, j, k, nf*NVAR + UMY)*rhoinv; ,
      //              amrex::Real uz = s(i, j, k, nf*NVAR + UMZ)*rhoinv;);
      // ei = s(i, j, k, nf*NVAR + UEDEN) * rhoinv - amrex::Real(0.5)*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));
      // for (int n = 0; n < NUM_SPECIES; ++n) {
      //     massfrac[n] = s(i, j, k, nf*NVAR + UFS+n) * rhoinv;
      // }
      // eos.EY2T(ei, massfrac, T);
      // eos.RTY2P(rho, T, massfrac, dp[nf-1]);
      // dp[nf-1] *= 2;            
      // mean_dp += 2*dp[nf-1];
    }
    mean_dp /= amrex::Real(NUM_FIELD);

    // u_i = u_i + 2.075*tke/eps*(mean_u_i - u_i)*dt + sqrt(2.1*eps)*dW - d(mean_p - p)/dx_i/rho*dt
    amrex::Real u_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
            u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
            
            u_sgs = 
            // u_sgs = 2.075 * tke/eps * (mean_u - u) * dt + sqrt(2.1*eps)*dW[nf-1][dir]
            // u_sgs = 0.005 * (mean_u - u) * dt
                      - 0.5 * (mean_dp - dp[nf-1]) / dx[dir] / s(i, j, k, nf*NVAR + URHO) * dt; //central difference for dpdx
// amrex::Print() << u_sgs << " " << tke << " " << eps << " " << dW[nf-1][dir] << std::endl;
            // amrex::Print() << u_sgs << " ";
            
            s(i, j, k, nf*NVAR + UMX+dir) += s(i, j, k, nf*NVAR + URHO) * u_sgs;

            // Enforce consistent e
            // s(i, j, k, nf*NVAR + UEDEN) += - 0.5 * s(i, j, k, nf*NVAR + URHO) * u * u 
            //                                + 0.5 * s(i, j, k, nf*NVAR + URHO) * (u+u_sgs) * (u+u_sgs);
        }
    }
    // amrex::Print() << "\n";
  });
}

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
spdf_iem_model (const amrex::Box& bx,
                amrex::Array4<amrex::Real>& s,
                const amrex::Real dt)
{   
    // amrex::ParallelFor(bx, 
    // [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    //     // Compute turbulent statistics
    //     amrex::Real k   = 1.0
    //     amrex::Real eps = 1.0;

    //     // u = u + k/eps*(mean_y - y)*dt
    //     for (int nf = 1; nf <= NUM_FIELD; ++nf) {
    //         for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    //             mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
    //             u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
                
    //             u += 2.075 * k/eps * (mean_u - u) * dt + sqrt(2.1*eps)*dW[nf-1][dir];
                
    //             s(i, j, k, nf*NVAR + UMX+dir) = s(i, j, k, URHO) * u;
    //         }
    //     }
    // });
}

#endif