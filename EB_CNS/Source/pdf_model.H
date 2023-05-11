#ifndef CNS_PDF_MODEL_H_
#define CNS_PDF_MODEL_H_

#include "CNS.H"
#include "index_macros.H"

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
vpdf_psgs_model (const amrex::Box& bx,
                 amrex::Array4<amrex::Real>& s,
                 const amrex::Real dt,
                 const amrex::Real(*)[] /*dW[AMREX_SPACEDIM][NUM_FIELD]*/,
                 const amrex::Real dx[AMREX_SPACEDIM])
{
  amrex::ParallelFor(bx, 
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    // Compute pressure (store pressure of each field in p[])
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);
      amrex::Real mean_dp = 0.0;
      amrex::Real dp[NUM_FIELD];
      amrex::Real mean_dup = 0.0;
      amrex::Real dup[NUM_FIELD];
      amrex::Real rho, rhoinv, ei, T, p, massfrac[NUM_SPECIES], up;
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
          massfrac[n] = s(iv-iv_dir, nf*NVAR + UFS+n)*rhoinv;
        }
        eos.REY2T(rho, ei, massfrac, T); 
        eos.RTY2P(rho, T, massfrac, p);
        if (isnan(p)) amrex::Print() << i << ": p_{i-1} is nan!!! T = " << T << ", ei = " << ei << 
                            ", rho = " << rho << " = 1/" << 1 / rho << ", ux = " << ux << " = " << s(i-1, j, k, nf*NVAR + UMX) / rho << std::endl;
        up = s(iv-iv_dir, nf*NVAR + UMX+dir)*rhoinv*p;

        // i+1
        rho = s(iv+iv_dir, nf*NVAR + URHO);
        rhoinv = 1.0 / rho;
        AMREX_D_TERM(ux = s(iv+iv_dir, nf*NVAR + UMX)*rhoinv; ,
                    uy = s(iv+iv_dir, nf*NVAR + UMY)*rhoinv; ,
                    uz = s(iv+iv_dir, nf*NVAR + UMZ)*rhoinv;);
        ei = s(iv+iv_dir, nf*NVAR + UEDEN) * rhoinv - amrex::Real(0.5)*(AMREX_D_TERM(ux*ux, + uy*uy, + uz*uz));
        for (int n = 0; n < NUM_SPECIES; ++n) {
            massfrac[n] = s(iv+iv_dir, nf*NVAR + UFS+n)*rhoinv;
        }
        eos.REY2T(rho, ei, massfrac, T);
        eos.RTY2P(rho, T, massfrac, dp[nf-1]);
        if (isnan(dp[nf-1])) amrex::Print() << i << ": p_{i+1} is nan!!! T = " << T << ", ei = " << ei << 
                                                    ", rho = " << rho << ", mx = " << s(i+1, j, k, nf*NVAR + UMX) << ", ux = " << ux << std::endl;
        dup[nf-1] = s(iv+iv_dir, nf*NVAR + UMX+dir)*rhoinv*dp[nf-1] - up;
        mean_dup += dup[nf-1];

        dp[nf-1] = dp[nf-1] - p; // p_{i+1} - p{i-1}
        mean_dp += dp[nf-1];
      }
      mean_dp /= amrex::Real(NUM_FIELD);
      mean_dup /= amrex::Real(NUM_FIELD);

      amrex::Real mean_u, u, u_sgs;
      for (int nf = 1; nf <= NUM_FIELD; ++nf) {        
        mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
        u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
        
        // u_sgs = -(mean_dp - dp[nf-1]) / 0.5*dx[dir] * dt; //central difference for dpdx
        // this should be mom_sgs instead
        
        // s(i, j, k, nf*NVAR + UMX+dir) += u_sgs;
        s(i, j, k, nf*NVAR + UMX+dir) -= (mean_dp - dp[nf-1]) / 0.5*dx[dir] * dt;

        // Enforce consistent e
        // s(i, j, k, nf*NVAR + UEDEN) += 0.5 * s(i, j, k, nf*NVAR + URHO) * ((u+u_sgs)*(u+u_sgs) - u*u);      
        // s(i, j, k, nf*NVAR + UEDEN) += s(i, j, k, nf*NVAR + URHO) * u * u_sgs; //vspdf21
        s(i, j, k, nf*NVAR + UEDEN) -= (mean_dup - dup[nf-1]) / 0.5*dx[dir] * dt; //vspdf22
      }
    }
  });  
}

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
vpdf_langevin_model (const amrex::Box& bx,
                     amrex::Array4<amrex::Real>& s,
                     const amrex::Real dt,
                     const amrex::Real dW[AMREX_SPACEDIM][NUM_FIELD],
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
    tke = 0.5 * tke / amrex::Real(NUM_FIELD-1); // k = 0.5*(u^2 + v^2 + w^2)
    // if (tke < 0.0) { 
    //   amrex::Print() << tke; 
    //   amrex::Error("TKE < 0!!!"); 
    // }
    amrex::Real deltainv = 1.0 / std::pow(AMREX_D_TERM(dx[0], * dx[1], * dx[2]), 1.0/amrex::Real(AMREX_SPACEDIM));
    amrex::Real eps = 0.53 * std::pow(tke, 1.5) * deltainv; //or 1.05? 0.53?
    // eps = amrex::max(eps, 1e-10);
    tke = amrex::max(tke, 1e-10); //prevent divide by 0

    // rhou_i = rhou_i + rho_i*Gij*(mean_u_i - u_i)*dt + rho_i*sqrt(C0*eps)*dW
    //                                                                    ^ need to /mean_rho?
    amrex::Real u_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        mean_u = s(i, j, k, UMX+dir) / s(i, j, k, URHO);
        u = s(i, j, k, nf*NVAR + UMX+dir) / s(i, j, k, nf*NVAR + URHO);
        // Gij = 2.075*eps/tke, C0 = 2.1
        u_sgs = 2.075 * eps/tke * (mean_u - u) * dt + sqrt(2.1*eps)*dW[dir][nf-1];

        s(i, j, k, nf*NVAR + UMX+dir) += s(i, j, k, nf*NVAR + URHO) * u_sgs;
      }
    }
    
    // Test dW
    // s(i, j, k, UFA+11) = dW[0][0];
  });
}

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
spdf_iem_model (const amrex::Box& bx,
                amrex::Array4<amrex::Real>& s,
                const amrex::Real dt,
                const amrex::Real(*)[] /*dW[AMREX_SPACEDIM][NUM_FIELD]*/,
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
    tke = 0.5 * tke / amrex::Real(NUM_FIELD-1); // k = 0.5*(u'^2 + v'^2 + w'^2)
    amrex::Real deltainv = 1.0 / std::pow(AMREX_D_TERM(dx[0], * dx[1], * dx[2]), 1.0/amrex::Real(AMREX_SPACEDIM));
    amrex::Real eps = 0.53 * std::pow(tke, 1.5) * deltainv; //or 1.05?
    // if (isnan(eps)) amrex::Print() << i << "," << j << ": eps is nan!!!";
    // eps = amrex::max(eps, 1e-10);
    tke = amrex::max(tke, 1e-10); //prevent divide by 0

    // rhoQ_i = rhoQ_i + rho_i*0.5*C_Q*eps/k*(mean_y - y)*dt
    amrex::Real mean_Q, Q, Q_sgs;
    for (int nf = 1; nf <= NUM_FIELD; ++nf) {
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        // Species (Q = Y), C_Y = 2
        mean_Q = s(i, j, k, UFS+ns) / s(i, j, k, URHO);
        Q = s(i, j, k, nf*NVAR + UFS+ns) / s(i, j, k, nf*NVAR + URHO);        
        Q_sgs = eps/tke * (mean_Q - Q) * dt;        
        s(i, j, k, nf*NVAR + UFS+ns) += s(i, j, k, nf*NVAR + URHO) * Q_sgs;
      }

      // Energy (Q = E), C_E = 2
      mean_Q = s(i, j, k, UEDEN) / s(i, j, k, URHO);
      Q = s(i, j, k, nf*NVAR + UEDEN) / s(i, j, k, nf*NVAR + URHO);        
      Q_sgs = eps/tke * (mean_Q - Q) * dt;        
      s(i, j, k, nf*NVAR + UEDEN) += s(i, j, k, nf*NVAR + URHO) * Q_sgs;
    }
    });
}

#endif