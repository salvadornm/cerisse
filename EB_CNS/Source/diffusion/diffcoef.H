// #ifndef CNS_DIFFCOEF_H_
// #define CNS_DIFFCOEF_H_

// #include "index_macros.H"
// #include "parm.H"

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_diffcoef_eb(int i, int j, int k,
//                 amrex::Array4<amrex::Real const> const& q,
//                 amrex::Array4<amrex::EBCellFlag const> const& flag,
//                 amrex::Array4<amrex::Real> const& coefs,
//                 Parm const& parm) noexcept
// {
//     using amrex::Real;

//     if (q(i,j,k,QTEMP) < 0.) {
//         amrex::Abort("Negative temperature sent to cns_diffcoef_eb");
//     }

//     bool cov = flag(i,j,k).isCovered();
//     coefs(i,j,k,CMU) = cov ? -1.e10 : parm.C_S * std::sqrt(q(i,j,k,QTEMP)) * q(i,j,k,QTEMP) / (q(i,j,k,QTEMP)+parm.T_S);
//     coefs(i,j,k,CXI)  = cov ? -1.e10 : Real(0.0);
//     coefs(i,j,k,CLAM) = cov ? -1.e10 : coefs(i,j,k,CMU)*parm.cp/parm.Pr;
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_constcoef_eb(int i, int j, int k,
//                  amrex::Array4<amrex::EBCellFlag const> const& flag,
//                  amrex::Array4<amrex::Real> const& coefs,
//                  Parm const& parm) noexcept
// {
//     using amrex::Real;

//     bool cov = flag(i,j,k).isCovered();
//     coefs(i,j,k,CMU) = cov ? -1.e01 : parm.const_visc_mu;
//     coefs(i,j,k,CXI)  = cov ? -1.e01 : parm.const_visc_ki;
//     coefs(i,j,k,CLAM) = cov ? -1.e01 : parm.const_lambda;
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_diffcoef (int i, int j, int k,
//               amrex::Array4<amrex::Real const> const& q,
//               amrex::Array4<amrex::Real> const& coefs,
//               Parm const& parm) noexcept
// {
//     using amrex::Real;

//     if (q(i,j,k,QTEMP) < 0.) {
// //        std::cout << "Computing diff coeffs at " << IntVect(i,j,k) <<
// //             " with  temp = " << q(i,j,k,UTEMP)  << std::endl;
//         amrex::Abort("Negative temperature sent to cns_diffcoef_eb");
//     }

//     coefs(i,j,k,CMU) = parm.C_S * std::sqrt(q(i,j,k,QTEMP)) * q(i,j,k,QTEMP) / (q(i,j,k,QTEMP)+parm.T_S);
//     coefs(i,j,k,CXI)  = Real(0.0);
//     coefs(i,j,k,CLAM) = coefs(i,j,k,CMU)*parm.cp/parm.Pr;
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void
// cns_constcoef (int i, int j, int k,
//               amrex::Array4<amrex::Real> const& coefs,
//               Parm const& parm) noexcept
// {
//     coefs(i,j,k,CMU) = parm.const_visc_mu;
//     coefs(i,j,k,CXI)  = parm.const_visc_ki;
//     coefs(i,j,k,CLAM) = parm.const_lambda;
// }

// #endif