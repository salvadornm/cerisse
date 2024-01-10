#ifndef CNS_NSCBC_H_
#define CNS_NSCBC_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include <PelePhysics.H>
#include <EOS.H>

#include "index_macros.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void one_side_derivative(
  amrex::IntVect iv, int idir, int isign,
  amrex::Real deltainv, amrex::Real& dp, amrex::Real& dun, amrex::Real& dut,
  amrex::Real& dutt, amrex::Real& drho, const amrex::Array4<const amrex::Real>& q)
{
  const auto iv_dir = amrex::IntVect::TheDimensionVector(idir) * isign;

  AMREX_D_TERM(int QUN = QU + idir;, int QUT = QU + (idir + 1) % (AMREX_SPACEDIM - 1);
               , int QUTT = QU + (idir + 2) % (AMREX_SPACEDIM - 1);)

  dp = amrex::Real(isign) *
       (-1.5 * q(iv, QPRES) + 2.0 * q(iv + iv_dir, QPRES) -
        0.5 * q(iv + 2 * iv_dir, QPRES)) *
       deltainv;
  AMREX_D_TERM(dun = amrex::Real(isign) *
                     (-1.5 * q(iv, QUN) + 2.0 * q(iv + iv_dir, QUN) -
                      0.5 * q(iv + 2 * iv_dir, QUN)) *
                     deltainv;
               , dut = amrex::Real(isign) *
                       (-1.5 * q(iv, QUT) + 2.0 * q(iv + iv_dir, QUT) -
                        0.5 * q(iv + 2 * iv_dir, QUT)) *
                       deltainv;
               , dutt = amrex::Real(isign) *
                        (-1.5 * q(iv, QUTT) + 2.0 * q(iv + iv_dir, QUTT) -
                         0.5 * q(iv + 2 * iv_dir, QUTT)) *
                        deltainv;)
  drho = amrex::Real(isign) *
         (-1.5 * q(iv, QRHO) + 2.0 * q(iv + iv_dir, QRHO) -
          0.5 * q(iv + 2 * iv_dir, QRHO)) *
         deltainv;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void central_derivative(
  amrex::IntVect iv, int idir, amrex::Real deltainv, amrex::Real& dp,
  amrex::Real& dun, amrex::Real& dut, amrex::Real& dutt, amrex::Real& drho,
  const amrex::Array4<const amrex::Real>& q)
{
  const auto iv_dir = amrex::IntVect::TheDimensionVector(idir);

  AMREX_D_TERM(int QUN = QU + idir;, int QUT = QU + (idir + 1) % (AMREX_SPACEDIM - 1);
               , int QUTT = QU + (idir + 2) % (AMREX_SPACEDIM - 1);)

  dp = 0.5 * (q(iv + iv_dir, QPRES) - q(iv - iv_dir, QPRES)) * deltainv;
  AMREX_D_TERM(dun = 0.5 * (q(iv + iv_dir, QUN) - q(iv - iv_dir, QUN)) * deltainv;
               , dut = 0.5 * (q(iv + iv_dir, QUT) - q(iv - iv_dir, QUT)) * deltainv;
               , dutt =
                   0.5 * (q(iv + iv_dir, QUTT) - q(iv - iv_dir, QUTT)) * deltainv;)
  drho = 0.5 * (q(iv + iv_dir, QRHO) - q(iv - iv_dir, QRHO)) * deltainv;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
extrap_prim(amrex::IntVect iv, amrex::IntVect iv_dir, int index, int isign,
            amrex::Real delta, amrex::Real dqdx, const amrex::Array4<amrex::Real>& q)
{  
  amrex::Real dq = dqdx * delta;
  // q(iv - iv_dir, index) = q(iv + iv_dir, index) - 2 * isign * dq; // original
  q(iv - iv_dir, index) = q(iv, index) - isign * dq;
  q(iv - 2 * iv_dir, index) = -2 * q(iv + iv_dir, index) - 3 * q(iv, index) +
                              6 * q(iv - iv_dir, index) + 6 * isign * dq;
  q(iv - 3 * iv_dir, index) = 3 * q(iv + iv_dir, index) + 10 * q(iv, index) -
                              18 * q(iv - iv_dir, index) +
                              6 * q(iv - 2 * iv_dir, index) - 12 * isign * dq;
  q(iv - 4 * iv_dir, index) =
    -2 * q(iv + iv_dir, index) - 13 * q(iv, index) + 24 * q(iv - iv_dir, index) -
    12 * q(iv - 2 * iv_dir, index) + 4 * q(iv - 3 * iv_dir, index) + 12 * isign * dq;
  q(iv - 5 * iv_dir, index) =
    5 * q(iv + iv_dir, index) + 38.5 * q(iv, index) - 75 * q(iv - iv_dir, index) +
    50 * q(iv - 2 * iv_dir, index) - 25 * q(iv - 3 * iv_dir, index) +
    7.5 * q(iv - 4 * iv_dir, index) - 30 * isign * dq;
  q(iv - 6 * iv_dir, index) =
    -6 * q(iv + iv_dir, index) - 60.9 * q(iv, index) + 126 * q(iv - iv_dir, index) -
    105 * q(iv - 2 * iv_dir, index) + 70 * q(iv - 3 * iv_dir, index) -
    31.5 * q(iv - 4 * iv_dir, index) + 8.4 * q(iv - 5 * iv_dir, index) +
    42 * isign * dq;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
nscbc_c2prim(int i, int j, int k, int ns, amrex::Array4<amrex::Real const> const& u,
              amrex::Array4<amrex::Real> const& q) noexcept
{
  amrex::Real rho = u(i, j, k, ns + URHO);
  amrex::Real rhoinv = amrex::Real(1.0) / rho;
  AMREX_D_TERM(amrex::Real ux = u(i, j, k, ns + UMX) * rhoinv;
               , amrex::Real uy = u(i, j, k, ns + UMY) * rhoinv;
               , amrex::Real uz = u(i, j, k, ns + UMZ) * rhoinv;);
  amrex::Real ei = u(i, j, k, ns + UEDEN) * rhoinv -
                   amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz));

  amrex::Real Y[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y[n] = u(i, j, k, ns + UFS + n) * rhoinv;
    q(i, j, k, QFS + n) = Y[n];
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real T = 0, p, cs = 0, gamma, wbar; 
  eos.REY2T(rho, ei, Y, T);
  eos.RTY2P(rho, T, Y, p);
  eos.RTY2Cs(rho, T, Y, cs);
  eos.RTY2G(rho, T, Y, gamma);
  eos.Y2WBAR(Y, wbar);

  q(i, j, k, QRHO) = rho;
  AMREX_D_TERM(q(i, j, k, QU) = ux;, q(i, j, k, QV) = uy;, q(i, j, k, QW) = uz;);
  q(i, j, k, QPRES) = p;
  q(i, j, k, QC) = cs;
  q(i, j, k, QG) = gamma;
  q(i, j, k, QEINT) = pele::physics::Constants::RU / wbar; // use the place for ei to store R, sorry for the confusion 
  q(i, j, k, QTEMP) = T;

  // // Monitor problem cell
  std::ostream& ss = std::cout;

  if ((T <= 0.0) || amrex::isnan(cs) || amrex::isnan(ux) || 
      (rho <= 0.0) || (p <= 0.0) || amrex::isnan(p) || amrex::isnan(rho)) {
      ss << "n@" << i << "," << j << "," << k << ": ns=" << ns << ": rho=" << rho
      << " T=" << T << " p=" << p << " cs=" << cs << " gamma=" << gamma
      << " u=" << ux << " rhoE=" << u(i,j,k,ns+UEDEN) << " ei=" << ei << '\n';
      for (int n = 0; n < NUM_SPECIES; n++) {
          ss << " rhoY" << n << "=" << u(i,j,k,ns+UFS+n) << " ";
      }
      ss << '\n';
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
nscbc_p2cons(int i, int j, int k, int ns, amrex::Array4<amrex::Real> const& u,
           amrex::Array4<amrex::Real const> const& q) noexcept
{
  amrex::Real rho = q(i, j, k, QRHO);
  AMREX_D_TERM(amrex::Real ux = q(i, j, k, QU);, amrex::Real uy = q(i, j, k, QV);
               , amrex::Real uz = q(i, j, k, QW);)
  amrex::Real p = q(i, j, k, QPRES);

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real ei, y[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) { y[n] = q(i, j, k, QFS + n); }
  eos.RYP2E(rho, y, p, ei);

  u(i, j, k, ns + URHO) = rho;
  AMREX_D_TERM(u(i, j, k, ns + UMX) = rho * ux;, u(i, j, k, ns + UMY) = rho * uy;
               , u(i, j, k, ns + UMZ) = rho * uz;);
  u(i, j, k, ns + UEDEN) =
    rho * (ei + amrex::Real(0.5) * (AMREX_D_TERM(ux * ux, +uy * uy, +uz * uz)));
  for (int n = 0; n < NUM_SPECIES; ++n) { u(i, j, k, ns + UFS + n) = rho * y[n]; }
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_transverse_terms(
//   int i, int j, int k, int idir, amrex::Real* T, amrex::Real dpdx, amrex::Real dudx,
//   amrex::Real dvdx, amrex::Real dwdx, amrex::Real drhodx, amrex::Real dpdy,
//   amrex::Real dudy, amrex::Real dvdy, amrex::Real dwdy, amrex::Real drhody,
//   amrex::Real dpdz, amrex::Real dudz, amrex::Real dvdz, amrex::Real dwdz,
//   amrex::Real drhodz, const amrex::Array4<amrex::Real>& q)
// {
//   const amrex::Real inv_rho = 1.0 / q(i, j, k, QRHO);
//   amrex::Real T1, T2, T3, T4, T5;
//   if (idir == 0) {
//     T1 = (q(i, j, k, QV) * (dpdy - q(i, j, k, QRHO) * q(i, j, k, QC) * dudy)) +
//          (q(i, j, k, QW) * (dpdz - q(i, j, k, QRHO) * q(i, j, k, QC) * dudz)) +
//          (q(i, j, k, QG) * q(i, j, k, QPRES) * (dvdy + dwdz));
//     T2 = (q(i, j, k, QV) * ((q(i, j, k, QC) * q(i, j, k, QC) * drhody) - dpdy)) +
//          (q(i, j, k, QW) * ((q(i, j, k, QC) * q(i, j, k, QC) * drhodz) - dpdz));
//     T3 = q(i, j, k, QV) * dvdy + q(i, j, k, QW) * dvdz + dpdy * inv_rho;
//     T4 = q(i, j, k, QV) * dwdy + q(i, j, k, QW) * dwdz + dpdz * inv_rho;
//     T5 = (q(i, j, k, QV) * (dpdy + q(i, j, k, QRHO) * q(i, j, k, QC) * dudy)) +
//          (q(i, j, k, QW) * (dpdz + q(i, j, k, QRHO) * q(i, j, k, QC) * dudz)) +
//          (q(i, j, k, QG) * q(i, j, k, QPRES) * (dvdy + dwdz));
//   } else if (idir == 1) {
//     T1 = (q(i, j, k, QU) * (dpdx - q(i, j, k, QRHO) * q(i, j, k, QC) * dvdx)) +
//          (q(i, j, k, QW) * (dpdz - q(i, j, k, QRHO) * q(i, j, k, QC) * dvdz)) +
//          (q(i, j, k, QG) * q(i, j, k, QPRES) * (dudx + dwdz));
//     T2 = q(i, j, k, QU) * dudx + q(i, j, k, QW) * dudz + dpdx * inv_rho;
//     T3 = (q(i, j, k, QU) * ((q(i, j, k, QC) * q(i, j, k, QC) * drhodx) - dpdx)) +
//          (q(i, j, k, QW) * ((q(i, j, k, QC) * q(i, j, k, QC) * drhodz) - dpdz));

//     T4 = q(i, j, k, QU) * dwdx + q(i, j, k, QW) * dwdz + dpdz * inv_rho;
//     T5 = (q(i, j, k, QU) * (dpdx + q(i, j, k, QRHO) * q(i, j, k, QC) * dvdx)) +
//          (q(i, j, k, QW) * (dpdz + q(i, j, k, QRHO) * q(i, j, k, QC) * dvdz)) +
//          (q(i, j, k, QG) * q(i, j, k, QPRES) * (dudx + dwdz));
//   } else if (idir == 2) {
//     T1 = (q(i, j, k, QU) * (dpdx - q(i, j, k, QRHO) * q(i, j, k, QC) * dwdx)) +
//          (q(i, j, k, QV) * (dpdy - q(i, j, k, QRHO) * q(i, j, k, QC) * dwdy)) +
//          (q(i, j, k, QG) * q(i, j, k, QPRES) * (dudx + dvdy));
//     T2 = q(i, j, k, QU) * dudx + q(i, j, k, QV) * dudy + dpdx * inv_rho;
//     T3 = q(i, j, k, QU) * dvdx + q(i, j, k, QV) * dvdy + dpdy * inv_rho;
//     T4 = (q(i, j, k, QU) * ((q(i, j, k, QC) * q(i, j, k, QC) * drhodx) - dpdx)) +
//          (q(i, j, k, QV) * ((q(i, j, k, QC) * q(i, j, k, QC) * drhody) - dpdy));
//     T5 = (q(i, j, k, QU) * (dpdx + q(i, j, k, QRHO) * q(i, j, k, QC) * dwdx)) +
//          (q(i, j, k, QV) * (dpdy + q(i, j, k, QRHO) * q(i, j, k, QC) * dwdy)) +
//          (q(i, j, k, QG) * q(i, j, k, QPRES) * (dudx + dvdy));
//   } else {
//     amrex::Abort("Problem of idir in compute_transverse_terms()");
//   }
//   T[0] = T1;
//   T[1] = T2;
//   T[2] = T3;
//   T[3] = T4;
//   T[4] = T5;
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void tangential_derivative(
//   int i, int j, int k, int idir, amrex::Real delta, amrex::Real& dp, amrex::Real& du,
//   amrex::Real& dv, amrex::Real& dw, amrex::Real& drho,
//   const amrex::Array4<amrex::Real>& q)
// {
//   // Warning, idir means the tangential direction, this is different from 2D
//   // (sorry)
//   if (idir == 0) {
//     // 2nd order Central
//     dp = (q(i + 1, j, k, QPRES) - q(i - 1, j, k, QPRES)) / (2.0 * delta);
//     du = (q(i + 1, j, k, QU) - q(i - 1, j, k, QU)) / (2.0 * delta);
//     dv = (q(i + 1, j, k, QV) - q(i - 1, j, k, QV)) / (2.0 * delta);
//     dw = (q(i + 1, j, k, QW) - q(i - 1, j, k, QW)) / (2.0 * delta);
//     drho = (q(i + 1, j, k, QRHO) - q(i - 1, j, k, QRHO)) / (2.0 * delta);
//   } else if (idir == 1) {
//     // 2nd order Central
//     dp = (q(i, j + 1, k, QPRES) - q(i, j - 1, k, QPRES)) / (2.0 * delta);
//     du = (q(i, j + 1, k, QU) - q(i, j - 1, k, QU)) / (2.0 * delta);
//     dv = (q(i, j + 1, k, QV) - q(i, j - 1, k, QV)) / (2.0 * delta);
//     dw = (q(i, j + 1, k, QW) - q(i, j - 1, k, QW)) / (2.0 * delta);
//     drho = (q(i, j + 1, k, QRHO) - q(i, j - 1, k, QRHO)) / (2.0 * delta);
//   } else if (idir == 2) {
//     // 2nd order Central
//     dp = (q(i, j, k + 1, QPRES) - q(i, j, k - 1, QPRES)) / (2.0 * delta);
//     du = (q(i, j, k + 1, QU) - q(i, j, k - 1, QU)) / (2.0 * delta);
//     dv = (q(i, j, k + 1, QV) - q(i, j, k - 1, QV)) / (2.0 * delta);
//     dw = (q(i, j, k + 1, QW) - q(i, j, k - 1, QW)) / (2.0 * delta);
//     drho = (q(i, j, k + 1, QRHO) - q(i, j, k - 1, QRHO)) / (2.0 * delta);
//   } else {
//     amrex::Abort("Problem of idir in impose_NSCBC_3d:tangential_derivative");
//   }
// }

// AMREX_GPU_DEVICE
// AMREX_FORCE_INLINE
// void update_ghost_cells(int i, int j, int k, int bc_type, int idir, int isign,
//                         amrex::Real delta, const int* domlo, const int* domhi,
//                         const amrex::Real* L, const amrex::Array4<amrex::Real>& uin,
//                         const amrex::Array4<amrex::Real>& q,
//                         const amrex::Array4<amrex::Real>& qaux)
// {
//   int idx_gc1, idx_gc2, idx_gc3, idx_gc4, idx_int1, idx_int2, idx_int3;
//   int idx_start, idx_end, local_index;
//   amrex::Real drho, du, dv, dw, dp, wall_sign;
//   const amrex::Real small = 1.e-8;
//   AMREX_ASSERT((idir == 0) || (idir == 1) || (idir == 2));
//   AMREX_ASSERT((isign == 1) || (isign == -1));
//   amrex::Real L1 = L[0];
//   amrex::Real L2 = L[1];
//   amrex::Real L3 = L[2];
//   amrex::Real L4 = L[3];
//   amrex::Real L5 = L[4];

//   // Compute new spatial derivative
//   if (idir == 0) {
//     local_index = i;
//     drho = (L2 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
//     du = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
//     dv = L3;
//     dw = L4;
//     dp = 0.5 * (L1 + L5);
//   } else if (idir == 1) {
//     local_index = j;
//     drho = (L3 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
//     du = L2;
//     dv = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
//     dw = L4;
//     dp = 0.5 * (L1 + L5);
//   } else if (idir == 2) {
//     local_index = k;
//     drho = (L4 + 0.5 * (L1 + L5)) / (qaux(i, j, k, QC) * qaux(i, j, k, QC));
//     du = L2;
//     dv = L3;
//     dw = (L5 - L1) / (2.0 * qaux(i, j, k, QC) * q(i, j, k, QRHO));
//     dp = 0.5 * (L1 + L5);
//   }

//   if (isign == 1) {
//     idx_gc1 = local_index - 1;
//     idx_gc2 = local_index - 2;
//     idx_gc3 = local_index - 3;
//     idx_gc4 = local_index - 4;
//     idx_int1 = local_index + 1;
//     idx_int2 = local_index + 2;
//     idx_int3 = local_index + 3;
//     idx_start = domlo[idir] - 1;
//     idx_end = domlo[idir] - 4;
//   } else if (isign == -1) {
//     idx_gc1 = local_index + 1;
//     idx_gc2 = local_index + 2;
//     idx_gc3 = local_index + 3;
//     idx_gc4 = local_index + 4;
//     idx_int1 = local_index - 1;
//     idx_int2 = local_index - 2;
//     idx_int3 = local_index - 3;
//     idx_start = domhi[idir] + 1;
//     idx_end = domhi[idir] + 4;
//   }

//   if (idir == 0) {
//     // Update ghost cells
//     // 2nd order
//     q(idx_gc1, j, k, QU) = q(idx_int1, j, k, QU) - 2.0 * delta * du * isign;
//     q(idx_gc1, j, k, QV) = q(idx_int1, j, k, QV) - 2.0 * delta * dv * isign;
//     q(idx_gc1, j, k, QW) = q(idx_int1, j, k, QW) - 2.0 * delta * dw * isign;
//     q(idx_gc1, j, k, QRHO) = q(idx_int1, j, k, QRHO) - 2.0 * delta * drho * isign;
//     q(idx_gc1, j, k, QPRES) = q(idx_int1, j, k, QPRES) - 2.0 * delta * dp * isign;

//     q(idx_gc2, j, k, QU) = -2.0 * q(idx_int1, j, k, QU) - 3.0 * q(i, j, k, QU) +
//                            6.0 * q(idx_gc1, j, k, QU) + 6.0 * delta * du * isign;
//     q(idx_gc2, j, k, QV) = -2.0 * q(idx_int1, j, k, QV) - 3.0 * q(i, j, k, QV) +
//                            6.0 * q(idx_gc1, j, k, QV) + 6.0 * delta * dv * isign;
//     q(idx_gc2, j, k, QW) = -2.0 * q(idx_int1, j, k, QW) - 3.0 * q(i, j, k, QW) +
//                            6.0 * q(idx_gc1, j, k, QW) + 6.0 * delta * dw * isign;
//     q(idx_gc2, j, k, QRHO) = -2.0 * q(idx_int1, j, k, QRHO) -
//                              3.0 * q(i, j, k, QRHO) + 6.0 * q(idx_gc1, j, k, QRHO) +
//                              6.0 * delta * drho * isign;
//     q(idx_gc2, j, k, QPRES) =
//       -2.0 * q(idx_int1, j, k, QPRES) - 3.0 * q(i, j, k, QPRES) +
//       6.0 * q(idx_gc1, j, k, QPRES) + 6.0 * delta * dp * isign;

//     q(idx_gc3, j, k, QU) = 3.0 * q(idx_int1, j, k, QU) + 10.0 * q(i, j, k, QU) -
//                            18.0 * q(idx_gc1, j, k, QU) + 6.0 * q(idx_gc2, j, k, QU) -
//                            12.0 * delta * du * isign;
//     q(idx_gc3, j, k, QV) = 3.0 * q(idx_int1, j, k, QV) + 10.0 * q(i, j, k, QV) -
//                            18.0 * q(idx_gc1, j, k, QV) + 6.0 * q(idx_gc2, j, k, QV) -
//                            12.0 * delta * dv * isign;
//     q(idx_gc3, j, k, QW) = 3.0 * q(idx_int1, j, k, QW) + 10.0 * q(i, j, k, QW) -
//                            18.0 * q(idx_gc1, j, k, QW) + 6.0 * q(idx_gc2, j, k, QW) -
//                            12.0 * delta * dw * isign;
//     q(idx_gc3, j, k, QRHO) =
//       3.0 * q(idx_int1, j, k, QRHO) + 10.0 * q(i, j, k, QRHO) -
//       18.0 * q(idx_gc1, j, k, QRHO) + 6.0 * q(idx_gc2, j, k, QRHO) -
//       12.0 * delta * drho * isign;
//     q(idx_gc3, j, k, QPRES) =
//       3.0 * q(idx_int1, j, k, QPRES) + 10.0 * q(i, j, k, QPRES) -
//       18.0 * q(idx_gc1, j, k, QPRES) + 6.0 * q(idx_gc2, j, k, QPRES) -
//       12.0 * delta * dp * isign;

//     q(idx_gc4, j, k, QU) = -2.0 * q(idx_int1, j, k, QU) - 13.0 * q(i, j, k, QU) +
//                            24.0 * q(idx_gc1, j, k, QU) -
//                            12.0 * q(idx_gc2, j, k, QU) + 4.0 * q(idx_gc3, j, k, QU) +
//                            12.0 * delta * du * isign;
//     q(idx_gc4, j, k, QV) = -2.0 * q(idx_int1, j, k, QV) - 13.0 * q(i, j, k, QV) +
//                            24.0 * q(idx_gc1, j, k, QV) -
//                            12.0 * q(idx_gc2, j, k, QV) + 4.0 * q(idx_gc3, j, k, QV) +
//                            12.0 * delta * dv * isign;
//     q(idx_gc4, j, k, QW) = -2.0 * q(idx_int1, j, k, QW) - 13.0 * q(i, j, k, QW) +
//                            24.0 * q(idx_gc1, j, k, QW) -
//                            12.0 * q(idx_gc2, j, k, QW) + 4.0 * q(idx_gc3, j, k, QW) +
//                            12.0 * delta * dw * isign;
//     q(idx_gc4, j, k, QRHO) =
//       -2.0 * q(idx_int1, j, k, QRHO) - 13.0 * q(i, j, k, QRHO) +
//       24.0 * q(idx_gc1, j, k, QRHO) - 12.0 * q(idx_gc2, j, k, QRHO) +
//       4.0 * q(idx_gc3, j, k, QRHO) + 12.0 * delta * drho * isign;
//     q(idx_gc4, j, k, QPRES) =
//       -2.0 * q(idx_int1, j, k, QPRES) - 13.0 * q(i, j, k, QPRES) +
//       24.0 * q(idx_gc1, j, k, QPRES) - 12.0 * q(idx_gc2, j, k, QPRES) +
//       4.0 * q(idx_gc3, j, k, QPRES) + 12.0 * delta * dp * isign;

//     // If BC is SlipWall or NoSlipWall
//     if ((bc_type == 4) || (bc_type == 5)) {
//       if (bc_type == 5) {
//         wall_sign = -1.0;
//       } else if (bc_type == 4) {
//         wall_sign = 1.0;
//       }

//       q(idx_gc1, j, k, QU) = -q(i, j, k, QU);
//       q(idx_gc2, j, k, QU) = -q(idx_int1, j, k, QU);
//       q(idx_gc3, j, k, QU) = -q(idx_int2, j, k, QU);
//       q(idx_gc4, j, k, QU) = -q(idx_int3, j, k, QU);

//       q(idx_gc1, j, k, QV) = wall_sign * q(i, j, k, QV);
//       q(idx_gc2, j, k, QV) = wall_sign * q(idx_int1, j, k, QV);
//       q(idx_gc3, j, k, QV) = wall_sign * q(idx_int2, j, k, QV);
//       q(idx_gc4, j, k, QV) = wall_sign * q(idx_int3, j, k, QV);

//       q(idx_gc1, j, k, QW) = wall_sign * q(i, j, k, QW);
//       q(idx_gc2, j, k, QW) = wall_sign * q(idx_int1, j, k, QW);
//       q(idx_gc3, j, k, QW) = wall_sign * q(idx_int2, j, k, QW);
//       q(idx_gc4, j, k, QW) = wall_sign * q(idx_int3, j, k, QW);

//       q(idx_gc1, j, k, QRHO) = q(i, j, k, QRHO);
//       q(idx_gc2, j, k, QRHO) = q(idx_int1, j, k, QRHO);
//       q(idx_gc3, j, k, QRHO) = q(idx_int2, j, k, QRHO);
//       q(idx_gc4, j, k, QRHO) = q(idx_int3, j, k, QRHO);

//       q(idx_gc1, j, k, QPRES) = q(i, j, k, QPRES);
//       q(idx_gc2, j, k, QPRES) = q(idx_int1, j, k, QPRES);
//       q(idx_gc3, j, k, QPRES) = q(idx_int2, j, k, QPRES);
//       q(idx_gc4, j, k, QPRES) = q(idx_int3, j, k, QPRES);
//     }

//   } else if (idir == 1) {
//     // Update ghost cells
//     // 2nd order
//     q(i, idx_gc1, k, QU) = q(i, idx_int1, k, QU) - 2.0 * delta * du * isign;
//     q(i, idx_gc1, k, QV) = q(i, idx_int1, k, QV) - 2.0 * delta * dv * isign;
//     q(i, idx_gc1, k, QW) = q(i, idx_int1, k, QW) - 2.0 * delta * dw * isign;
//     q(i, idx_gc1, k, QRHO) = q(i, idx_int1, k, QRHO) - 2.0 * delta * drho * isign;
//     q(i, idx_gc1, k, QPRES) = q(i, idx_int1, k, QPRES) - 2.0 * delta * dp * isign;

//     q(i, idx_gc2, k, QU) = -2.0 * q(i, idx_int1, k, QU) - 3.0 * q(i, j, k, QU) +
//                            6.0 * q(i, idx_gc1, k, QU) + 6.0 * delta * du * isign;
//     q(i, idx_gc2, k, QV) = -2.0 * q(i, idx_int1, k, QV) - 3.0 * q(i, j, k, QV) +
//                            6.0 * q(i, idx_gc1, k, QV) + 6.0 * delta * dv * isign;
//     q(i, idx_gc2, k, QW) = -2.0 * q(i, idx_int1, k, QW) - 3.0 * q(i, j, k, QW) +
//                            6.0 * q(i, idx_gc1, k, QW) + 6.0 * delta * dw * isign;
//     q(i, idx_gc2, k, QRHO) = -2.0 * q(i, idx_int1, k, QRHO) -
//                              3.0 * q(i, j, k, QRHO) + 6.0 * q(i, idx_gc1, k, QRHO) +
//                              6.0 * delta * drho * isign;
//     q(i, idx_gc2, k, QPRES) =
//       -2.0 * q(i, idx_int1, k, QPRES) - 3.0 * q(i, j, k, QPRES) +
//       6.0 * q(i, idx_gc1, k, QPRES) + 6.0 * delta * dp * isign;

//     q(i, idx_gc3, k, QU) = 3.0 * q(i, idx_int1, k, QU) + 10.0 * q(i, j, k, QU) -
//                            18.0 * q(i, idx_gc1, k, QU) + 6.0 * q(i, idx_gc2, k, QU) -
//                            12.0 * delta * du * isign;
//     q(i, idx_gc3, k, QV) = 3.0 * q(i, idx_int1, k, QV) + 10.0 * q(i, j, k, QV) -
//                            18.0 * q(i, idx_gc1, k, QV) + 6.0 * q(i, idx_gc2, k, QV) -
//                            12.0 * delta * dv * isign;
//     q(i, idx_gc3, k, QW) = 3.0 * q(i, idx_int1, k, QW) + 10.0 * q(i, j, k, QW) -
//                            18.0 * q(i, idx_gc1, k, QW) + 6.0 * q(i, idx_gc2, k, QW) -
//                            12.0 * delta * dw * isign;
//     q(i, idx_gc3, k, QRHO) =
//       3.0 * q(i, idx_int1, k, QRHO) + 10.0 * q(i, j, k, QRHO) -
//       18.0 * q(i, idx_gc1, k, QRHO) + 6.0 * q(i, idx_gc2, k, QRHO) -
//       12.0 * delta * drho * isign;
//     q(i, idx_gc3, k, QPRES) =
//       3.0 * q(i, idx_int1, k, QPRES) + 10.0 * q(i, j, k, QPRES) -
//       18.0 * q(i, idx_gc1, k, QPRES) + 6.0 * q(i, idx_gc2, k, QPRES) -
//       12.0 * delta * dp * isign;

//     q(i, idx_gc4, k, QU) = -2.0 * q(i, idx_int1, k, QU) - 13.0 * q(i, j, k, QU) +
//                            24.0 * q(i, idx_gc1, k, QU) -
//                            12.0 * q(i, idx_gc2, k, QU) + 4.0 * q(i, idx_gc3, k, QU) +
//                            12.0 * delta * du * isign;
//     q(i, idx_gc4, k, QV) = -2.0 * q(i, idx_int1, k, QV) - 13.0 * q(i, j, k, QV) +
//                            24.0 * q(i, idx_gc1, k, QV) -
//                            12.0 * q(i, idx_gc2, k, QV) + 4.0 * q(i, idx_gc3, k, QV) +
//                            12.0 * delta * dv * isign;
//     q(i, idx_gc4, k, QW) = -2.0 * q(i, idx_int1, k, QW) - 13.0 * q(i, j, k, QW) +
//                            24.0 * q(i, idx_gc1, k, QW) -
//                            12.0 * q(i, idx_gc2, k, QW) + 4.0 * q(i, idx_gc3, k, QW) +
//                            12.0 * delta * dw * isign;
//     q(i, idx_gc4, k, QRHO) =
//       -2.0 * q(i, idx_int1, k, QRHO) - 13.0 * q(i, j, k, QRHO) +
//       24.0 * q(i, idx_gc1, k, QRHO) - 12.0 * q(i, idx_gc2, k, QRHO) +
//       4.0 * q(i, idx_gc3, k, QRHO) + 12.0 * delta * drho * isign;
//     q(i, idx_gc4, k, QPRES) =
//       -2.0 * q(i, idx_int1, k, QPRES) - 13.0 * q(i, j, k, QPRES) +
//       24.0 * q(i, idx_gc1, k, QPRES) - 12.0 * q(i, idx_gc2, k, QPRES) +
//       4.0 * q(i, idx_gc3, k, QPRES) + 12.0 * delta * dp * isign;

//     if ((bc_type == 4) || (bc_type == 5)) {
//       if (bc_type == 5) {
//         wall_sign = -1.0;
//       } else if (bc_type == 4) {
//         wall_sign = 1.0;
//       }

//       q(i, idx_gc1, k, QU) = wall_sign * q(i, j, k, QU);
//       q(i, idx_gc2, k, QU) = wall_sign * q(i, idx_int1, k, QU);
//       q(i, idx_gc3, k, QU) = wall_sign * q(i, idx_int2, k, QU);
//       q(i, idx_gc4, k, QU) = wall_sign * q(i, idx_int3, k, QU);

//       q(i, idx_gc1, k, QV) = -q(i, j, k, QV);
//       q(i, idx_gc2, k, QV) = -q(i, idx_int1, k, QV);
//       q(i, idx_gc3, k, QV) = -q(i, idx_int2, k, QV);
//       q(i, idx_gc4, k, QV) = -q(i, idx_int3, k, QV);

//       q(i, idx_gc1, k, QW) = wall_sign * q(i, j, k, QW);
//       q(i, idx_gc2, k, QW) = wall_sign * q(i, idx_int1, k, QW);
//       q(i, idx_gc3, k, QW) = wall_sign * q(i, idx_int2, k, QW);
//       q(i, idx_gc4, k, QW) = wall_sign * q(i, idx_int3, k, QW);

//       q(i, idx_gc1, k, QRHO) = q(i, j, k, QRHO);
//       q(i, idx_gc2, k, QRHO) = q(i, idx_int1, k, QRHO);
//       q(i, idx_gc3, k, QRHO) = q(i, idx_int2, k, QRHO);
//       q(i, idx_gc4, k, QRHO) = q(i, idx_int3, k, QRHO);

//       q(i, idx_gc1, k, QPRES) = q(i, j, k, QPRES);
//       q(i, idx_gc2, k, QPRES) = q(i, idx_int1, k, QPRES);
//       q(i, idx_gc3, k, QPRES) = q(i, idx_int2, k, QPRES);
//       q(i, idx_gc4, k, QPRES) = q(i, idx_int3, k, QPRES);
//     }
//   } else if (idir == 2) {
//     // Update ghost cells
//     // 2nd order
//     q(i, j, idx_gc1, QU) = q(i, j, idx_int1, QU) - 2.0 * delta * du * isign;
//     q(i, j, idx_gc1, QV) = q(i, j, idx_int1, QV) - 2.0 * delta * dv * isign;
//     q(i, j, idx_gc1, QW) = q(i, j, idx_int1, QW) - 2.0 * delta * dw * isign;
//     q(i, j, idx_gc1, QRHO) = q(i, j, idx_int1, QRHO) - 2.0 * delta * drho * isign;
//     q(i, j, idx_gc1, QPRES) = q(i, j, idx_int1, QPRES) - 2.0 * delta * dp * isign;

//     q(i, j, idx_gc2, QU) = -2.0 * q(i, j, idx_int1, QU) - 3.0 * q(i, j, k, QU) +
//                            6.0 * q(i, j, idx_gc1, QU) + 6.0 * delta * du * isign;
//     q(i, j, idx_gc2, QV) = -2.0 * q(i, j, idx_int1, QV) - 3.0 * q(i, j, k, QV) +
//                            6.0 * q(i, j, idx_gc1, QV) + 6.0 * delta * dv * isign;
//     q(i, j, idx_gc2, QW) = -2.0 * q(i, j, idx_int1, QW) - 3.0 * q(i, j, k, QW) +
//                            6.0 * q(i, j, idx_gc1, QW) + 6.0 * delta * dw * isign;
//     q(i, j, idx_gc2, QRHO) = -2.0 * q(i, j, idx_int1, QRHO) -
//                              3.0 * q(i, j, k, QRHO) + 6.0 * q(i, j, idx_gc1, QRHO) +
//                              6.0 * delta * drho * isign;
//     q(i, j, idx_gc2, QPRES) =
//       -2.0 * q(i, j, idx_int1, QPRES) - 3.0 * q(i, j, k, QPRES) +
//       6.0 * q(i, j, idx_gc1, QPRES) + 6.0 * delta * dp * isign;

//     q(i, j, idx_gc3, QU) = 3.0 * q(i, j, idx_int1, QU) + 10.0 * q(i, j, k, QU) -
//                            18.0 * q(i, j, idx_gc1, QU) + 6.0 * q(i, j, idx_gc2, QU) -
//                            12.0 * delta * du * isign;
//     q(i, j, idx_gc3, QV) = 3.0 * q(i, j, idx_int1, QV) + 10.0 * q(i, j, k, QV) -
//                            18.0 * q(i, j, idx_gc1, QV) + 6.0 * q(i, j, idx_gc2, QV) -
//                            12.0 * delta * dv * isign;
//     q(i, j, idx_gc3, QW) = 3.0 * q(i, j, idx_int1, QW) + 10.0 * q(i, j, k, QW) -
//                            18.0 * q(i, j, idx_gc1, QW) + 6.0 * q(i, j, idx_gc2, QW) -
//                            12.0 * delta * dw * isign;
//     q(i, j, idx_gc3, QRHO) =
//       3.0 * q(i, j, idx_int1, QRHO) + 10.0 * q(i, j, k, QRHO) -
//       18.0 * q(i, j, idx_gc1, QRHO) + 6.0 * q(i, j, idx_gc2, QRHO) -
//       12.0 * delta * drho * isign;
//     q(i, j, idx_gc3, QPRES) =
//       3.0 * q(i, j, idx_int1, QPRES) + 10.0 * q(i, j, k, QPRES) -
//       18.0 * q(i, j, idx_gc1, QPRES) + 6.0 * q(i, j, idx_gc2, QPRES) -
//       12.0 * delta * dp * isign;

//     q(i, j, idx_gc4, QU) = -2.0 * q(i, j, idx_int1, QU) - 13.0 * q(i, j, k, QU) +
//                            24.0 * q(i, j, idx_gc1, QU) -
//                            12.0 * q(i, j, idx_gc2, QU) + 4.0 * q(i, j, idx_gc3, QU) +
//                            12.0 * delta * du * isign;
//     q(i, j, idx_gc4, QV) = -2.0 * q(i, j, idx_int1, QV) - 13.0 * q(i, j, k, QV) +
//                            24.0 * q(i, j, idx_gc1, QV) -
//                            12.0 * q(i, j, idx_gc2, QV) + 4.0 * q(i, j, idx_gc3, QV) +
//                            12.0 * delta * dv * isign;
//     q(i, j, idx_gc4, QW) = -2.0 * q(i, j, idx_int1, QW) - 13.0 * q(i, j, k, QW) +
//                            24.0 * q(i, j, idx_gc1, QW) -
//                            12.0 * q(i, j, idx_gc2, QW) + 4.0 * q(i, j, idx_gc3, QW) +
//                            12.0 * delta * dw * isign;
//     q(i, j, idx_gc4, QRHO) =
//       -2.0 * q(i, j, idx_int1, QRHO) - 13.0 * q(i, j, k, QRHO) +
//       24.0 * q(i, j, idx_gc1, QRHO) - 12.0 * q(i, j, idx_gc2, QRHO) +
//       4.0 * q(i, j, idx_gc3, QRHO) + 12.0 * delta * drho * isign;
//     q(i, j, idx_gc4, QPRES) =
//       -2.0 * q(i, j, idx_int1, QPRES) - 13.0 * q(i, j, k, QPRES) +
//       24.0 * q(i, j, idx_gc1, QPRES) - 12.0 * q(i, j, idx_gc2, QPRES) +
//       4.0 * q(i, j, idx_gc3, QPRES) + 12.0 * delta * dp * isign;

//     if ((bc_type == 4) || (bc_type == 5)) {
//       if (bc_type == 5) {
//         wall_sign = -1.0;
//       } else if (bc_type == 4) {
//         wall_sign = 1.0;
//       }

//       q(i, j, idx_gc1, QU) = wall_sign * q(i, j, k, QU);
//       q(i, j, idx_gc2, QU) = wall_sign * q(i, j, idx_int1, QU);
//       q(i, j, idx_gc3, QU) = wall_sign * q(i, j, idx_int2, QU);
//       q(i, j, idx_gc4, QU) = wall_sign * q(i, j, idx_int3, QU);

//       q(i, j, idx_gc1, QV) = wall_sign * q(i, j, k, QV);
//       q(i, j, idx_gc2, QV) = wall_sign * q(i, j, idx_int1, QV);
//       q(i, j, idx_gc3, QV) = wall_sign * q(i, j, idx_int2, QV);
//       q(i, j, idx_gc4, QV) = wall_sign * q(i, j, idx_int3, QV);

//       q(i, j, idx_gc1, QW) = -q(i, j, k, QW);
//       q(i, j, idx_gc2, QW) = -q(i, j, idx_int1, QW);
//       q(i, j, idx_gc3, QW) = -q(i, j, idx_int2, QW);
//       q(i, j, idx_gc4, QW) = -q(i, j, idx_int3, QW);

//       q(i, j, idx_gc1, QRHO) = q(i, j, k, QRHO);
//       q(i, j, idx_gc2, QRHO) = q(i, j, idx_int1, QRHO);
//       q(i, j, idx_gc3, QRHO) = q(i, j, idx_int2, QRHO);
//       q(i, j, idx_gc4, QRHO) = q(i, j, idx_int3, QRHO);

//       q(i, j, idx_gc1, QPRES) = q(i, j, k, QPRES);
//       q(i, j, idx_gc2, QPRES) = q(i, j, idx_int1, QPRES);
//       q(i, j, idx_gc3, QPRES) = q(i, j, idx_int2, QPRES);
//       q(i, j, idx_gc4, QPRES) = q(i, j, idx_int3, QPRES);
//     }
//   }
//   auto eos = pele::physics::PhysicsType::eos();
//   // Recompute missing values thanks to EOS
//   int ii = i;
//   int jj = j;
//   int kk = k;
//   for (int hop = idx_start; hop < idx_end; hop = hop - isign) {
//     if (idir == 0) {
//       ii = hop;
//     } else if (idir == 1) {
//       jj = hop;
//     } else if (idir == 2) {
//       kk = hop;
//     }
//     amrex::Real eos_state_p = q(ii, jj, kk, QPRES);
//     amrex::Real eos_state_rho = q(ii, jj, kk, QRHO);
//     amrex::Real eos_state_massfrac[NUM_SPECIES];
//     for (int n = 0; n < NUM_SPECIES; n++) {
//       eos_state_massfrac[n] = q(ii, jj, kk, QFS + n);
//     }
//     amrex::Real eos_state_aux[NUM_AUX];
//     for (int n = 0; n < NUM_AUX; n++) { eos_state_aux[n] = q(ii, jj, kk, QFX + n); }

//     amrex::Real eos_state_T, eos_state_e, eos_state_gamma, eos_state_cs;
//     amrex::Real eos_state_dpdr_e, eos_state_dpde, wbar;
//     eos.Y2WBAR(eos_state_massfrac, wbar);
//     eos.RYP2E(eos_state_rho, eos_state_massfrac, eos_state_p, eos_state_e);
//     eos.EY2T(eos_state_e, eos_state_massfrac, eos_state_T);
//     eos.TY2G(eos_state_T, eos_state_massfrac, eos_state_gamma);
//     eos.RTY2Cs(eos_state_rho, eos_state_T, eos_state_massfrac, eos_state_cs);
//     eos.RTY2dpde_dpdre(eos_state_rho, eos_state_T, eos_state_massfrac,
//                        eos_state_dpde, eos_state_dpdr_e);
//     q(ii, jj, kk, QTEMP) = eos_state_T;
//     q(ii, jj, kk, QREINT) = eos_state_e * q(ii, jj, kk, QRHO);
//     q(ii, jj, kk, QGAME) = q(ii, jj, kk, QPRES) / q(ii, jj, kk, QREINT) + 1.0;

//     qaux(ii, jj, kk, QDPDR) = eos_state_dpdr_e;
//     qaux(ii, jj, kk, QDPDE) = eos_state_dpde;
//     qaux(ii, jj, kk, QGAMC) = eos_state_gamma;
//     qaux(ii, jj, kk, QC) = eos_state_cs;
//     qaux(ii, jj, kk, QCSML) = amrex::max<amrex::Real>(
//       constants::small_num(), constants::small_num() * eos_state_cs);
//     qaux(ii, jj, kk, QRSPEC) = pele::physics::Constants::RU / wbar;

//     // Here the update of the conservative variables uin seems to have only an
//     // impact on the application of artificial viscosity difmag.
//     uin(ii, jj, kk, URHO) = eos_state_rho;
//     uin(ii, jj, kk, UMX) = q(ii, jj, kk, QU) * eos_state_rho;
//     uin(ii, jj, kk, UMY) = q(ii, jj, kk, QV) * eos_state_rho;
//     uin(ii, jj, kk, UMZ) = q(ii, jj, kk, QW) * eos_state_rho;
//     uin(ii, jj, kk, UEINT) = eos_state_rho * eos_state_e;
//     uin(ii, jj, kk, UEDEN) =
//       eos_state_rho * (eos_state_e + 0.5 * (q(ii, jj, kk, QU) * q(ii, jj, kk, QU) +
//                                             q(ii, jj, kk, QV) * q(ii, jj, kk, QV) +
//                                             q(ii, jj, kk, QW) * q(ii, jj, kk, QW)));
//     uin(ii, jj, kk, UTEMP) = eos_state_T;
//     for (int n = 0; n < NUM_SPECIES; n++) {
//       uin(ii, jj, kk, UFS + n) = eos_state_rho * eos_state_massfrac[n];
//     }
//   }
// }

// AMREX_GPU_DEVICE
// AMREX_FORCE_INLINE
// void compute_waves(int i, int j, int k, int idir, int isign, int bc_type,
//                    const amrex::Real* problen, const amrex::Real bc_params[6],
//                    const amrex::Real bc_target[5], const amrex::Real T[5],
//                    amrex::Real* L, amrex::Real dp, amrex::Real du, amrex::Real dv,
//                    amrex::Real dw, amrex::Real drho,
//                    const amrex::Array4<amrex::Real>& q,
//                    const amrex::Array4<amrex::Real>& qaux)
// {
//   // Note: for convenience, all waves are called L because we have just 1
//   // direction here Waves M and N will be employed for corners

//   AMREX_ASSERT((idir == 0) || (idir == 1) || (idir == 2));

//   amrex::Real mach_local =
//     std::sqrt(q(i, j, k, QU) * q(i, j, k, QU) + q(i, j, k, QV) * q(i, j, k, QV) +
//               q(i, j, k, QW) * q(i, j, k, QW)) /
//     qaux(i, j, k, QC);

//   // Recasting targets values and numerical parameters
//   amrex::Real TARGET_VX = bc_target[0];
//   amrex::Real TARGET_VY = bc_target[1];
//   amrex::Real TARGET_VZ = bc_target[2];
//   amrex::Real TARGET_TEMPERATURE = bc_target[3];
//   amrex::Real TARGET_PRESSURE = bc_target[4];

//   amrex::Real relax_T = bc_params[0];
//   amrex::Real relax_U = bc_params[1];
//   amrex::Real relax_V = bc_params[2];
//   amrex::Real relax_W = bc_params[3];
//   // Here we have the abilities to set beta=local Mach number
//   // it may works better for outflow BCs
//   amrex::Real beta = bc_params[4];
//   if (bc_params[4] < 0.0) { beta = mach_local; }
//   amrex::Real sigma_out = bc_params[5];

//   amrex::Real L1, L2, L3, L4, L5, Kout;
//   const amrex::Real T1 = T[0];
//   const amrex::Real T2 = T[1];
//   const amrex::Real T3 = T[2];
//   const amrex::Real T4 = T[3];
//   const amrex::Real T5 = T[4];
//   // Computing known numerical LODI waves
//   if (idir == 0) {
//     // Numerical LODI waves along X
//     L1 = (q(i, j, k, QU) - qaux(i, j, k, QC)) *
//          (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * du);
//     L2 = q(i, j, k, QU) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
//     L3 = q(i, j, k, QU) * dv;
//     L4 = q(i, j, k, QU) * dw;
//     L5 = (q(i, j, k, QU) + qaux(i, j, k, QC)) *
//          (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * du);
//   } else if (idir == 1) {
//     // Numerical LODI waves along Y
//     L1 = (q(i, j, k, QV) - qaux(i, j, k, QC)) *
//          (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dv);
//     L2 = q(i, j, k, QV) * du;
//     L3 = q(i, j, k, QV) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
//     L4 = q(i, j, k, QV) * dw;
//     L5 = (q(i, j, k, QV) + qaux(i, j, k, QC)) *
//          (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dv);
//   } else if (idir == 2) {
//     // Numerical LODI waves along Z
//     L1 = (q(i, j, k, QW) - qaux(i, j, k, QC)) *
//          (dp - (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dw);
//     L2 = q(i, j, k, QW) * du;
//     L3 = q(i, j, k, QW) * dv;
//     L4 = q(i, j, k, QW) * (((qaux(i, j, k, QC) * qaux(i, j, k, QC)) * drho) - dp);
//     L5 = (q(i, j, k, QW) + qaux(i, j, k, QC)) *
//          (dp + (q(i, j, k, QRHO) * qaux(i, j, k, QC)) * dw);
//   }

//   // Computing missing LODI waves from BC model
//   if (bc_type == 7) {
//     // Inflow
//     if (idir == 0) {
//       if (isign == 1) {
//         L5 = relax_U *
//                ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
//                 (1.0 - mach_local * mach_local) / problen[idir]) *
//                (q(i, j, k, QU) - TARGET_VX) -
//              ((1.0 - beta) * T5);
//       } else if (isign == -1) {
//         L1 = relax_U *
//                ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
//                 (1.0 - mach_local * mach_local) / problen[idir]) *
//                (q(i, j, k, QU) - TARGET_VX) -
//              ((1.0 - beta) * T1);
//       }
//       L2 = relax_T *
//              (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
//               problen[idir]) *
//              (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
//            ((1.0 - beta) * T2);
//       L3 = relax_V * (qaux(i, j, k, QC) / problen[idir]) *
//              (q(i, j, k, QV) - TARGET_VY) -
//            ((1.0 - beta) * T3);
//       L4 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
//              (q(i, j, k, QW) - TARGET_VZ) -
//            ((1.0 - beta) * T4);
//     } else if (idir == 1) {
//       if (isign == 1) {
//         L5 = relax_V *
//                ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
//                 (1.0 - mach_local * mach_local) / problen[idir]) *
//                (q(i, j, k, QV) - TARGET_VY) -
//              ((1.0 - beta) * T5);
//       } else if (isign == -1) {
//         L1 = relax_V *
//                ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
//                 (1.0 - mach_local * mach_local) / problen[idir]) *
//                (q(i, j, k, QV) - TARGET_VY) -
//              ((1.0 - beta) * T1);
//       }
//       L2 = relax_U * (qaux(i, j, k, QC) / problen[idir]) *
//              (q(i, j, k, QU) - TARGET_VX) -
//            ((1.0 - beta) * T2);
//       L3 = relax_T *
//              (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
//               problen[idir]) *
//              (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
//            ((1.0 - beta) * T3);
//       L4 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
//              (q(i, j, k, QW) - TARGET_VZ) -
//            ((1.0 - beta) * T4);
//     } else if (idir == 2) {
//       if (isign == 1) {
//         L5 = relax_W *
//                ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
//                 (1.0 - mach_local * mach_local) / problen[idir]) *
//                (q(i, j, k, QW) - TARGET_VZ) -
//              ((1.0 - beta) * T5);
//       } else if (isign == -1) {
//         L1 = relax_W *
//                ((q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QC)) *
//                 (1.0 - mach_local * mach_local) / problen[idir]) *
//                (q(i, j, k, QW) - TARGET_VZ) -
//              ((1.0 - beta) * T1);
//       }
//       L2 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
//              (q(i, j, k, QU) - TARGET_VX) -
//            ((1.0 - beta) * T2);
//       L3 = relax_W * (qaux(i, j, k, QC) / problen[idir]) *
//              (q(i, j, k, QV) - TARGET_VY) -
//            ((1.0 - beta) * T3);
//       L4 = relax_T *
//              (q(i, j, k, QRHO) * qaux(i, j, k, QC) * qaux(i, j, k, QRSPEC) /
//               problen[idir]) *
//              (q(i, j, k, QTEMP) - TARGET_TEMPERATURE) -
//            ((1.0 - beta) * T4);
//     } else {
//       amrex::Abort("Error:: Wait, is this the fourth dimension?");
//     }
//   } else if ((bc_type == 4) || (bc_type == 5)) {
//     // Slipwall and NoSlipWall
//     // Values long Y will be computed by mirror functions below
//     // but we set waves values to 0 to avoid undefined variables
//     L1 = 0.0;
//     L2 = 0.0;
//     L3 = 0.0;
//     L4 = 0.0;
//     L5 = 0.0;
//   } else if (bc_type == 8) {
//     // Outflow
//     Kout = sigma_out * (1.0 - (mach_local * mach_local)) *
//            (qaux(i, j, k, QC) / problen[idir]);
//     if (isign == 1) {
//       L5 = (Kout * (q(i, j, k, QPRES) - TARGET_PRESSURE)) - ((1.0 - beta) * T5);
//     } else if (isign == -1) {
//       L1 = (Kout * (q(i, j, k, QPRES) - TARGET_PRESSURE)) - ((1.0 - beta) * T1);
//     }
//   } else {
//     amrex::Abort("Error:: This BC is not yet implemented for x dir in "
//                  "characteristic form");
//   }

//   // Shaping the waves to be at the good dimension
//   if (idir == 0) {
//     L1 = L1 / (q(i, j, k, QU) - qaux(i, j, k, QC));
//     L5 = L5 / (q(i, j, k, QU) + qaux(i, j, k, QC));
//     if (q(i, j, k, QU) == 0.0) {
//       L2 = 0.0;
//       L3 = 0.0;
//       L4 = 0.0;
//     } else {
//       L2 = L2 / q(i, j, k, QU);
//       L3 = L3 / q(i, j, k, QU);
//       L4 = L4 / q(i, j, k, QU);
//     }
//   } else if (idir == 1) {
//     L1 = L1 / (q(i, j, k, QV) - qaux(i, j, k, QC));
//     L5 = L5 / (q(i, j, k, QV) + qaux(i, j, k, QC));
//     if (q(i, j, k, QV) == 0.0) {
//       L2 = 0.0;
//       L3 = 0.0;
//       L4 = 0.0;
//     } else {
//       L2 = L2 / q(i, j, k, QV);
//       L3 = L3 / q(i, j, k, QV);
//       L4 = L4 / q(i, j, k, QV);
//     }
//   } else if (idir == 2) {
//     L1 = L1 / (q(i, j, k, QW) - qaux(i, j, k, QC));
//     L5 = L5 / (q(i, j, k, QW) + qaux(i, j, k, QC));
//     if (q(i, j, k, QW) == 0.0) {
//       L2 = 0.0;
//       L3 = 0.0;
//       L4 = 0.0;
//     } else {
//       L2 = L2 / q(i, j, k, QW);
//       L3 = L3 / q(i, j, k, QW);
//       L4 = L4 / q(i, j, k, QW);
//     }
//   }
//   L[0] = L1;
//   L[1] = L2;
//   L[2] = L3;
//   L[3] = L4;
//   L[4] = L5;
// }

#endif