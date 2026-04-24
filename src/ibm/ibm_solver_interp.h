#ifndef IBM_SOLVER_INTERP_H_
#define IBM_SOLVER_INTERP_H_
// ============================================================================
// ibm_solver_interp.h — Interpolation, extrapolation, and coordinate-transform helpers
//
// This file contains in-class definitions for private member
// functions of ibm_solver_t.  It is included inside the class body in ibm_solver.h
// and must NOT be included independently.
//
// Contents:
//   1. valid_mirror              — Count fluid cells in an interpolation stencil
//   2. search_optimal_image_point— Best first image-point placement (multiple attempts)
//   3. search_image_point        — Subsequent image-point placement (single step)
//   4. computeIPweights          — Bi-/tri-linear interpolation weights
//   5. interpolateIMs            — Interpolate primitives at image points
//   6. extrapolate               — Lagrangian extrapolation to ghost point
//   7. global2local / local2global — Velocity coordinate transforms
//   8. check_interpolation_stencil — Validate stencil containment in box
// ============================================================================
// ============================================================================
// 1. valid_mirror
// ============================================================================
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
int valid_mirror(
    int i, int j, int k,
    const Array4<const uint8_t>& ibMarkers)
{
    int fluid_count = 0;
#if (AMREX_SPACEDIM == 2)
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
          int ii = i + di;
          int jj = j + dj;
          if (ibMarkers(ii, jj, 0, 0) == 0) {
              ++fluid_count;
          }
      }
    }
    return fluid_count;
#else
    for (int di = 0; di <= 1; ++di) {
      for (int dj = 0; dj <= 1; ++dj) {
        for (int dk = 0; dk <= 1; ++dk) {
          int ii = i + di;
          int jj = j + dj;
          int kk = k + dk;
          if (ibMarkers(ii, jj, kk, 0) == 0) {
            ++fluid_count;
          }
        }
      }
    }
    return fluid_count;
#endif
}
// ============================================================================
// 2. search_optimal_image_point
// ============================================================================
template <int eorder_t, int iorder_t, typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void search_optimal_image_point(
    const Point& cp_start,
    const LocalFrame& localframe,
    int lev,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_lev,
    Real di,
    const Box& bxg,
    const Array4<uint8_t const>& ibMarkers,
    const IPDATA& ipData,
    int f_idx,
    Array2D<Real, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
    Array2D< int, 0, eorder_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
    Array1D<Real, 0, eorder_t - 1>& disIM,
    Array1D< int, 0, eorder_t - 1>& imp_ninterp)
{
    int best_fluid = 0;
    Array1D<Real, 0, AMREX_SPACEDIM - 1> candi_xyz;
    Array1D< int, 0, AMREX_SPACEDIM - 1> candi_ijk;
    constexpr int N_ATTEMPTS = (GP_OR_SURF ? N_ATTEMPTS_GP : N_ATTEMPTS_SURF);
    // Image-point placement distance multipliers (multiples of di = alpha * cell_diagonal).
    // The first attempt places the image point at 1.0*di; if the stencil there is
    // insufficient (too many solid neighbours), we fall back to 1.5*di, 2.0*di, etc.
    // Local constexpr arrays are required for device-side address validity.
    constexpr Real IMP_FACTOR_GP_local[]   = {1.0, 1.5, 2.0};
    constexpr Real IMP_FACTOR_SURF_local[] = {1.0, 1.5, 2.0, 2.5, 3.0};
    const Real* IMP_FACTOR = (GP_OR_SURF ? IMP_FACTOR_GP_local : IMP_FACTOR_SURF_local);
    
    for (int attempt = 0; attempt < N_ATTEMPTS; ++attempt) {
    
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          candi_xyz(d) = cp_start[d] + IMP_FACTOR[attempt] * di * localframe.normal[d];
          candi_ijk(d) = int(amrex::Math::floor(
              (candi_xyz(d) - prob_lo[d]) / dx_lev[d] - 0.5
          ));
      }
    
#if (AMREX_SPACEDIM == 2)
      bool in_box = check_interpolation_stencil<IPDATA>(candi_ijk(0), candi_ijk(1), 0, 
                                                bxg, lev,
                                                ipData, f_idx,
                                                (attempt == 0) ? CheckMode::Abort : CheckMode::Silent);
      int n_fluid = (in_box) ? valid_mirror(candi_ijk(0), candi_ijk(1), 0, ibMarkers) : -1;
#else
      bool in_box = check_interpolation_stencil<IPDATA>(candi_ijk(0), candi_ijk(1), candi_ijk(2),
                                                bxg, lev,
                                                ipData, f_idx,
                                                (attempt == 0) ? CheckMode::Abort : CheckMode::Silent);
      int n_fluid = (in_box) ? valid_mirror(candi_ijk(0), candi_ijk(1), candi_ijk(2), ibMarkers) : -1;
#endif
      if (n_fluid > best_fluid) {
        best_fluid = n_fluid;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            imp_xyz(0, d) = candi_xyz(d);
            imp_ijk(0, d) = candi_ijk(d);
        }
        imp_ninterp(0) = n_fluid; 
        disIM(0) = IMP_FACTOR[attempt] * di;
      }
      constexpr int IDEAL_NINTERP = ipow(iorder_t + 1, AMREX_SPACEDIM);
      if (best_fluid == IDEAL_NINTERP) {
        break;
      }
    } // end loop on attempt
    constexpr int INTERP_THRESHOLD = (GP_OR_SURF ? INTERP_THRESHOLD_GP : INTERP_THRESHOLD_SURF);
    if (best_fluid < INTERP_THRESHOLD) {
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("Not enough valid interpolation points for first image point! "
                            "lev=%d f_idx=%d best_fluid=%d threshold=%d\n",
                            lev, f_idx, best_fluid, INTERP_THRESHOLD);
        AMREX_ASSERT(false);
#else
        int current_geom = -1;
        int current_elem = ipData.elemIdx[f_idx];
        const char* point_label;
        Real p_x = 0.0, p_y = 0.0, p_z = 0.0;
        if constexpr (GP_OR_SURF == 1) {
            current_geom = ipData.geomIdx[f_idx];
            point_label = "Ghost Point";
            
            int gp_i = ipData.gp_ijk[f_idx](0);
            int gp_j = ipData.gp_ijk[f_idx](1);
            p_x = prob_lo[0] + (0.5_rt + gp_i) * dx_lev[0];
            p_y = prob_lo[1] + (0.5_rt + gp_j) * dx_lev[1];
#if (AMREX_SPACEDIM == 3)
            int gp_k = ipData.gp_ijk[f_idx](2);
            p_z = prob_lo[2] + (0.5_rt + gp_k) * dx_lev[2];
#endif
        } else {
            point_label = "Surface Point";
            
            p_x = cp_start[0];
            p_y = cp_start[1];
#if (AMREX_SPACEDIM == 3)
            p_z = cp_start[2];
#endif
        }
        
#if (AMREX_SPACEDIM == 3)
        std::printf("Not enough valid interpolation points found for the first image point!\n"
                    "  Level: %d\n"
                    "  %s: (%f, %f, %f)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Best Fluid Points Found: %d (Threshold: %d)\n",
                    lev, point_label, p_x, p_y, p_z, 
                    current_geom, current_elem,
                    best_fluid, INTERP_THRESHOLD);
#else
        std::printf("Not enough valid interpolation points found for the first image point!\n"
                    "  Level: %d\n"
                    "  %s: (%f, %f)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n"
                    "  Best Fluid Points Found: %d (Threshold: %d)\n",
                    lev, point_label, p_x, p_y, 
                    current_geom, current_elem,
                    best_fluid, INTERP_THRESHOLD);
#endif
        std::fflush(stdout);
#endif // AMREX_DEVICE_COMPILE
    }
}
// ============================================================================
// 3. search_image_point
// ============================================================================
template <int order_t, int iorder_t, typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void search_image_point(
    int jj,
    const Point& cp_start,
    const LocalFrame& localframe,
    int lev,
    const GpuArray<Real, AMREX_SPACEDIM>& prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_lev,
    Real di,
    const Box& bxg,
    const Array4<uint8_t const>& ibMarkers,
    const IPDATA& ipData,
    int f_idx,
    Array2D<Real, 0, order_t - 1, 0, AMREX_SPACEDIM - 1>& imp_xyz,
    Array2D< int, 0, order_t - 1, 0, AMREX_SPACEDIM - 1>& imp_ijk,
    Array1D<Real, 0, order_t - 1>& disIM,
    Array1D< int, 0, order_t - 1>& imp_ninterp)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        imp_xyz(jj, d) = cp_start[d] + di * localframe.normal[d];
        imp_ijk(jj, d) = int(amrex::Math::floor(
            (imp_xyz(jj, d) - prob_lo[d]) / dx_lev[d] - 0.5
        ));
    }
#if (AMREX_SPACEDIM == 2)
    bool in_box = check_interpolation_stencil<IPDATA>(imp_ijk(jj, 0), imp_ijk(jj, 1), 0, 
                                              bxg, lev,
                                              ipData, f_idx, 
                                              CheckMode::Silent);
    int fluid = (in_box) ? valid_mirror(imp_ijk(jj, 0), imp_ijk(jj, 1), 0, ibMarkers) : -1;
#else
    bool in_box = check_interpolation_stencil<IPDATA>(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2),
                                              bxg, lev,
                                              ipData, f_idx,
                                              CheckMode::Silent);   
    int fluid = (in_box) ? valid_mirror(imp_ijk(jj, 0), imp_ijk(jj, 1), imp_ijk(jj, 2), ibMarkers) : -1;
#endif
    disIM(jj) = (jj > 0) ? disIM(jj - 1) + di : di;
    imp_ninterp(jj) = fluid;
}
// ============================================================================
// 4. computeIPweights
// ============================================================================
template <int eorder_t, int iorder_t, typename IPDATA, int N_InterP = ipow(iorder_t + 1, AMREX_SPACEDIM), int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
void computeIPweights(
    Array2D<Real,0,eorder_t-1,0,N_InterP-1>&                     weights,
    Array3D< int,0,eorder_t-1,0,N_InterP-1,0,AMREX_SPACEDIM-1>&  ip_ijk,
    Array2D<Real,0,eorder_t-1,0,AMREX_SPACEDIM-1>&               imp_xyz,
    Array2D< int,0,eorder_t-1,0,AMREX_SPACEDIM-1>&               imp_ijk,
    Array1D< int,0,eorder_t-1>&                                  imp_ninterp,
    const GpuArray<Real, AMREX_SPACEDIM>&                        prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM>&                        dxyz,
    const Array4<uint8_t const>&                                 ibFab)
{
    constexpr int INTERP_THRESHOLD = (GP_OR_SURF ? INTERP_THRESHOLD_GP : INTERP_THRESHOLD_SURF);
    for (int iim = 0; iim < eorder_t; ++iim) {
      if (imp_ninterp(iim) < INTERP_THRESHOLD) {
        for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) = Real(0.0);
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = -99;
          }
        }
        continue;
      }
      int base_ijk[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {base_ijk[d] = imp_ijk(iim, d);}
      Real frac[AMREX_SPACEDIM];
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real lo = prob_lo[d] + Real(base_ijk[d] + 0.5_rt) * dxyz[d];
        frac[d] = (imp_xyz(iim, d) - lo) / dxyz[d];
      }
      int  sumfluid   = 0;
      Real sumweights = Real(0.0);
      for (int corner = 0; corner < N_InterP; ++corner) {
        int  ijk[AMREX_SPACEDIM];
        Real w = Real(1.0);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int bit = (corner >> d) & 1;
            ijk[d] = base_ijk[d] + bit;
            const Real fd = frac[d];
            w *= (bit ? fd : (Real(1.0) - fd));
        }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ip_ijk(iim, corner, d) = ijk[d];
        }
        int ii = ijk[0];
        int jj = ijk[1];
#if (AMREX_SPACEDIM == 3)
        int kk = ijk[2];
#else
        int kk = 0;
#endif
        int fluid = !ibFab(ii, jj, kk, 0);
        weights(iim, corner) = w * Real(fluid);
        sumfluid   += fluid;
        sumweights += weights(iim, corner);
      }
      AMREX_ASSERT_WITH_MESSAGE(
          sumweights > Real(0.0),
          "computeIPweights: sum of raw weights is zero (unexpected numerical error).");
      AMREX_ASSERT_WITH_MESSAGE(
          sumfluid == imp_ninterp(iim),
          "computeIPweights: mismatch in fluid stencil count.");
      // Runtime guard: in release builds ASSERT may be stripped; protect against FPE.
      Real inv_sum = (sumweights > Real(1.0e-30))
                   ? Real(1.0) / sumweights
                   : Real(0.0);
      Real check_sum = Real(0.0);
      for (int corner = 0; corner < N_InterP; ++corner) {
          weights(iim, corner) *= inv_sum;
          check_sum += weights(iim, corner);
      }
      AMREX_ASSERT_WITH_MESSAGE(
          amrex::Math::abs(check_sum - Real(1.0)) < Real(1.0e-9),
          "Interpolation point weights do not sum to 1.0");
    } // end loop over image points
}
// ============================================================================
// 5. interpolateIMs
// ============================================================================
template <int eorder_t, int iorder_t, int N_InterP = ipow(iorder_t + 1, AMREX_SPACEDIM)>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static void interpolateIMs(
    const Array3D< int, 0, eorder_t - 1, 0, N_InterP - 1, 0, AMREX_SPACEDIM-1>&  imp_ip_ijk,
    const Array2D<Real, 0, eorder_t - 1, 0, N_InterP - 1>&                       imp_ipweights,
    const Array4<Real>&                                                          prims,
    Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM-1>&                           primsNormal) noexcept
{
    for (int iim = 0; iim < eorder_t; ++iim) {
        for (int iip = 0; iip < N_InterP; ++iip) {
            const Real w = imp_ipweights(iim, iip);
            if (w == 0.0) continue;
            const int ii = imp_ip_ijk(iim, iip, 0);
            const int jj = imp_ip_ijk(iim, iip, 1);
        #if (AMREX_SPACEDIM == 3)
            const int kk = imp_ip_ijk(iim, iip, 2);
        #else
            const int kk = 0;
        #endif
            for (int n = 0; n < cls_t::NPRIM; ++n) {
                primsNormal(iim + 2, n) += prims(ii, jj, kk, n) * w;
            }
        }
    }
}
// ============================================================================
// 6. extrapolate
// ============================================================================
template <int eorder_t>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static void extrapolate(
    Array2D<Real, 0, eorder_t + 1, 0, cls_t::NPRIM - 1>& prims, 
    const Array1D< int, 0, eorder_t - 1>& imp_ninterp,
    const Real disGP, const Array1D<Real, 0, eorder_t - 1>& disIM)
{
    // Determine effective order based on INTERP_THRESHOLD
    int eff_order = eorder_t;
    for (int k = 0; k < eorder_t; ++k) {
        if (imp_ninterp(k) < INTERP_THRESHOLD_GP) {
            eff_order = k;
            break;
        }
    }
    // only extrapolate up to QLS (Last Species), skipping aux vars like QC, QG, QEINT.
    for (int n = 0; n <= cls_t::QLS; ++n) {
        
        if (eff_order >= 2) {
            // Quadratic Lagrange interpolation using surface value (slot 1)
            // and two image-point values (slots 2 and 3).
            // Note: this branch is only reachable when eorder_t >= 2.
            if constexpr (eorder_t >= 2) {
                Real u0 = prims(1, n);
                Real u1 = prims(2, n);
                Real u2 = prims(3, n);

                Real x1 = disIM(0);
                Real x2 = disIM(1);
                AMREX_ASSERT(x1 > 0 && x2 > 0);

                // Guard: if image points are nearly coincident, the quadratic
                // Lagrange denominator (x1-x2) vanishes. Fall back to linear.
                constexpr Real eps_dist = Real(1.0e-12);
                if (amrex::Math::abs(x1 - x2) < eps_dist * amrex::max(x1, x2)) {
                    // linear fallback using surface + first image point
                    Real slope = (u1 - u0) / x1;
                    prims(0, n) = u0 - slope * disGP;
                } else {
                    Real x = -disGP;
                    Real L0 = (x - x1) * (x - x2) / (x1 * x2);
                    Real L1 = x * (x - x2) / (x1 * (x1 - x2));
                    Real L2 = x * (x - x1) / (x2 * (x2 - x1));
                    prims(0, n) = u0 * L0 + u1 * L1 + u2 * L2;
                }
            }
        }
        else if (eff_order == 1) {
            Real val_surf = prims(1, n);
            Real val_im1  = prims(2, n);
            Real d_im1    = disIM(0);
            Real slope = (val_im1 - val_surf) / d_im1;
            prims(0, n) = val_surf - slope * disGP;
        }
        else {
            prims(0, n) = prims(1, n);
        }
    }
}
// ============================================================================
// 7. global2local / local2global
// ============================================================================
template <int eorder_t>
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
static void global2local(
    int iip,
    Array2D<Real,0,eorder_t+1,0,cls_t::NPRIM-1>& primsNormal,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2)
{
    const Real ux = primsNormal(iip, cls_t::QU);
    const Real uy = primsNormal(iip, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    const Real uz = primsNormal(iip, cls_t::QW);
#else
    const Real uz = Real(0.0);
#endif
    primsNormal(iip, cls_t::QU) =
        ux * norm(0) + uy * norm(1)
#if (AMREX_SPACEDIM == 3)
        + uz * norm(2)
#endif
        ;
    primsNormal(iip, cls_t::QV) =
        ux * tan1(0) + uy * tan1(1)
#if (AMREX_SPACEDIM == 3)
        + uz * tan1(2)
#endif
        ;
#if (AMREX_SPACEDIM == 3)
    primsNormal(iip, cls_t::QW) =
        ux * tan2(0) + uy * tan2(1) + uz * tan2(2);
#endif
}
template <int eorder_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static void local2global(
    int jj,
    Array2D<Real,0,eorder_t+1,0,cls_t::NPRIM-1>& primsNormal,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan1,
    const Array1D<Real,0,AMREX_SPACEDIM-1>& tan2)
{
    const Real un  = primsNormal(jj, cls_t::QU);
    const Real ut1 = primsNormal(jj, cls_t::QV);
#if (AMREX_SPACEDIM == 3)
    const Real ut2 = primsNormal(jj, cls_t::QW);
#else
    const Real ut2 = Real(0.0);
#endif
    primsNormal(jj, cls_t::QU) =
        un  * norm(0) + ut1 * tan1(0)
#if (AMREX_SPACEDIM == 3)
      + ut2 * tan2(0)
#endif
      ;
    primsNormal(jj, cls_t::QV) =
        un  * norm(1) + ut1 * tan1(1)
#if (AMREX_SPACEDIM == 3)
      + ut2 * tan2(1)
#endif
      ;
#if (AMREX_SPACEDIM == 3)
    primsNormal(jj, cls_t::QW) =
        un * norm(2) + ut1 * tan1(2) + ut2 * tan2(2);
#endif
}
// ============================================================================
// 8. check_interpolation_stencil
// ============================================================================
template <typename IPDATA, int GP_OR_SURF = is_gpData_t<IPDATA>::value ? 1 : 0>
static AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
bool check_interpolation_stencil(
    int i, int j, int k, 
    const amrex::Box& bx, 
    int lev,
    const IPDATA& ipData, 
    int f_idx,
    CheckMode mode)
{
#if (AMREX_SPACEDIM == 2)
    bool is_valid = bx.contains(amrex::IntVect(i, j)) && 
                    bx.contains(amrex::IntVect(i+1, j+1));
#else
    bool is_valid = bx.contains(amrex::IntVect(i, j, k)) && 
                    bx.contains(amrex::IntVect(i+1, j+1, k+1));
#endif
    if (!is_valid) {
        if (mode == CheckMode::Silent) {
            return false;
        }
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("Interpolation stencil out of box bounds! "
                            "lev=%d base=(%d,%d,%d) f_idx=%d\n",
                            lev, i, j, k, f_idx);
        AMREX_ASSERT(false);
#else
        int current_geom = -1;
        int current_elem = ipData.elemIdx[f_idx];
        if constexpr (GP_OR_SURF == 1) {
             current_geom = ipData.geomIdx[f_idx];
        }
#if (AMREX_SPACEDIM == 3)
        std::printf("Interpolation stencil out of box bounds!\n"
                    "  Level: %d\n"
                    "  Stencil Base: (%d, %d, %d)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n",
                    lev, i, j, k,
                    current_geom, current_elem);
#else
        std::printf("Interpolation stencil out of box bounds!\n"
                    "  Level: %d\n"
                    "  Stencil Base: (%d, %d)\n"
                    "  Geometry Index: %d\n"
                    "  Element Index:  %d\n",
                    lev, i, j,
                    current_geom, current_elem);
#endif
        std::fflush(stdout);
        if (mode == CheckMode::Warn) {
            amrex::Warning("Interpolation stencil out of box bounds!");
        } else if (mode == CheckMode::Abort) {
            amrex::Abort("Interpolation stencil out of box bounds!");
        }
#endif // AMREX_DEVICE_COMPILE
    }
    return is_valid;
}
#endif // IBM_SOLVER_INTERP_H_
