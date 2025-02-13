#ifndef Weno_H
#define Weno_H

#include <AMReX_FArrayBox.H>

#include "Closures.h"

#define POWER2(x) ((x) * (x))
#define POWER6(x) ((x) * (x) * (x) * (x) * (x) * (x))

// Choices are WenoZ5, Teno5. (Teno6 is slightly unstable)
namespace ReconScheme {
struct WenoZ5 {
  static constexpr int ng = 3;

  /**
   * @brief From a 2*ng x NCONS array, extract the left stencil.
   * @param[in] n  Index of the component.
   * @param[in] fp Flux array.
   * @param[out] s Stencil array [i+1, i, i-1, i-2, i-3].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void left_stencil(
      int n, amrex::Real const fp[2 * ng][NCONS], amrex::Real s[5]) noexcept {
    for (int m = 1; m < 2 * ng; ++m) s[m - 1] = fp[2 * ng - 1 - m][n];
  }

  /**
   * @brief From a 2*ng x NCONS array, extract the right stencil.
   * @param[in] n  Index of the component.
   * @param[in] fm Flux array.
   * @param[out] s Stencil array [i-2, i-1, i, i+1, i+2].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void right_stencil(
      int n, amrex::Real const fm[2 * ng][NCONS], amrex::Real s[5]) noexcept {
    for (int m = 1; m < 2 * ng; ++m) s[m - 1] = fm[m][n];
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void smoothness_indicator(
      const amrex::Real s[5], amrex::Real beta[3]) noexcept {
    // constexpr Real w13o12 = 13.0 / 12.0;
    // constexpr Real w1o4 = 0.25;
    // beta[2] = w13o12 * POWER2(s[4] - 2.0 * s[3] + s[2]) +
    //           w1o4 * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
    // beta[1] =
    //     w13o12 * POWER2(s[3] - 2.0 * s[2] + s[1]) + w1o4 * POWER2(s[3] - s[1]);
    // beta[0] = w13o12 * POWER2(s[2] - 2.0 * s[1] + s[0]) +
    //           w1o4 * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);

    beta[2] = Real(13. / 12.) * POWER2(s[4] - 2.0 * s[3] + s[2]) +
              Real(0.25) * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
    beta[1] = Real(13. / 12.) * POWER2(s[3] - 2.0 * s[2] + s[1]) +
              Real(0.25) * POWER2(s[3] - s[1]);
    beta[0] = Real(13. / 12.) * POWER2(s[2] - 2.0 * s[1] + s[0]) +
              Real(0.25) * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void linear_polynomial_recon(
      const amrex::Real s[5], amrex::Real vr[3]) noexcept {
    vr[2] = 11.0 * s[2] - 7.0 * s[3] + 2.0 * s[4];
    vr[1] = -s[3] + 5.0 * s[2] + 2.0 * s[1];
    vr[0] = 2.0 * s[2] + 5.0 * s[1] - s[0];
  }

  /**
   * \brief (FV)WENO-Z5. Ref https://doi.org/10.1016/j.jcp.2010.11.028.
   * \param[in] s Stencil values at cells [i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon(
      const amrex::Real s[5]) noexcept {
    using amrex::Real;
    constexpr Real eps = 1e-40;
    Real vr[3], beta[3], tmp;

    smoothness_indicator(s, beta);
    tmp = std::abs(beta[2] - beta[0]);
    beta[2] = 1.0 + tmp / (eps + beta[2]);
    beta[1] = (1.0 + tmp / (eps + beta[1])) * 6.0;
    beta[0] = (1.0 + tmp / (eps + beta[0])) * 3.0;
    tmp = 1.0 / (beta[2] + beta[1] + beta[0]);

    linear_polynomial_recon(s, vr);

    return tmp / Real(6.0) *
           (beta[2] * vr[2] + beta[1] * vr[1] + beta[0] * vr[0]);
  }
};  // struct WenoZ5

struct Teno5 : public WenoZ5 {
  /**
   * \brief (FV)TENO-5. Ref https://doi.org/10.1016/j.jcp.2015.10.037.
   * \param[in] s Stencil values at cells [i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon(
      const amrex::Real s[5]) noexcept {
    using amrex::Real;

    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    constexpr Real cutoff = 1e-4;
    Real vr[3], beta[3], tmp;
    
    smoothness_indicator(s, beta);
    tmp = std::abs(std::abs(beta[2] - beta[0]) -
                   (beta[2] + 4.0 * beta[1] + beta[0]) / 6.0);
    beta[2] = POWER6(1.0 + tmp / (eps + beta[2]));
    beta[1] = POWER6(1.0 + tmp / (eps + beta[1]));
    beta[0] = POWER6(1.0 + tmp / (eps + beta[0]));
    tmp = 1.0 / (beta[2] + beta[1] + beta[0]);
    beta[2] = beta[2] * tmp < cutoff ? 0.0 : 1.0;
    beta[1] = beta[1] * tmp < cutoff ? 0.0 : 6.0;
    beta[0] = beta[0] * tmp < cutoff ? 0.0 : 3.0;
    tmp = 1.0 / (beta[2] + beta[1] + beta[0]);

    linear_polynomial_recon(s, vr);

    return tmp / Real(6.0) *
           (beta[2] * vr[2] + beta[1] * vr[1] + beta[0] * vr[0]);
  }
};  // struct Teno5

struct Teno6 {
  static constexpr int ng = 3;

  /**
   * @brief From a 2*ng x NCONS array, extract the left stencil.
   * @param[in] n  Index of the component.
   * @param[in] fp Flux array.
   * @param[out] s Stencil array [i+2, i+1, i, i-1, i-2, i-3].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void left_stencil(
      int n, amrex::Real const fp[2 * ng][NCONS], amrex::Real s[6]) noexcept {
    for (int m = 0; m < 2 * ng; ++m) s[m] = fp[2 * ng - 1 - m][n];
  }

  /**
   * @brief From a 2*ng x NCONS array, extract the left stencil.
   * @param[in] n  Index of the component.
   * @param[in] fm Flux array.
   * @param[out] s Stencil array [i-3, i-2, i-1, i, i+1, i+2].
   */
  template <size_t NCONS>
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void right_stencil(
      int n, amrex::Real const fm[2 * ng][NCONS], amrex::Real s[6]) noexcept {
    for (int m = 0; m < 2 * ng; ++m) s[m] = fm[m][n];
  }

  /**
   * \brief (FV)TENO-6. Ref https://doi.org/10.1016/j.jcp.2015.10.037.
   * \param[in] s Stencil values at cells [i-3, i-2, i-1, i, i+1, i+2].
   * \return Reconstructed value at i-1/2.
   */
  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real recon(
      const amrex::Real s[6]) noexcept {
    using amrex::Real;

    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    constexpr Real cutoff = 1e-5;
    Real vr[4], beta[4], tmp;

    constexpr Real w13o12 = 13.0 / 12.0;
    constexpr Real w1o4 = 0.25;
    beta[3] =
        Real(1. / 36.) *
            POWER2(-11.0 * s[3] + 18.0 * s[2] - 9.0 * s[1] + 2.0 * s[0]) +
        w13o12 * POWER2(2.0 * s[3] - 5.0 * s[2] + 4.0 * s[1] - s[0]) +
        Real(781. / 720.) * POWER2(-s[3] + 3.0 * s[2] - 3.0 * s[1] + s[0]);
    beta[2] = w13o12 * POWER2(s[4] - 2.0 * s[3] + s[2]) +
              w1o4 * POWER2(s[4] - 4.0 * s[3] + 3.0 * s[2]);
    beta[1] =
        w13o12 * POWER2(s[3] - 2.0 * s[2] + s[1]) + w1o4 * POWER2(s[3] - s[1]);
    beta[0] = w13o12 * POWER2(s[2] - 2.0 * s[1] + s[0]) +
              w1o4 * POWER2(3.0 * s[2] - 4.0 * s[1] + s[0]);

    tmp = std::abs(beta[3] - (beta[2] + 4.0 * beta[1] + beta[0]) / 6.0);
    beta[3] = POWER6(1.0 + tmp / (eps + beta[3]));
    beta[2] = POWER6(1.0 + tmp / (eps + beta[2]));
    beta[1] = POWER6(1.0 + tmp / (eps + beta[1]));
    beta[0] = POWER6(1.0 + tmp / (eps + beta[0]));
    tmp = 1.0 / (beta[3] + beta[2] + beta[1] + beta[0]);
    beta[3] = beta[3] * tmp < cutoff ? 0.0 : 4.0;
    beta[2] = beta[2] * tmp < cutoff ? 0.0 : 1.0;
    beta[1] = beta[1] * tmp < cutoff ? 0.0 : 9.0;
    beta[0] = beta[0] * tmp < cutoff ? 0.0 : 6.0;
    tmp = 1.0 / (beta[3] + beta[2] + beta[1] + beta[0]);

    vr[3] = 0.5 * (s[0] - 5.0 * s[1] + 13.0 * s[2] + 3.0 * s[3]);
    vr[2] = 11.0 * s[3] - 7.0 * s[4] + 2.0 * s[5];
    vr[1] = -s[4] + 5.0 * s[3] + 2.0 * s[2];
    vr[0] = 2.0 * s[3] + 5.0 * s[2] - s[1];

    return tmp / Real(6.0) *
           (beta[3] * vr[3] + beta[2] * vr[2] + beta[1] * vr[1] +
            beta[0] * vr[0]);
  }
};  // struct Teno6
};  // namespace ReconScheme

// TODO: this class name is not very accurate. llf? high_order_fd?
template <typename Scheme, typename cls_t>
class weno_t {
  static constexpr int ng = Scheme::ng;  // number of ghost cells on each side
  static_assert(ng <= cls_t::NGHOST);    // ensure enough ghost cells

 public:
  AMREX_GPU_HOST_DEVICE
  weno_t() {}

  AMREX_GPU_HOST_DEVICE
  ~weno_t() {}

#if (AMREX_USE_GPIBM || CNS_USE_EB )
  void inline eflux_ibm(const amrex::Geometry& geom, const amrex::MFIter& mfi,
                        const amrex::Array4<const amrex::Real>& prims_in,
                        std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                        const amrex::Array4<amrex::Real>& rhs, const cls_t* cls,
                        const amrex::Array4<const bool>& ibMarkers)
#else
  void inline eflux(const amrex::Geometry& geom, const amrex::MFIter& mfi,
                    const amrex::Array4<const amrex::Real>& prims,
                    std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const amrex::Array4<amrex::Real>& rhs, const cls_t* cls)
#endif
  {
    using amrex::Array4, amrex::Box, amrex::Dim3, amrex::IntVect, amrex::Real;

    const Box& bx = mfi.tilebox();
    const auto dxinv = geom.InvCellSizeArray();

    // clear rhs
    ParallelFor(bx, cls_t::NCONS,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                  rhs(i, j, k, n) = 0.0;
                });

    // for each direction
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {

      const Box& flxbx = amrex::surroundingNodes(bx, dir);

      auto const& flx = flxt[dir]->array(); // snm

      ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        IntVect iv(AMREX_D_DECL(i, j, k));
        IntVect ivd(IntVect::TheDimensionVector(dir));

#if (AMREX_USE_GPIBM || CNS_USE_EB )
        Real prims_ptr[2 * ng * cls_t::NPRIM];
        Dim3 prims_begin = (iv - ng * ivd).dim3();
        Dim3 prims_end = (iv + (ng - 1) * ivd).dim3();
        prims_end.x += 1;
        prims_end.y += 1;
        prims_end.z += 1;
        Array4<Real> prims(prims_ptr, prims_begin, prims_end, cls_t::NPRIM);
        if (!fill_solid_prims(iv, ivd, dir, prims_in, prims, ibMarkers)) {
          return;  // skip solid cells
        }
#endif

        const Real alpha = cls->max_char_speed(iv, dir, ng, prims);
        const auto roe_avg = cls->roe_avg_state(iv, dir, prims);

        Real cons[cls_t::NCONS], f[cls_t::NCONS], fp[2 * ng][cls_t::NCONS],
            fm[2 * ng][cls_t::NCONS];
        for (int m = 0; m < 2 * ng; ++m) {
          // LLF splitting into left- and right-running fluxes
          cls->prims2flux(iv + (m - ng) * ivd, dir, prims, f);
          cls->prims2cons(iv + (m - ng) * ivd, prims, cons);

          for (int n = 0; n < cls_t::NCONS; ++n) {
            fp[m][n] = 0.5 * (cons[n] + f[n] / alpha);
            fm[m][n] = 0.5 * (cons[n] - f[n] / alpha);
          }

          // Convert into characteristic variables
          cls->cons2char(roe_avg, fp[m]);
          cls->cons2char(roe_avg, fm[m]);
        }

        // Reconstruct with upwind stencils
        Real fpL[cls_t::NCONS], fmR[cls_t::NCONS], s[2 * ng];
        for (int n = 0; n < cls_t::NCONS; ++n) {
          Scheme::left_stencil(n, fp, s);
          fpL[n] = Scheme::recon(s);

          Scheme::right_stencil(n, fm, s);
          fmR[n] = Scheme::recon(s);
        }

        // Convert back to conservative variables
        cls->char2cons(roe_avg, fpL);
        cls->char2cons(roe_avg, fmR);

        for (int n = 0; n < cls_t::NCONS; ++n) {
          flx(iv, n) = alpha * (fpL[n] - fmR[n]);
        }
      });

      // add flux derivative to rhs
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    IntVect iv(AMREX_D_DECL(i, j, k));
                    IntVect ivp = iv + IntVect::TheDimensionVector(dir);
                    rhs(i, j, k, n) += dxinv[dir] * (flx(iv, n) - flx(ivp, n));
                  });
    }  // end of for each direction
  }

#if (AMREX_USE_GPIBM || CNS_USE_EB )
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool fill_solid_prims(
      amrex::IntVect iv, amrex::IntVect ivd, int cdir,
      amrex::Array4<const amrex::Real> const &prims_in,
      amrex::Array4<amrex::Real> const &prims,
      amrex::Array4<const bool> const &ib_mask) {
    // Find the first solid point, which will be the ghost point
    int gl = ng, gr = ng;  // ghost point position on left and right
    for (int m = 0; m < ng; ++m) {
      if (ib_mask(iv + m * ivd,0)) {  // snm
        gr = std::min(gr, m);
      }
      if (ib_mask(iv - (m + 1) * ivd,0)) { //snm
        gl = std::min(gl, m);
      }
    }
    if (gl == 0 && gr == 0) return false;  // skip solid cells


     //printf("  gl gr = %d %d \n ",gl,gr);

    // // Option1: use first order recon up to ng cells away from GP
    // int glr = std::min(gl, gr);
    // if (glr < ng) {
    //   glr = 0;
    // }
    // for (int m = 0; m < ng; ++m) {
    //   if (m <= glr) {
    //     for (int n = 0; n < cls_t::NPRIM; ++n) {
    //       prims(iv + m * ivd, n) = prims_in(iv + m * ivd, n);
    //       prims(iv - (m + 1) * ivd, n) = prims_in(iv - (m + 1) * ivd, n);
    //     }
    //   } else {
    //     for (int n = 0; n < cls_t::NPRIM; ++n) {
    //       prims(iv + m * ivd, n) = prims_in(iv + glr * ivd, n);
    //       prims(iv - (m + 1) * ivd, n) = prims_in(iv - (glr + 1) * ivd, n);
    //     }
    //   }
    // }

    // 2. fill the other solid points with the ghost point values
    // if (gl == 0 || gr == 0) {
    //   gl = 0;
    //   gr = 0; // use 1st order if next to GP
    // }
    // for (int m = 0; m < ng; ++m) {
    //   if (m <= gr) {
    //     for (int n = 0; n < cls_t::NPRIM; ++n) {
    //       prims(iv + m * ivd, n) = prims_in(iv + m * ivd, n);
    //     }
    //   } else {
    //     for (int n = 0; n < cls_t::NPRIM; ++n) {
    //       prims(iv + m * ivd, n) = prims_in(iv + gr * ivd, n);
    //     }
    //   }
    // }
    // for (int m = 0; m < ng; ++m) {
    //   if (m <= gl) {
    //     for (int n = 0; n < cls_t::NPRIM; ++n) {
    //       prims(iv - (m + 1) * ivd, n) = prims_in(iv - (m + 1) * ivd, n);
    //     }
    //   } else {
    //     for (int n = 0; n < cls_t::NPRIM; ++n) {
    //       prims(iv - (m + 1) * ivd, n) = prims_in(iv - (gl + 1) * ivd, n);
    //     }
    //   }
    // }

    // 3. fill using 1st order BC (adiabatic no-slip)
    for (int m = 0; m < 2 * ng; ++m) {

      //  printf("  line 390 m = %d ng=%d \n ",m,ng);
      //  printf("iv = %d %d %d \n ",iv);
      //  printf("ivd = %d %d %d \n ",ivd);
      //  printf("iv - (m -ng)*ivd = %d %d %d \n ",iv - (m -ng)*ivd);

      // prims  <--- prims_in   at iv (i + h,j,k)    why looping?  
      for (int n = 0; n < cls_t::NPRIM; ++n) {        
        prims(iv + (m - ng) * ivd, n) = prims_in(iv + (m - ng) * ivd, n);
      }
    }


    for (int m = 0; m < ng; ++m) {
   
      if (m >= gl) {
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(iv - (m + 1) * ivd, n) = prims(iv - (2 * gl - m) * ivd, n);
        }
        // prims(iv - (m + 1) * ivd, cls_t::QU + cdir) *= -1; // this is non-permeable
        for (int c = 0; c < amrex::SpaceDim; ++c) 
	    prims(iv - (m + 1) * ivd, cls_t::QU + c) *= -1; // this is no-slip non-permeable
      }
      if (m >= gr) {
        for (int n = 0; n < cls_t::NPRIM; ++n) {
          prims(iv + m * ivd, n) = prims(iv + (2 * gr - m - 1) * ivd, n);
        }
        // prims(iv + m * ivd, cls_t::QU + cdir) *= -1; // this is non-permeable
	for (int c = 0; c < amrex::SpaceDim; ++c)
            prims(iv + m * ivd, cls_t::QU + c) *= -1; // this is no-slip non-permeable
      }
    }

    return true;
  }
#endif
};

#endif
