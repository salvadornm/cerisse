#ifndef Weno_H
#define Weno_H

#include <AMReX_FArrayBox.H>
// #include <Index.h>

///
/// \brief Template class for Riemann solvers
///
/// \param iSplit Flux vector splitting (0 - Global lax-friedrichs) 
/// \param isChar Characteristic transform 
///
/// ```
/// {rst}
///
/// :math:`f_{i+1/2}= `
/// ```
///
template <int iSplit, bool isChar, typename cls_t>
class weno_t {
  public: 
  AMREX_GPU_HOST_DEVICE
  weno_t() {}
  ~weno_t() {}

  /// \brief Template class for Riemann solvers
  ///
  /// \param ivars interpolation vars
  ///
  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& ivars,
                    const Array4<Real>& rhs,
                    const cls_t& cls) {
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(NGHOST);
    const Box& bxn = mfi.grownnodaltilebox(-1, 0);  // 0,N+1 all directions

    FArrayBox varsf(bxg, NCONS, The_Async_Arena());
    FArrayBox lambdaf(bxg, 1, The_Async_Arena());
    FArrayBox tempf(bxg, NCONS, The_Async_Arena());

    Array4<Real> const& vars= varsf.array();
    Array4<Real> const& lambda= lambdaf.array();
    Array4<Real> const& temp= tempf.array();

    // zero rhs
    ParallelFor(bxg, NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int
    n) noexcept {rhs(i, j, k, n)=0.0;});

    // for each direction
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

      // prims to fluxes
      // compute eigenvalues
      ParallelFor(bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cls.prims2fluxes(i, j, k, prims, ivars, vdir);
        Real cs = cls.sos(prims(i,j,k,cls.QT));
        Real udir = prims(i,j,k,cls.QU) * vdir[0] + prims(i,j,k,cls.QV) * vdir[1] + prims(i,j,k,cls.QW) * vdir[2];
        lambda(i, j, k, 0) = std::abs(max(udir + cs, udir - cs, udir));
      });


      // split flux
      for (int is=-1; is<1; is+=2) {

      }

      // compute interpolation vars
      if constexpr (isChar) {
        Array2D<Real,0,NCONS-1,0,NCONS-1> jac;
        ParallelFor(bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          // R = 
          // ivars(i,j,k,) = R*var;
        });
      }

      // interpolate split-flux at i-1/2, j-1/2, k-1/2
      ParallelFor(bxn, NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    // flux_dir(i, j, k, n, prims, ivars, lambda, temp);

                    // weno5js()
                  });

      // // transform interpolation vars to fluxes at i-1/2, j-1/2 and k-1/2
      // if constexpr (isChar) {
      //   printf("WENO characteristic transform to do");
      //   exit(0);
      // }

      // // add flux derivative to rhs = -(fi+1 - fi)/dx = (fi - fi+1)/dx
      // ParallelFor(
      //     bx, NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      //       rhs(i, j, k, n) +=
      //           dxinv[dir] * (temp(i, j, k, n) -
      //                         temp(i + vdir[0], j + vdir[1], k + vdir[2], n));
      //     });
    }

  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(int& i, int& j, int& k, int& n, const Array4<Real>& prims, const Array4<Real>& vars, const Array4<Real>& lambda, Array4<Real>& temp) {
    // Real maxeigen = max(lambda(i - 3, j, k, 0), lambda(i - 2, j, k, 0),
    //                     lambda(i - 1, j, k, 0), lambda(i, j, k, 0),
    //                     lambda(i + 1, j, k, 0), lambda(i + 2, j, k, 0));

    // df[0] = maxeigen * cons(i - 3, j, k, n);
    // df[1] = maxeigen * cons(i - 2, j, k, n);
    // df[2] = maxeigen * cons(i - 1, j, k, n);
    // df[3] = maxeigen * cons(i, j, k, n);
    // df[4] = maxeigen * cons(i + 1, j, k, n);
    // df[5] = maxeigen * cons(i + 2, j, k, n);

    // // f(i-1/2)^+
    // sten[0] = Real(0.5) * (pfx(i - 3, j, k, n) + df[0]);
    // sten[1] = Real(0.5) * (pfx(i - 2, j, k, n) + df[1]);
    // sten[2] = Real(0.5) * (pfx(i - 1, j, k, n) + df[2]);
    // sten[3] = Real(0.5) * (pfx(i, j, k, n) + df[3]);
    // sten[4] = Real(0.5) * (pfx(i + 1, j, k, n) + df[4]);
    // Real ivar = weno5js(sten);

    // // f(i-1/2)^- (note, we are flipping the stencil)
    // sten[4] = Real(0.5) * (pfx(i - 2, j, k, n) - df[1]);
    // sten[3] = Real(0.5) * (pfx(i - 1, j, k, n) - df[2]);
    // sten[2] = Real(0.5) * (pfx(i, j, k, n) - df[3]);
    // sten[1] = Real(0.5) * (pfx(i + 1, j, k, n) - df[4]);
    // sten[0] = Real(0.5) * (pfx(i + 2, j, k, n) - df[5]);
    // ivar += weno5js(sten);

    // temp(i, j, k, n) = flx;
  }


  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real weno5js(const GpuArray<Real, 5>& s) {
    constexpr Real eps = 1e-40;
    constexpr int q = 2;  // tunable positive integer
    Real vl[3], beta[3], alpha[3];
    beta[0] = (13.0_rt / 12.0_rt) * (s[0] - 2.0_rt * s[1] + s[2]) *
                  (s[0] - 2.0_rt * s[1] + s[2]) +
              0.25 * (s[0] - 4.0_rt * s[1] + 3.0_rt * s[2]) *
                  (s[0] - 4.0_rt * s[1] + 3.0_rt * s[2]);
    beta[1] =
        (13.0_rt / 12.0_rt) * (s[1] - 2.0_rt * s[2] + s[3]) *
            (s[1] - 2.0_rt * s[2] + s[3]) +
        0.25 * (s[1] - s[3]) * (s[1] - s[3]) * (s[1] - s[3]) * (s[1] - s[3]);
    beta[2] = (13.0_rt / 12.0_rt) * (s[2] - 2.0_rt * s[3] + s[4]) *
                  (s[2] - 2.0_rt * s[3] + s[4]) +
              0.25_rt * (3.0_rt * s[2] - 4.0_rt * s[3] + s[4]) *
                  (3.0_rt * s[2] - 4.0_rt * s[3] + s[4]);
    alpha[0] = 0.1_rt / pow(eps + beta[0], q);
    alpha[1] = 0.6_rt / pow(eps + beta[1], q);
    alpha[2] = 0.3_rt / pow(eps + beta[2], q);
    Real sum = 1.0_rt / (alpha[0] + alpha[1] + alpha[2]);
    vl[0] = 2.0_rt * s[0] - 7.0_rt * s[1] + 11.0_rt * s[2];
    vl[1] = -s[1] + 5.0_rt * s[2] + 2.0_rt * s[3];
    vl[2] = 2.0_rt * s[2] + 5.0_rt * s[3] - s[4];
    Real fr = (Real(1.0) / 6.0_rt) * sum *
              (alpha[0] * vl[0] + alpha[1] * vl[1] + alpha[2] * vl[2]);
    return fr;
  };


};




// // Computes Euler fluxes with global lax-friedrichs splitting at
// // i-1/2,j-1/2,k-1/2
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void numericalflux_globallaxsplit(
//     int i, int j, int k, int n, const auto& cons, const auto& pfx,
//     const auto& pfy, const auto& pfz, const auto& lambda, const auto& nfx,
//     const auto& nfy, const auto& nfz) {
//   GpuArray<Real, 5> sten;
//   GpuArray<Real, 6> df = {0};

//   // x direction /////////////////////////////////////////////////////////
//   Real maxeigen = max(lambda(i - 3, j, k, 0), lambda(i - 2, j, k, 0),
//                       lambda(i - 1, j, k, 0), lambda(i, j, k, 0),
//                       lambda(i + 1, j, k, 0), lambda(i + 2, j, k, 0));

//   df[0] = maxeigen * cons(i - 3, j, k, n);
//   df[1] = maxeigen * cons(i - 2, j, k, n);
//   df[2] = maxeigen * cons(i - 1, j, k, n);
//   df[3] = maxeigen * cons(i, j, k, n);
//   df[4] = maxeigen * cons(i + 1, j, k, n);
//   df[5] = maxeigen * cons(i + 2, j, k, n);

//   // f(i-1/2)^+
//   sten[0] = Real(0.5) * (pfx(i - 3, j, k, n) + df[0]);
//   sten[1] = Real(0.5) * (pfx(i - 2, j, k, n) + df[1]);
//   sten[2] = Real(0.5) * (pfx(i - 1, j, k, n) + df[2]);
//   sten[3] = Real(0.5) * (pfx(i, j, k, n) + df[3]);
//   sten[4] = Real(0.5) * (pfx(i + 1, j, k, n) + df[4]);
//   Real flx = weno5js(sten);

//   // f(i-1/2)^- (note, we are flipping the stencil)
//   sten[4] = Real(0.5) * (pfx(i - 2, j, k, n) - df[1]);
//   sten[3] = Real(0.5) * (pfx(i - 1, j, k, n) - df[2]);
//   sten[2] = Real(0.5) * (pfx(i, j, k, n) - df[3]);
//   sten[1] = Real(0.5) * (pfx(i + 1, j, k, n) - df[4]);
//   sten[0] = Real(0.5) * (pfx(i + 2, j, k, n) - df[5]);
//   flx += weno5js(sten);
//   nfx(i, j, k, n) = flx;

//   // y direction /////////////////////////////////////////////////////////
//   maxeigen = max(lambda(i, j - 3, k, 1), lambda(i, j - 2, k, 1),
//                  lambda(i, j - 1, k, 1), lambda(i, j, k, 1),
//                  lambda(i, j + 1, k, 1), lambda(i, j + 2, k, 1));
//   df[0] = maxeigen * cons(i, j - 3, k, n);
//   df[1] = maxeigen * cons(i, j - 2, k, n);
//   df[2] = maxeigen * cons(i, j - 1, k, n);
//   df[3] = maxeigen * cons(i, j, k, n);
//   df[4] = maxeigen * cons(i, j + 1, k, n);
//   df[5] = maxeigen * cons(i, j + 2, k, n);

//   // f(j-1/2)^+
//   sten[0] = Real(0.5) * (pfy(i, j - 3, k, n) + df[0]);
//   sten[1] = Real(0.5) * (pfy(i, j - 2, k, n) + df[1]);
//   sten[2] = Real(0.5) * (pfy(i, j - 1, k, n) + df[2]);
//   sten[3] = Real(0.5) * (pfy(i, j, k, n) + df[3]);
//   sten[4] = Real(0.5) * (pfy(i, j + 1, k, n) + df[4]);
//   flx = weno5js(sten);

//   // f(j-1/2)^- (note, we are flipping the stencil)
//   sten[4] = Real(0.5) * (pfy(i, j - 2, k, n) - df[1]);
//   sten[3] = Real(0.5) * (pfy(i, j - 1, k, n) - df[2]);
//   sten[2] = Real(0.5) * (pfy(i, j, k, n) - df[3]);
//   sten[1] = Real(0.5) * (pfy(i, j + 1, k, n) - df[4]);
//   sten[0] = Real(0.5) * (pfy(i, j + 2, k, n) - df[5]);
//   flx += weno5js(sten);
//   nfy(i, j, k, n) = flx;

//   // z direction /////////////////////////////////////////////////////////
//   maxeigen = max(lambda(i, j, k - 3, 2), lambda(i, j, k - 2, 2),
//                  lambda(i, j, k - 1, 2), lambda(i, j, k, 2),
//                  lambda(i, j, k + 1, 2), lambda(i, j, k + 2, 2));
//   df[0] = maxeigen * cons(i, j, k - 3, n);
//   df[1] = maxeigen * cons(i, j, k - 2, n);
//   df[2] = maxeigen * cons(i, j, k - 1, n);
//   df[3] = maxeigen * cons(i, j, k, n);
//   df[4] = maxeigen * cons(i, j, k + 1, n);
//   df[5] = maxeigen * cons(i, j, k + 2, n);

//   // f(k-1/2)^+
//   sten[0] = Real(0.5) * (pfz(i, j, k - 3, n) + df[0]);
//   sten[1] = Real(0.5) * (pfz(i, j, k - 2, n) + df[1]);
//   sten[2] = Real(0.5) * (pfz(i, j, k - 1, n) + df[2]);
//   sten[3] = Real(0.5) * (pfz(i, j, k, n) + df[3]);
//   sten[4] = Real(0.5) * (pfz(i, j, k + 1, n) + df[4]);
//   flx = weno5js(sten);

//   // f(k-1/2)^- (note, we are flipping the stencil)
//   sten[4] = Real(0.5) * (pfz(i, j, k - 2, n) - df[1]);
//   sten[3] = Real(0.5) * (pfz(i, j, k - 1, n) - df[2]);
//   sten[2] = Real(0.5) * (pfz(i, j, k, n) - df[3]);
//   sten[1] = Real(0.5) * (pfz(i, j, k + 1, n) - df[4]);
//   sten[0] = Real(0.5) * (pfz(i, j, k + 2, n) - df[5]);
//   flx += weno5js(sten);
//   nfz(i, j, k, n) = flx;
// }

#endif