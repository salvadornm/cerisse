#ifndef KEEP_H_
#define KEEP_H_

#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_GpuLaunch.H>
#include <array>
#include <cstdint>
#include <cmath>

// EOS/species-aware KEEP split flux from CentralKEEP_SPEC.h, with the
// parameterized interface and marker convention used by Skew.h.
// param: order (2/4/6), dissipation, C2skew and C4skew.
// With dissipation enabled, kinetic-energy/entropy preservation is no longer
// expected. Arithmetic species averages are the same baseline as KEEP_SPEC;
// they are not a general entropy-conservative multicomponent flux.
template <typename param, typename cls_t>
class keep_t {
 public:
  static constexpr int order = param::order;
  static constexpr int halfsten = order / 2;
  static_assert(order == 2 || order == 4 || order == 6,
                "KEEP order must be 2, 4, or 6");
  static_assert(cls_t::NGHOST >= halfsten,
                "KEEP stencil exceeds available ghost cells");
  static_assert(!param::dissipation || cls_t::NGHOST >= 2,
                "KEEP dissipation sensor requires at least two ghost cells");

  amrex::Real Cshock = param::C2skew;
  amrex::Real Cdamp = param::C4skew;

  AMREX_GPU_HOST_DEVICE keep_t() = default;
  AMREX_GPU_HOST_DEVICE ~keep_t() = default;

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void eflux_ibm(const amrex::Geometry&, const amrex::MFIter& mfi,
                const amrex::Array4<amrex::Real>& prims,
                const std::array<amrex::FArrayBox*, AMREX_SPACEDIM>& flxt,
                const amrex::Array4<amrex::Real>& cons, const cls_t* cls,
                const amrex::Array4<uint8_t>& marker) {
#else
  void eflux(const amrex::Geometry&, const amrex::MFIter& mfi,
             const amrex::Array4<amrex::Real>& prims,
             const std::array<amrex::FArrayBox*, AMREX_SPACEDIM>& flxt,
             const amrex::Array4<amrex::Real>& cons, const cls_t* cls) {
#endif
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      // Only faces normal to dir; do not evaluate unused transverse nodal rows.
      const auto faces = amrex::surroundingNodes(mfi.tilebox(), dir);
      const auto flx = flxt[dir]->array();
      amrex::ParallelFor(faces,
        [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if (AMREX_USE_GPIBM || CNS_USE_EB)
          flux_dir_ibm(i,j,k,dir,prims,cons,flx,cls,marker);
#else
          flux_dir(i,j,k,dir,prims,cons,flx,cls);
#endif
        });
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void flux_dir(int i, int j, int k, int dir,
                const amrex::Array4<amrex::Real>& q,
                const amrex::Array4<amrex::Real>& cons,
                const amrex::Array4<amrex::Real>& flx, const cls_t* cls) const {
    keep_flux(i,j,k,dir,order,q,flx,cls);
    if constexpr (param::dissipation) {
      add_dissipation(i,j,k,dir,q,cons,flx,false,amrex::Array4<uint8_t>{});
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void flux_dir_ibm(int i, int j, int k, int dir,
                    const amrex::Array4<amrex::Real>& q,
                    const amrex::Array4<amrex::Real>& cons,
                    const amrex::Array4<amrex::Real>& flx, const cls_t* cls,
                    const amrex::Array4<uint8_t>& marker) const {
    const int di = int(dir == 0), dj = int(dir == 1), dk = int(dir == 2);
    // Component 0: solid; component 1: IBM ghost point or EB cut cell.
    if (marker(i,j,k,0) && marker(i-di,j-dj,k-dk,0)) {
      clear_flux(i,j,k,flx);
      return;
    }
    bool near_wall = false;
    // Check the whole high-order footprint, not just the adjacent cells:
    // a distant covered cell must never enter a fourth/sixth-order flux.
    constexpr int radius = param::dissipation && halfsten < 2 ? 2 : halfsten;
    for (int l = -radius; l < radius; ++l) {
      near_wall = near_wall || marker(i+l*di,j+l*dj,k+l*dk,0) ||
                               marker(i+l*di,j+l*dj,k+l*dk,1);
    }
    keep_flux(i,j,k,dir,near_wall ? 2 : order,q,flx,cls);
    if constexpr (param::dissipation) {
      add_dissipation(i,j,k,dir,q,cons,flx,near_wall,marker);
    }
  }

 private:
  using Real = amrex::Real;
  using amrex_array_t = amrex::Array4<Real>;

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void clear_flux(int i, int j, int k, const amrex_array_t& flx) const {
    for (int n=0; n<cls_t::NCONS; ++n) flx(i,j,k,n) = Real(0.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real keep_coefficient(int face_order, int l) const {
    if (face_order == 2) return Real(1.0);
    if (face_order == 4) return l == 1 ? Real(4.0/3.0) : Real(-1.0/6.0);
    return l == 1 ? Real(1.5) : (l == 2 ? Real(-0.3) : Real(1.0/30.0));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void keep_flux(int i, int j, int k, int dir, int face_order,
                 const amrex_array_t& q, const amrex_array_t& flx,
                 const cls_t* cls) const {
    using amrex::GpuArray;
    clear_flux(i, j, k, flx);

    const GpuArray<int, 3> vdir = {
        int(dir == 0), int(dir == 1), int(dir == 2)};

    for (int l = 1; l <= face_order / 2; ++l) {
      for (int m = 0; m <= l - 1; ++m) {
        const int i1 = i + m*vdir[0];
        const int j1 = j + m*vdir[1];
        const int k1 = k + m*vdir[2];
        const int i2 = i + (m-l)*vdir[0];
        const int j2 = j + (m-l)*vdir[1];
        const int k2 = k + (m-l)*vdir[2];

        const Real rho1 = q(i1,j1,k1,cls_t::QRHO);
        const Real rho2 = q(i2,j2,k2,cls_t::QRHO);
        const Real un1 = q(i1,j1,k1,cls_t::QU + dir);
        const Real un2 = q(i2,j2,k2,cls_t::QU + dir);
        const Real ux1 = q(i1,j1,k1,cls_t::QU);
        const Real ux2 = q(i2,j2,k2,cls_t::QU);
        const Real uy1 = q(i1,j1,k1,cls_t::QV);
        const Real uy2 = q(i2,j2,k2,cls_t::QV);
        const Real uz1 = q(i1,j1,k1,cls_t::QW);
        const Real uz2 = q(i2,j2,k2,cls_t::QW);
        const Real p1 = q(i1,j1,k1,cls_t::QPRES);
        const Real p2 = q(i2,j2,k2,cls_t::QPRES);

        Real Y1[NUM_SPECIES];
        Real Y2[NUM_SPECIES];
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          Y1[ns] = q(i1,j1,k1,cls_t::QFS + ns);
          Y2[ns] = q(i2,j2,k2,cls_t::QFS + ns);
        }

        Real e1 = Real(0.0);
        Real e2 = Real(0.0);
        cls->TYR2E(q(i1,j1,k1,cls_t::QT), Y1, rho1, e1);
        cls->TYR2E(q(i2,j2,k2,cls_t::QT), Y2, rho2, e2);

        const Real mass_flux = fg_quad(rho1, rho2, un1, un2);
        const Real kinetic_cross = Real(0.5)*
            (ux1*ux2 + uy1*uy2 + uz1*uz2);
        const Real weight = keep_coefficient(face_order, l);

        flx(i,j,k,cls_t::UMX) += weight*
            (fgh_cubic(rho1,rho2,un1,un2,ux1,ux2) +
             Real(vdir[0])*f_div(p1,p2));
        flx(i,j,k,cls_t::UMY) += weight*
            (fgh_cubic(rho1,rho2,un1,un2,uy1,uy2) +
             Real(vdir[1])*f_div(p1,p2));
        flx(i,j,k,cls_t::UMZ) += weight*
            (fgh_cubic(rho1,rho2,un1,un2,uz1,uz2) +
             Real(vdir[2])*f_div(p1,p2));
        flx(i,j,k,cls_t::UET) += weight*
            (mass_flux*kinetic_cross +
             fgh_cubic(rho1,rho2,un1,un2,e1,e2) +
             fg_div(p1,p2,un1,un2));

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          flx(i,j,k,cls_t::UFS + ns) +=
              weight*mass_flux*species_average(Y1[ns],Y2[ns]);
        }
      }
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real f_div(Real a, Real b) const { return Real(0.5)*(a+b); }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fg_div(Real a, Real b, Real c, Real d) const {
    return Real(0.5)*(a*d+b*c);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fg_quad(Real a, Real b, Real c, Real d) const {
    return Real(0.25)*(a+b)*(c+d);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fgh_cubic(Real a, Real b, Real c, Real d, Real e, Real f) const {
    return Real(0.125)*(a+b)*(c+d)*(e+f);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real species_average(Real yl, Real yr) const { return Real(0.5)*(yl+yr); }

  // JST-style shock capturing and high-order damping, as in Skew.h.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real sensor(Real left, Real center, Real right) const {
    const Real denom = right + Real(2.0)*center + left +
        std::abs(right-center) + std::abs(center-left) + Real(1.e-40);
    return std::abs(Real(2.0)*(right-Real(2.0)*center+left)/denom);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void add_dissipation(int i, int j, int k, int dir,
                       const amrex_array_t& q, const amrex_array_t& cons,
                       const amrex_array_t& flx, bool near_wall,
                       const amrex::Array4<uint8_t>& marker) const {
    const int di = int(dir == 0), dj = int(dir == 1), dk = int(dir == 2);
    auto sample = [=] AMREX_GPU_DEVICE(int l, int n) noexcept -> Real {
      return q(i+l*di,j+l*dj,k+l*dk,n);
    };
    // Adjacent states are already supplied by the EB/IBM boundary treatment.
    // For a wall sensor, use only triplets whose extra cell is fluid or an
    // IBM ghost point. If neither exists, use full two-point shock damping.
    bool left_ok = true, right_ok = true;
    if (near_wall) {
      left_ok = !marker(i-2*di,j-2*dj,k-2*dk,0) ||
                 marker(i-2*di,j-2*dj,k-2*dk,1);
      right_ok = !marker(i+di,j+dj,k+dk,0) || marker(i+di,j+dj,k+dk,1);
    }
    Real numerator = Real(0.0), denominator = Real(1.e-16);
    const int variables[2] = {cls_t::QRHO, cls_t::QPRES};
    for (int n : variables) {
      Real s = Real(0.0);
      if (left_ok) s = sensor(sample(-2,n),sample(-1,n),sample(0,n));
      if (right_ok) s = amrex::max(s,sensor(sample(-1,n),sample(0,n),sample(1,n)));
      if (!left_ok && !right_ok) s = Real(1.0);
      numerator += s*s;
      denominator += s;
    }
    const Real s = numerator/denominator;
    const Real radius = amrex::max(std::abs(sample(-1,cls_t::QU+dir))+sample(-1,cls_t::QC),
                                   std::abs(sample(0,cls_t::QU+dir))+sample(0,cls_t::QC));
    const Real eps2 = Cshock*radius*s;
    const Real eps4 = amrex::max(Real(0.0),Cdamp*radius-eps2);
    if (near_wall) {
      for (int n=0; n<cls_t::NCONS; ++n) {
        flx(i,j,k,n) -= eps2*(cons(i,j,k,n)-cons(i-di,j-dj,k-dk,n));
      }
      return;
    }

    for (int l=0; l<order; ++l) {
      Real shock, damp;
      if constexpr (order == 2) {
        shock = l == 0 ? Real(-1.0) : Real(1.0);
        damp = Real(0.0);
      } else if constexpr (order == 4) {
        const Real shocks[4] = {Real(-1.0/6.0),Real(-5.0/6.0),Real(7.0/6.0),Real(-1.0/6.0)};
        const Real damps[4] = {Real(-1.0),Real(3.0),Real(-3.0),Real(1.0)};
        shock = shocks[l]; damp = damps[l];
      } else {
        const Real shocks[6] = {Real(0.0),Real(-1.0/6.0),Real(-5.0/6.0),Real(7.0/6.0),Real(-1.0/6.0),Real(0.0)};
        const Real damps[6] = {Real(1.0),Real(-5.0),Real(10.0),Real(-10.0),Real(5.0),Real(-1.0)};
        shock = shocks[l]; damp = damps[l];
      }
      const int shift = l-halfsten;
      for (int n=0; n<cls_t::NCONS; ++n) {
        flx(i,j,k,n) += (eps4*damp-eps2*shock)*cons(i+shift*di,j+shift*dj,k+shift*dk,n);
      }
    }
  }
};
#endif
