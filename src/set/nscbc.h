#ifndef NSCBC_H_
#define NSCBC_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>

namespace nscbc {

static constexpr int LMINUS = 0;
static constexpr int LENT   = 1;
static constexpr int LTAN1  = 2;
static constexpr int LTAN2  = 3;
static constexpr int LPLUS  = 4;
static constexpr int LSP    = 5;

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int qvel (int dir)
{
    return cls_t::QU + dir; // assumes QU,QV,QW contiguous
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void tangent_dirs (int dir, int& t1, int& t2)
{
#if AMREX_SPACEDIM == 1
    t1 = -1; t2 = -1;
#elif AMREX_SPACEDIM == 2
    t1 = 1 - dir; t2 = -1;
#else
    t1 = (dir + 1) % 3;
    t2 = (dir + 2) % 3;
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void one_sided_deriv_prim (
    amrex::IntVect const& iv,
    int dir,
    int side_sign,          // low side: +1, high side: -1
    amrex::Real dxinv,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real& drho,
    amrex::Real& dun,
    amrex::Real& dut1,
    amrex::Real& dut2,
    amrex::Real& dT,
    amrex::Real* dY,
    bool second_order = true)
{
    using amrex::Real;

    const auto e = amrex::IntVect::TheDimensionVector(dir) * side_sign;
    const Real sdx = Real(side_sign) * dxinv;

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);

    auto D1 = [&] AMREX_GPU_DEVICE (int n) noexcept -> Real {
        return sdx * (-q(iv,n) + q(iv+e,n));
    };

    auto D2 = [&] AMREX_GPU_DEVICE (int n) noexcept -> Real {
        return sdx * (
            -Real(1.5)*q(iv,n)
            +Real(2.0)*q(iv+e,n)
            -Real(0.5)*q(iv+e+e,n));
    };

    auto D = [&] AMREX_GPU_DEVICE (int n) noexcept -> Real {
        return second_order ? D2(n) : D1(n);
    };

    drho = D(cls_t::QRHO);
    dun  = D(qvel<cls_t>(dir));
    dut1 = (t1 >= 0) ? D(qvel<cls_t>(t1)) : Real(0.0);
    dut2 = (t2 >= 0) ? D(qvel<cls_t>(t2)) : Real(0.0);
    dT   = D(cls_t::QT);

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        dY[ns] = D(cls_t::QFS + ns);
    }
#else
    amrex::ignore_unused(dY);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_L_lodi (
    amrex::IntVect const& iv,
    int dir,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real drho,
    amrex::Real dun,
    amrex::Real dut1,
    amrex::Real dut2,
    amrex::Real dT,
    amrex::Real const* dY,
    amrex::Real* L)
{
    using amrex::Real;

    const Real rho = q(iv, cls_t::QRHO);
    const Real T   = q(iv, cls_t::QT);
    const Real p   = q(iv, cls_t::QPRES);
    const Real c   = q(iv, cls_t::QC);
    const Real un  = q(iv, qvel<cls_t>(dir));

    // gamma = rho c^2 / p, consistent with ideal-gas p = rho R T.
    const Real gamma = rho*c*c / p;

    const Real lam_m = un - c;
    const Real lam_0 = un;
    const Real lam_p = un + c;

    L[LMINUS] = lam_m * (
          drho / (Real(2.0)*gamma)
        - rho*dun / (Real(2.0)*c)
        + rho*dT  / (Real(2.0)*gamma*T));

    L[LENT] = lam_0 * (
        -(gamma - Real(1.0))*T*drho/rho + dT/gamma);

    L[LTAN1] = lam_0 * dut1;
    L[LTAN2] = lam_0 * dut2;

    L[LPLUS] = lam_p * (
          drho / (Real(2.0)*gamma)
        + rho*dun / (Real(2.0)*c)
        + rho*dT  / (Real(2.0)*gamma*T));

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        L[LSP + ns] = lam_0 * dY[ns];
    }
#else
    amrex::ignore_unused(dY);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow (
    amrex::IntVect const& iv,
    int side_sign, // low side: +1, high side: -1
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real* L)
{
    amrex::ignore_unused(iv, q);

    // Low boundary: incoming acoustic is L+.
    // High boundary: incoming acoustic is L-.
    if (side_sign > 0) {
        L[LPLUS] = amrex::Real(0.0);
    } else {
        L[LMINUS] = amrex::Real(0.0);
    }
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void primitive_rhs_from_L (
    amrex::IntVect const& iv,
    int dir,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real const* L,
    amrex::Real& drhodt,
    amrex::Real& dudt,
    amrex::Real& dvdt,
    amrex::Real& dwdt,
    amrex::Real& dTdt,
    amrex::Real* dYdt)
{
    using amrex::Real;

    const Real rho = q(iv, cls_t::QRHO);
    const Real T   = q(iv, cls_t::QT);
    const Real p   = q(iv, cls_t::QPRES);
    const Real c   = q(iv, cls_t::QC);
    const Real gamma = rho*c*c / p;

    drhodt = -L[LMINUS] - L[LPLUS] + rho*L[LENT]/T;
    dTdt   = -(gamma - Real(1.0))*T*(L[LMINUS] + L[LPLUS])/rho - L[LENT];

    dudt = Real(0.0);
    dvdt = Real(0.0);
    dwdt = Real(0.0);

    const Real dundt = -c*(L[LPLUS] - L[LMINUS])/rho;

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);

    Real velrhs[3] = {Real(0.0), Real(0.0), Real(0.0)};
    velrhs[dir] = dundt;
#if AMREX_SPACEDIM >= 2
    if (t1 >= 0) velrhs[t1] = -L[LTAN1];
#endif
#if AMREX_SPACEDIM == 3
    if (t2 >= 0) velrhs[t2] = -L[LTAN2];
#endif

    dudt = velrhs[0];
#if AMREX_SPACEDIM >= 2
    dvdt = velrhs[1];
#endif
#if AMREX_SPACEDIM == 3
    dwdt = velrhs[2];
#endif

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        dYdt[ns] = -L[LSP + ns];
    }
#else
    amrex::ignore_unused(dYdt);
#endif
}

template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void add_lodi_rhs_to_cons (
    amrex::IntVect const& iv,
    int dir,
    int side_sign,
    amrex::Real dxinv,
    cls_t const*,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& rhs,
    bool second_order = true)
{
    using amrex::Real;

    Real drho, dun, dut1, dut2, dT;
    Real dY[NUM_SPECIES] = {Real(0.0)};

    one_sided_deriv_prim<cls_t>(
        iv, dir, side_sign, dxinv, q,
        drho, dun, dut1, dut2, dT, dY, second_order);

    Real L[5 + NUM_SPECIES] = {Real(0.0)};

    compute_L_lodi<cls_t>(
        iv, dir, q, drho, dun, dut1, dut2, dT, dY, L);

    apply_lodi_outflow<cls_t>(iv, side_sign, q, L);

    Real drhodt, dudt, dvdt, dwdt, dTdt;
    Real dYdt[NUM_SPECIES] = {Real(0.0)};

    primitive_rhs_from_L<cls_t>(
        iv, dir, q, L,
        drhodt, dudt, dvdt, dwdt, dTdt, dYdt);

    const Real rho = q(iv, cls_t::QRHO);
    const Real u   = q(iv, cls_t::QU);
#if AMREX_SPACEDIM >= 2
    const Real v = q(iv, cls_t::QV);
#else
    const Real v = Real(0.0);
#endif
#if AMREX_SPACEDIM == 3
    const Real w = q(iv, cls_t::QW);
#else
    const Real w = Real(0.0);
#endif

    rhs(iv, cls_t::URHO) += drhodt;
    rhs(iv, cls_t::UMX ) += u*drhodt + rho*dudt;
#if AMREX_SPACEDIM >= 2
    rhs(iv, cls_t::UMY ) += v*drhodt + rho*dvdt;
#endif
#if AMREX_SPACEDIM == 3
    rhs(iv, cls_t::UMZ ) += w*drhodt + rho*dwdt;
#endif

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        const Real Y = q(iv, cls_t::QFS + ns);
        rhs(iv, cls_t::UFS + ns) += Y*drhodt + rho*dYdt[ns];
    }
#endif

    // Ideal-gas LODI energy closure for first implementation.
    const Real p = q(iv, cls_t::QPRES);
    const Real T = q(iv, cls_t::QT);
    const Real c = q(iv, cls_t::QC);
    const Real gamma = rho*c*c / p;
    const Real Rmix = p/(rho*T);
    const Real cv = Rmix/(gamma - Real(1.0));
    const Real e = cv*T;

    const Real dke =
        u*dudt
#if AMREX_SPACEDIM >= 2
        + v*dvdt
#endif
#if AMREX_SPACEDIM == 3
        + w*dwdt
#endif
        ;

    const Real ke = Real(0.5)*(u*u + v*v + w*w);

    rhs(iv, cls_t::UET) += (e + ke)*drhodt + rho*(cv*dTdt + dke);
}

} // namespace nscbc

#endif
