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

static constexpr int NLWAVES = NUM_SPECIES + 5;

static constexpr amrex::Real r1_2 =  amrex::Real(1.0) / amrex::Real(2.0); 
static constexpr amrex::Real r1_3 =  amrex::Real(1.0) / amrex::Real(3.0); 


static constexpr bool use_transverse_terms = false;    //later can be passed as a parameter
static constexpr amrex::Real transverse_scale = 1.0;  // use it as parameter


// Default parameter object for LODI-only NSCBCs.
// The solver may also pass any user-defined nscbc_parm_t with the same fields.
// (defiend as well in CNS.h coudl be cleaned!!)
struct NSCBCParm {
    amrex::Real Lchar   = 1.0;
    amrex::Real Mmax    = 0.1;
    amrex::Real Ptarget = 101325.0;
    amrex::Real sigma   = 0.28;
    amrex::Real utarget = 0.0;
    amrex::Real vtarget = 0.0;
    amrex::Real wtarget = 0.0;
    amrex::Real Ttarget = 300.0;
    amrex::Real eta     = 1.0;        
};


template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void check_nscbc (Vector<int> &nslo, Vector<int> &nshi, nscbc_parm_t const& nscbc_parm)
{
    amrex::Print() << " Using NSBC:  \n";
    
    // check if NSBC boundaries are well posed
    int type = 0;
    for (int n=0;n< AMREX_SPACEDIM;n++){
        for (int bcside=0; bcside< 2; bcside++) {
            if (bcside==0) {
                //type = cls_t::nscbc_type_lo[n]; 
                type = nslo[n];   
                amrex::Print() << " dim = " << n <<  "  nscbc_type_lo " <<  type << "\n";
            }
            else {
                //type = cls_t::nscbc_type_hi[n];  
                type = nshi[n]; 
                amrex::Print() << " dim = " << n <<  "  nscbc_type_hi " <<  type << "\n";  
            }
            // classify BC
            if (type ==1 ){
                amrex::Print() << " non-reflecting inflow with target relaxation  \n";
                amrex::Print() << "      utarget "   <<  nscbc_parm.utarget << "\n";
                amrex::Print() << "      vtarget "   <<  nscbc_parm.vtarget << "\n";
                amrex::Print() << "      wtarget "   <<  nscbc_parm.wtarget << "\n";
                amrex::Print() << "      Ttarget "   <<  nscbc_parm.Ttarget << "\n";
                amrex::Print() << "      eta     "   <<  nscbc_parm.eta     << "\n";
            }
            else if (type == 2) {
                amrex::Print() << " purely non-reflecting outflow \n";
            }
            else if (type == 3) {
                amrex::Print() << " non-reflecting outflow  targeting pressure  \n";        
                amrex::Print() << "      p_target "   <<  nscbc_parm.Ptarget << "\n";  
                amrex::Print() << "      sigma "      <<  nscbc_parm.sigma << "\n";
                amrex::Print() << "      Lchar "      <<  nscbc_parm.Lchar  << "\n";                         
                amrex::Print() << "      Mmax  "      <<  nscbc_parm.Mmax  << "\n";                                                       
            }
            else if (type == 0) {                          
                amrex::Print() << " no NSBC  \n";        
            }
        }    
    }
    //

}
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
//-----------------------------------------------------------------
// Normal derivative at a physical boundary face.
//
// qbc is face-centred in direction dir and stores Q_BC(y,z) (or the
// corresponding tangential coordinates). q is the cell-centred primitive
// state. iv_inner is the first cell inside the domain.
//
// side_sign = +1 at a low boundary and -1 at a high boundary.
// Thus iv_inner + side_sign*e_dir is always the second interior cell.
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void one_sided_deriv_prim (
    amrex::IntVect const& iv_face,
    amrex::IntVect const& iv_inner,
    int dir,
    int side_sign,
    amrex::Real dxinv,
    amrex::Array4<const amrex::Real> const& qbc,
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

    const auto ein = amrex::IntVect::TheDimensionVector(dir) * side_sign;

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);

    auto D = [&] AMREX_GPU_DEVICE (int n) noexcept -> Real {
        const Real qb = qbc(iv_face, n);
        const Real q1 = q(iv_inner, n);

        if (second_order) {
            const Real q2 = q(iv_inner + ein, n);
            // Low side : (-8 qb + 9 q1 - q2)/(3 dx)            
            return Real(side_sign) * r1_3 * dxinv * (-Real(8.0)*qb + Real(9.0)*q1 - q2);                 
        }

        // Distance from the boundary face to the first cell centre is dx/2.
        return Real(side_sign) * Real(2.0)*dxinv * (q1 - qb);
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
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void one_sided_deriv_ghost_prim(
    IntVect const& iv,
    int dir,
    int side_sign,
    Real dxinv,
    Array4<const Real> const& q,
    Real& drho,
    Real& dun,
    Real& dut1,
    Real& dut2,
    Real& dT,
    Real* dY,
    int order)
{
    const IntVect ein =
        side_sign * IntVect::TheDimensionVector(dir);

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);

    auto D = [=] AMREX_GPU_DEVICE (int n) noexcept -> Real
    {
        const Real q0 = q(iv, n);
        const Real q1 = q(iv + ein, n);

        if (order <= 1) {
            return Real(side_sign) * dxinv * (q1 - q0);
        }

        const Real q2 = q(iv + 2*ein, n);
        if (order == 2) {
            return Real(side_sign) * Real(0.5) * dxinv * (-Real(3.0)*q0 +Real(4.0)*q1 -q2);
        }

        const Real q3 = q(iv + 3*ein, n);
        if (order == 3) {
            return Real(side_sign) * dxinv / Real(6.0) * (-Real(11.0)*q0 +Real(18.0)*q1 -Real(9.0)*q2 +Real(2.0)*q3);
        }

        const Real q4 = q(iv + 4*ein, n);
        return Real(side_sign) * dxinv / Real(12.0) * (-Real(25.0)*q0 +Real(48.0)*q1 -Real(36.0)*q2 +Real(16.0)*q3 -Real(3.0)*q4);
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
//-----------------------------------------------------------------
// Compute the normal LODI waves from the local face primitive state.
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_L_lodi (
    amrex::IntVect const& iv_face,
    int dir,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real drho,
    amrex::Real dun,
    amrex::Real dut1,
    amrex::Real dut2,
    amrex::Real dT,
    amrex::Real const* dY,
    amrex::Real* L)
{
    using amrex::Real;

    const Real rho   = qbc(iv_face, cls_t::QRHO);
    const Real T     = qbc(iv_face, cls_t::QT);
    const Real c     = qbc(iv_face, cls_t::QC);
    const Real un    = qbc(iv_face, qvel<cls_t>(dir));
    const Real gamma = qbc(iv_face, cls_t::QG);

    const Real lam_m = un - c;
    const Real lam_0 = un;
    const Real lam_p = un + c;
    const Real o_2gamma = Real(0.5) / gamma;

    L[LMINUS] = lam_m * o_2gamma *  (drho - gamma*rho*dun/c + rho*dT/T);

    // L[LENT] = lam_0 *
    //     (-(gamma - Real(1.0))*T*drho/rho + dT/gamma);

    L[LENT] = lam_0 / gamma * (dT - (gamma - Real(1.0))*T*drho/rho);  //NEw      

    L[LTAN1] = lam_0 * dut1;
    L[LTAN2] = lam_0 * dut2;

    L[LPLUS] = lam_p * o_2gamma *  (drho + gamma*rho*dun/c + rho*dT/T);

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        L[LSP + ns] = lam_0 * dY[ns];
    }
#else
    amrex::ignore_unused(dY);
#endif
}

//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow (
    amrex::IntVect const& iv_face,
    int side_sign,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real* L)
{
    amrex::ignore_unused(iv_face, qbc);

    // Low boundary: L+ is incoming. High boundary: L- is incoming.
    if (side_sign > 0) {
        L[LPLUS] = amrex::Real(0.0);
    } else {
        L[LMINUS] = amrex::Real(0.0);
    }
}
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow_transverse(
    int side_sign,
    amrex::Real const* Tchar,
    amrex::Real* L)
{
    using amrex::Real;

    // Low boundary: L+ + T+ is incoming.
    if (side_sign > 0) {
        L[LPLUS]  = -transverse_scale*Tchar[LPLUS];
    }
    // High boundary: L- + T- is incoming.
    else {
        L[LMINUS] = -transverse_scale*Tchar[LMINUS];
    }
}
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow_pressure_relaxation (
    amrex::IntVect const& iv_face,
    int side_sign,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real* L,
    nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;

    const Real p = qbc(iv_face, cls_t::QPRES);
    const Real c = qbc(iv_face, cls_t::QC);
    const Real dp = p - nscbc_parm.Ptarget;

    const Real relax = nscbc_parm.sigma
        * (Real(1.0) - nscbc_parm.Mmax*nscbc_parm.Mmax)
        / (Real(2.0) * c * nscbc_parm.Lchar)
        * dp;

    if (side_sign > 0) {
        L[LPLUS] = relax;
    } else {
        L[LMINUS] = relax;
    }
}
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow_pressure_relaxation_transverse(
    amrex::IntVect const& iv_face,
    int side_sign,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real const* Tchar,
    amrex::Real* L,
    nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;

    const Real p  = qbc(iv_face, cls_t::QPRES);
    const Real c  = qbc(iv_face, cls_t::QC);
    const Real dp = p - nscbc_parm.Ptarget;

    const Real relax = nscbc_parm.sigma * ( Real(1.0) - nscbc_parm.Mmax*nscbc_parm.Mmax)
        / ( Real(2.0) * c * nscbc_parm.Lchar )* dp;

    if (side_sign > 0) {
        L[LPLUS]  = relax - transverse_scale*Tchar[LPLUS];
    } else {
        L[LMINUS] = relax - transverse_scale*Tchar[LMINUS];
    }
}
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real target_velocity_component (
    int dir,
    amrex::Real utarget,
    amrex::Real vtarget,
    amrex::Real wtarget)
{
    if (dir == 0) return utarget;
#if AMREX_SPACEDIM >= 2
    if (dir == 1) return vtarget;
#endif
#if AMREX_SPACEDIM == 3
    if (dir == 2) return wtarget;
#endif
    amrex::ignore_unused(wtarget);
    return amrex::Real(0.0);
}
//-----------------------------------------------------------------
// LODI-only subsonic non-reflecting inflow with target relaxation.
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_inflow_relaxation (
    amrex::IntVect const& iv_face,
    int dir,
    int side_sign,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real* L,
    nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;

    const Real rho = qbc(iv_face, cls_t::QRHO);
    const Real c   = qbc(iv_face, cls_t::QC);
    const Real T   = qbc(iv_face, cls_t::QT);
    const Real eta = nscbc_parm.eta;
    amrex::ignore_unused(T);

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);

    const Real u_target[3] = {
        nscbc_parm.utarget,
#if AMREX_SPACEDIM >= 2
        nscbc_parm.vtarget,
#else
        Real(0.0),
#endif
#if AMREX_SPACEDIM == 3
        nscbc_parm.wtarget
#else
        Real(0.0)
#endif
    };

    L[LENT] = eta * (qbc(iv_face, cls_t::QT) - nscbc_parm.Ttarget);

#if AMREX_SPACEDIM >= 2
    if (t1 >= 0) {
        L[LTAN1] = eta * (qbc(iv_face, qvel<cls_t>(t1)) - u_target[t1]);
    }
#endif
#if AMREX_SPACEDIM == 3
    if (t2 >= 0) {
        L[LTAN2] = eta * (qbc(iv_face, qvel<cls_t>(t2)) - u_target[t2]);
    }
#endif

    const Real un      = qbc(iv_face, qvel<cls_t>(dir));
    const Real un_targ = u_target[dir];
    const Real acoustic_delta = rho/c * eta * (un - un_targ);

    if (side_sign > 0) {
        L[LPLUS]  = L[LMINUS] + acoustic_delta;
    } else {
        L[LMINUS] = L[LPLUS] - acoustic_delta;
    }

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        L[LSP + ns] = Real(0.0);
    }
#endif
}

//-----------------------------------------------------------------
// Convert the L waves into a primitive-variable RHS at the face.
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void primitive_rhs_from_L (
    amrex::IntVect const& iv_face,
    int dir,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real const* L,
    amrex::Real& drhodt,
    amrex::Real& dudt,
    amrex::Real& dvdt,
    amrex::Real& dwdt,
    amrex::Real& dTdt,
    amrex::Real* dYdt)
{
    using amrex::Real;

    const Real rho   = qbc(iv_face, cls_t::QRHO);
    const Real T     = qbc(iv_face, cls_t::QT);
    const Real c     = qbc(iv_face, cls_t::QC);
    const Real gamma = qbc(iv_face, cls_t::QG);

    drhodt = -L[LMINUS] - L[LPLUS] + rho*L[LENT]/T;
    dTdt   = -(gamma - Real(1.0))*T*
             (L[LMINUS] + L[LPLUS])/rho - L[LENT];

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
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_transverse_characteristics(
    amrex::IntVect const& iv_face,
    amrex::IntVect const& iv_inner,
    int dir,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real* Tchar)
{
    using amrex::Real;

    // Characteristic basis evaluated from the evolved boundary state.
    const Real rho   = qbc(iv_face, cls_t::QRHO);
    const Real T     = qbc(iv_face, cls_t::QT);
    const Real c     = qbc(iv_face, cls_t::QC);
    const Real gamma = qbc(iv_face, cls_t::QG);

    // Transverse convection evaluated from the first interior plane.
    Real vel[3] = {
        q(iv_inner, cls_t::QU),
#if AMREX_SPACEDIM >= 2
        q(iv_inner, cls_t::QV),
#else
        Real(0.0),
#endif
#if AMREX_SPACEDIM == 3
        q(iv_inner, cls_t::QW)
#else
        Real(0.0)
#endif
    };

    auto Dc =
        [&] AMREX_GPU_DEVICE(int tdir, int n) noexcept -> Real
    {
        const auto e =
            amrex::IntVect::TheDimensionVector(tdir);

        return Real(0.5)*dxinv[tdir]*
            (
                q(iv_inner + e, n)
              - q(iv_inner - e, n)
            );
    };

    Real Trho = Real(0.0);
    Real TT   = Real(0.0);

    Real Tvel[3] = {
        Real(0.0),
        Real(0.0),
        Real(0.0)
    };

#if NUM_SPECIES > 1
    Real TY[NUM_SPECIES] = {Real(0.0)};
#endif

    for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {

        if (tdir == dir) continue;

        const Real ut = vel[tdir];

        Trho +=
              ut*Dc(tdir, cls_t::QRHO)
            + rho*Dc(tdir, qvel<cls_t>(tdir));

        TT +=
              ut*Dc(tdir, cls_t::QT)
            + (gamma - Real(1.0))*T
              * Dc(tdir, qvel<cls_t>(tdir));

        for (int m = 0; m < AMREX_SPACEDIM; ++m) {
            Tvel[m] +=
                ut*Dc(tdir, qvel<cls_t>(m));
        }

        Tvel[tdir] +=
            Dc(tdir, cls_t::QPRES)/rho;

#if NUM_SPECIES > 1
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            TY[ns] +=
                ut*Dc(tdir, cls_t::QFS + ns);
        }
#endif
    }

    const Real Tun = Tvel[dir];
    const Real inv_2gamma = Real(0.5)/gamma;

    Tchar[LMINUS] =
        inv_2gamma*
        (
            Trho
          - gamma*rho*Tun/c
          + rho*TT/T
        );

    Tchar[LENT] =
        (
            TT
          - (gamma - Real(1.0))*T*Trho/rho
        )/gamma;

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);

    Tchar[LTAN1] =
        (t1 >= 0) ? Tvel[t1] : Real(0.0);

    Tchar[LTAN2] =
        (t2 >= 0) ? Tvel[t2] : Real(0.0);

    Tchar[LPLUS] =
        inv_2gamma*
        (
            Trho
          + gamma*rho*Tun/c
          + rho*TT/T
        );

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        Tchar[LSP + ns] = TY[ns];
    }
#endif
}
//-----------------------------------------------------------------
// Compute the conservative RHS of the face-centred U_BC state.
// qbc and rhs_bc are face-centred; q is the cell-centred interior primitive
// field. All local boundary quantities vary independently with iv_face.
//
//-----------------------------------------------------------------
//
//
//-----------------------------------------------------------------


//------------------------------------------------------------------------------
// 
//
//-------//-----------------------------------------------------------------
// Compute the conservative NSCBC RHS at one cell-centred ghost cell.
//
// iv
//     Cell-centred index of the ghost-shell cell being advanced.
//
// dir
//     Direction normal to the physical boundary.
//
// side_sign
//     +1 at a low physical boundary.
//     -1 at a high physical boundary.
//
//     Therefore:
//
//         iv + side_sign*e_dir
//
//     always points inward.
//
// q
//     Cell-centred primitive field for the complete RK-stage state.
//     The persistent ghost shell must already have overwritten the
//     physical ghost cells of the conservative stage state before q
//     is constructed.
//
// rhs_ghost
//     Conservative ghost-shell RHS.
//
// derivative_order
//     1: first-order inward derivative.
//     2: second-order inward derivative.
//
// use_transverse_here
//     Should be true only for face-interior ghost cells. Until the
//     Lodato edge/corner treatment is implemented, pass false for
//     cells lying outside the domain in two or more directions.
//-----------------------------------------------------------------
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void add_lodi_ghost_rhs_to_cons(
    amrex::IntVect const& iv,
    int dir,
    int side_sign,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
    cls_t const* cls,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& rhs_ghost,
    int derivative_order,
    int nscbc_type,
    nscbc_parm_t const& nscbc_parm)    
{
    using amrex::Real;

    amrex::ignore_unused(cls);

    // ==============================================================
    // 1. Normal primitive-variable derivatives
    // ==============================================================

    Real drho = Real(0.0);
    Real dun  = Real(0.0);
    Real dut1 = Real(0.0);
    Real dut2 = Real(0.0);
    Real dT   = Real(0.0);

    Real dY[NUM_SPECIES] = {Real(0.0)};

    one_sided_deriv_ghost_prim<cls_t>(iv,dir,side_sign,dxinv[dir],q,drho,dun,dut1,dut2,dT,dY,derivative_order);

    // ==============================================================
    // 2. Normal characteristic-wave amplitudes
    //
    // compute_L_lodi only requires a local index and a primitive
    // Array4, so the cell-centred primitive field can be passed
    // directly in place of the previous face-centred qbc.
    // ==============================================================

    Real L[NLWAVES] = {Real(0.0)};

    compute_L_lodi<cls_t>(iv,dir,q,drho,dun,dut1,dut2,dT,dY,L);

    // ==============================================================
    // 3. Transverse characteristic terms
    //
    // The existing routine evaluates its characteristic basis at
    // qbc(iv_face) and its centred tangential derivatives around
    // q(iv_inner). Passing iv for both indices and q for both fields
    // gives the required cell-centred evaluation.
    //
    // Do not enable this at edges or corners until their coupled
    // characteristic treatment is implemented.
    // ==============================================================

    // const bool use_transverse_here = use_transverse_terms && face_interior;
    const bool use_transverse_here = false; // temp

    const bool include_transverse = use_transverse_here && use_transverse_terms && (nscbc_type == 2 || nscbc_type == 3);
    Real Tchar[NLWAVES] = {Real(0.0)};

    if (include_transverse) {
        compute_transverse_characteristics<cls_t>(
            iv,       // local characteristic state
            iv,       // centre for tangential derivatives
            dir,
            dxinv,
            q,        // local primitive state
            q,        // neighbouring primitive states
            Tchar);
    }

    // ==============================================================
    // 4. Prescribe the incoming characteristic waves
    // ==============================================================

    if (nscbc_type == 1) {

        // Subsonic non-reflecting inflow with target relaxation.
        // The current implementation is LODI-only.
        apply_lodi_inflow_relaxation<cls_t>( iv, dir, side_sign, q, L, nscbc_parm);

    } else if (nscbc_type == 2) {

        // Purely non-reflecting outflow.
        if (include_transverse) {
            apply_lodi_outflow_transverse<cls_t>(side_sign,Tchar,L);
        } else {
            apply_lodi_outflow<cls_t>(iv,side_sign, q,L);
        }

    } else if (nscbc_type == 3) {

        // Non-reflecting outflow with pressure relaxation.
        if (include_transverse) {
            apply_lodi_outflow_pressure_relaxation_transverse<cls_t>(iv, side_sign, q, Tchar, L, nscbc_parm);
        } else {
            apply_lodi_outflow_pressure_relaxation<cls_t>(iv, side_sign, q, L, nscbc_parm);
        }
    } else {
        // This branch should be excluded by the host-side setup.
        // Retain the device-side failure for defensive checking.
        amrex::Abort( "NSCBC: nscbc_type must be 1, 2, or 3");
    }

    // ==============================================================
    // 5. Form the complete characteristic balance A = L + T
    //
    // For an incoming acoustic wave, the routines above have already
    // included -T in the imposed value. Adding T here therefore gives
    // the desired imposed incoming balance.
    // ==============================================================

    Real A[NLWAVES] = {Real(0.0)};

    for (int n = 0; n < NLWAVES; ++n) {
        A[n] = L[n];
        if (include_transverse) {
            A[n] += transverse_scale*Tchar[n];
        }
    }

    // ==============================================================
    // 6. Reconstruct primitive-variable time derivatives
    // ==============================================================

    Real drhodt = Real(0.0);
    Real dudt   = Real(0.0);
    Real dvdt   = Real(0.0);
    Real dwdt   = Real(0.0);
    Real dTdt   = Real(0.0);

    Real dYdt[NUM_SPECIES] = {Real(0.0)};

    primitive_rhs_from_L<cls_t>(iv,dir,q,A,drhodt,dudt,dvdt,dwdt,dTdt, dYdt);

    // ==============================================================
    // 7. Convert primitive RHS to conservative RHS
    // ==============================================================

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

    // Density
    rhs_ghost(iv, cls_t::URHO) = drhodt;

    // Momentum
    rhs_ghost(iv, cls_t::UMX) = u*drhodt + rho*dudt;

#if AMREX_SPACEDIM >= 2
    rhs_ghost(iv, cls_t::UMY) = v*drhodt + rho*dvdt;
#endif

#if AMREX_SPACEDIM == 3
    rhs_ghost(iv, cls_t::UMZ) = w*drhodt + rho*dwdt;
#endif

    // Species densities
#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        const Real Y = q(iv, cls_t::QFS + ns);
        rhs_ghost(iv, cls_t::UFS + ns) = Y*drhodt + rho*dYdt[ns];
    }
#endif

    // ==============================================================
    // 8. Total-energy RHS
    //
    // This retains the calorically-perfect-gas reconstruction from
    // the current implementation:
    //
    //     rho E = rho(e + ke)
    //
    //     d(rho E)/dt
    //       = (e + ke) drho/dt
    //       + rho [cv dT/dt + d(ke)/dt].
    //
    // This should eventually be replaced by an EOS-specific helper
    // for variable-cp or reacting mixtures.
    // ==============================================================

    const Real p = q(iv, cls_t::QPRES);
    const Real T = q(iv, cls_t::QT);

    const Real gamma =q(iv, cls_t::QG);
    const Real e = p / (rho*(gamma - Real(1.0)));
    const Real cv = e/T;

    const Real ke = Real(0.5)*(u*u + v*v + w*w);

    const Real dkedt = u*dudt + v*dvdt + w*dwdt;

    rhs_ghost(iv, cls_t::UET) = (e + ke)*drhodt + rho*(cv*dTdt + dkedt);
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void physical_flux_from_qbc(
    amrex::IntVect const& iv_face,
    int dir,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Array4<amrex::Real> const& flux,
    cls_t const& cls) noexcept
{
    using amrex::Real;

    const Real rho = qbc(iv_face, cls_t::QRHO);
    const Real p   = qbc(iv_face, cls_t::QPRES);

    const Real u = qbc(iv_face, cls_t::QU);
#if AMREX_SPACEDIM >= 2
    const Real v = qbc(iv_face, cls_t::QV);
#else
    const Real v = Real(0.0);
#endif

#if AMREX_SPACEDIM == 3
    const Real w = qbc(iv_face, cls_t::QW);
#else
    const Real w = Real(0.0);
#endif

    Real vel[3] = {u, v, w};
    const Real un = vel[dir];

    Real Y[NUM_SPECIES] = {Real(0.0)};

#if NUM_SPECIES > 1
    Real sumY = Real(0.0);

    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        Y[ns] = amrex::max(
            qbc(iv_face, cls_t::QFS + ns),
            Real(0.0));
        sumY += Y[ns];
    }

    if (sumY > Real(0.0)) {
        const Real inv_sum = Real(1.0)/sumY;

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            Y[ns] *= inv_sum;
        }
    }
#else
    Y[0] = Real(1.0);
#endif

    Real etot;
    cls.RYP2E(rho, Y, p, etot);

    etot += Real(0.5)*(u*u + v*v + w*w);

    flux(iv_face, cls_t::URHO) = rho*un;

    flux(iv_face, cls_t::UMX) = rho*un*u;
    flux(iv_face, cls_t::UMY) = rho*un*v;
    flux(iv_face, cls_t::UMZ) = rho*un*w;

    flux(iv_face, cls_t::UMX + dir) += p;

    flux(iv_face, cls_t::UET) =
        un*(rho*etot + p);

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        flux(iv_face, cls_t::UFS + ns) =
            rho*un*Y[ns];
    }
#endif
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void physical_flux_from_cell_primitive(
    amrex::IntVect const& iv_inner,
    amrex::IntVect const& iv_face,
    int dir,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& flux,
    cls_t const& cls) noexcept
{
    using amrex::Real;

    const Real rho = q(iv_inner, cls_t::QRHO);
    const Real p   = q(iv_inner, cls_t::QPRES);

    const Real u = q(iv_inner, cls_t::QU);

#if AMREX_SPACEDIM >= 2
    const Real v = q(iv_inner, cls_t::QV);
#else
    const Real v = Real(0.0);
#endif

#if AMREX_SPACEDIM == 3
    const Real w = q(iv_inner, cls_t::QW);
#else
    const Real w = Real(0.0);
#endif

    Real vel[3] = {u, v, w};
    const Real un = vel[dir];

    Real Y[NUM_SPECIES] = {Real(0.0)};

#if NUM_SPECIES > 1
    Real sumY = Real(0.0);

    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        Y[ns] = amrex::max(
            q(iv_inner, cls_t::QFS + ns),
            Real(0.0));

        sumY += Y[ns];
    }

    if (sumY > Real(0.0)) {
        const Real inv_sum = Real(1.0) / sumY;

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            Y[ns] *= inv_sum;
        }
    } else {
        Y[0] = Real(1.0);
    }
#else
    Y[0] = Real(1.0);
#endif

    Real etot;

    /*
     * RYP2E returns specific internal energy.
     */
    cls.RYP2E(rho, Y, p, etot);

    etot += Real(0.5) * (u*u + v*v + w*w);

    flux(iv_face, cls_t::URHO) = rho * un;

    flux(iv_face, cls_t::UMX) = rho * un * u;
    flux(iv_face, cls_t::UMY) = rho * un * v;
    flux(iv_face, cls_t::UMZ) = rho * un * w;

    flux(iv_face, cls_t::UMX + dir) += p;

    flux(iv_face, cls_t::UET) =
        un * (rho * etot + p);

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        flux(iv_face, cls_t::UFS + ns) =
            rho * un * Y[ns];
    }
#endif
}
//---------------------------------

} // namespace nscbc

#endif
