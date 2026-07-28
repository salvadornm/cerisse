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

//Real pt_factor = 0.0; 

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
            // Low side : (-3 qb + 4 q1 - q2)/(2 dx)
            // High side: ( 3 qb - 4 q1 + q2)/(2 dx)
            return Real(side_sign) * dxinv *
                   (-Real(1.5)*qb + Real(2.0)*q1 - Real(0.5)*q2);
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

    L[LMINUS] = lam_m * o_2gamma *
        (drho - gamma*rho*dun/c + rho*dT/T);

    L[LENT] = lam_0 *
        (-(gamma - Real(1.0))*T*drho/rho + dT/gamma);

    L[LTAN1] = lam_0 * dut1;
    L[LTAN2] = lam_0 * dut2;

    L[LPLUS] = lam_p * o_2gamma *
        (drho + gamma*rho*dun/c + rho*dT/T);

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        L[LSP + ns] = lam_0 * dY[ns];
    }
#else
    amrex::ignore_unused(dY);
#endif
}

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
// LODI pressure-relaxed subsonic outflow.
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow_pressure_relaxation (
    amrex::IntVect const& iv_face,
    int dir,
    int side_sign,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real* L,
    nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;
    amrex::ignore_unused(dir);

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

    L[LENT] = eta *
        (qbc(iv_face, cls_t::QT) - nscbc_parm.Ttarget);

#if AMREX_SPACEDIM >= 2
    if (t1 >= 0) {
        L[LTAN1] = eta *
            (qbc(iv_face, qvel<cls_t>(t1)) - u_target[t1]);
    }
#endif
#if AMREX_SPACEDIM == 3
    if (t2 >= 0) {
        L[LTAN2] = eta *
            (qbc(iv_face, qvel<cls_t>(t2)) - u_target[t2]);
    }
#endif

    const Real un      = qbc(iv_face, qvel<cls_t>(dir));
    const Real un_targ = u_target[dir];
    const Real acoustic_delta = rho/c * eta * (un - un_targ);

    if (side_sign > 0) {
        L[LPLUS] = L[LMINUS] + acoustic_delta;
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

//-----------------------------------------------------------------
// Lodato/Rathore transverse correction evaluated on the face-centred Q_BC
// field. This requires one valid tangential neighbour on each side; therefore
// allocate/fill at least one tangential ghost face before enabling it.
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void transverse_rhs_lodato (
    amrex::IntVect const& iv_face,
    int dir,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Real& drhodt,
    amrex::Real& dudt,
    amrex::Real& dvdt,
    amrex::Real& dwdt,
    amrex::Real& dTdt,
    amrex::Real* dYdt,
    amrex::Real transverse_scale = amrex::Real(1.0))
{
    using amrex::Real;

    const Real rho   = qbc(iv_face, cls_t::QRHO);
    const Real T     = qbc(iv_face, cls_t::QT);
    const Real gamma = qbc(iv_face, cls_t::QG);

    Real vel[3] = {Real(0.0), Real(0.0), Real(0.0)};
    vel[0] = qbc(iv_face, cls_t::QU);
#if AMREX_SPACEDIM >= 2
    vel[1] = qbc(iv_face, cls_t::QV);
#endif
#if AMREX_SPACEDIM == 3
    vel[2] = qbc(iv_face, cls_t::QW);
#endif

    auto Dc = [&] AMREX_GPU_DEVICE (int tdir, int n) noexcept -> Real {
        const auto e = amrex::IntVect::TheDimensionVector(tdir);
        return Real(0.5) * dxinv[tdir] *
               (qbc(iv_face + e, n) - qbc(iv_face - e, n));
    };

    Real trho = Real(0.0);
    Real tT   = Real(0.0);
    Real tu[3] = {Real(0.0), Real(0.0), Real(0.0)};
    Real div_ut = Real(0.0);
#if NUM_SPECIES > 1
    Real tY[NUM_SPECIES] = {Real(0.0)};
#endif

    for (int tdir = 0; tdir < AMREX_SPACEDIM; ++tdir) {
        if (tdir == dir) continue;

        const Real ut = vel[tdir];
        const Real dpdx = Dc(tdir, cls_t::QPRES);

        trho += ut*Dc(tdir, cls_t::QRHO)
               + rho*Dc(tdir, qvel<cls_t>(tdir));
        tT += ut*Dc(tdir, cls_t::QT);
        div_ut += Dc(tdir, qvel<cls_t>(tdir));

        tu[0] += ut*Dc(tdir, cls_t::QU);
#if AMREX_SPACEDIM >= 2
        tu[1] += ut*Dc(tdir, cls_t::QV);
#endif
#if AMREX_SPACEDIM == 3
        tu[2] += ut*Dc(tdir, cls_t::QW);
#endif
        tu[tdir] += dpdx/rho;

#if NUM_SPECIES > 1
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            tY[ns] += ut*Dc(tdir, cls_t::QFS + ns);
        }
#endif
    }

    drhodt -= transverse_scale*trho;
    dudt   -= transverse_scale*tu[0];
#if AMREX_SPACEDIM >= 2
    dvdt   -= transverse_scale*tu[1];
#endif
#if AMREX_SPACEDIM == 3
    dwdt   -= transverse_scale*tu[2];
#endif
    dTdt   -= transverse_scale*
              (tT + (gamma - Real(1.0))*T*div_ut);

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        dYdt[ns] -= transverse_scale*tY[ns];
    }
#else
    amrex::ignore_unused(dYdt);
#endif
}

//-----------------------------------------------------------------
// Compute the conservative RHS of the face-centred U_BC state.
// qbc and rhs_bc are face-centred; q is the cell-centred interior primitive
// field. All local boundary quantities vary independently with iv_face.
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void add_lodi_rhs_to_cons (
    amrex::IntVect const& iv_face,
    amrex::IntVect const& iv_inner,
    int dir,
    int side_sign,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
    cls_t const*,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& rhs_bc,
    bool second_order,
    int nscbc_type,
    nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;

    Real drho, dun, dut1, dut2, dT;
    Real dY[NUM_SPECIES] = {Real(0.0)};

    one_sided_deriv_prim<cls_t>(
        iv_face, iv_inner, dir, side_sign, dxinv[dir], qbc, q,
        drho, dun, dut1, dut2, dT, dY, second_order);

    Real L[5 + NUM_SPECIES] = {Real(0.0)};

    compute_L_lodi<cls_t>(
        iv_face, dir, qbc,
        drho, dun, dut1, dut2, dT, dY, L);

    if (nscbc_type == 1) {
        apply_lodi_inflow_relaxation<cls_t>(
            iv_face, dir, side_sign, qbc, L, nscbc_parm);
    } else if (nscbc_type == 2) {
        apply_lodi_outflow<cls_t>(iv_face, side_sign, qbc, L);
    } else if (nscbc_type == 3) {
        apply_lodi_outflow_pressure_relaxation<cls_t>(
            iv_face, dir, side_sign, qbc, L, nscbc_parm);
    } else {
        amrex::Error("NSCBC: nscbc_type must be 1, 2, or 3");
    }

    Real drhodt, dudt, dvdt, dwdt, dTdt;
    Real dYdt[NUM_SPECIES] = {Real(0.0)};

    primitive_rhs_from_L<cls_t>(
        iv_face, dir, qbc, L,
        drhodt, dudt, dvdt, dwdt, dTdt, dYdt);

    // if (use_transverse_terms) {
    //     transverse_rhs_lodato<cls_t>(
    //         iv_face, dir, dxinv, qbc,
    //         drhodt, dudt, dvdt, dwdt, dTdt, dYdt,
    //         transverse_scale);
    // }

    const Real rho = qbc(iv_face, cls_t::QRHO);
    const Real u   = qbc(iv_face, cls_t::QU);
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

    rhs_bc(iv_face, cls_t::URHO) = drhodt;
    rhs_bc(iv_face, cls_t::UMX ) = u*drhodt + rho*dudt;
#if AMREX_SPACEDIM >= 2
    rhs_bc(iv_face, cls_t::UMY ) = v*drhodt + rho*dvdt;
#endif
#if AMREX_SPACEDIM == 3
    rhs_bc(iv_face, cls_t::UMZ ) = w*drhodt + rho*dwdt;
#endif

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        const Real Y = qbc(iv_face, cls_t::QFS + ns);
        rhs_bc(iv_face, cls_t::UFS + ns) = Y*drhodt + rho*dYdt[ns];
    }
#endif

    // Calorically-perfect-gas energy closure retained from the original
    // implementation. Replace this later with a closure/EOS-specific routine.
    const Real P     = qbc(iv_face, cls_t::QPRES);
    const Real T     = qbc(iv_face, cls_t::QT);
    const Real gamma = qbc(iv_face, cls_t::QG);
    const Real e     = P / (rho*(gamma - Real(1.0)));
    const Real cv    = e/T;

    const Real dke = u*dudt
#if AMREX_SPACEDIM >= 2
        + v*dvdt
#endif
#if AMREX_SPACEDIM == 3
        + w*dwdt
#endif
        ;

    const Real ke = Real(0.5)*(u*u + v*v + w*w);

    rhs_bc(iv_face, cls_t::UET) = (e + ke)*drhodt + rho*(cv*dTdt + dke);
}
//------------------------------------------------------------------------------
// 
//
//------------------------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hllc_primitive_flux(
    amrex::IntVect const& iv_face,
    amrex::IntVect const& iv_inner,
    int dir,
    int side_sign,
    amrex::Array4<const amrex::Real> const& qbc,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& flux,
    cls_t const& cls) noexcept
{
    using amrex::Real;

    // Momentum-component mapping for a face normal to dir.
    const int UM1 = cls_t::UMX + dir;

    int t1, t2;
    tangent_dirs<cls_t>(dir, t1, t2);
    const int UM2 = (t1 >= 0) ? cls_t::UMX + t1 : cls_t::UMY;
    const int UM3 = (t2 >= 0) ? cls_t::UMX + t2 : cls_t::UMZ;

    const int QU1 = qvel<cls_t>(dir);
    const int QU2 = (t1 >= 0) ? qvel<cls_t>(t1) : cls_t::QV;
    const int QU3 = (t2 >= 0) ? qvel<cls_t>(t2) : cls_t::QW;

    /*  side_sign:
     *   +1: low boundary
     *       left  state = QBC
     *       right state = interior
     *   -1: high boundary
     *       left  state = interior
     *       right state = QBC
     */
    const bool low_side = (side_sign > 0);

    Real rl, ul, pl, ut1l, ut2l, cl;
    Real rr, ur, pr, ut1r, ut2r, cr;

    Real yl[NUM_SPECIES] = {Real(0.0)};
    Real yr[NUM_SPECIES] = {Real(0.0)};

    if (low_side) {
        rl   = qbc(iv_face, cls_t::QRHO);
        ul   = qbc(iv_face, QU1);
        pl   = qbc(iv_face, cls_t::QPRES);
        ut1l = (t1 >= 0) ? qbc(iv_face, QU2) : Real(0.0);
        ut2l = (t2 >= 0) ? qbc(iv_face, QU3) : Real(0.0);
        cl   = qbc(iv_face, cls_t::QC);

        rr   = q(iv_inner, cls_t::QRHO);
        ur   = q(iv_inner, QU1);
        pr   = q(iv_inner, cls_t::QPRES);
        ut1r = (t1 >= 0) ? q(iv_inner, QU2) : Real(0.0);
        ut2r = (t2 >= 0) ? q(iv_inner, QU3) : Real(0.0);
        cr   = q(iv_inner, cls_t::QC);

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            yl[ns] = qbc(iv_face, cls_t::QFS + ns);
            yr[ns] = q(iv_inner, cls_t::QFS + ns);
        }

    } else {

        rl   = q(iv_inner, cls_t::QRHO);
        ul   = q(iv_inner, QU1);
        pl   = q(iv_inner, cls_t::QPRES);
        ut1l = (t1 >= 0) ? q(iv_inner, QU2) : Real(0.0);
        ut2l = (t2 >= 0) ? q(iv_inner, QU3) : Real(0.0);
        cl   = q(iv_inner, cls_t::QC);

        rr   = qbc(iv_face, cls_t::QRHO);
        ur   = qbc(iv_face, QU1);
        pr   = qbc(iv_face, cls_t::QPRES);
        ut1r = (t1 >= 0) ? qbc(iv_face, QU2) : Real(0.0);
        ut2r = (t2 >= 0) ? qbc(iv_face, QU3) : Real(0.0);
        cr   = qbc(iv_face, cls_t::QC);

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            yl[ns] = q(iv_inner, cls_t::QFS + ns);
            yr[ns] = qbc(iv_face, cls_t::QFS + ns);
        }
    }

#if NUM_SPECIES > 1
    auto normalize_Y = [] AMREX_GPU_DEVICE (Real* Y) noexcept
    {
        Real sumY = Real(0.0);
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            Y[ns] = amrex::max(Y[ns], Real(0.0));
            sumY += Y[ns];
        }
        if (sumY > Real(0.0)) {
            const Real inv_sum = Real(1.0)/sumY;
            for (int ns = 0; ns < NUM_SPECIES; ++ns) { Y[ns] *= inv_sum;}
        } else {
            Y[0] = Real(1.0);
            for (int ns = 1; ns < NUM_SPECIES; ++ns) { Y[ns] = Real(0.0);}
        }
    };
    normalize_Y(yl);
    normalize_Y(yr);
#endif

    Real el, er;

    cls.RYP2E(rl, yl, pl, el);
    cls.RYP2E(rr, yr, pr, er);

    el += Real(0.5) *(ul*ul + ut1l*ut1l + ut2l*ut2l);
    er += Real(0.5) * (ur*ur + ut1r*ut1r + ut2r*ut2r);

    // HLLC signal speeds
    Real sl = amrex::min(ul - cl, ur - cr);
    Real sr = amrex::max(ul + cl, ur + cr);

    const Real density_ratio =std::sqrt(rr/rl);

    const Real uroe = (ul + density_ratio*ur) / (Real(1.0) + density_ratio);
    const Real croe = (cl + density_ratio*cr) / (Real(1.0) + density_ratio);

    sl = amrex::min(sl, uroe - croe);
    sr = amrex::max(sr, uroe + croe);

    Real frho;
    Real fmn;
    Real fmt1;
    Real fmt2;
    Real fenergy;
    Real fspecies[NUM_SPECIES] = {Real(0.0)};

    // HLLC RS
    if (sl >= Real(0.0)) {
        frho    = rl*ul;
        fmn     = rl*ul*ul + pl;
        fmt1    = rl*ul*ut1l;
        fmt2    = rl*ul*ut2l;
        fenergy = ul*(rl*el + pl);

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            fspecies[ns] = rl*ul*yl[ns];
        }

    } else if (sr <= Real(0.0)) {

        frho    = rr*ur;
        fmn     = rr*ur*ur + pr;
        fmt1    = rr*ur*ut1r;
        fmt2    = rr*ur*ut2r;
        fenergy = ur*(rr*er + pr);

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            fspecies[ns] = rr*ur*yr[ns];
        }

    } else {

        const Real denominator = rl*(sl - ul) - rr*(sr - ur);

        const Real sstar = (pr - pl + rl*ul*(sl - ul) - rr*ur*(sr - ur)) / denominator;

        if (sstar >= Real(0.0)) {

            const Real rstar = rl*(sl - ul)/(sl - sstar);

            const Real estar = el + (sstar - ul) * (sstar + pl/(rl*(sl - ul)));

            frho = rl*ul      + sl*(rstar - rl);
            fmn  = rl*ul*ul   + pl +sl*(rstar*sstar - rl*ul);
            fmt1 = rl*ul*ut1l + sl*(rstar*ut1l - rl*ut1l);

            fmt2 = rl*ul*ut2l + sl*(rstar*ut2l - rl*ut2l);

            fenergy = ul*(rl*el + pl) + sl*(rstar*estar - rl*el);

            for (int ns = 0; ns < NUM_SPECIES; ++ns) {
                fspecies[ns] = rl*ul*yl[ns] + sl*(rstar*yl[ns] - rl*yl[ns]);
            }

        } else {

            const Real rstar = rr*(sr - ur)/(sr - sstar);

            const Real estar = er + (sstar - ur) * (sstar + pr/(rr*(sr - ur)));

            frho = rr*ur + sr*(rstar - rr);

            fmn = rr*ur*ur + pr + sr*(rstar*sstar - rr*ur);

            fmt1 = rr*ur*ut1r + sr*(rstar*ut1r - rr*ut1r);

            fmt2 = rr*ur*ut2r + sr*(rstar*ut2r - rr*ut2r);

            fenergy = ur*(rr*er + pr) + sr*(rstar*estar - rr*er);

            for (int ns = 0; ns < NUM_SPECIES; ++ns) {
                fspecies[ns] = rr*ur*yr[ns] + sr*(rstar*yr[ns] - rr*yr[ns]);
            }
        }
    }

    // Store the flux in global Cartesian conservative ordering.
    flux(iv_face, cls_t::URHO) = frho;
    flux(iv_face, UM1)         = fmn;

#if AMREX_SPACEDIM >= 2
    flux(iv_face, UM2) = fmt1;
#else
    amrex::ignore_unused(fmt1);
#endif

#if AMREX_SPACEDIM == 3
    flux(iv_face, UM3) = fmt2;
#else
    /*
     * Cerisse retains UMZ in 2D. For a 2D x/y calculation this component
     * should still be the advected z-momentum flux.
     */
    flux(iv_face, cls_t::UMZ) = low_side ? frho*qbc(iv_face, cls_t::QW) : frho*q(iv_inner, cls_t::QW);
#endif

    flux(iv_face, cls_t::UET) = fenergy;

    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        flux(iv_face, cls_t::UFS + ns) = fspecies[ns];
    }
}



//--

} // namespace nscbc

#endif
