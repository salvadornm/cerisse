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
// compute one-sided derivatives
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

    // derivatives of deisty, velocity and temperature (and species)
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
// compute L-waves
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

//-----------------------------------------------------------------
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
//-----------------------------------------------------------------
// LODI-only pressure relaxation
// high-x outflow: LMINUS incoming
// low-x  outflow: LPLUS  incoming
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_outflow_pressure_relaxation (
    amrex::IntVect const& iv, int dir,
    int side_sign,          // low: +1, high: -1
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real* L, nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;
    amrex::ignore_unused(dir);

    const Real p = q(iv, cls_t::QPRES);
    const Real c = q(iv, cls_t::QC);

    const Real dp = p - nscbc_parm.Ptarget;

    // Lodato/Poinsot-style pressure relaxation.
    // Same parameter-passing style as apply_lodi_inflow_relaxation().
    // Thesis convention used here:
    //   L_in = sigma * (1-Mmax^2)/(2*c*Lchar) * (p-Ptarget)
    const Real relax = nscbc_parm.sigma
                       * (Real(1.0) - nscbc_parm.Mmax*nscbc_parm.Mmax)
                       / (Real(2.0) * c * nscbc_parm.Lchar)
                       * dp;

    if (side_sign > 0) {
        // low boundary outflow: incoming acoustic is L+
        L[LPLUS] = relax;
    } else {
        // high boundary outflow: incoming acoustic is L-
        L[LMINUS] = relax;
    }
}
//-----------------------------------------------------------------
template <typename cls_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real target_velocity_component (int dir,
    amrex::Real utarget, amrex::Real vtarget, amrex::Real wtarget)
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
// LODI-only subsonic non-reflective inflow with target relaxation.
// Rathore thesis §4.1.3: all waves except the upstream-travelling
// acoustic wave are incoming at a subsonic inflow.  In this LODI
// implementation there are no transverse corrections, so the incoming
// entropy, tangential, normal-acoustic and species waves are relaxed
// toward target T, velocity and composition.
//
// low-side inflow  with un > 0: outgoing acoustic is L-, incoming is L+.
// high-side inflow with un < 0: outgoing acoustic is L+, incoming is L-.
// The normal acoustic relation is written in the scaling used by
// primitive_rhs_from_L so that the target contribution gives
// d(un)/dt = -eta*(un-utarget), apart from the outgoing-wave coupling.
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void apply_lodi_inflow_relaxation (
    amrex::IntVect const& iv, int dir,
    int side_sign,          // low: +1, high: -1
    amrex::Array4<const amrex::Real> const& q,
    amrex::Real* L, nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;

    const Real rho = q(iv, cls_t::QRHO);
    const Real c   = q(iv, cls_t::QC);
    const Real T   = q(iv, cls_t::QT);
    const Real eta = nscbc_parm.eta;

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

    // Temperature/entropy wave.
    L[LENT] = eta * (T - nscbc_parm.Ttarget);

    // Tangential velocity waves.
#if AMREX_SPACEDIM >= 2
    if (t1 >= 0) {
        L[LTAN1] = eta * (q(iv, qvel<cls_t>(t1)) - u_target[t1]);
    }
#endif
#if AMREX_SPACEDIM == 3
    if (t2 >= 0) {
        L[LTAN2] = eta * (q(iv, qvel<cls_t>(t2)) - u_target[t2]);
    }
#endif

    // Normal velocity acoustic relaxation, preserving the outgoing acoustic wave.
    const Real un      = q(iv, qvel<cls_t>(dir));
    const Real un_targ = u_target[dir];
    const Real acoustic_delta = rho/c * eta * (un - un_targ);

    if (side_sign > 0) {
        // low boundary inflow: incoming acoustic is L+; L- exits the domain
        L[LPLUS] = L[LMINUS] + acoustic_delta;
    } else {
        // high boundary inflow: incoming acoustic is L-; L+ exits the domain
        L[LMINUS] = L[LPLUS] - acoustic_delta;
    }

#if NUM_SPECIES > 1
    // No species target array has been added to NSCBCParm yet.  Keep species
    // locally non-reflective until Ytarget support is added.
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        L[LSP + ns] = Real(0.0);
    }
#endif
}
//-----------------------------------------------------------------
// compute  drho/dt , du/dt, .. etc from L-waves
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

    const Real rho   = q(iv, cls_t::QRHO);
    const Real T     = q(iv, cls_t::QT);
    const Real p     = q(iv, cls_t::QPRES);
    const Real c     = q(iv, cls_t::QC);
    const Real gamma = q(iv, cls_t::QG); 

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
//-----------------------------------------------------------------
// add Lwave sto RHS in the Ghost points
template <typename cls_t, typename nscbc_parm_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void add_lodi_rhs_to_cons (
    amrex::IntVect const& iv,
    int dir,
    int side_sign,
    amrex::Real dxinv,
    cls_t const*,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<amrex::Real> const& rhs,
    bool second_order,
    int nscbc_type, nscbc_parm_t const& nscbc_parm)
{
    using amrex::Real;

    Real drho, dun, dut1, dut2, dT;
    Real dY[NUM_SPECIES] = {Real(0.0)};

    one_sided_deriv_prim<cls_t>(
        iv, dir, side_sign, dxinv, q,
        drho, dun, dut1, dut2, dT, dY, second_order);

    Real L[5 + NUM_SPECIES] = {Real(0.0)};

    compute_L_lodi<cls_t>(iv, dir, q, drho, dun, dut1, dut2, dT, dY, L);

    // select NSCBC type:
    //   1: subsonic non-reflecting inflow with target relaxation
    //   2: purely non-reflecting outflow
    //   3: subsonic outflow with pressure relaxation
    if (nscbc_type==1){
        apply_lodi_inflow_relaxation<cls_t>(iv, dir, side_sign, q, L, nscbc_parm);
    }
    else if (nscbc_type==2){ 
        apply_lodi_outflow<cls_t>(iv, side_sign, q, L);
    }
    else if (nscbc_type==3) {          
        apply_lodi_outflow_pressure_relaxation<cls_t>( iv, dir, side_sign, q, L,nscbc_parm);             
    }    
    else {
        amrex::Error("NSCBC:: wrong bc_type use 1, 2 or 3");
    }

    ///

    Real drhodt, dudt, dvdt, dwdt, dTdt;
    Real dYdt[NUM_SPECIES] = {Real(0.0)};

    primitive_rhs_from_L<cls_t>( iv, dir, q, L, drhodt, dudt, dvdt, dwdt, dTdt, dYdt);
 
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

    rhs(iv, cls_t::URHO) = drhodt;
    rhs(iv, cls_t::UMX ) = u*drhodt + rho*dudt;
#if AMREX_SPACEDIM >= 2
    rhs(iv, cls_t::UMY ) = v*drhodt + rho*dvdt;
#endif
#if AMREX_SPACEDIM == 3
    rhs(iv, cls_t::UMZ ) = w*drhodt + rho*dwdt;
#endif

#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        const Real Y = q(iv, cls_t::QFS + ns);
        rhs(iv, cls_t::UFS + ns) = Y*drhodt + rho*dYdt[ns];
    }
#endif

    // Ideal-gas LODI energy closure for first implementation.
    const Real P = q(iv, cls_t::QPRES);
    const Real T = q(iv, cls_t::QT);
    const Real c = q(iv, cls_t::QC);
    const Real gamma = q(iv, cls_t::QG); 
    const Real e = P / (rho*(gamma - Real(1.0)) );  //e perfect gas
    const Real cv = e/T;

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

    rhs(iv, cls_t::UET) = (e + ke)*drhodt + rho*(cv*dTdt + dke);
}

} // namespace nscbc

#endif
