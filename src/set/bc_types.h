#ifndef BCTYPES_H_
#define BCTYPES_H_

#include <AMReX_FArrayBox.H>


// shortcuts for different manual bc to be used within bcnormal in prob.h

// to activate add: 

template <typename cls_t>
class manual_bc_t
{
  public:

    manual_bc_t() {}

    ~manual_bc_t() {}

  // .......................................................................//
  // \brief  fixed mass flow rate at given temperature anc composition
  //         pressure will adapt in ghost point 
  // \param  nx,ny,nz : normal face  pointing into domain
  // \param  rhoUfix  : fixed mass flow rate per unit area [kg / m2 s] 
  // \param  Tfix     : fixed Temperature  [K]
  // \param  Yfix     : fixed composition (mass fraction array)
  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_inlet_fixmassflow(
  const Real nx, const Real ny, Real nz, const cls_t* cls,
  const Real rhoUfix, const Real Tfix, const Real* Yfix,
  const Real Uinner[cls_t::NCONS], Real* Ughost )
  {
    // calculate primitive array
    Real Qinner[cls_t::NPRIM]={0.0};      
    cls->cons2prims_point(Uinner,Qinner);

    // pghost = pinner
    const Real P = Qinner[cls_t::QPRES];            

    // ensure mass flow rate is constant - no acoustic reflection at inlet      
    Ughost[cls_t::UMX] = nx*rhoUfix;
    Ughost[cls_t::UMY] = ny*rhoUfix;
    Ughost[cls_t::UMZ] = nz*rhoUfix;
        
    // adjust density ghost to match fix temperature with fix P (and Y)
#if NUM_SPECIES > 1       
    Real rho=0.0; cls->PYT2R(P,Yfix,Tfix,rho); // recalculate rho in case of multiple species
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Ughost[cls_t::UFS+n] = rho*Yfix[n];		
    } 
#else
    Real rho  = Qinner[cls_t::QRHO] * Qinner[cls_t::QT]/Tfix;  
    Ughost[cls_t::URHO] = rho;
#endif     
    // internal specific energy ghost point
    Real e_ext    = 0.0; cls->RYP2E(rho, Yfix, P, e_ext);      
    // kinetic energy ghost point
    Real rhoe_kin = 0.5*( Ughost[cls_t::UMX]*Ughost[cls_t::UMX] + 
                          Ughost[cls_t::UMY]*Ughost[cls_t::UMY] +
                          Ughost[cls_t::UMZ]*Ughost[cls_t::UMZ])/rho; 
    Ughost[cls_t::UET] = rho* e_ext + rhoe_kin;    
  }
  // .......................................................................//
  // \brief  fixed pressure
  // \param  nx,ny,nz: normal face  pointing into domain
  // \param  P0      : fix pressure [Pa]

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_fixP(
    const Real /*nx*/, const Real /*ny*/, Real /*nz*/, const cls_t* cls,
    const Real P0,
    const Real Uinner[cls_t::NCONS], Real* Ughost )
    {
    // calculate primitive array
    Real Qinner[cls_t::NPRIM]={0.0};      
    cls->cons2prims_point(Uinner,Qinner);

    // pghost = pfix
    const Real P = P0;      

    // T and species 0-gradient
    const Real T = Qinner[cls_t::QT];
    Real Y[NUM_SPECIES] = {1.0};
#if NUM_SPECIES > 1      
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = Qinner[cls_t::QFS+n];		
    } 
#endif
    // compute density and internal energy
    Real rho = 0.0; Real e_ext=0.0;
    cls->PYT2R(P,Y,T,rho); cls->RYP2E(rho,Y,P,e_ext);
  
    // assign rho
#if NUM_SPECIES > 1      
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Ughost[cls_t::UFS+n] = rho*Y[n];		
    } 
#else
    Ughost[cls_t::URHO] = rho;
#endif  
    //  assign momentum
    Ughost[cls_t::UMX] = rho*Qinner[cls_t::QU];
    Ughost[cls_t::UMY] = rho*Qinner[cls_t::QV];
    Ughost[cls_t::UMZ] = rho*Qinner[cls_t::QW];
  
    // kinetic energy ghost point
    Real rhoe_kin = 0.5*( Ughost[cls_t::UMX]*Ughost[cls_t::UMX] + 
                        Ughost[cls_t::UMY]*Ughost[cls_t::UMY] +
                        Ughost[cls_t::UMZ]*Ughost[cls_t::UMZ])/rho; 
    Ughost[cls_t::UET] = rho* e_ext + rhoe_kin;    
  }
  // .......................................................................//
  // \brief  subsonic outflow fixed pressure at outlet
  // \param  nx,ny,nz: normal face  pointing into domain
  // \param  P0      : outlet pressure
  // poor's man NSBC based on Whitfield et al. AIAA-84-1552(1984)
  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void bc_subsonic_outflow_fixP(
    const Real nx, const Real ny, Real nz, const cls_t* cls,
    const Real P0,
    const Real Uinner[cls_t::NCONS], Real* Ughost )
    {
	
    // calculate primitive array
    Real Qinner[cls_t::NPRIM]={0.0};      
    cls->cons2prims_point(Uinner,Qinner);

    // pghost = poutlet
    const Real P    = P0;   
    // inner P and rho
    const Real Pi   = Qinner[cls_t::QPRES];  // pressure
    const Real one_over_rho0 = 1.0/Qinner[cls_t::QRHO];   // rho
    const Real one_over_c0   = 1.0/Qinner[cls_t::QC];     // sound speed

    const Real rho = Qinner[cls_t::QRHO]+ (P -Pi)*one_over_c0*one_over_c0;

    // Pout > Pinside : if nx> 0  du>0   , if nx < 0  du<0
    const Real u = Qinner[cls_t::QU] - (P -Pi)*one_over_c0*nx*one_over_rho0;
    const Real v = Qinner[cls_t::QV] - (P -Pi)*one_over_c0*ny*one_over_rho0;
    const Real w = Qinner[cls_t::QW] - (P -Pi)*one_over_c0*nz*one_over_rho0;

    // convert back to conservative vars

    //  assign momentum
    Ughost[cls_t::UMX] = rho*u;
    Ughost[cls_t::UMY] = rho*v;
    Ughost[cls_t::UMZ] = rho*w;
  
    // assign rho
    Real Y[NUM_SPECIES] = {1.0};
#if NUM_SPECIES > 1      
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = Qinner[cls_t::QFS+n];
      Ughost[cls_t::UFS+n] = rho*Y[n];		      
    } 
#else
    Ughost[cls_t::URHO] = rho;
#endif  

    // compute and assign energy
    const Real e_kin = 0.5*(u*u+ v*v + w*w);
    Real e_ext=0.0; cls->RYP2E(rho,Y,P,e_ext);
  
    Ughost[cls_t::UET] = rho* e_ext + rho*e_kin;   


  }

};

#endif
