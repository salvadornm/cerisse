#ifndef TFM_H_
#define TFM_H_

#include <CNSconstants.h>

//////////////////////////////// TFM TEMPLATE /////////////////////////////////
/**
 * \brief Template for calculating ATF/TFM related subroutines
 *  includes, sensor and wrinkling functions
 *  requires that an LES template exists ??
 */
template <typename param, typename idx_t>
class TFM_t {

 
  // input parameters
  amrex::Real F0=1.0;            // Thickening factor   for ATF
  amrex::Real SLO  = 2.0;        // laminar flame speed for wrinking function 
                                 // (default H2/air premixed flame at 1 atm and 298K)  ~ 2 m/s
  amrex::Real deltaLO  = 0.0001; // laminar flame thickness H2 stoichiomnetric   ~ 0.1 mm
  amrex::Real Ret = 10;          // Turbulent Reynolds number for wrinkling

  // other parameters
  int sensor_id = idx_t::QT;    // index of variable used for sensor (e.g. progress variable)  (default QT)
  amrex::Real c1= 2200.0;        // limits to define progress variable for sensor (e.g. c=(T-Tu)/(Tb-Tu) )  
  amrex::Real c0= 298.0;         // (default Tb=2200 and Tu=298)  Tb=c1  Tu=c0 
    
  // some constants for sensor function
  static constexpr amrex::Real cms = 0.28 ;       // wrinkling constant  for efficiency


  // derived
  amrex::Real beta =  1.0;


  public:
  AMREX_GPU_HOST_DEVICE
  TFM_t() {
     amrex::ParmParse pp("atf");
      if (!pp.query("thickening_factor", F0)) {
        amrex::Print() << " ATF: using no thickening factor (default 1) \n ";   
      }
      if (!pp.query("laminar_flame_speed", SLO)) {
        amrex::Print() << " ATF: using default laminar flame speed (2 m/s) \n ";     
      }
      if (!pp.query("Ret", Ret)) {
        amrex::Print() << " ATF: using default Ret = 10 \n ";     
      }
      // if (!pp.(query"cms", cms)) {
      //   amrex::Print() << " ATF: default cms  = 0.28 \n ";   
      // }
      if (!pp.query("laminar_flame_thickness", deltaLO)) {
        amrex::Print() << " ATF: using default laminar flame thickness (0.1 mm) \n ";     
      }

      beta = max(2.0*std::log(2.0)/(3.0*cms*(std::sqrt(Ret)-1.0 +1.e-8)),0.0); // > 0

  }

  AMREX_GPU_HOST_DEVICE
  ~TFM_t() {}

  /**
   * \brief calculates flame sensor based on a primitive variable for ATF/TFM
   * \param[in]  i,j,k        cell index 
   * \param[in]  q            primitive variables array 
   * \param[in]  id_SENSOR    index of sensor variable in q
   * \param[out] flamesensor  1: inside flame (max)  0:outside  (is clipped to 0 and 1)
   *  based on Durand and Polifke (2007) symmetric sensor formulation
  */  
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real flame_sensor(
  const int i, const int j, const int k, const Array4<const Real>& q) const
  {
    amrex::Real c = q(i, j, k, sensor_id);  c= (c - c0)/(c1-c0);
    // clip c between 0 and 1
    c = amrex::max(amrex::min(c, 1.0),0.0);    
    return( 16.0* c*c*(1.0-c)*(1.0-c));
  }  
  /**
   * \brief calculates dynamic thickening factor for ATF/TFM based on flame sensor
   * \param[in]  omega        flamesensor
   * \param[out] thickening   F0 (max)  inside flame  and 1: outside 
   * 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real thickening( amrex::Real omega) const
  {
    return( 1.0 + (F0-1.0)*omega); // only thicken in flame region (omega>0) 
  }

  /**
   * \brief calculates efficiency factor for ATF/TFM based on flame sensor
   * \param[in]  usgs       sub-grid velocity (LES should provide estimate)
   * \param[in]  Delta      filter width
   * \param[out] efficiency   efficiency factor > 1 
   * 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real efficiency(amrex::Real usgs, amrex::Real Delta) const
  {    
    Real w  = wrinkling(usgs, Delta, deltaLO);    // wrinkling factor laminar flame
    Real w1 = wrinkling(usgs, Delta, F0*deltaLO); // wrinkling factor thickened flame  (will be samaller than w)

    return(w/w1); // 
  }

  /**
   * \brief calculates wrinkling factor 
   * \param[in]  usgs         sub-grid velocity (LES should provide estimate)
   * \param[in]  Delta        filter width
   * \param[out] wrinkling    wrinkling factor > 1
   * 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real wrinkling (amrex::Real usgs, amrex::Real Delta, amrex::Real deltaL) const
  {
    Real Unorm = usgs/SLO; Real Lnorm = Delta/deltaL;
    Real wrink = 1.0 + beta*Unorm*gamma(Unorm, Lnorm); // wrinkling factor

    return( wrink); 
  }

  /**
   * \brief fitted Charlette function for wrinkling factor
   * \param[in]  U             non-dimensional sub-grid velocity (usgs/SLO)
   * \param[in]  L             non-dimensional filter width (Delta/deltaLO)
   * \param[out] Charlette     Charlette function for wrinkling factor > 0 to up to 5 for large U and L
   * 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real gamma (amrex::Real U, amrex::Real L) const
  {
    Real arg =  -1.2/std::pow(U,0.3);
    Real Charlette = (-0.15*std::exp(-0.15*L) - 0.25*std::exp(U) + 0.85*std::exp(arg))*std::pow(L,two_third );
    return(Charlette); 
  }


////////////////////////////////////////////////////////////////
};



////////////////////////////////////////////////////////////////////////////////
#endif

