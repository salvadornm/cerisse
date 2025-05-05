// This file contains on-the-fly statistics functions.

#include <AMReX_ParmParse.H>

#include "CNS.h"

using namespace amrex;

//----------------------------------------------------------------------------------
// function to initialise statistic arrays to 0, is called per level
void CNS::setupStats() {

  // amrex::Print( ) << " Initialise STATS " << std::endl; 
  // amrex::Print( ) << " NSTATS = " << PROB::ProbClosures::NSTAT << std::endl; 
  // amrex:: Print( ) << " level = " << level << std::endl;
    
  MultiFab &S_stats = get_new_data(Stats_Type);
  auto const &data_stats = S_stats.arrays();
  amrex::ParallelFor(
    S_stats, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      for (int n=0; n < PROB::ProbClosures::NSTAT; n++)  {
        data_stats[box_no](i, j, k, n) =  0; 
      }                 
    });
}
//----------------------------------------------------------------------------------
// functon to compute and store mean values <U> and mean square <UU>, called per level
void CNS::computeStats(){

  // create multifabs
  MultiFab &S_stats = get_new_data(Stats_Type);
  MultiFab &S       = get_new_data(State_Type);
  
  Real dt = parent->dtLevel(level);
  
  // time_stats += dt; // obsolete
  Real time_aux = time_stat_level[level];
  
  // amrex::Print( ) << " Compute  STATS " << std::endl; 
  // amrex:: Print( ) << " level = " << level << std::endl;
  // amrex::Print( ) << " dt= " << dt <<  std::endl; 
  // amrex::Print( ) << " time_stat= " << time_stat_level[level] <<  std::endl; 
  // amrex::Print( ) << " time_aux= " << time_aux <<  std::endl; 
  // amrex::Print( ) << " time_old= "  << time_stat_level[level]-dt <<  std::endl; 
  
  // closures
  const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
                        
  // storing stats 
  for (MFIter mfi(S, false); mfi.isValid(); ++mfi) {
    // create box
    const Box& bx = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(cls_h.NGHOST);
    // extract arrays
    Array4<Real> const& state_cons  = S.array(mfi);   
    Array4<Real> const& data_stats  = S_stats.array(mfi);   
    // create primitive array (including ghost)
    FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
    Array4<Real> const& prims= primf.array();
    // from state variable to primitives (Thermodynamic closure) also update ghost
    cls_h.cons2prims(mfi, state_cons, prims); 
    amrex::ParallelFor(bx, [=](int i, int j, int k) {       
      Real rho = prims(i,j,k,cls_h.QRHO);
      //Real o_rho = 1.0/rho;
      // reset array  (so is sum (u*dt)
      for (int n=0; n< cls_h.NSTAT; n++){
        data_stats(i, j, k, n) *= (time_aux - dt);    
      }                

      // store velocity mean and square  correlations
      if (cls_h.record_velocity > 0) { 
        Real vel[3]={0,0,0};         
        vel[0] = prims(i,j,k,cls_h.QU);vel[1] = prims(i,j,k,cls_h.QV);vel[2] = prims(i,j,k,cls_h.QW);
      
        for (int n=0; n< AMREX_SPACEDIM; n++){
          data_stats(i, j, k, n) += vel[n]*dt;            
        }                
        for (int n=0; n< AMREX_SPACEDIM ; n++) {
          data_stats(i, j, k, AMREX_SPACEDIM+n) += vel[n]*vel[n]*dt;            
        }       
#if (AMREX_SPACEDIM >= 2)          
        data_stats(i, j, k, 2*AMREX_SPACEDIM) += vel[0]*vel[1]*dt;
#endif
#if (AMREX_SPACEDIM == 3)
        data_stats(i, j, k, 2*AMREX_SPACEDIM+1) += vel[0]*vel[2]*dt;
        data_stats(i, j, k, 2*AMREX_SPACEDIM+2) += vel[1]*vel[2]*dt;
#endif
      }
    
      // store P,T,rho mean and square
      if (cls_h.record_PTrho > 0 ){
                       
        Real P = prims(i, j, k, cls_h.QPRES);
        Real T = prims(i, j, k, cls_h.QT);            
        data_stats(i, j, k, INDEX_THERM)   += P*dt;
        data_stats(i, j, k, INDEX_THERM+1) += T*dt;
        data_stats(i, j, k, INDEX_THERM+2) += rho*dt;
        data_stats(i, j, k, INDEX_THERM+3) += P*P*dt;
        data_stats(i, j, k, INDEX_THERM+4) += T*T*dt;
        data_stats(i, j, k, INDEX_THERM+5) += rho*rho*dt;    
      }
    
      // store mean values  
      for (int n=0; n< cls_h.NSTAT; n++){
        data_stats(i, j, k, n) /= time_aux;            
      }           
    }); 

  } 
  
}

//----------------------------------------------------------------------------------
// debugging function to access stat array
void CNS::debugStats(int is, int js, int ks) {
  
  MultiFab &S_stats = get_new_data(Stats_Type);
  auto const &data_stats = S_stats.arrays();

  printf(" *** DEBUG point >>>  i=%d j=%d k=%d \n",is,js,ks);
  amrex::ParallelFor(
    S_stats, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      
      if ((i==is) && (j==js))     {  
        printf(" box_no=%d \n ",box_no );
        for (int n=0; n< PROB::ProbClosures::NSTAT; n++) {      
        printf(" n=%d  DATAstat= %f \n ",n,data_stats[box_no](i, j, k, n) );
        }
      }         
    });    

  printf("  ***  \n ") ; 

}
